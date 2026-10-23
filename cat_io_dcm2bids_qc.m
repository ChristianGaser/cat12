function [Pr,QM,Pqc] = cat_io_dcm2bids_qc(P, type, opts, Pprotocols)
%cat_io_dcm2bids_qc. Basic image quality control of cat_io_dcm2bids.
%  Quality measures of the image P (of the BIDS datatype type): 
%    BSM .. between scan motion in 4D data (average correction of the realignment)
%    WSM .. within scan motion in 4D data (image variance between scans)
%    ISR .. inhomogeneity to signal ratio
%    NSR .. noise to signal ratio
%    RES .. RMS of the voxel resolution
%  4D data is realigned (rp_*.txt with the parameters of each volume, DWI 
%  with separate realignments of the b0 and of the diffusion-weighted 
%  volumes), and the realigned volumes are resliced in memory (see 
%  realignedVolumes in cat_io_dcm2bids_helper), i.e. without a resliced 
%  copy (prefix r of former versions, removed). The measures are rated (marks 1-6 for
%  the [best worst] measures, see qualityRating) by the QC definition of the
%  protocol file Pprotocols (qc<protocol>.json), otherwise by the default 
%  definition of the datatype (opts.qcdefault, see cat_io_dcm2bids_defaults),
%  and saved as catDCM2BIDSqc_<image>.mat/json (with the used definition). 
%  The ratings are printed in the scan line of the command line output. 
%  The denoised average of 4D data (of the first 8 volumes) is written next
%  to the image (<image>_denoisedavg.nii, zipped for opts.gzipi, see
%  cat_io_dcm2bids for the export as derivative and cat_io_dcm2bids_render
%  for the rendering), whereas the noise maps (noise map, noise field, and
%  slice motion) are only used for the measures, as they were not reliable
%  enough as maps.
%
%  [Pr,QM,Pqc] = cat_io_dcm2bids_qc(P,type,opts,Pprotocols)
%
%  P          .. image
%  type       .. BIDS datatype (anat, func, dwi, fmap)
%  opts       .. options of cat_io_dcm2bids (gzipi, hrrealign)
%  Pprotocols .. protocol file of the scan (empty: print the measures)
%  Pr         .. the image P (the realignment is in rp_*.txt)
%  QM         .. quality measures
%  Pqc        .. QC file (mat)
%
%  See also cat_io_dcm2bids.

  FNQC = {'BSM' 'WSM' 'ISR' 'NSR' 'RES'}; 

  opts.MarkColor  = cat_io_colormaps('marks+',40); 
  col2mark = @(val) opts.MarkColor(min(size(opts.MarkColor,1)-3,max(1,floor( val/9.5 * ...
    size(opts.MarkColor,1)))),:); 

  if opts.gzipi
    Pqc = spm_file(spm_file(strrep(P,'.nii.gz','.nii'),'ext',''),'prefix','catDCM2BIDSqc_','ext','mat');
  else
    Pqc = spm_file(strrep(P,'.nii.gz','.nii'),'prefix','catDCM2BIDSqc_','ext','mat');
  end
  Pqcj = spm_file(Pqc,'ext','json');

  % DWI realigned by former versions (one realignment with the b0 and the
  % diffusion-weighted volumes as sessions, see the realignment below,
  % recognized by the moved header sidecar rp_*.mat) are processed again
  Pbase = regexprep(P,'\.nii(\.gz)?$','');
  if strcmp(type,'dwi') && exist(spm_file([Pbase '.mat'],'prefix','rp_'),'file')
    for Pold = {spm_file([Pbase '.mat'],'prefix','rp_'), spm_file([Pbase '.txt'],'prefix','rp_'), Pqc, Pqcj, ...
        spm_file([Pbase '.nii'],'prefix','r'), spm_file([Pbase '.nii.gz'],'prefix','r')}
      if exist(Pold{1},'file'), delete(Pold{1}); end
    end
  end

  % resliced copies of former versions (the realigned volumes are resliced
  % in memory, see realignedVolumes in cat_io_dcm2bids_helper)
  for Pold = {spm_file([Pbase '.nii'],'prefix','r'), spm_file([Pbase '.nii.gz'],'prefix','r')}
    if exist(Pold{1},'file'), delete(Pold{1}); end
  end

  % echo images with missing slices are measured by the merged echoes of
  % their series (complete but with mixed echo times, see conversionInfo in
  % cat_io_dcm2bids_helper), also if the QC is older than the merged echoes
  [~,~,Pmerged] = cat_io_dcm2bids_helper('conversionInfo',P);

  if exist(Pqc,'file') && ~opts.rerun
    Pr = P;
    load(Pqc,'QM'); QM.SQR = nan;
    % test if the structure fits
    QM0 = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);
    run = ~strcmp( char(join(sort(fieldnames(QM0)))) ,char(join(sort(fieldnames(QM)))) );
    if ~isempty(Pmerged), d = [dir(Pqc) dir(Pmerged)]; run = run || d(1).datenum < d(2).datenum; end
  else
    run = 1;
  end

  if run
    % gunzip 
    try
      P = cat_io_dcm2bids_helper('prepNii',P,opts,0);
    catch
      Pr = P; 
      QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);  
      cat_io_cprintf([0.5 0 0],'QC-failed\n');
      return
    end
    V = spm_vol(P); Pimg = P; % (image of the measures)
    if ~isempty(Pmerged) && isscalar(V)
      Pmerged = cat_io_dcm2bids_helper('prepNii',Pmerged,opts,0);
      V = spm_vol(char(Pmerged)); Pimg = char(Pmerged);
    end
    vx_vol = sqrt(sum(V(1).mat(1:3,1:3).^2));
    
    % realignment parameters of 4D data (rp_*.txt)
    Pr  = P; 
    Prp = spm_file(P,'prefix','rp_','ext','.txt');
   
    if V(1).dim(3) < 5
      % spectroscopy overview image
      QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);  
      cat_io_cprintf('blue','spectroscopy preview?\n'); 
      return
    end

%%
    if 0
      %%
      clear matlabbatch
      if 1
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj = struct('vol', {{P}}, 'resample', {{P}});
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg  = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm      = {fullfile(spm('dir'),'tpm','TPM.nii')};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg   = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.r        = [0 0 0.1 0.01 0.04];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm     = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp     = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb       = [-90 -90 -70; 90 120 110];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox      = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp   = 4;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix   = 'w';
      end

      spm_jobman('run',matlabbatch); % just for debugging
      %evalc('spm_jobman(''run'',matlabbatch);');  
    else
      
    end


    try
      Y = single(spm_read_vols(V));
    catch e
      if cat_io_contains(e.message,'Unknown datatype.')
        cat_io_cprintf('err','Cannot read Nifti.\n');
        Pr = P; 
        QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);  
        cat_io_cprintf([0.5 0 0],'QC-failed\n');
        return
      end
      %rmdir(fileparts(P),'Recursive',true);
      %error('cat_io_dcm2bids:runQC','Remove "%s" to reimport. Please rerun!\n',P); 
    end
    sig75 = nan(1,size(Y,4)); 
    for di = 1:size(Y,4)
      Ydi = Y(:,:,:,di); 
      sig75(di) = cat_stat_prctile(Ydi(:),75); 
      clear Ydi; 
    end 
    if strcmp(type,'dwi') 
      isepi = sig75 > mean(sig75); 
    end
  
    if numel(V) > 1 && (~exist(Prp,'file') || opts.rerun) % (only 4D data)
      %% realignment (SPM, register to mean) of all volumes (rp_*.txt with 
      % the parameters of all volumes relative to the first one, without 
      % changes of the header, i.e. without reslicing, as the realigned 
      % volumes are resliced in memory into the grid of the first volume, 
      % see realignedVolumes): DWI with separate realignments of the b0 volumes and of the
      % diffusion-weighted volumes (the first volume of each set unmoved), as
      % one realignment with both sets as sessions registers their first
      % volumes across the contrasts (least squares), which shifted the
      % diffusion-weighted volumes by several mm relative to the b0 volumes
      V = spm_vol(P);
      if strcmp(type,'dwi'), ses = {find(isepi), find(~isepi)}; else, ses = {1:numel(V)}; end
      if opts.hrrealign
        ef = struct('quality',0.95,'sep',1.5,'fwhm',1);
      else
        ef = struct('quality',0.8,'sep',4,'fwhm',2);
      end
      ef.rtm = 1; ef.wrap = [0 0 0]; ef.interp = 2; ef.weight = ''; ef.graphics = 0;
      for si = find(cellfun(@numel,ses) > 1)
        evalc('Vs = spm_realign(V(ses{si}),ef);'); V(ses{si}) = Vs;
      end
      rp = zeros(numel(V),6);
      for vi = 1:numel(V), q = spm_imatrix(V(vi).mat / V(1).mat); rp(vi,:) = q(1:6); end
      save(spm_file(P,'prefix','rp_','ext','.txt'),'rp','-ascii');
    end
  

    %% main QC estimation 
    %  - special options: 0-none, 1-run for QC, 2-keep
    %  - this needs further work to run quick or use data permanently 
    %  
    opts.sliceMotionCor = 1; 
    opts.biasCor        = 1; 
    opts.denoise        = 1; 
    opts.tlim           = min(size(Y,4),8); 
    Pmaps               = {}; 
    % realigned volumes (resliced in memory into the grid of the first one)
    Vr = cat_io_dcm2bids_helper('realignedVolumes',Pimg,Prp,opts);
    for run = 1%:1 + strcmp(type,'dwi')
      if strcmp(type,'dwi')
        % sub-set-wise correction
        if run == 1
          Yr = realignedData(Vr,find(~isepi));
        else
          Yr = realignedData(Vr,find(isepi));
        end
      else
        Yr = realignedData(Vr,1:numel(Vr));
      end
 
      % basic segmentation of object/background
      % Ym  .. mean image
      % Ybb .. boundary map
      % Yb  .. brain/head/object mask - non-noise area
      % Yw  .. bias map
      % Yg  .. gradient/edge/noise map
      % g0  .. threshold in Yg
      % s1  .. signal intensity in Ym (for GM-WM intensity)
      Ym  = real(single(cat_stat_nanmean(Yr,4)));
      bbs = max(1,floor(min(4,size(Ym)/4))); % boundary width (integer, also for few slices)
      Ybb = true(size(Ym)); Ybb(bbs(1):end-bbs(1),bbs(2):end-bbs(2),bbs(3):end-bbs(3)) = false; 
      if  nnz(Ym(:)<0) ./ numel(Ym)  < .3 %~strcmp(type,'fmap')  &&  (
        %% typical image with mostly positive values and high intensity object
        Yg  = cat_vol_grad(Ym) ./ Ym; % this function has issues with negative non-noise structures
        g0  = cat_stat_prctile(Yg(:),10) * 2;
        Yb  = ~cat_vol_morph(Yg > g0*2,'ldc',2); 
        g0  = cat_stat_prctile(Yg(~Yb(:)),10) * 2;
        s0  = cat_stat_prctile(Ym(Yg(:) < g0 & Ym(:) > cat_stat_prctile(Ym(:),80) & Ym(:) < cat_stat_prctile(Ym(:),95)),80); 
        Yo  = Yg < g0  &  Ym > s0*.4  &  Ym < s0*1.5;
        Yo2 = cat_vol_morph(Yg<g0/2 & Yo,'lo'); % WM (within the tissue, e.g. not the background of EPIs)
        % first correction 
        Yw  = cat_vol_smooth3X(cat_vol_approx( abs(Yr(:,:,:,1)) .* Yo),16); 
        Yw  = Yw .* wmScale(Ym,Yw,Yo2);
        % second correction (only with enough tissue of the expected intensity,
        % as cat_vol_approx of an empty image gives zeros, i.e., NaN measures)
        Ysc = Yo & Ym./Yw > .9 & Ym./Yw < 1.5; 
        if nnz(Ysc) > 0.1 * nnz(Yo)
          Yw  = cat_vol_smooth3X(cat_vol_approx( abs(Yr(:,:,:,1)) .* Ysc),4);
          Yw  = Yw .* wmScale(Ym,Yw,Yo2);
        end
        clear Ysc
        Yb  = cat_vol_morph(Yg .* (Ym./Yw) > g0,'ldc',2); 
        Ynr = cat_vol_localstat(Ym./Yw,Yb,1,4); nr = cat_stat_nanmean(Ynr(Yb(:)))*2; 
        Ywm = cat_vol_morph(Yg<g0*2 & (Ym./Yw)>.8-nr & (Ym./Yw)<1.2+nr,'ldo',0);
        s1  = cat_stat_prctile(Ym(Ywm(:)),90); 
      else
        %% typical field-map with positive and negative values
        %  - here the bias is the information (so no correction) surrounded 
        %    by heavy noise that defines the signal intensity 
        Yg  = cat_vol_localstat(Ym,true(size(Ym)),2,4) ./ ...
              cat_vol_localstat(cat_vol_smooth3X(abs(Ym),1),true(size(Ym)),1,4);
        g0  = cat_stat_prctile(Yg(:),10) * 2;
        Yb  = cat_vol_morph( cat_vol_morph(Yg < g0 & ~Ybb,'ldo',1), 'ldc', 4); 
        s1  = cat_stat_prctile(Ym(~Yb(:)),90);
        Yw  = ones(size(Yg)) * s1;
      end  
      clear Yg g0 Ybb 
      


      %% 4D data data evaluation
      if ~strcmp(type,'anat') %strcmp(type,'dwi') &&  strcmp(type,'func')
        if strcmp(type,'dwi') 
          % in case of dwi, we need to split between the EPI images and
          % direction weighted scans
          if run==1, epiids = find(~isepi); else, epiids = find(isepi); end
        else
          epiids = 1:opts.tlim;
        end

        %%
        Ymn = zeros(V(1).dim,'single');
        Yr  = realignedData(Vr,epiids);
        Yd  = zeros([size(Yr,1) size(Yr,2) size(Yr,3) min(size(Yr,4),opts.tlim)],'single');
        WSM = nan(1,size(Yr,4));
        for vi = 1:min(size(Yr,4),opts.tlim)
        % for each time-point / direction apply the general bias correction
          Ya  = Yr(:,:,:,vi) ./ Yw;
    

          % If the data was realigned, we can correct for slice-wise motion
          % artifacts and interpret this as within-slice motion (WSM).
          if opts.sliceMotionCor 
            Yw1 = Yw;
            for zi = 1:size(Yr,3)
              if zi == 1
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi+1:zi+1),3); 
              elseif zi == size(Yr,3)
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi-1:2:zi-1),3); 
              else
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi-1:2:zi+1),3); 
              end
              Ytmp = cat_vol_smooth3X( repmat(Ytmp,1,1,3), mean(4./vx_vol(1:2)) );
              Yw1(:,:,zi) = Ytmp(:,:,1);
            end
            Yas = Ya - Yw1 .* abs(Yw1).^.25;
            Yas(Ya==0) = 0; % apply defacing
            WSM(vi) = cat_stat_nanmean( (Yas(:) - Ya(:)).^2 ).^.5;
            if opts.sliceMotionCor > 2, Ya = Yas; end
          end
    

          % Denoising of a single slice to quantify the amount of noise 
          % in the difference image
          Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end
          Ymn = Ymn + Yas / min(size(Yr,4),opts.tlim); % (the corrected volumes are not written)
          %if run==1 % why only run 1?
          Yd(:,:,:,vi) = sqrt( (Ya-Yas).^2 * 2 ); % Rician noise
          %end
        end
        
        % temporary corrected volumes of former versions
        Pc = spm_file(P,'prefix','c'); 
        if exist(Pc,'file'), delete(Pc); end

        if run==1 % why only run 1?
          % noise map and its smooth noise field (NSR), which are not written
          % as they were not reliable enough as maps
          Yn   = mean(Yd(:,:,:,1:min(size(Yr,4),opts.tlim)), 4) .* (Ya~=0);
          Yns  = cat_vol_approx(cat_vol_median3(Yn)) .* (Ya~=0);

          %% write output
          Vn = V(1); Vn.dt(1) = 16;
          Vn.fname = spm_file(V(1).fname,'suffix','_denoisedavg');
          Vn.descrip = sprintf(...
            'denoise, bias-adjusted average image of %d input volums < %s', ...
            min(size(Yr,4),opts.tlim), V(1).descrip);
          spm_write_vol(Vn,Ymn);
          Pmaps{end+1} = Vn.fname; %#ok<AGROW>
        end
      end
    end
    

    % do measurements
    %Ym  = Ym ./ Yw; % bias corrected
    Ys  = cat_stat_nanstd(Yr,4) ./ Yw; 
    try
      Yss = cat_vol_approx(cat_vol_median3(Ys)); 
    catch
      Yss = cat_vol_approx(smooth3(Ys)); 
    end

    % get motion parameters
    Pm  = spm_file(P,'prefix','rp_','ext','.txt');
    if exist(Pm,'file'), rp = load(Pm); else, rp = NaN; end
  
    % final measures
    QM.BSM  = cat_stat_nanmean(cat_stat_nanstd(rp,1).^2).^.5; % average motion (between scan movement)
    QM.ISR  = cat_stat_nanstd(Yw(Yb(:))) ./ s1;   % homogeneity to signal rating
    if exist('Yns','var') && mean(Yns(:))~=0 % (strcmp(type,'dwi') || strcmp(type,'func')) && 
      QM.NSR = cat_stat_nanmean(Yns(Yb(:)));      % noise to signal rating based on the denoising
      QM.WSM = cat_stat_nanmean(WSM(1:min(numel(WSM),opts.tlim)).^2).^.5;      % within slice motion 
    

    elseif   strcmp(type,'anat')  &&  ( nnz(Ym(:)<0) ./ numel(Ym)  < .3 ) 
      [Ya,Ybr]  = cat_vol_resize({mean(Yr,4)./Yw,single(Yb)},'reduceV',vx_vol,2,32,'meanm');
      
      % measure variance in background and foreground
      Yg     = abs(cat_vol_grad(Ya+min(Ya(:)))) ./ abs(Ya+min(Ya(:))); 
      Ytis   = Yg<cat_stat_prctile(Yg(Ybr(:)>.5),50) & Ya>.5 & Ya<1.5 & Ybr>.5; 
      Ytis   = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ytis,'l',1),'lc',1),'de',1);
      Ytis   = Yg<cat_stat_prctile(Yg(Ytis(:)>.5),50) & Ya>.5 & Ya<1.5 & Ybr>.5; 
      Ytis   = cat_vol_morph(Ytis,'ldo',1);
      [Ygr,Ybgr,Ytisr] = cat_vol_resize({Ya,~Ybr,Ytis},'reduceV',1,2,32,'meanm');
      NSRbg  = cat_vol_localstat(Ygr,Ybgr>.9,2,4);  NSRbg  = cat_stat_nanmedian(NSRbg(Ybgr(:)>.9)); 
      NSRtis = cat_vol_localstat(Ygr,Ytisr>.9,2,4); NSRtis = cat_stat_nanmedian(NSRtis(Ytisr(:)>.9)); 
      QM.NSR = min([NSRbg,NSRtis]);       % noise to signal rating based on the approximated 
  
      %% estimate denoising difference (without noise maps, as they were not
      %  reliable enough as maps)
      Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end

      % average
      QM.NSR = QM.NSR; % , cat_stat_nanmean((Ya(:) - Yas(:)).^2).^.5 );       % noise to signal rating based on the approximated 
      QM.WSM = cat_stat_nanmean((Ya(:) - Yas(:)).^2).^.5; 
    else
      QM.NSR = cat_stat_nanmean(Yss(Yb(:)));       % noise to signal rating based on the approximated bias field
      if exist('WSM','var') 
        QM.WSM = cat_stat_nanmean(WSM(1:opts.tlim).^2).^.5;
      else
        QM.WSM = NaN; 
      end
    end
    QM.vx_vol = vx_vol; 
    QM.RES    = cat_stat_nanmean(QM.vx_vol.^2).^.5;   % RESolution rating
    QM.SQR    = nan; 
  end

  [QR,qcdef] = qualityRating(QM,Pprotocols,type,opts);
  printvals = isempty(fieldnames(rmfield(QR,'SQR'))) || isnan(QR.SQR); % (measures without any rating)
  for fni = 1:numel(FNQC)
    if (~printvals && isnan(QR.(FNQC{fni}))) || ...
       ( printvals && isnan(QM.(FNQC{fni})))
      fprintf('    -')
    else
      if printvals % original values
        fprintf('%5.2f',QM.(FNQC{fni}));
      else % ratings 
        cat_io_cprintf(col2mark(QR.(FNQC{fni})),'%5.1f',QR.(FNQC{fni}));
      end
    end
  end
  if isnan(QR.SQR)
    fprintf('    -')
  else
    cat_io_cprintf(col2mark(QR.SQR),'%5.1f',QR.SQR);
  end
  fprintf(' \n');

  save(Pqc,'QM');
  J = struct('qualitymeasures',QM,'qualityrating',QR,'qcdefinition',qcdef);
  if ~isempty(Pmerged) % (image of the measures)
    J.image = sprintf('%s (merged echoes with missing slices, mixed echo times)',spm_file(regexprep(char(Pmerged),'\.gz$',''),'filename'));
  end
  cat_io_json(Pqcj,J);

  Pr = cat_io_dcm2bids_helper('prepNiigz',P,opts); % repack (the image itself, no resliced copy)
  if ~isempty(Pmerged), cat_io_dcm2bids_helper('prepNiigz',Pmerged,opts); end
  if exist('Pmaps','var'), cat_io_dcm2bids_helper('prepNiigz',Pmaps,opts); end

end
% =========================================================================
function Y = realignedData(V,ids)
%realignedData. Volumes ids of the volumes V as 4D array (single), resliced
%  in memory with the realignment of 4D data (see realignedVolumes in 
%  cat_io_dcm2bids_helper).
  Y = zeros([V(1).dim numel(ids)],'single');
  for k = 1:numel(ids), Y(:,:,:,k) = cat_io_dcm2bids_helper('realignedVolume',V,ids(k)); end
end
% =========================================================================
function f = wmScale(Ym,Yw,Ywm)
%wmScale. Scaling of the bias field Yw to the intensity of the mean image 
%  Ym in the WM mask Ywm (1 if not finite, e.g. for an empty mask).
  f = mean(Ym(Ywm(:)) ./ Yw(Ywm(:))); 
  if ~isfinite(f) || f <= 0, f = 1; end
end
% =========================================================================
function [QR,qcdef] = qualityRating(QM,Pprotocols,type,opts)
%qualityRating. Ratings (marks) of the quality measures QM by the QC
%  definition of the protocol (qc<protocol>.json), otherwise by the default
%  definition of the BIDS datatype type (opts.qcdefault.<type>, or .other),
%  with the [best worst] measure of the marks 1 and 6 (linear, clipped to
%  0.5-10.5 as in cat_stat_marks) or a single value that has to fit exactly
%  (mark 1, otherwise 10.5). NaN skips a measure. SQR is the RMS of all 
%  ratings. qcdef names the used definition.
  Pqcprotocols = '';
  if ~isempty(Pprotocols), Pqcprotocols = spm_file(Pprotocols,'prefix','qc','ext','.json'); end
  if ~isempty(Pqcprotocols) && exist(Pqcprotocols,'file')
    QCP = cat_io_json(Pqcprotocols); qcdef = spm_file(Pqcprotocols,'filename');
  elseif isfield(opts,'qcdefault') && isfield(opts.qcdefault,char(type))
    QCP = opts.qcdefault.(char(type)); qcdef = sprintf('default (opts.qcdefault.%s)',char(type));
  elseif isfield(opts,'qcdefault') && isfield(opts.qcdefault,'other')
    QCP = opts.qcdefault.other; qcdef = 'default (opts.qcdefault.other)';
  else
    QCP = struct(); qcdef = '';
  end

  % prepare output
  FN = fieldnames(QM); 
  for fni = 1:numel(FN)
    QR.(FN{fni}) = nan;
  end

  % quality rating
  FN = fieldnames(QCP); 
  if isfield(QCP,'RES'), QCP.vx_vol = QCP.RES; FN = [FN; intersect(fieldnames(QM),'vx_vol')]; end
  for fni = 1:numel(FN)
    if isfield(QM,FN{fni}) && all(~(isnan(QCP.(FN{fni})))) && all(~isnan(QM.(FN{fni})))
      if numel(QCP.(FN{fni}))==2 && abs(diff(QCP.(FN{fni}))) > 0.001
        for i=1:numel(QM.(FN{fni}))
          QR.(FN{fni})(i) = max(0.5,min(10.5, (QM.(FN{fni})(i) - QCP.(FN{fni})(1) ) / ...
            ( QCP.(FN{fni})(2) - QCP.(FN{fni})(1) ) * 5 + 1));
        end
        QR.(FN{fni}) = mean( QR.(FN{fni}).^2 )^.5;  
      else
        QR.(FN{fni}) = 10.5 - 9.5*(abs(QM.(FN{fni}) - QCP.(FN{fni})(1))<.001);
      end
    else
      QR.(FN{fni}) = nan;
    end
  end
  
  % averaging
  if isfield(QR,'vx_vol'), QR0 = rmfield(QR,'vx_vol'); end
  if isempty(FN)
    QR.SQR = nan;
  else
    fc = 2;
    QR.SQR = min(10.5,max(0.5, cat_stat_nanmean( cell2mat(struct2cell(QR0)).^fc ).^(1/fc)));
  end
end
