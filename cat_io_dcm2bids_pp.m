function varargout = cat_io_dcm2bids_pp(action,varargin)
%cat_io_dcm2bids_pp. Image processing of cat_io_dcm2bids.
%  Anonymization by defacing (anonymize), basic preprocessing of anatomical
%  scans by the SPM segmentation (segmentanat), and the affine registration
%  of each session to MNI space with the orientation rating of the session
%  and the coverage rating of its scans (writeSessionAffines, e.g. for the
%  rendering, see cat_io_dcm2bids_render, and cat_io_dcm2bids_overlap).
%
%  varargout = cat_io_dcm2bids_pp(action,varargin)
%
%  Actions (see the help of the local functions):
%    Pout = cat_io_dcm2bids_pp('anonymize',Pin,opts,datatype)
%    [Pm,Pc0,Pwc1,Pseg] = cat_io_dcm2bids_pp('segmentanat',Pin,datatype,opts)
%    S = cat_io_dcm2bids_pp('writeSessionAffines',Pdbdirpath,Pdbnii,datatype,suffix,BIDSpathd,sub,ses,opts,BIDSfile,Pprot)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'anonymize','segmentanat','writeSessionAffines'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Pout = anonymize(Pin,opts,datatype) 
  
  if (~strcmp( datatype, 'anat') && opts.anonymize < 2) || opts.anonymize == 0
    Pout = Pin; 
    return
  end

  % prepare gz-input and setup prefix for output
  Pout = spm_file(Pin,'prefix','anon_');
  
  % be lazy if the file already exist
  if exist(Pout,'file') && ~opts.rerun; return; end

  % gunzip raw input file 
  Pin = cat_io_dcm2bids_helper('prepNii',Pin,opts,0);

  % defacing
  try
    V = spm_vol(Pin); 
    if numel(V)>1
      % create average to run defacing on this one
      Y  = spm_read_vols(V);
      Ym = cat_stat_nanmean(Y,4);
      
      %%%%%%%%%% this might bias the realignment!
      Vm = V(1); Vm.fname = spm_file(Vm.fname,'prefix','anon_');
      spm_write_vol(Vm,Ym);
      Pmsk = spm_deface( struct( 'images' , Vm.fname )); 
      Ymsk = spm_read_vols(spm_vol(Pmsk)) > 0; 
      delete(Vm.fname); delete(Pmsk);

      % apply masking
      Va   = V; 
      for vi = 1:numel(Va)
        Va(vi).fname = spm_file(Va(vi).fname,'prefix','anon_');
        Y(:,:,:,vi) = Y(:,:,:,vi) .* Ymsk;
        spm_write_vol(Va(vi),Y(:,:,:,vi));
      end

    else
      spm_deface( struct( 'images' , Pin ));
    end
  catch
    copyfile(Pin,Pout);
    %%% mark as failed?
  end

  % zip the output and remove the unzipped copy of the raw input file
  Pout = cat_io_dcm2bids_helper('prepNiigz',Pout,opts);
  cat_io_dcm2bids_helper('prepNiigz',Pin,opts); 

end
% =========================================================================
function [Pm,Pc0,Pwc1,Pseg] = segmentanat(Pin,datatype,opts)
  
  %% segment on anonymized data!
  Pm = ''; Pc0 = ''; Pwc1 = ''; Pseg = ''; 

  if ~opts.preprocessing, return; end

  if (strcmp( datatype, 'anat') || opts.preprocessing > 1) 
    if opts.denoise, prefn = 'sanlm_'; else, prefn = ''; end

    Pm   = spm_file(Pin,'prefix',['m' prefn]); 
    Pc0  = spm_file(Pin,'prefix',['c0' prefn]); 
    Pwc1 = spm_file(Pin,'prefix',['wc1' prefn]); 
    Pmat = spm_file(strrep(Pin,'.nii.gz','.nii'),'prefix',prefn,'suffix','_seg8','ext','mat'); 
    Pseg = spm_file(strrep(strrep(Pin,'.nii.gz','.nii'),[filesep 'anon_'], filesep), ...
      'prefix','catDCM2BIDSsegus_','ext','json');

    if ~exist( Pc0, 'file') || opts.rerun
      % denoising
      if opts.denoise && ~exist( spm_file(Pin,'prefix',prefn), 'file')
        cat_vol_sanlm(struct('data', {{Pin}},'verb',0,'prefix',prefn));
      end
  

      % run SPM segmentation  
      if opts.denoise, Pin2 = spm_file(Pin,'prefix',prefn); else, Pin2 = Pin; end
      SPMsegment( Pin2 ,opts);


      % CAT QC 
      if opts.runqc > 1 
        qcversion = 'cat_vol_qa201901x';
        if opts.denoise, prefix = ['c0' prefn]; else, prefix = 'c0'; end
        Pin2 = cat_io_dcm2bids_helper('prepNii',spm_file({Pin},'prefix',prefix),opts,0);
        cat_vol_qa('p0',Pin2,Pin2,Pin2,'','',...
          struct('prefix',[qcversion '_'],'version',qcversion,'rerun',opts.rerun,'verb',0) );
      end
    
   
      % evaluate segmentation  
      seg8 = load(Pmat);
      if strcmp(spm_file(Pc0,'ext'),'gz') && exist(Pc0,'file')
        try
          evalc('V = spm_vol( Pc0 );'); 
        catch
          SPMsegment( Pin2 ,opts);
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        end
      elseif exist(spm_file(Pc0,'ext',''),'file')
        try 
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        catch
          SPMsegment( Pin2 ,opts);
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        end
      end
      Y      = spm_read_vols(V); 
      vx_vol = sqrt(sum(V(1).mat(1:3,1:3).^2));

      % tissue volumes
      spmus.TIV   = nnz(Y(:)>0.5) * prod(vx_vol) / 1000;
      spmus.aGMV  = nnz(round(Y(:))==2) * prod(vx_vol) / 1000; 
      spmus.aWMV  = nnz(round(Y(:))==3) * prod(vx_vol) / 1000; 
      spmus.aCSFV = nnz(round(Y(:))==1) * prod(vx_vol) / 1000; 
      spmus.rGMV  = spmus.aGMV  ./ spmus.TIV; 
      spmus.rWMV  = spmus.aWMV  ./ spmus.TIV; 
      spmus.rCSFV = spmus.aCSFV ./ spmus.TIV; 
      
      % tissue intensities
      spmus.iGM   = seg8.mn(seg8.lkp==1) * seg8.mg(seg8.lkp==1);
      spmus.iWM   = max(seg8.mn(seg8.lkp==2));
      spmus.iCSF  = min(seg8.mn(seg8.lkp==3));

      % QC like parameters
      %   ll  = log-likelihood
      %   NCR = noise-to-contrast-ratio as minimum brain tissue variance 
      %         divided by the average tissue contrast
      spmus.qc.TPMll  = seg8.ll;
      spmus.qc.NCR    = min( shiftdim(seg8.vr(seg8.lkp(:)<4).^.5) ) ./ ...
                        mean( [ abs(spmus.iGM-spmus.iWM) abs(spmus.iGM-spmus.iCSF) ...
                                abs(spmus.iWM-spmus.iCSF)] * 2 * 3); 
      cat_io_json(Pseg,spmus);
      
    
      % zip the segmentation outputs and remove all unzipped copies (also 
      % of the input and denoised image, and the unzipped c0 of the CAT QC)
      if opts.gzipi
        prefixes = {'','c0','wc0','wc1','mwc1','l0','wl0','m','y_'};
        prefns   = unique({'',prefn}); 
        for pri = 1:numel(prefixes)
          for pfi = 1:numel(prefns)
            file = spm_file(strrep(Pin,'.nii.gz','.nii'),'prefix',[prefixes{pri} prefns{pfi}]); 
            if exist(file,'file'), cat_io_dcm2bids_helper('prepNiigz',file,opts); end
          end
        end
      end
    end
  end
end
% =========================================================================
function matlabbatch = SPMsegment(Pfiles,opts)
%SPMsegment. Run SPM segmentation for anatomical data

  if exist( spm_file(Pfiles, 'prefix', 'l0'), 'file')
    return
  end

  Pfiles = cellstr(cat_io_dcm2bids_helper('prepNii',Pfiles,opts,0));

  % TPM setting 
  if exist(fullfile(spm('dir'),'tpm','mni0R1p5_TPM7blr.nii'),'file')
    Ptpm = fullfile(spm('dir'),'tpm','mni0R1p5_TPM7blr.nii'); ngaus = [1 1 1 2 1 1 3];
  else
    Ptpm = fullfile(spm('dir'),'tpm','TPM.nii'); ngaus = [1 1 2 3 4 2];
  end
  Vtpm = spm_vol(Ptpm);

  % SPM segmentation 
  mi = 1; 
  matlabbatch{mi}.spm.spatial.preproc.channel.vols     = Pfiles; 
  matlabbatch{mi}.spm.spatial.preproc.channel.biasreg  = 0.001;
  % in general a bit more is better 
  matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm = 45;  % default = 60 
  matlabbatch{mi}.spm.spatial.preproc.channel.write    = [0 1];
  for ci = 1:numel(Vtpm)
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).tpm    = {sprintf('%s,%d',Ptpm,ci)};
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).ngaus  = ngaus(ci);
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).native = [ci<numel(Vtpm) 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).warped = [ci<numel(Vtpm) ci<2]; % unmod mod
  end
  % MRF remove fine anatomical details and it is better to live with random noise/artifacts 
  matlabbatch{mi}.spm.spatial.preproc.warp.mrf     = 0.1; % default = 1 
  matlabbatch{mi}.spm.spatial.preproc.warp.cleanup = 1;
  matlabbatch{mi}.spm.spatial.preproc.warp.reg     = [0 0.0001 0.05 0.005 0.02];
  % we are now in MNI space and this performs better
  matlabbatch{mi}.spm.spatial.preproc.warp.affreg  = 'subj'; % default = 'mni' 
  matlabbatch{mi}.spm.spatial.preproc.warp.fwhm    = 0; 
  matlabbatch{mi}.spm.spatial.preproc.warp.samp    = 6;      % default = 3 
  matlabbatch{mi}.spm.spatial.preproc.warp.write   = [0 1];  % backward forward


  % create label map for quick review
  for wi = 0:1 % subjects/template space
    for li = 0:1 
      mi = mi + 1; 
      for ci = 1:numel(Vtpm)-1
        if li==1, label='l'; else, label='c'; end
        if wi
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.images{ci}(1) = ...
            cfg_dep(sprintf('Segment: wc%d Images',ci), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{ci}, '.','wc', '()',{':'}));
        else
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.images{ci}(1) = ...
            cfg_dep(sprintf('Segment: c%d Images',ci), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
        end
      end
      if wi
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.prefix = ['\f\f\fw' label '0'];
      else
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.prefix = ['\f\f' label '0'];
      end
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.suffix          = '';
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.outdir          = {''};
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.BIDSdir         = '';
      if li == 0
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = 'i1*2 + i2*3 + i3*1';
      else 
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = 'round(i1)*1'; 
        for cii = 2:numel(Vtpm)-1
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = [ ...
            matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression , ...
            sprintf('+round(i%d)*%d',cii,cii)]; 
        end
      end
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.var             = struct('name', {}, 'value', {});
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.dmtx    = 0;
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.mask    = 0; % masking is not working here
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.interp  = 1 - li;
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.dtype   = 2; % 2
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.coreg   = 0;
    end
  end


  % cleanup 
  %  - remove classes ([w]c1-c#) that were useful for the label map but are
  %    not further required
  mi = mi + 1; 
  for fi = 1:numel(Pfiles)
    for ci = 1:numel(Vtpm)-1
      if fi==1 && ci == 1
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
          cfg_dep(sprintf('Segment: c%d Images',ci), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
      else
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(end+1) = ...
          cfg_dep(sprintf('Segment: c%d Images',ci), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
      end                
    end
    for ci = 2:numel(Vtpm)-1 % keep GM 
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(end+1) = ...
        cfg_dep( sprintf('Segment: wc%d Images',ci), ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','tiss', '()',{ci}, '.','wc', '()',{':'}));
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
  end

  % run batch
  %spm_jobman('run',matlabbatch); % just for debugging
  evalc('spm_jobman(''run'',matlabbatch);');  
end
% =========================================================================
function S = writeSessionAffines(Pdbdirpath,Pdbnii,datatype,suffix,BIDSpathd,sub,ses,opts,BIDSfile,Pprot)
%writeSessionAffines. Affine registration to MNI space for each session 
%  with the orientation rating of the session and the coverage rating of 
%  its scans.
%  The scans of a session share the same scanner space and therefore the 
%  same registration. The Affine of the segmentation (seg8.mat) with the 
%  best log-likelihood is used, otherwise spm_maff8 estimates it for the 
%  (best) anatomical scan. The Affine and its rigid part are saved as 
%    derivatives/sub-*/ses-*/catDCM2BIDSaffine_sub-*_ses-*.json
%  that maps world (scanner) coordinates to MNI coordinates. 
%  The orientation of the session is measured by the RMSE of the translation
%  (OTM, offset of the image origin from the MNI origin in mm) and of the 
%  rotation (ORM, in degree) of the registration with their ratings (OTR, 
%  ORR, marks 1-6 by opts.orientation) and their RMS as session orientation
%  rating (SOR), which are saved in this JSON (qualitymeasures/qualityrating).
%  The RMSE fits better than the largest translation/rotation, as a 
%  rotation around one axis is registered faster than the same rotation 
%  around all axes (whereas the registration itself was robust up to 100 mm
%  and 60 degree in tests). The coverage of each scan (BBM/BBR, see 
%  cat_io_dcm2bids_overlap) with the mask of its protocol Pprot 
%  (<protocol>_msk.nii if available) is added to its QC JSON in the 
%  derivatives (catDCM2BIDSqc_<BIDSfile>.json). The ratings are printed in 
%  one line per session. 
%  S is a structure array with the database session directory (dbses), the
%  Affine (empty without registration), and the SOR of each session. 

  S = struct('dbses',{},'Affine',{},'SOR',{}); 
  if ~exist('BIDSfile','var'), BIDSfile = cell(size(Pdbnii)); end
  if ~exist('Pprot','var'),    Pprot    = cell(size(Pdbnii)); end
  valid = find(~cellfun('isempty',BIDSpathd) & ~cellfun('isempty',Pdbdirpath)); 
  if isempty(valid), return; end
  MarkColor = cat_io_colormaps('marks+',40); 
  col2mark  = @(val) MarkColor(min(size(MarkColor,1)-3,max(1,floor( val/9.5 * size(MarkColor,1)))),:); 
  mark      = @(x,s) max(0.5,min(10.5,(x - s(1)) / (s(2) - s(1)) * 5 + 1)); % as qualityRating in cat_io_dcm2bids_qc
  
  % group by database session directory (sub-*/ses-*/snr-*)
  dbses     = cellfun(@fileparts, Pdbdirpath(valid), 'UniformOutput', false); 
  [udb,~,u] = unique(dbses); 
  for si = 1:numel(udb)
    ids = valid(u==si); 
    S(si).dbses = udb{si}; S(si).Affine = []; S(si).SOR = nan; 
    sesline = sprintf('session %s',spm_str_manip(spm_file(udb{si},'basename'),'l50')); 
    [Affine,method,src,ll] = getSessionAffine(udb{si}, Pdbnii(ids), datatype(ids), suffix(ids), opts); 
    if isempty(Affine)
      cat_io_cprintf([.5 .5 .5],'%60s : no registration (no anatomical scan)\n',sesline); 
      continue
    end
    imat  = spm_imatrix(Affine); 
    Rigid = spm_matrix(imat(1:6)); 

    % orientation measures (RMSE of the translation and rotation) and ratings
    QM = struct('OTM',sqrt(mean(imat(1:3).^2)),'ORM',sqrt(mean((imat(4:6)*180/pi).^2))); 
    QR = struct('OTR',mark(QM.OTM,opts.orientation.OTM),'ORR',mark(QM.ORM,opts.orientation.ORM)); 
    QR.SOR = sqrt(mean([QR.OTR QR.ORR].^2)); 
    S(si).Affine = Affine; S(si).SOR = QR.SOR; 

    % save in the derivatives session directory of each BIDS (sub)directory
    Psesd = cellfun(@fileparts, BIDSpathd(ids), 'UniformOutput', false); 
    [uses,ui] = unique(Psesd); 
    for di = 1:numel(uses)
      id = ids(ui(di)); 
      if isempty(ses{id}), sesname = sub{id}; else, sesname = [sub{id} '_' ses{id}]; end
      if ~exist(uses{di},'dir'), mkdir(uses{di}); end
      cat_io_json(fullfile(uses{di},['catDCM2BIDSaffine_' sesname '.json']), ...
        struct('Affine',Affine,'Rigid',Rigid,'method',method,'source',src,'ll',ll, ...
        'qualitymeasures',QM,'qualityrating',QR, ...
        'info','Affine and rigid transformation from world (scanner) to MNI coordinates.')); 
    end

    % coverage of each scan (added to its QC JSON in the derivatives)
    BBM = nan(size(ids)); BBR = BBM; 
    for ii = 1:numel(ids)
      id = ids(ii); 
      if isempty(Pdbnii{id}) || ~exist(Pdbnii{id},'file') || isempty(BIDSfile{id}), continue; end
      if isempty(Pprot{id}), Pmsk = ''; else, Pmsk = spm_file(Pprot{id},'suffix','_msk','ext','.nii'); end
      try
        [BBM(ii),BBR(ii)] = cat_io_dcm2bids_overlap('coverage',Pdbnii{id},Affine,Pmsk,opts); 
      catch e
        cat_io_cprintf('err','  coverage failed for "%s": %s\n',Pdbnii{id},e.message); 
        continue
      end
      Pqcj = fullfile(BIDSpathd{id},['catDCM2BIDSqc_' spm_file(BIDSfile{id},'basename') '.json']); 
      if exist(Pqcj,'file'), J = cat_io_json(Pqcj); else, J = struct(); end
      J.qualitymeasures.BBM = BBM(ii); J.qualityrating.BBR = BBR(ii); 
      if ~exist(BIDSpathd{id},'dir'), mkdir(BIDSpathd{id}); end
      cat_io_json(Pqcj,J); 
    end

    % session line: orientation measures and ratings, worst coverage
    fprintf('%60s : ',sesline); 
    cat_io_cprintf([0 0 0.5],'%-61s',sprintf('MNI %s: origin %0.1f mm, rotation %0.1f%s (RMSE)', ...
      method,QM.OTM,QM.ORM,char(176))); 
    for fn = {'OTR','ORR','SOR'}
      fprintf(' %s',fn{1}); cat_io_cprintf(col2mark(QR.(fn{1})),'%4.1f',QR.(fn{1})); 
    end
    if any(~isnan(BBR))
      [~,wi] = max(BBR); 
      fprintf('  coverage BBR'); cat_io_cprintf(col2mark(BBR(wi)),'%4.1f',BBR(wi)); 
      fprintf(' (%0.1f%% missing: %s)',BBM(wi),regexprep(spm_file(BIDSfile{ids(wi)},'basename'),'^sub-[^_]+_(ses-[^_]+_)?','')); 
    end
    fprintf('\n'); 
  end
end
% =========================================================================
function [Affine,method,src,ll] = getSessionAffine(Pdbses,Pnii,datatype,suffix,opts)
%getSessionAffine. Affine of the session: (1) seg8.mat with the best 
%  log-likelihood, (2) cached or new spm_maff8 registration of the best 
%  anatomical scan (T1w preferred), or (3) empty. 
  Affine = []; method = ''; src = ''; ll = nan; 

  % (1) segmentation
  Pseg8 = cat_vol_findfiles(Pdbses,'*_seg8.mat'); 
  for pi = 1:numel(Pseg8)
    try
      seg8 = load(Pseg8{pi},'Affine','ll'); 
      if isnan(ll) || seg8.ll > ll
        Affine = seg8.Affine; ll = seg8.ll; method = 'seg8'; 
        src    = Pseg8{pi}(numel(Pdbses)+2:end); 
      end
    end
  end
  if ~isempty(Affine), return; end

  % (2) spm_maff8 (cached in the database session directory)
  Pcache = fullfile(Pdbses,'catDCM2BIDSmaff8.mat'); 
  if exist(Pcache,'file') && ~opts.rerun
    load(Pcache,'Affine','method','src','ll'); 
    return
  end
  anat = find(strcmp(datatype,'anat') & ~cellfun('isempty',Pnii)); 
  if isempty(anat), return; end
  t1w  = anat(strcmpi(suffix(anat),'T1w')); 
  if ~isempty(t1w), anat = t1w; end
  P    = Pnii{anat(1)}; 
  try
    Pu = cat_io_dcm2bids_helper('prepNii',P,opts,0); 
    V  = spm_vol(char(Pu)); V = V(1); 
    [Affine,ll] = runMaff8(V); 
    method = 'maff8'; src = P(numel(Pdbses)+2:end); 
    save(Pcache,'Affine','method','src','ll'); 
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts); 
  catch e
    cat_io_cprintf('err','  spm_maff8 failed for "%s": %s\n',P,e.message); 
    Affine = []; 
  end
end
% =========================================================================
function [Affine,ll] = runMaff8(V)
%runMaff8. Affine registration to the TPM similar to SPM's segmentation 
%  (spm_preproc_run/run_affine): a coarse registration from the image 
%  origin and from the image center, and a refinement of the better one.
  persistent tpm
  if isempty(tpm)
    if exist(fullfile(spm('dir'),'tpm','mni0R1p5_TPM7blr.nii'),'file')
      Ptpm = fullfile(spm('dir'),'tpm','mni0R1p5_TPM7blr.nii'); 
    else
      Ptpm = fullfile(spm('dir'),'tpm','TPM.nii'); 
    end
    tpm = spm_load_priors8(Ptpm); 
  end
  samp = 6; fwhm = 0; affreg = 'mni'; 

  % coarse registration with the origin at the center of the field of view
  V1 = V; M = V1.mat; c = (V1.dim+1)/2; V1.mat(1:3,4) = -M(1:3,1:3)*c(:); 
  [Affine1,ll1] = spm_maff8(V1,8,(fwhm+1)*16,tpm,[],affreg); 
  Affine1       = Affine1*(V1.mat/M); 
  % coarse registration with the origin from the header
  [Affine2,ll2] = spm_maff8(V,8,(fwhm+1)*16,tpm,[],affreg); 
  if ll1 > ll2, Affine = Affine1; else, Affine = Affine2; end

  % refinement
  Affine      = spm_maff8(V,samp,(fwhm+1)*16,tpm,Affine,affreg); 
  [Affine,ll] = spm_maff8(V,samp, fwhm,      tpm,Affine,affreg); 
end
