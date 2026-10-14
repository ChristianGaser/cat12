function varargout = cat_io_dcm2bids_pp_anat(action,varargin)
%cat_io_dcm2bids_pp_anat. Anatomical and shared processing of cat_io_dcm2bids.
%  Anonymization by defacing (anonymize), basic preprocessing of anatomical
%  scans by the SPM segmentation (segmentanat), the affine registration of
%  each session to MNI space with the orientation rating of the session 
%  and the coverage rating of its scans (writeSessionAffines, e.g. for the
%  rendering, see cat_io_dcm2bids_render, and cat_io_dcm2bids_overlap), and
%  the test of the phase-encoding reversal of scans with opposite phase-
%  encoding direction, e.g. AP/PA field maps (checkPEreversal), and helper
%  functions of the DWI and rs-fMRI preprocessing (see cat_io_dcm2bids_pp_dwi
%  and cat_io_dcm2bids_pp_rsfmri): the segmentation of the session
%  (sessionAnat), its labels in the space of other scans (sampleAnat), the
%  rigid coregistration (coregImages), the distortion correction by SCOPE
%  (scopeVDM, unwarpMatrix), and the parameters of the preprocessing level
%  opts.pplevel (levelParameters, levelUpdate).
%
%  varargout = cat_io_dcm2bids_pp_anat(action,varargin)
%
%  Actions (see the help of the local functions):
%    Pout = cat_io_dcm2bids_pp_anat('anonymize',Pin,opts,datatype)
%    p = cat_io_dcm2bids_pp_anat('levelParameters',level)
%    run = cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts)
%    [Pm,Pc0,Pwc1,Pseg] = cat_io_dcm2bids_pp_anat('segmentanat',Pin,datatype,opts)
%    S = cat_io_dcm2bids_pp_anat('writeSessionAffines',Pdbdirpath,Pdbnii,datatype,suffix,BIDSpathd,sub,ses,opts,BIDSfile,Pprot)
%    R = cat_io_dcm2bids_pp_anat('checkPEreversal',Pdbnii,Vjson,datatype,BIDSpathd,BIDSfile,opts)
%    A = cat_io_dcm2bids_pp_anat('sessionAnat',Pdbnii,datatype,suffix)
%    L = cat_io_dcm2bids_pp_anat('sampleAnat',V,M,A,Patlas,opts)
%    M = cat_io_dcm2bids_pp_anat('coregImages',Ysrc,Vsrc,Yref,Vref)
%    [vdm,info] = cat_io_dcm2bids_pp_anat('scopeVDM',Ppair,vols,Vt,Yt,flip,TRT,opts)
%    W = cat_io_dcm2bids_pp_anat('unwarpMatrix',dim,vdm,ind)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_pp_dwi, cat_io_dcm2bids_pp_rsfmri.

  switch action
    case {'anonymize','levelParameters','levelUpdate','segmentanat','writeSessionAffines','checkPEreversal', ...
          'sessionAnat','sampleAnat','coregImages','scopeVDM','unwarpMatrix'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp_anat:unknownAction','Unknown action "%s".',action);
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
      Pmsk = char(spm_deface( struct( 'images' , Vm.fname )));
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
  catch e
    cat_io_cprintf('err','  Defacing failed (%s), the image is not defaced.\n',e.message);
    copyfile(Pin,regexprep(Pout,'\.gz$','')); % (zipped below)
    %%% mark as failed?
  end

  % zip the output and remove the unzipped copy of the raw input file
  Pout = cat_io_dcm2bids_helper('prepNiigz',Pout,opts);
  cat_io_dcm2bids_helper('prepNiigz',Pin,opts);

end
% =========================================================================
function p = levelParameters(level)
%levelParameters. Parameters of the preprocessing level (opts.pplevel), a
%  continuous level of the accuracy and processing time of the DTI and 
%  rs-fMRI measures (the segmentation of the anatomical scans is the same 
%  for all levels):
%    1 - full quality (reference): SPM default parameters of the distortion
%        correction (SCOPE) at the native resolution and all scans
%    2 - optimized: fast settings and a resolution limit of SCOPE 
%    3 - quick review: lower resolution of SCOPE and only the first scan of
%        each type of a session (e.g. one T1w, one DWI, and one rs-fMRI)
%    4 - preview: lowest resolution of SCOPE (rough maps)
%  The parameters are given for the levels 1 to 4 and interpolated between
%  them (log-linear, levels below 1 are level 1, levels above 4 are level 4):
%    scope.res    .. resolution limit of SCOPE in mm (see scopeVDM, 1.5 mm
%                    is the native resolution of typical DWI)
%    scope.fwhm   .. spatial scales of SCOPE in mm (SPM default [8 4 2 1 0] 
%                    for level 1, otherwise the scales [8 4 2 1] of at least 
%                    half the resolution, as the number of conjugate gradient
%                    iterations, i.e. the time, depends little on the 
%                    resolution)
%    scope.interp .. B-spline degree of SCOPE (SPM default 2 for level 1, 
%                    otherwise linear)
%    onescan      .. only the first scan of each type of a session
%
%  Evaluation with the Berlin pilot session of the DZPG (T1w, 2 DWI, and 
%  rs-fMRI scans, see 'pplevels' of cat_DCM2BIDS_selftest): processing 
%  time of the segmentation, DTI, and rs-fMRI, and the differences to 
%  level 1 of the global measures (WM FA, tSNR, and segregation) and the 
%  correlation of the FA map and the connectivity matrix with level 1: 
%
%    level            time   WM FA  r(FA map)  tSNR  segregation  r(conn)
%    1 full quality  26 min  0.352      1      38.3     0.142         1
%    1.5             17 min  +0.000   0.95     -0.2    +0.003       0.99
%    2 optimized      9 min  +0.000   0.93     +0.0    +0.003       0.99
%    2.5              7 min  +0.000   0.85     -0.1    +0.009       0.99
%    3 quick review   4 min  +0.000   0.82     +0.1    +0.019       0.99
%    3.5            2.5 min  +0.000   0.78     -0.3    +0.015       0.97
%    4 preview      2.5 min  +0.000   0.79     +0.1    +0.027       0.98
%
%  I.e., the global measures are the same for all levels (the segregation 
%  varies most), whereas the maps differ slightly. The time of the levels 1
%  and 2 is dominated by SCOPE (about 95% of the DTI and rs-fMRI time), 
%  whereas the segmentation (about 2 min) and further steps (about 0.5 min
%  per session) limit the speed-up of higher levels (no further speed-up 
%  above level 3.5). 
%
%  Not used (tested with this session): 
%   - Subsets of the volumes save almost no time (reading all volumes takes
%     only seconds, the realignment is part of the QC), but bias the 
%     measures: 13 of 64 DWI directions (farthest-point sampling, better 
%     than the first or random directions) increase the WM FA by 5%; the 
%     first 20% of the rs-fMRI (1.6 min) increase the within-network 
%     connectivity by 50% (connectivity matrix r = 0.64 with all volumes);
%     every 5th volume keeps the connectivity pattern (r = 0.90) but lowers
%     the segregation by 26% and the FD is wrong (motion of 4 s). 
%   - A larger sampling distance of the segmentation saves little time 
%     (102 s for 3 mm, 68 s for 6 mm) but changes the volumes (TIV -4%, CSF
%     -11% for 6 mm), and a lower resolution of the T1w (e.g. 1.5 mm, 30 s)
%     changes the GM and WM volumes by about 20% (partial volume effects).
  lev = [1   2   3   4];                                      % levels
  res = [1.5 2.5 4   6];                                      % SCOPE resolution limit in mm
  L   = min(4,max(1,level));
  p.level      = level;
  p.scope.res  = exp(interp1(lev,log(res),L));
  if L < 1.5
    p.scope.fwhm = [8 4 2 1 0]; p.scope.interp = 2;         % SPM default
  else                                                      % fast: linear, scales of at least half the resolution
    f = [8 4 2 1]; p.scope.fwhm = f(f >= p.scope.res / 2); p.scope.interp = 1;
  end
  p.onescan    = L >= 3;
end
% =========================================================================
function run = levelUpdate(Pjson,opts)
%levelUpdate. True if the result Pjson of the preprocessing (with the
%  preprocessing level in method.level) has to be (re)processed, i.e. if it
%  does not exist, for opts.rerun, or if it was processed with a higher
%  level (lower accuracy) than opts.pplevel (results of a lower level are
%  kept).
  run = ~exist(Pjson,'file') || opts.rerun;
  if ~run
    J   = cat_io_json(Pjson);
    run = ~isfield(J,'method') || ~isfield(J.method,'level') || J.method.level > opts.pplevel;
  end
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

    % (re)processing for another sampling distance of the segmentation (e.g.
    % of older versions), that is the SPM default for all preprocessing 
    % levels, as the volumes depend on it (e.g. TIV -4% and CSF -11% for 6 mm,
    % whereas the time is only 1.5 times lower)
    samp = 3;
    run  = ~exist(Pc0,'file') || ~exist(Pseg,'file') || opts.rerun;
    if ~run
      J   = cat_io_json(Pseg);
      run = ~isfield(J,'method') || ~isfield(J.method,'samp') || J.method.samp ~= samp;
    end
    if run
      stime = tic;

      % denoising
      if opts.denoise && ~exist( spm_file(Pin,'prefix',prefn), 'file')
        cat_vol_sanlm(struct('data', {{Pin}},'verb',0,'prefix',prefn));
      end


      % run SPM segmentation
      if opts.denoise, Pin2 = spm_file(Pin,'prefix',prefn); else, Pin2 = Pin; end
      SPMsegment( Pin2 ,opts, samp);


      % CAT QC 
      if opts.runqc > 1 
        qcversion = 'cat_vol_qa201901x';
        if opts.denoise, prefix = ['c0' prefn]; else, prefix = 'c0'; end
        Pin2 = cat_io_dcm2bids_helper('prepNii',spm_file({Pin},'prefix',prefix),opts,0);
        cat_vol_qa('p0',Pin2,Pin2,Pin2,'','',...
          struct('prefix',[qcversion '_'],'version',qcversion,'rerun',opts.rerun,'verb',0) );
      end
    
   
      % evaluate segmentation, i.e. the new (unzipped) label map, as the 
      % zipped one of a previous segmentation is only replaced at the end
      seg8 = load(Pmat);
      Pc0u = regexprep(Pc0,'\.gz$','');
      if exist(Pc0u,'file'), Pc0r = Pc0u; else, Pc0r = Pc0; end %#ok<NASGU> (evalc)
      try
        evalc('V = spm_vol( Pc0r );'); 
      catch
        SPMsegment( Pin2 ,opts, samp);
        evalc('V = spm_vol( Pc0u );'); 
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
      spmus.method    = struct('samp',samp,'duration',round(toc(stime)));
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
function matlabbatch = SPMsegment(Pfiles,opts,samp)
%SPMsegment. Run SPM segmentation for anatomical data with the sampling 
%  distance samp (see segmentanat).

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
  matlabbatch{mi}.spm.spatial.preproc.warp.samp    = samp;   % default = 3
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
function R = checkPEreversal(Pdbnii,Vjson,datatype,BIDSpathd,BIDSfile,opts)
%checkPEreversal. Phase-encoding reversal measure (PRM) of the scans of a 
%  session with an opposite phase-encoding direction (e.g. AP/PA spin-echo
%  field maps, the PA refscan of the DWI, or PA runs of the rs-fMRI): the 
%  correlation of the first volume of a scan with the first volume of its 
%  partner, i.e., the closest scan (series number) of the same datatype, 
%  image size, and echo time with the opposite phase-encoding direction 
%  (in the sidecar or by the name, e.g. AP and PA, as a scan that was not 
%  reversed has the same direction in the sidecar). 
%  The reversed distortions result in low correlations (0.4-0.7 in the DZPG
%  pilots), whereas scans that were not reversed by the protocol (e.g. by 
%  the combination of a 180 degree in-plane rotation and the inverted 
%  phase-encoding polarity) are almost identical (r > 0.9). The PRM and its
%  rating PRR (marks 1-6 by opts.reversal.PRM) are added to the QC JSON of
%  the scan in the derivatives (catDCM2BIDSqc_<BIDSfile>.json), and pairs 
%  with a PRR above 3 are printed. Single-band references (SBRef) are only 
%  compared with each other. 
%  R is a structure array with the scan and partner index, PRM, and PRR. 

  R    = struct('id',{},'partner',{},'PRM',{},'PRR',{}); 
  mark = @(x,s) max(0.5,min(10.5,(x - s(1)) / (s(2) - s(1)) * 5 + 1)); % as qualityRating in cat_io_dcm2bids_qc
  n    = numel(Pdbnii); 
  pe   = repmat({''},1,n); te = nan(1,n); snr = nan(1,n); dims = cell(1,n); ok = false(1,n); sb = ok; 
  nd   = repmat({''},1,n); % direction by name (AP, PA, LR, RL)
  for i = 1:n
    V = Vjson{i}; 
    if ~isstruct(V) || ~isfield(V,'PhaseEncodingDirection') || ~isfield(V,'Dimensions') || ...
       isempty(Pdbnii{i}) || ~exist(Pdbnii{i},'file') || isempty(BIDSfile{i}), continue; end
    if isfield(V,'ImageType') && any(ismember(upper(cellstr(V.ImageType)),{'P','PHASE'})), continue; end % magnitude only
    pe{i}   = V.PhaseEncodingDirection; 
    dims{i} = V.Dimensions(1:min(3,end)); 
    if isfield(V,'EchoTime'),     te(i)  = V.EchoTime; end
    if isfield(V,'SeriesNumber'), snr(i) = V.SeriesNumber; end
    if isfield(V,'SeriesDescription'), sb(i) = ~isempty(regexpi(V.SeriesDescription,'sbref','once')); end
    if isfield(V,'ProtocolName')
      t = regexpi(V.ProtocolName,'(?:^|[_ -])(AP|PA|LR|RL)(?:$|[_ -])','tokens','once'); 
      if ~isempty(t), nd{i} = upper(t{1}); end
    end
    ok(i)   = ~isempty(regexp(pe{i},'^[ijk]-?$','once')); 
  end
  Y  = cell(1,n); % first volumes (read once)
  ax = cellfun(@(x) x(1:min(1,end)),pe,'UniformOutput',false); % phase-encoding axis
  opp = @(a,b) ~isempty(a) && ~isempty(b) && any(strcmp([a b],{'APPA','PAAP','LRRL','RLLR'})); 
  for i = find(ok)
    % partner: same datatype, image size, echo time, single-band reference 
    % (or not), and phase-encoding axis with opposite direction (sidecar or name)
    p = find(ok & sb == sb(i) & strcmp(datatype,datatype{i}) & cellfun(@(d) isequal(d,dims{i}),dims) & ...
      abs(te - te(i)) <= 0.01*te(i) & strcmp(ax,ax{i}) & ...
      (cellfun(@numel,pe) ~= numel(pe{i}) | cellfun(@(x) opp(x,nd{i}),nd))); 
    if isempty(p), continue; end
    [~,pk] = min(abs(snr(p) - snr(i))); p = p(pk); 
    try
      for j = [i p]
        if isempty(Y{j}), Y{j} = cat_io_dcm2bids_helper('niiVolume',Pdbnii{j},1); end
      end
      % correlation within the brain and other high intensity voxels (upper 30%)
      yi = sort(Y{i}(:)); yp = sort(Y{p}(:)); 
      m  = Y{i} > yi(round(0.7*numel(yi))) | Y{p} > yp(round(0.7*numel(yp))); 
      c  = corrcoef(Y{i}(m),Y{p}(m)); 
    catch e
      cat_io_cprintf('err','  phase-encoding reversal test failed for "%s": %s\n',Pdbnii{i},e.message); 
      continue
    end
    R(end+1) = struct('id',i,'partner',p,'PRM',c(1,2),'PRR',mark(c(1,2),opts.reversal.PRM)); %#ok<AGROW>
    Pqcj = fullfile(BIDSpathd{i},['catDCM2BIDSqc_' spm_file(BIDSfile{i},'basename') '.json']); 
    if exist(Pqcj,'file'), J = cat_io_json(Pqcj); else, J = struct(); end
    J.qualitymeasures.PRM = R(end).PRM; J.qualityrating.PRR = R(end).PRR; 
    if ~exist(BIDSpathd{i},'dir'), mkdir(BIDSpathd{i}); end
    cat_io_json(Pqcj,J); 
    if R(end).PRR > 3 && ~any([R(1:end-1).id] == p & [R(1:end-1).partner] == i) % each pair once
      nam = @(k) regexprep(spm_file(BIDSfile{k},'basename'),'^sub-[^_]+_(ses-[^_]+_)?',''); 
      cat_io_cprintf('err','%60s : %s (%s) and %s (%s) are similar (r=%0.2f, PRR %0.1f)\n', ...
        'phase-encoding not reversed?',nam(i),pe{i},nam(p),pe{p},R(end).PRM,R(end).PRR); 
    end
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
% =========================================================================
function A = sessionAnat(Pdbnii,datatype,suffix)
%sessionAnat. Segmentation of the best T1w scan of a session (otherwise of
%  another anatomical scan) by the log-likelihood of the SPM segmentation 
%  (seg8.mat, see segmentanat): bias-corrected image Pm (for the 
%  coregistration), label map Pc0 (CSF 1, GM 2, WM 3), and deformation 
%  field Py to MNI space. A is empty without segmentation.
  A    = struct([]);
  anat = find(strcmp(datatype,'anat') & ~cellfun('isempty',Pdbnii));
  t1w  = anat(strcmpi(suffix(anat),'T1w'));
  for ids = {t1w, setdiff(anat,t1w)}
    ll = -inf;
    for i = reshape(ids{1},1,[])
      Pseg8 = cat_vol_findfiles(fileparts(Pdbnii{i}),'*_seg8.mat');
      for si = 1:numel(Pseg8)
        [pth,nam] = fileparts(Pseg8{si}); nam = regexprep(nam,'_seg8$','');
        F = struct('Pm',niiFile(pth,['m' nam]),'Pc0',niiFile(pth,['c0' nam]),'Py',niiFile(pth,['y_' nam]));
        if any(cellfun('isempty',struct2cell(F))), continue; end
        try
          seg8 = load(Pseg8{si},'ll');
        catch
          continue
        end
        if seg8.ll > ll, A = F; ll = seg8.ll; end
      end
    end
    if ~isempty(A), return; end
  end
end
% =========================================================================
function P = niiFile(pth,nam)
%niiFile. NIfTI nam (.nii.gz or .nii) in the directory pth ('' if missing).
  P = '';
  for ext = {'.nii.gz','.nii'}
    if exist(fullfile(pth,[nam ext{1}]),'file'), P = fullfile(pth,[nam ext{1}]); return; end
  end
end
% =========================================================================
function L = sampleAnat(V,M,A,Patlas,opts)
%sampleAnat. Anatomical labels in the grid V (dim and mat, e.g. of the
%  undistorted mean of an EPI) by the transformation M from the world space
%  of V to the world space of the anatomical scan A (see sessionAnat and 
%  coregImages): the tissue label c0 (CSF 1, GM 2, WM 3, linear 
%  interpolation) and the labels of the atlases Patlas (cell of atlases in 
%  MNI space, e.g. of CAT, nearest neighbor) by the deformation of the 
%  segmentation, i.e., without normalization of the images. The atlas 
%  labels are only given in the brain (c0 > 0.5).
  [x,y,z] = ndgrid(1:V.dim(1),1:V.dim(2),1:V.dim(3));
  W   = M * V.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))]; clear x y z
  P   = cat_io_dcm2bids_helper('prepNii',{A.Pc0; A.Py},opts,0);
  mni = ones(4,size(W,2));
  try
    Vc   = spm_vol(P{1}); X = Vc.mat \ W;
    L.c0 = reshape(spm_sample_vol(Vc,X(1,:),X(2,:),X(3,:),1),V.dim);
    N    = nifti(P{2}); X = N.mat \ W;
    for d = 1:3, mni(d,:) = spm_sample_vol(double(N.dat(:,:,:,1,d)),X(1,:),X(2,:),X(3,:),1); end
  catch e
    cat_io_dcm2bids_helper('prepNiigz',P,opts); rethrow(e);
  end
  cat_io_dcm2bids_helper('prepNiigz',P,opts);
  L.atlas = cell(size(Patlas));
  for ai = 1:numel(Patlas)
    Va = spm_vol(Patlas{ai}); X = Va.mat \ mni;
    L.atlas{ai} = reshape(round(spm_sample_vol(Va,X(1,:),X(2,:),X(3,:),0)),V.dim);
    L.atlas{ai}(L.c0 < 0.5) = 0;
  end
end
% =========================================================================
function M = coregImages(Ysrc,Vsrc,Yref,Vref)
%coregImages. Rigid coregistration (spm_coreg, normalized mutual 
%  information) of the image Ysrc (array with the mat of Vsrc) to the 
%  reference image Yref (array with the mat of Vref) with the 
%  transformation M from the world space of the source to the world space
%  of the reference.
  T = tempname; mkdir(T);
  spm_write_vol(struct('fname',fullfile(T,'src.nii'),'dim',size(Ysrc),'mat',Vsrc.mat,'dt',[16 0],'pinfo',[1 0 0]'),Ysrc);
  spm_write_vol(struct('fname',fullfile(T,'ref.nii'),'dim',size(Yref),'mat',Vref.mat,'dt',[16 0],'pinfo',[1 0 0]'),Yref);
  try
    evalc('x = spm_coreg(spm_vol(fullfile(T,''ref.nii'')),spm_vol(fullfile(T,''src.nii'')),struct(''cost_fun'',''nmi'',''graphics'',0));');
  catch e
    rmdir(T,'s'); rethrow(e);
  end
  rmdir(T,'s');
  M = inv(spm_matrix(x));
end
% =========================================================================
function [vdm,info] = scopeVDM(Ppair,vols,Vt,Yt,flip,TRT,opts)
%scopeVDM. Voxel displacement map (VDM, in voxels along the phase-encoding
%  axis j) of a target image (grid Vt with the distorted mean image Yt) by
%  SCOPE (susceptibility correction using opposite phase-encoding, spm_scope
%  of SPM25, a re-implementation of FSL topup) of the mean images of the 
%  volumes vols (all if empty) of the database images Ppair (e.g. the AP/PA
%  spin-echo field maps or the b0 volumes of DWI scans with opposite phase-
%  encoding direction, the realigned images of the QC if available). 
%  The VDM of the pair is estimated with the settings and the resolution
%  limit of the preprocessing level (see levelParameters, i.e. the SPM
%  default settings at the native resolution for level 1 and faster
%  settings at lower resolutions for higher levels) and cached in the
%  database (catDCM2BIDSscope_*.nii next to the first image, with the
%  resolution and settings in the name). It is mapped to the target by the rigid
%  coregistration of the pair image with the phase-encoding direction of 
%  the target (the first image, or the second one for flip) to the target 
%  and scaled by the ratio of the total readout times TRT = [target pair] 
%  (and inverted for flip), i.e., the undistorted target is sampled at 
%  Yt(i,j + vdm(i,j,k),k) (see unwarpMatrix). info describes the correction.
  vdm  = [];
  if ~exist('spm_scope','file'), addpath(fullfile(spm('dir'),'toolbox','Spatial')); end
  if ~exist('spm_noise_estimate','file'), addpath(fullfile(spm('dir'),'toolbox','Longitudinal')); end
  pname = @(P) spm_file(fileparts(P),'basename'); % series directory, e.g. snr-0007-SpinEchoFieldMap_AP

  % mean images of the pair 
  m = cell(1,2); Vp = cell(1,2);
  for k = 1:2
    P = spm_file(Ppair{k},'prefix','r'); if ~exist(P,'file'), P = Ppair{k}; end
    keep = exist(regexprep(P,'\.gz$',''),'file') > 0; % unzipped by the caller (e.g. the target itself)
    Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
    try
      V = spm_vol(Pu{1}); vi = vols{k}; if isempty(vi), vi = 1:numel(V); end
      m{k} = 0; for i = vi(:)', m{k} = m{k} + spm_read_vols(V(i)); end
      m{k} = m{k} / numel(vi); Vp{k} = V(1);
    catch e
      if ~keep, cat_io_dcm2bids_helper('prepNiigz',Pu,opts); end
      rethrow(e);
    end
    if ~keep, cat_io_dcm2bids_helper('prepNiigz',Pu,opts); end
  end
  if ~isequal(Vp{1}.dim,Vp{2}.dim), info = 'none (pair of different size)'; return; end

  % VDM of the pair (cached), at the resolution limit (isotropic voxels of
  % the same field of view), in voxels of the first image
  p    = levelParameters(opts.pplevel);
  vx   = sqrt(sum(Vp{1}.mat(1:3,1:3).^2));
  res  = max([vx p.scope.res]);
  stag = sprintf('%s_b%d',sprintf('%g',p.scope.fwhm),p.scope.interp); % settings, e.g. 8421_b1
  snr2 = regexp(pname(Ppair{2}),'^snr-(\d+)','tokens','once'); if isempty(snr2), snr2 = {'x'}; end
  Pvdm = fullfile(fileparts(Ppair{1}),sprintf('catDCM2BIDSscope_vdm_snr-%s_%0.1fmm_%s.nii',snr2{1},res,stag));
  if ~exist(Pvdm,'file') && ~exist([Pvdm '.gz'],'file') || opts.rerun
    T  = tempname; mkdir(T); 
    sc = res ./ vx; dimN = ceil(Vp{1}.dim ./ sc);
    MN = Vp{1}.mat * [diag(sc) (0.5 - sc(:)/2); 0 0 0 1];
    [x,y,z] = ndgrid(1:dimN(1),1:dimN(2),1:dimN(3)); X = Vp{1}.mat \ MN * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))];
    for k = 1:2 % (spm_scope realigns its input images, i.e. only temporary copies)
      Y = reshape(spm_sample_vol(m{k},X(1,:),X(2,:),X(3,:),1),dimN);
      spm_write_vol(struct('fname',fullfile(T,sprintf('pe%d.nii',k)),'dim',dimN,'mat',MN,'dt',[16 0],'pinfo',[1 0 0]'),Y);
    end
    try
      runScope(T,p.scope.fwhm,p.scope.interp);
      Vv = spm_vol(fullfile(T,'vdm5_pe1.nii'));
    catch e
      rmdir(T,'s'); rethrow(e);
    end
    [x,y,z] = ndgrid(1:Vp{1}.dim(1),1:Vp{1}.dim(2),1:Vp{1}.dim(3)); X = Vv.mat \ Vp{1}.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))];
    v = reshape(spm_sample_vol(Vv,X(1,:),X(2,:),X(3,:),1),Vp{1}.dim) * sc(2);
    rmdir(T,'s');
    if nnz(isfinite(v)) < 0.5 * numel(v) || ~any(abs(v(isfinite(v))) > 0.01) % failed estimation (no cache)
      info = 'none (SCOPE failed)'; return
    end
    v(~isfinite(v)) = 0;
    Vo = struct('fname',Pvdm,'dim',Vp{1}.dim,'mat',Vp{1}.mat,'dt',[16 0],'pinfo',[1 0 0]', ...
      'descrip',sprintf('SCOPE VDM of %s and %s (%0.1f mm, %s)',pname(Ppair{1}),pname(Ppair{2}),res,stag));
    spm_write_vol(Vo,v);
    cat_io_dcm2bids_helper('prepNiigz',{Pvdm},opts);
  else
    Pu = cat_io_dcm2bids_helper('prepNii',{niiFile(fileparts(Pvdm),spm_file(Pvdm,'basename'))},opts,0);
    v  = spm_read_vols(spm_vol(Pu{1}));
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
  end

  % mapping to the target (the phase-encoding axes have to be parallel)
  pe = 1 + flip;
  M  = coregImages(Yt,Vt,m{pe},Vp{pe});
  aj = M(1:3,1:3) * Vt.mat(1:3,2); ap = Vp{1}.mat(1:3,2);
  if abs(aj' * ap) / (norm(aj) * norm(ap)) < cosd(10), info = 'none (other phase-encoding axis)'; return; end
  [x,y,z] = ndgrid(1:Vt.dim(1),1:Vt.dim(2),1:Vt.dim(3)); X = Vp{1}.mat \ M * Vt.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))];
  r   = TRT(1) / TRT(2); if ~isfinite(r), r = 1; end    % ratio of the total readout times
  vdm = reshape(spm_sample_vol(v,X(1,:),X(2,:),X(3,:),1),Vt.dim) * r * (1 - 2*flip);
  if nnz(isfinite(vdm)) < 0.5 * numel(vdm), vdm = []; info = 'none (mapping to the scan failed)'; return; end
  vdm(~isfinite(vdm)) = 0;
  info = sprintf('SCOPE (%s and %s, %0.1f mm, %s)',pname(Ppair{1}),pname(Ppair{2}),res,stag);
end
% =========================================================================
function runScope(T,fwhm,interp) %#ok<INUSD> (evalc)
%runScope. SCOPE (spm_scope) of the images pe1.nii and pe2.nii of the 
%  directory T (VDM in vdm5_pe1.nii) with the spatial scales fwhm and the 
%  B-spline degree interp
%   (1) with the default boundary condition of spm_diffeo (0) that spm_scope
%       requires, whereas the SPM segmentation sets 1 (spm_preproc8), which 
%       results in wrong VDMs (e.g. 2 mm median error and WM FA as without 
%       correction), and 
%   (2) without visible figures, as spm_scope draws each iteration into the
%       SPM Graphics window (created if missing and raised if visible) and
%       its title into the current figure, i.e. here into a hidden figure 
%       (existing Graphics windows get another tag until the end).
  bc = spm_diffeo('boundary'); spm_diffeo('boundary',0);
  Fg = findall(allchild(0),'Flat','Tag','Graphics'); set(Fg,'Tag','Graphics (hidden by cat_io_dcm2bids_pp_anat)');
  F0 = get(0,'CurrentFigure');
  Fh = figure('Visible','off','Tag','Graphics','IntegerHandle','off','HandleVisibility','on');
  c  = onCleanup(@() restoreScope(bc,Fh,Fg,F0)); % (also on errors)
  evalc('spm_scope(fullfile(T,''pe1.nii''),fullfile(T,''pe2.nii''),fwhm,[0 10 100],interp,1,''vdm5_'',T);');
end
% =========================================================================
function restoreScope(bc,Fh,Fg,F0)
%restoreScope. Boundary condition of spm_diffeo, figures, and current 
%  figure before runScope.
  spm_diffeo('boundary',bc);
  if isvalid(Fh), delete(Fh); end
  set(Fg(isvalid(Fg)),'Tag','Graphics');
  if ~isempty(F0) && isvalid(F0), set(0,'CurrentFigure',F0); end
end
% =========================================================================
function W = unwarpMatrix(dim,vdm,ind)
%unwarpMatrix. Sparse matrix W that samples a distorted volume Y at the 
%  undistorted positions of the voxels ind by the voxel displacement map 
%  vdm along the phase-encoding axis j (linear interpolation, see 
%  scopeVDM, no correction for an empty vdm), i.e., the undistorted values
%  of the voxels ind are W * Y(:).
  ind = ind(:); n = numel(ind);
  [i,j,k] = ind2sub(dim,ind);
  if isempty(vdm), p = j; else, p = j + double(vdm(ind)); end
  p  = min(dim(2),max(1,p)); j0 = min(dim(2)-1,floor(p)); w = p - j0;
  W  = sparse([1:n 1:n]',[sub2ind(dim,i,j0,k); sub2ind(dim,i,j0+1,k)],[1-w; w],n,prod(dim));
end
