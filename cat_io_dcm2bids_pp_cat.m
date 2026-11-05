function varargout = cat_io_dcm2bids_pp_cat(action,varargin)
%cat_io_dcm2bids_pp_cat. CAT preprocessing of the T1w scans of cat_io_dcm2bids.
%  Optional final preprocessing of the exported BIDS directory by the CAT
%  segmentation with surfaces (opts.cat.run, not by default), independent
%  of the basic preprocessing of cat_io_dcm2bids (SPM segmentation, DTI,
%  and rs-fMRI). The results are written in own derivatives directories
%  of the main CAT version (e.g. derivatives/CAT26), as
%  CAT changes over time and results of different versions should not be
%  mixed, i.e. a new CAT version starts a new derivatives directory and
%  the results of former versions are kept (e.g. for pooling studies of
%  one version):
%    (1) segmentation of one T1w scan per session (the one with the best
%        quality rating SQR of its QC, see cat_io_dcm2bids_qc) with the
%        default CAT settings and surfaces (cat.estwrite), and for subjects
%        with several sessions the longitudinal pipeline of CAT (cat.long,
%        derivatives/CAT26_long, see opts.cat.long)
%    (2) further surface measures of the central surfaces (gyrification GI,
%        sulcal depth SD, and surface area, cat.stools.surfextract) and the
%        resampling and smoothing of the thickness (opts.cat.fwhm(1)) and
%        of these measures (opts.cat.fwhm(2)) to the 32k mesh with merged
%        hemispheres (cat.stools.surfresamp)
%  Both steps are SPM batches that are saved in the derivatives directory
%  (catDCM2BIDS_batch_<step>_<date>.mat, i.e. they can be loaded into the
%  SPM batch editor, e.g. to rerun them with other settings) and that only
%  process scans without results (i.e. new sessions of further runs).
%
%  varargout = cat_io_dcm2bids_pp_cat(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_pp_cat('processBIDS',Pbids,opts,nproc)
%    S = cat_io_dcm2bids_pp_cat('sessionScans',Pbids)
%    matlabbatch = cat_io_dcm2bids_pp_cat('segmentBatch',P,long,Pbids,nproc)
%    matlabbatch = cat_io_dcm2bids_pp_cat('surfaceBatch',Pcentral,opts,nproc)
%    name = cat_io_dcm2bids_pp_cat('derivativesName',long)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_pp_anat, cat_long_main.

  switch action
    case {'processBIDS','sessionScans','segmentBatch','surfaceBatch','derivativesName'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp_cat:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function processBIDS(Pbids,opts,nproc)
%processBIDS. CAT segmentation (cross-sectional, and longitudinal for
%  subjects with several sessions, see opts.cat.long) of one T1w scan per
%  session of the BIDS directory Pbids (see sessionScans) without results
%  in the derivatives directory of the CAT version (see derivativesName),
%  and the further surface measures and the resampling of all central
%  surfaces without resampled thickness (see surfaceBatch), with nproc
%  parallel CAT processes (0 - in this process).
  fprintf('\nDCM2BIDS - CAT preprocessing (%s).\n',derivativesName(0));
  S = sessionScans(Pbids);
  if isempty(S), fprintf('  No T1w scans.\n'); return; end
  long = arrayfun(@(s) numel(s.scans) > 1,S);
  if opts.cat.long == 0, long(:) = false; end
  % (cross-sectional and longitudinal processing of the subjects)
  if opts.cat.long == 2, modes = {{S,0}}; else, modes = {{S(~long),0}}; end
  if opts.cat.long > 0 && any(long), modes{end+1} = {S(long),1}; end

  % (1) segmentation of the scans without results
  for mi = 1:numel(modes)
    [Sm,lm] = modes{mi}{:}; Pder = fullfile(Pbids,'derivatives',derivativesName(lm));
    if lm
      P = arrayfun(@(s) s.scans(:),Sm,'UniformOutput',false);
      P = P(cellfun(@(p) ~all(cellfun(@(f) hasResults(Pbids,Pder,f),p)),P));
    else
      P = [Sm.scans]'; P = P(~cellfun(@(f) hasResults(Pbids,Pder,f),P));
    end
    if isempty(P), continue; end
    fprintf('  CAT segmentation of %d %s (%s)\n',numel(P),plural(lm),Pder);
    matlabbatch = segmentBatch(P,lm,Pbids,nproc);
    runBatch(matlabbatch,Pder,'segment');
  end

  % (2) surface measures and resampling of the central surfaces without
  % resampled thickness (also of former runs that were interrupted)
  for lm = unique([0 any(long) * (opts.cat.long > 0)])
    Pder = fullfile(Pbids,'derivatives',derivativesName(lm));
    Pcentral = centralSurfaces(Pder);
    if isempty(Pcentral), continue; end
    fprintf('  CAT surface measures of %d scans (%s)\n',numel(Pcentral),Pder);
    matlabbatch = surfaceBatch(Pcentral,opts,nproc);
    runBatch(matlabbatch,Pder,'surface');
  end
  fprintf('DCM2BIDS - CAT preprocessing done.\n');
end
% =========================================================================
function S = sessionScans(Pbids)
%sessionScans. T1w scans of the BIDS directory Pbids, one per session (the
%  one with the best quality rating SQR of its QC JSON in the derivatives
%  of cat_io_dcm2bids, otherwise the first one), as structure of the
%  subjects with the sessions and scans (cell, sessions in alphabetical
%  order, i.e. by date for the date labels of cat_io_dcm2bids).
  S = struct('sub',{},'ses',{},'scans',{});
  D = dir(fullfile(Pbids,'sub-*','**','anat','*_T1w.nii*'));
  if isempty(D), return; end
  P = sort(fullfile({D.folder},{D.name}))';
  P = P(~contains(P,[filesep 'derivatives' filesep]));
  sub = regexp(P,'(sub-[^\\/_]+)','tokens','once'); sub = [sub{:}]';
  ses = regexp(P,'[\\/](ses-[^\\/_]+)[\\/]','tokens','once'); ses(cellfun('isempty',ses)) = {{''}}; ses = [ses{:}]';
  for su = unique(sub)'
    S(end+1).sub = su{1}; %#ok<AGROW>
    for se = unique(ses(strcmp(sub,su{1})))'
      Ps = P(strcmp(sub,su{1}) & strcmp(ses,se{1}));
      if numel(Ps) > 1
        q = cellfun(@(f) qcRating(Pbids,f),Ps); [~,bi] = min(q); Ps = Ps(bi);
      end
      S(end).ses{end+1} = se{1}; S(end).scans{end+1} = Ps{1};
    end
  end
end
% =========================================================================
function q = qcRating(Pbids,P)
%qcRating. Quality rating SQR of the scan P of the BIDS directory Pbids by
%  its QC JSON in the derivatives of cat_io_dcm2bids (Inf without).
  q   = inf;
  rel = strrep(fileparts(P),[Pbids filesep],'');
  Pqc = fullfile(Pbids,'derivatives','catDCM2BIDS',rel,['catDCM2BIDSqc_' regexprep(spm_file(P,'filename'),'\.nii(\.gz)?$','') '.json']);
  if exist(Pqc,'file')
    try %#ok<TRYNC>
      J = cat_io_json(Pqc);
      if isfield(J,'qualityrating') && isfield(J.qualityrating,'SQR') && isfinite(J.qualityrating.SQR), q = J.qualityrating.SQR; end
    end
  end
end
% =========================================================================
function matlabbatch = segmentBatch(P,long,Pbids,nproc)
%segmentBatch. SPM batch of the CAT segmentation with surfaces (default
%  settings of CAT) of the scans P (cell) or of the longitudinal pipeline
%  for long (P with the scans of each subject as cell), with the output in
%  the derivatives directory of the CAT version (see derivativesName,
%  relative to the BIDS root Pbids of the scans) and nproc parallel
%  processes.
  if long
    matlabbatch{1}.spm.tools.cat.long.datalong.subjects = reshape(P,1,[]);
    matlabbatch{1}.spm.tools.cat.long.nproc             = nproc;
    matlabbatch{1}.spm.tools.cat.long.output.surface    = 1;
    matlabbatch{1}.spm.tools.cat.long.output.BIDS.BIDSyes.BIDSfolder = fullfile('derivatives',derivativesName(1));
  else
    matlabbatch{1}.spm.tools.cat.estwrite.data          = reshape(P,[],1);
    matlabbatch{1}.spm.tools.cat.estwrite.nproc         = nproc;
    matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSyes.BIDSfolder = fullfile('derivatives',derivativesName(0));
  end
  if ~exist(fullfile(Pbids,'dataset_description.json'),'file')
    cat_io_cprintf('warn','  No dataset_description.json in %s (BIDS root of the CAT output).\n',Pbids);
  end
end
% =========================================================================
function matlabbatch = surfaceBatch(Pcentral,opts,nproc)
%surfaceBatch. SPM batch of the further surface measures of the left
%  central surfaces Pcentral (cell, both hemispheres are processed):
%  gyrification GI, sulcal depth SD, and surface area (default measures of
%  cat.stools.surfextract), and the resampling and smoothing to the 32k
%  mesh with merged hemispheres of the thickness (FWHM opts.cat.fwhm(1))
%  and of these measures (FWHM opts.cat.fwhm(2), area by the areal
%  interpolation of CAT) by the dependencies of the batch.
  Pcentral = reshape(Pcentral,[],1);
  Pthick   = cellfun(@(p) strrep(p,'lh.central.','lh.thickness.'),Pcentral,'UniformOutput',false);
  Pthick   = regexprep(Pthick,'\.gii$','');
  matlabbatch{1}.spm.tools.cat.stools.surfextract.data_surf = Pcentral;
  matlabbatch{1}.spm.tools.cat.stools.surfextract.nproc     = nproc;
  matlabbatch{1}.spm.tools.cat.stools.surfextract.lazy      = 1;
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.data_surf  = Pthick;
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.mesh32k    = 1;
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.fwhm_surf  = opts.cat.fwhm(1);
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.lazy       = 1;
  matlabbatch{2}.spm.tools.cat.stools.surfresamp.nproc      = nproc;
  % (dependencies of the measures of the initialized surfextract module)
  jid = cfg_util('initjob',matlabbatch(1)); [~,~,~,~,dep] = cfg_util('showjob',jid); cfg_util('deljob',jid);
  dep = dep{1}; sn = cell(1,numel(dep)); for di = 1:numel(dep), sn{di} = dep(di).sname; end
  for name = {'MNI gyrification','sulcal depth','surface area'}
    di = find(endsWith(sn,['Left ' name{1}]),1); if isempty(di), continue; end
    matlabbatch{end+1}.spm.tools.cat.stools.surfresamp = struct('data_surf',dep(di),'merge_hemi',1,'mesh32k',1, ...
      'fwhm_surf',opts.cat.fwhm(2),'lazy',1,'nproc',nproc); %#ok<AGROW>
  end
end
% =========================================================================
function name = derivativesName(long)
%derivativesName. Name of the derivatives directory of the main CAT 
%  version (e.g. CAT26 of CAT26.0.2, i.e. without minor version and 
%  revision), with the suffix _long for the longitudinal pipeline.
  name = regexp(cat_version,'^CAT\d+','match','once'); if isempty(name), name = cat_version; end
  if long, name = [name '_long']; end
end
% =========================================================================
function ok = hasResults(Pbids,Pder,P)
%hasResults. True if the scan P of the BIDS directory Pbids has the XML
%  report of the CAT segmentation (cat_<name>.xml) in the derivatives
%  directory Pder (in the subdirectory of the scan or in its subfolders).
  rel = strrep(fileparts(P),[Pbids filesep],'');
  nam = regexprep(spm_file(P,'filename'),'\.nii(\.gz)?$','');
  ok  = ~isempty(dir(fullfile(Pder,rel,'**',['cat_*' nam '.xml'])));
end
% =========================================================================
function Pcentral = centralSurfaces(Pder)
%centralSurfaces. Left central surfaces of the derivatives directory Pder
%  without resampled thickness (s<fwhm>.mesh.thickness.resampled_32k.*,
%  e.g. new scans or interrupted runs), without the average surfaces of
%  the longitudinal pipeline.
  Pcentral = {};
  if ~exist(Pder,'dir'), return; end
  D = dir(fullfile(Pder,'sub-*','**','lh.central.*.gii'));
  for di = 1:numel(D)
    nam = regexprep(D(di).name,'^lh\.central\.|\.gii$','');
    if startsWith(nam,'avg_'), continue; end
    R = dir(fullfile(D(di).folder,['s*.mesh.thickness.resampled_32k.' nam '.gii']));
    if isempty(R), Pcentral{end+1,1} = fullfile(D(di).folder,D(di).name); end %#ok<AGROW>
  end
end
% =========================================================================
function runBatch(matlabbatch,Pder,step)
%runBatch. Saves the SPM batch of the step in the derivatives directory
%  Pder (catDCM2BIDS_batch_<step>_<date>.mat, e.g. for the batch editor)
%  and runs it.
  if ~exist(Pder,'dir'), mkdir(Pder); end
  Pmb = fullfile(Pder,sprintf('catDCM2BIDS_batch_%s_%s.mat',step,char(datetime('now','Format','yyyyMMdd-HHmmss'))));
  save(Pmb,'matlabbatch');
  try
    spm_jobman('run',matlabbatch);
  catch e
    cat_io_cprintf('err','  CAT %s failed: %s (batch %s)\n',step,e.message,Pmb);
  end
end
% =========================================================================
function s = plural(long)
%plural. Unit of the CAT segmentation (scans or longitudinal subjects).
  if long, s = 'subjects (longitudinal)'; else, s = 'scans'; end
end
