function varargout = cat_io_dcm2bids_pp_dwi(action,varargin)
%cat_io_dcm2bids_pp_dwi. Basic preprocessing of diffusion-weighted images.
%  Diffusion tensor imaging (DTI) of each DWI scan of a session with at
%  least 12 diffusion-weighted volumes with b <= 1300 s/mm2 (DTI shell, 
%  otherwise the lowest shell, e.g. of single-shell data with b1500, see
%  dtiShell),
%  with the parameters of the preprocessing level opts.pplevel (see 
%  levelParameters in cat_io_dcm2bids_pp_anat, i.e. the resolution and settings
%  of SCOPE, and only the first scan of each b-value scheme for higher 
%  levels, e.g. the standard and the advanced protocol, whereas subsets of
%  the directions were tested but not used, as they bias the FA):
%    (1) realigned images of the QC (otherwise the original images),
%        denoised by sanlm for opts.denoise (each volume)
%    (2) distortion correction by the voxel displacement map of SCOPE of
%        the b0 volumes of the scan and of a DWI scan with opposite phase-
%        encoding direction (e.g. the PA refscan of the AP DWI, see scopeVDM
%        in cat_io_dcm2bids_pp_anat), i.e., the volumes are sampled at their
%        undistorted positions without resampling of the 4D data
%    (3) log-linear least squares fit of the tensor (b0 and DTI shell) in
%        the brain (T1w segmentation of the session by the rigid and 
%        non-linear coregistration of the b0 image to the T1w with the
%        Dice of the brain masks as test, see coregAnat and nonlinCoreg in
%        cat_io_dcm2bids_pp_anat, eroded by one voxel, otherwise b0
%        threshold, only voxels with signal in all volumes, i.e. not the
%        edge slices of the field of view that are missing in some volumes
%        after the realignment) with the
%        fractional anisotropy (FA), mean (MD), axial (AD), and radial (RD)
%        diffusivity (in mm2/s), the direction-encoded color map (DEC, the
%        absolute principal direction in the head (MNI) orientation 
%        weighted by the FA and a soft brain mask, i.e. red left-right, 
%        green posterior-anterior, and blue inferior-superior, see 
%        estimateDTI), and the fit residual of each volume (outlier volumes
%        with more than 3 SD of the residuals)
%    (4) global values: medians of FA, MD, AD, and RD in the white matter
%        (WM of the T1w segmentation with FA > 0.2), the peak width of the
%        MD in the white matter (95th - 5th percentile, similar to the PSMD
%        marker of small vessel disease but without skeleton), and the
%        median MD of the brain
%    (5) regional values: median FA and MD of the white matter regions of
%        the Mori atlas of CAT (by the deformation of the T1w segmentation)
%  The results are saved in the database (catDCM2BIDSdwi_<scan>.json with
%  the coregistration to the T1w, and the maps
%  catDCM2BIDSdwi_<FA|DEC>_<scan>.nii in the undistorted space of the scan
%  (also MD, AD, and RD for opts.dtimaps > 1), DEC with the 3 volumes red,
%  green, and blue, and the FA in MNI space wFA by the deformation of the
%  T1w segmentation, e.g. for voxel-based analyses, see warpToMNI, and the 
%  DEC map in the space of the T1w DECt1w by the non-linear coregistration,
%  e.g. as overlay of the T1w, see toT1w) and copied to the BIDS 
%  derivatives (the maps only for opts.preprocessing > 1, otherwise only in
%  the derivatives index of the database, rendered with 
%  their scan besides wFA and DECt1w, see cat_io_dcm2bids_render). The
%  segmentation is the one of the session (T1w, otherwise UNIT1, T2w, or
%  FLAIR, see anatScan in cat_io_dcm2bids_pp_anat). The gradient
%  directions are not rotated by the realignment (small effect for typical
%  motion of a few degree).
%  Mixed DTI (opts.dtimix, see processMix): sessions with several DWI scans
%  (e.g. the standard and the advanced protocol of the DZPG) get an 
%  additional DTI of the b0 and DTI shell volumes of all scans in the grid
%  of the scan with the most DTI volumes (further scans by the rigid 
%  coregistration of their FA maps and their own distortion correction, 
%  gradient directions rotated by the coregistration, S0 of each scan, e.g.
%  for different echo times), saved as catDCM2BIDSdwimix_<scan>.json and
%  catDCM2BIDSdwimix_<map>_<scan>.nii of this scan. The mixed values depend
%  on the b-values of the protocols (e.g. b850 of the advanced and b1200 of
%  the standard DZPG protocol: 3% lower WM FA and 3% higher MD in a pilot,
%  as the DTI depends on the b-value), i.e. the report uses the best scan.
%
%  varargout = cat_io_dcm2bids_pp_dwi(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_pp_dwi('processSession',Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
%    (with the mixed DTI of several DWI scans of the session, see processMix)
%
%  See also cat_io_dcm2bids_pp_anat, cat_io_dcm2bids_pp_rsfmri.

  switch action
    case {'processSession'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp_dwi:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function processSession(Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
%processSession. DTI of the DWI scans of a session (see the help of this
%  file) with the database images Pdbnii, sidecars Vjson, BIDS datatype and
%  suffix, derivatives directory BIDSpathd, and BIDS file of all scans of
%  the session.
  dwi = find(strcmp(datatype,'dwi') & strcmp(suffix,'dwi') & ~cellfun('isempty',Pdbnii) & ...
    ~cellfun('isempty',BIDSfile));
  if isempty(dwi), return; end
  b = cell(size(Pdbnii));
  for i = dwi
    Pb = [regexprep(Pdbnii{i},'\.nii(\.gz)?$','') '.bval'];
    if exist(Pb,'file'), b{i} = reshape(load(Pb),1,[]); end
  end
  dti = dwi(cellfun(@(x) isfinite(dtiShell(x)), b(dwi)));
  if isempty(dti), return; end
  A = cat_io_dcm2bids_pp_anat('sessionAnat',Pdbnii,datatype,suffix);
  if ~isempty(A) % BIDS file of the segmented scan (e.g. the T1w as underlay of the rendering), also of another session
    A.bids = ''; if ~isempty(BIDSfile{A.i}), A.bids = spm_file(BIDSfile{A.i},'basename'); end
    A.other = ~strcmp(fileparts(fileparts(Pdbnii{A.i})),fileparts(fileparts(Pdbnii{dti(1)}))); 
  end
  p = cat_io_dcm2bids_pp_anat('levelParameters',opts.pplevel);
  % for higher levels only the first scan of each b-value scheme (shells
  % rounded to 500 s/mm2, e.g. the standard and the advanced protocol of the
  % DZPG, but not a repeated scan of a protocol), others only if processed before
  key  = cellfun(@(x) mat2str(unique(500 * round(x(:)' / 500))),b(dti),'UniformOutput',false);
  [~,first] = unique(key,'stable'); comp = ~p.onescan | ismember(1:numel(dti),first);
  for ii = 1:numel(dti)
    i = dti(ii);
    try
      processScan(i,Pdbnii,Vjson,b,dwi,A,BIDSpathd{i},BIDSfile{i},opts,comp(ii));
    catch e
      cat_io_cprintf('err','%60s : DTI failed: %s (%s, line %d)\n',scanName(BIDSfile{i}), ...
        e.message,e.stack(1).name,e.stack(1).line);
    end
  end
  % mixed DTI of several DWI scans of the session (e.g. the standard and
  % the advanced protocol of the DZPG)
  if opts.dtimix && numel(dti) > 1
    try
      processMix(dti,Pdbnii,Vjson,b,dwi,A,BIDSpathd,BIDSfile,opts);
    catch e
      cat_io_cprintf('err','%60s : mixed DTI failed: %s (%s, line %d)\n',scanName(BIDSfile{dti(1)}), ...
        e.message,e.stack(1).name,e.stack(1).line);
    end
  end
end
% =========================================================================
function processScan(i,Pdbnii,Vjson,b,dwi,A,BIDSpathd,BIDSfile,opts,compute)
%processScan. DTI of the DWI scan i of the session (see processSession),
%  that is skipped if the results of the same or a lower preprocessing
%  level exist (besides opts.rerun, see levelUpdate in 
%  cat_io_dcm2bids_pp_anat). Without compute, only existing results are 
%  exported (results of older versions are processed again).
  [pth,nam] = fileparts(regexprep(Pdbnii{i},'\.gz$',''));
  pre   = 'catDCM2BIDSdwi_'; Pjson = fullfile(pth,[pre nam '.json']);
  [maps,req] = mapNames(A,opts);
  % results of older versions or without the maps of the options (e.g. the
  % wFA, or MD, AD, and RD for opts.dtimaps = 2) are processed again (also
  % without compute), where wFA needs the deformation of the segmentation
  % (and DECt1w also the non-linear coregistration)
  missing = exist(Pjson,'file') && ~all(cellfun(@(m) hasMap(pth,pre,m,nam),req));
  % (and results without T1w, if a T1w is available now, e.g. of another session)
  missing = missing || (exist(Pjson,'file') && ~isempty(A) && ~cat_io_dcm2bids_pp_anat('hasT1w',Pjson));
  compute = compute || (exist(Pjson,'file') && cat_io_dcm2bids_pp_anat('levelOutdated',Pjson)) || missing;
  if ~compute && ~exist(Pjson,'file'), return; end
  if compute && (missing || cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts))
    stime = tic;
    J = estimateScans(i,Pdbnii,Vjson,b,dwi,A,opts,pth,pre,nam);
    J.method.duration = round(toc(stime));
    cat_io_json(Pjson,J);
  else
    J = cat_io_json(Pjson);
  end
  exportResults(Pjson,pth,pre,nam,maps,BIDSpathd,BIDSfile,opts);
  G = J.global;
  fprintf('%60s : DTI WM FA %0.3f, MD %0.2f, AD %0.2f, RD %0.2f (10^-3 mm2/s), outlier volumes %d, %s\n', ...
    scanName(BIDSfile),G.FA_WM,G.MD_WM*1000,G.AD_WM*1000,G.RD_WM*1000,J.qc.outlier_volumes, ...
    J.method.distortion_correction);
end
% =========================================================================
function processMix(dti,Pdbnii,Vjson,b,dwi,A,BIDSpathd,BIDSfile,opts)
%processMix. Mixed DTI of the DWI scans dti of a session with DTI results
%  (at least two, e.g. the standard and the advanced protocol of the DZPG):
%  one tensor fit of the b0 and DTI shell volumes of all scans in the grid
%  of the scan with the most DTI volumes (see estimateDTI), saved as
%    catDCM2BIDSdwimix_<scan>.json and catDCM2BIDSdwimix_<map>_<scan>.nii
%  of this scan (and in its BIDS derivatives), processed again if one of
%  the DTI results of the scans is newer (or for a higher level or an older
%  version, see levelUpdate in cat_io_dcm2bids_pp_anat).
  pj  = @(k) fullfile(fileparts(Pdbnii{k}),['catDCM2BIDSdwi_' spm_file(regexprep(Pdbnii{k},'\.gz$',''),'basename') '.json']);
  dti = dti(arrayfun(@(k) exist(pj(k),'file') > 0,dti));
  if numel(dti) < 2, return; end
  nd  = cellfun(@(x) sum(x > 500 & x <= dtiShell(x)),b(dti)); [~,o] = sort(nd,'descend'); sc = dti(o);
  i   = sc(1); [pth,nam] = fileparts(regexprep(Pdbnii{i},'\.gz$',''));
  pre = 'catDCM2BIDSdwimix_'; Pjson = fullfile(pth,[pre nam '.json']);
  [maps,req] = mapNames(A,opts);
  fdate = @(f) subsref(dir(f),substruct('.','datenum'));
  if ~exist(Pjson,'file') || any(arrayfun(@(k) fdate(pj(k)) > fdate(Pjson),sc)) || ...
      ~all(cellfun(@(m) hasMap(pth,pre,m,nam),req)) || cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts)
    stime = tic;
    J = estimateScans(sc,Pdbnii,Vjson,b,dwi,A,opts,pth,pre,nam);
    J.method.duration = round(toc(stime));
    J.method.mix.scan = scanName(BIDSfile{i});
    for k = 1:numel(J.method.mix.further), J.method.mix.further(k).scan = scanName(BIDSfile{sc(k+1)}); end
    cat_io_json(Pjson,J);
  else
    J = cat_io_json(Pjson);
  end
  exportResults(Pjson,pth,pre,nam,maps,BIDSpathd{i},BIDSfile{i},opts);
  G = J.global;
  fprintf('%60s : mixed DTI of %d scans (%d volumes) WM FA %0.3f, MD %0.2f (10^-3 mm2/s), outlier volumes %d\n', ...
    scanName(BIDSfile{i}),numel(sc),J.method.volumes,G.FA_WM,G.MD_WM*1000,J.qc.outlier_volumes);
end
% =========================================================================
function J = estimateScans(sc,Pdbnii,Vjson,b,dwi,A,opts,pth,pre,nam)
%estimateScans. DTI of the scans sc (see estimateDTI) with the realignment
%  of the QC (realigned volumes resliced in memory, see realignedVolumes in
%  cat_io_dcm2bids_helper, otherwise the original volumes) with the maps saved as
%  <pre><map>_<nam>.nii (e.g. catDCM2BIDSdwi_FA_<scan>.nii) in pth.
  P  = cell(1,numel(sc)); Prp = P;
  for k = 1:numel(sc)
    [P{k},Prp{k}] = cat_io_dcm2bids_pp_anat('qcFiles',Pdbnii{sc(k)},opts,'dwi'); % input and realignment of the QC (of the defaced copy for level 2)
  end
  Pu = cat_io_dcm2bids_helper('prepNii',P,opts,0);
  try
    [J,Ymap] = estimateDTI(sc,Pu,Prp,Pdbnii,Vjson,b,dwi,A,opts);
  catch e
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts); rethrow(e);
  end
  cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
  for mi = 1:numel(Ymap) % (maps of the scan and wFA in MNI space, DEC with the volumes red, green, and blue)
    V = Ymap(mi).V; V.fname = fullfile(pth,[pre Ymap(mi).name '_' nam '.nii']);
    V.descrip = sprintf('DTI %s (cat_io_dcm2bids_pp_dwi)',Ymap(mi).name);
    if exist(V.fname,'file'), delete(V.fname); end
    for ci = 1:size(Ymap(mi).Y,4), V.n = [ci 1]; spm_write_vol(V,Ymap(mi).Y(:,:,:,ci)); end
    cat_io_dcm2bids_helper('prepNiigz',{V.fname},opts);
  end
end
% =========================================================================
function exportResults(Pjson,pth,pre,nam,maps,BIDSpathd,BIDSfile,opts)
%exportResults. Copy of the results Pjson and of the maps <pre><map>_<nam>
%  (for opts.preprocessing > 1, otherwise only registered in the 
%  derivatives index of the database, see registerDerivative in 
%  cat_io_dcm2bids_helper) to the BIDS derivatives BIDSpathd with the
%  BIDS name of the scan (e.g. catDCM2BIDSdwi_FA_<BIDS name>).
  bnam = spm_file(BIDSfile,'basename');
  if ~exist(BIDSpathd,'dir'), mkdir(BIDSpathd); end
  copyfile(Pjson,fullfile(BIDSpathd,[pre bnam '.json']));
  for mi = 1:numel(maps)
    Pmap = fullfile(pth,[pre maps{mi} '_' nam '.nii']);
    if ~exist(Pmap,'file'), Pmap = fullfile(pth,[pre maps{mi} '_' nam '.nii.gz']); end
    if exist(Pmap,'file')
      cat_io_dcm2bids_helper('registerDerivative',Pmap,fullfile(BIDSpathd,[pre maps{mi} '_' bnam]));
      if opts.preprocessing > 1
        cat_io_dcm2bids_helper('exportNii',Pmap,fullfile(BIDSpathd,[pre maps{mi} '_' bnam]),opts);
      end
    end
  end
end
% =========================================================================
function [maps,req] = mapNames(A,opts)
%mapNames. Maps of the DTI (MD, AD, and RD for opts.dtimaps > 1) and the
%  required maps (wFA only with the segmentation A, DECt1w only with the
%  non-linear coregistration).
  maps = {'FA','DEC','wFA','DECt1w'}; if opts.dtimaps > 1, maps = [maps {'MD','AD','RD'}]; end
  req  = setdiff(maps,{'DECt1w'},'stable'); if isempty(A), req = setdiff(req,{'wFA'},'stable'); end
end
% =========================================================================
function ok = hasMap(pth,pre,m,nam)
%hasMap. True if the map <pre><m>_<nam>.nii(.gz) exists in pth.
  ok = exist(fullfile(pth,[pre m '_' nam '.nii']),'file') || exist(fullfile(pth,[pre m '_' nam '.nii.gz']),'file');
end
% =========================================================================
function [J,Ymap] = estimateDTI(sc,P,Prp,Pdbnii,Vjson,b,dwi,A,opts)
%estimateDTI. DTI of the (unzipped) images P (cell) of the DWI scans sc of
%  the session with the realignments Prp of the QC (cell of rp_*.txt, the
%  realigned volumes are resliced in memory, see realignedVolumes in 
%  cat_io_dcm2bids_helper), i.e. of one scan or of several scans (mixed DTI, see
%  processMix) with the volumes of all scans in the grid of the first scan,
%  with the result structure J (JSON, with the coregistration to the T1w)
%  and the maps Ymap (structure array with the name, image Y, and grid V:
%  FA, MD, AD, and RD for opts.dtimaps > 1, DEC with 4 dimensions, and wFA
%  in MNI space, i.e. not modulated, as the FA is no volume measure). The
%  gradient directions (frame of dcm2niix, i.e. of FSL: voxel axes with a
%  flipped first axis for images with a positive determinant) are
%  transformed into world (RAS) coordinates of the first scan by the
%  rotation of the image (and of the coregistration of further scans), and
%  the principal direction of the tensor into the head (MNI) orientation 
%  by the rotation of the registration of the T1w (see headRotation in
%  cat_io_dcm2bids_pp_anat), i.e. the colors of the DEC map are independent
%  of the slice orientation and of the head position. The DEC map is
%  weighted by a soft brain mask (eroded by one voxel and smoothed by 1 mm).
%  Further scans are sampled at the positions of the brain voxels of the
%  first scan by the rigid coregistration of their undistorted b0 images
%  (spm_coreg) and their own distortion correction (trilinear, see
%  sampleMatrix), and the tensor fit has an S0 of each scan (e.g. for 
%  different echo times).
  i  = sc(1); Vi = cat_io_dcm2bids_helper('realignedVolumes',P{1},Prp{1},opts); dim = Vi(1).dim;
  [bi,g,sel,Yb0] = scanVolumes(i,Vi,Pdbnii,b);

  % (2) distortion correction (sampling at the undistorted positions)
  [vdm,dc] = distortion(i,Pdbnii,Vjson,b,dwi,Vi(1),Yb0,opts);
  Wall = cat_io_dcm2bids_pp_anat('unwarpMatrix',dim,vdm,(1:prod(dim))');
  Yb0u = reshape(Wall * Yb0(:),dim);

  % (3) brain, white matter, and atlas (T1w segmentation of the session)
  M = eye(4); Qc = struct(); Rh = eye(3); Dn = []; 
  if ~isempty(A) % (rigid and non-linear coregistration to the T1w)
    [M,Qc,Dn] = cat_io_dcm2bids_pp_anat('coregAnat',Yb0u,Vi(1),A,opts);
    L     = cat_io_dcm2bids_pp_anat('sampleAnat',Vi(1),M,A,{atlasFile('mori')},opts,Dn);
    brain = cat_vol_morph(L.c0 > 0.5,'e',1); wmT1 = L.c0 > 2.5; atlas = L.atlas{1}; mask = 'T1w segmentation (eroded by one voxel, with signal in all volumes)';
    Rh    = cat_io_dcm2bids_pp_anat('headRotation',A,M);
  else
    brain = Yb0u > 0.2 * cat_stat_prctile(Yb0u(:),99); wmT1 = true(dim); atlas = zeros(dim); mask = 'b0 threshold (with signal in all volumes)';
  end
  soft = zeros(dim); spm_smooth(double(brain),soft,1 ./ sqrt(sum(Vi(1).mat(1:3,1:3).^2))); % soft mask (1 mm)
  W = Wall(brain(:),:); clear Wall
  [x1,x2,x3] = ind2sub(dim,find(brain)); XYZ = Vi(1).mat * [x1'; x2'; x3'; ones(1,numel(x1))]; clear x1 x2 x3

  % volumes of the DTI of each scan (denoised for opts.denoise) with the
  % gradient directions in world coordinates of the first scan, and the 
  % tensor fit of the voxels with signal in all volumes (e.g. not the edge
  % slices of the field of view that are missing in some volumes after the
  % realignment, with high FA in the direction of the slices, or outside 
  % of the field of view of a further scan)
  Y = cell(1,numel(sc)); bs = []; gs = zeros(3,0); is = []; bad = zeros(nnz(brain),1);
  mix = struct('scan',{},'volumes',{},'distortion_correction',{},'coregistration',{},'translation',{},'rotation',{});
  for k = 1:numel(sc)
    if k == 1
      Vk = Vi; bk = bi; gk = g; selk = sel; Wk = W; Rk = eye(3);
    else
      Vk = cat_io_dcm2bids_helper('realignedVolumes',P{k},Prp{k},opts); [bk,gk,selk,Yb0k] = scanVolumes(sc(k),Vk,Pdbnii,b); dk = Vk(1).dim;
      [vdmk,dck] = distortion(sc(k),Pdbnii,Vjson,b,dwi,Vk(1),Yb0k,opts);
      Yb0k = reshape(cat_io_dcm2bids_pp_anat('unwarpMatrix',dk,vdmk,(1:prod(dk))') * Yb0k(:),dk);
      % world of scan k to world of the first scan: coregistration of the FA
      % maps of the DTI of the scans (position of the diffusion-weighted 
      % volumes, e.g. 0.5 mm more accurate than the b0 images in the DZPG 
      % pilots), otherwise of the undistorted b0 images
      [Ff,Vff] = faMap(Pdbnii{i},opts); [Fk,Vfk] = faMap(Pdbnii{sc(k)},opts);
      if ~isempty(Ff) && ~isempty(Fk)
        Mk = cat_io_dcm2bids_pp_anat('coregImages',Fk,Vfk,Ff,Vff); cm = 'FA maps (spm_coreg, nmi)';
      else
        Mk = cat_io_dcm2bids_pp_anat('coregImages',Yb0k,Vk(1),Yb0u,Vi(1)); cm = 'undistorted b0 images (spm_coreg, nmi)';
      end
      U    = Vk(1).mat \ (Mk \ XYZ);                                        % undistorted voxel positions in scan k
      [Wk,inside] = sampleMatrix(dk,vdmk,U(1:3,:)); bad = bad + ~inside(:);
      Rk   = Mk(1:3,1:3); Rk = Rk ./ sqrt(sum(Rk.^2,1));
      q    = spm_imatrix(Mk);
      mix(end+1) = struct('scan',sprintf('%d',k),'volumes',numel(selk),'distortion_correction',dck, ... % (BIDS name, see processMix)
        'coregistration',cm,'translation',round(q(1:3),1),'rotation',round(q(4:6) * 180/pi,1)); %#ok<AGROW>
    end
    Rv = Vk(1).mat(1:3,1:3); Rv = Rv ./ sqrt(sum(Rv.^2,1));                 % voxel axes in world coordinates
    Fx = diag([1 - 2*(det(Rv) > 0) 1 1]);                                   % FSL convention of the bvec of dcm2niix
    Y{k} = zeros(nnz(brain),numel(selk),'single');
    for si = 1:numel(selk)
      Yv  = single(cat_io_dcm2bids_helper('realignedVolume',Vk,selk(si)));
      bad = bad + Wk * double(~(Yv(:) > 0));                              % (also NaN, before the denoising)
      if opts.denoise, cat_sanlm(Yv,3,1,1); end
      Y{k}(:,si) = Wk * double(Yv(:));
    end
    bs = [bs bk(selk)]; gs = [gs Rk * Rv * Fx * gk(:,selk)]; is = [is k * ones(1,numel(selk))]; %#ok<AGROW>
  end
  Y  = [Y{:}];
  vb = find(brain); brain(vb(bad > 0)) = false; Y = Y(bad == 0,:);
  [x,res]          = tensorFit(Y,bs,gs,is);
  [FA,MD,AD,RD,V1] = tensorMaps(x);
  V1  = Rh * V1; V1 = V1 ./ max(eps,sqrt(sum(V1.^2,1)));   % head (MNI) orientation
  DEC = abs(V1) .* FA .* soft(brain)';                     % red, green, blue (soft brain mask)

  % (4) global values in the white matter and (5) Mori WM regions
  wm = wmT1(brain)' & FA > 0.2;
  G  = struct('FA_WM',median(FA(wm)),'MD_WM',median(MD(wm)),'AD_WM',median(AD(wm)),'RD_WM',median(RD(wm)), ...
    'MD_WM_peakwidth',diff(cat_stat_prctile(MD(wm),[5 95])),'MD_brain',median(MD),'WM_voxels',nnz(wm));
  R  = struct();
  [id,abbr] = atlasRegions('mori','Vwm');
  lab = atlas(brain)';
  for ri = 1:numel(id)
    v = lab == id(ri) & wmT1(brain)';
    if nnz(v) >= 10, R.(abbr{ri}) = struct('FA',median(FA(v)),'MD',median(MD(v)),'voxels',nnz(v)); end
  end
  mr = median(res); sr = 1.4826 * median(abs(res - mr));

  p = cat_io_dcm2bids_pp_anat('levelParameters',opts.pplevel);
  J.method = struct('level',opts.pplevel,'version',p.version,'DTI', ...
    sprintf('log-linear least squares fit (b0 and b <= %g s/mm2)',round(dtiShell(b{sc(1)}))), ...
    'volumes',numel(bs),'bvalues',unique(50 * round(bs / 50)),'denoising',denoising(opts), ...
    'distortion_correction',dc,'brain',mask,'atlas','mori (CAT)','gradient_rotation','none', ...
    'DEC','absolute principal direction in the head (MNI) orientation times FA and soft brain mask (red x, green y, blue z)', ...
    'coregistration',Qc);
  if numel(sc) > 1
    J.method.DTI = 'log-linear least squares fit (b0 and the DTI shell of each scan) of the volumes of all scans with the S0 of each scan';
    J.method.gradient_rotation = 'rotation of the coregistration of further scans (not of the realignment)';
    J.method.mix = struct('scan','1','volumes',nnz(is == 1),'further',mix);
  end
  if ~isempty(A), J.method.T1w = A.bids; J.method.T1w_other_session = isfield(A,'other') && A.other; end
  J.coreg    = M; % world space of the scan to the world space of the T1w
  J.global   = G;
  J.regional = R;
  J.qc       = struct('outlier_volumes',sum(res > mr + 3*sr),'residual',reshape(res,1,[]), ...
    'voxels_without_signal',nnz(bad));                     % (brain voxels without signal in some volumes)
  Vo = struct('dim',dim,'mat',Vi(1).mat,'dt',[16 0],'pinfo',[1 0 0]');
  if opts.dtimaps > 1, nm = {'FA','MD','AD','RD'}; Ym = {FA,MD,AD,RD}; else, nm = {'FA'}; Ym = {FA}; end
  Ymap = struct('name',nm,'Y',Ym,'V',Vo);
  for mi = 1:numel(Ymap), Y = nan(dim,'single'); Y(brain) = Ymap(mi).Y; Ymap(mi).Y = Y; end
  Yd = zeros([dim 3],'single');
  for ci = 1:3, Y = zeros(dim,'single'); Y(brain) = DEC(ci,:); Yd(:,:,:,ci) = Y; end
  Ymap(end+1) = struct('name','DEC','Y',Yd,'V',Vo);
  if ~isempty(A) % FA in MNI space (deformation of the T1w), e.g. for voxel-based analyses
    Py = cat_io_dcm2bids_helper('prepNii',{A.Py},opts,0);
    try
      % FA and its mask (normalized interpolation at the edge of the brain,
      % NaN outside of the brain)
      Yf = Ymap(1).Y; [Yw,Vw] = cat_io_dcm2bids_pp_anat('warpToMNI',cat(4,Yf,single(isfinite(Yf))),Vo,M,Py{1},Dn);
      w  = Yw(:,:,:,2); Yw = Yw(:,:,:,1) ./ max(eps,w); Yw(w < 0.5) = nan;
      Ymap(end+1) = struct('name','wFA','Y',Yw,'V',Vw);
      % DEC map in the space of the T1w (non-linear coregistration, e.g. as
      % overlay of the T1w in the rendering)
      if ~isempty(Dn)
        Ymap(end+1) = struct('name','DECt1w','Y',cat_io_dcm2bids_pp_anat('toT1w',Yd,Dn), ...
          'V',struct('dim',Dn.Vf.dim,'mat',Dn.Vf.mat,'dt',[16 0],'pinfo',[1 0 0]'));
      end
    catch e
      cat_io_cprintf('err','  FA in MNI space failed: %s\n',e.message);
    end
    cat_io_dcm2bids_helper('prepNiigz',Py,opts);
  end
end
% =========================================================================
function [Y,V] = faMap(P,opts)
%faMap. FA map of the DTI of the database image P (NaN as 0, empty if not
%  available, e.g. for the coregistration of the scans of the mixed DTI).
  Y = []; V = [];
  [pth,nam] = fileparts(regexprep(P,'\.gz$',''));
  Pf = fullfile(pth,['catDCM2BIDSdwi_FA_' nam '.nii']); if ~exist(Pf,'file'), Pf = [Pf '.gz']; end
  if ~exist(Pf,'file'), return; end
  Pu = cat_io_dcm2bids_helper('prepNii',{Pf},opts,0);
  try, V = spm_vol(Pu{1}); Y = spm_read_vols(V); Y(~isfinite(Y)) = 0; catch, Y = []; V = []; end
  cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
end
% =========================================================================
function [bi,g,sel,Yb0] = scanVolumes(i,Vi,Pdbnii,b)
%scanVolumes. b-values bi, gradient directions g (3 x volumes, see the 
%  bvec of dcm2niix), the volumes sel of the DTI (b0 and the DTI shell, see
%  dtiShell), and the mean b0 image Yb0 of the DWI scan i with the volumes Vi.
  bi = b{i};
  if numel(bi) ~= numel(Vi)
    error('cat_io_dcm2bids_pp_dwi:bval','%d b-values for %d volumes.',numel(bi),numel(Vi));
  end
  g   = load([regexprep(Pdbnii{i},'\.nii(\.gz)?$','') '.bvec']);
  ib0 = find(bi < 100);
  sel = find(bi <= dtiShell(bi));                           % b0 and DTI shell
  if isempty(ib0), error('cat_io_dcm2bids_pp_dwi:b0','No b0 volume.'); end
  Yb0 = 0; for vi = ib0, Yb0 = Yb0 + cat_io_dcm2bids_helper('realignedVolume',Vi,vi); end; Yb0 = Yb0 / numel(ib0);
end
% =========================================================================
function bmax = dtiShell(b)
%dtiShell. Upper b-value of the DTI shell of the b-values b of a scan, i.e. 
%  the DTI uses the b0 volumes and the volumes with 500 < b <= bmax: 1300 
%  s/mm2 for at least 12 volumes with b <= 1300 (e.g. b1000 or b1200 of 
%  multi-shell data), otherwise the lowest shell with at least 12 volumes
%  (b-values within 10% of its lowest value, e.g. single-shell data with
%  b1500), NaN without DTI shell.
  bmax = nan; b = b(:)';
  if sum(b > 500 & b <= 1300) >= 12, bmax = 1300; return; end
  for bk = sort(b(b > 500))
    if sum(b >= bk & b <= 1.1 * bk) >= 12, bmax = 1.1 * bk; return; end
  end
end
% =========================================================================
function [W,inside] = sampleMatrix(dim,vdm,U)
%sampleMatrix. Sparse matrix W (positions x prod(dim)) of the trilinear 
%  interpolation of a distorted volume (grid dim) at the undistorted voxel
%  positions U (3 x positions, e.g. of the brain voxels of another scan),
%  i.e. at U shifted by the voxel displacement map vdm along the phase-
%  encoding axis j (interpolated, no shift for an empty vdm, see 
%  unwarpMatrix in cat_io_dcm2bids_pp_anat), with the positions inside of
%  the volume (only these with weights).
  n = size(U,2);
  if ~isempty(vdm)
    d = spm_sample_vol(double(vdm),U(1,:),U(2,:),U(3,:),1); d(~isfinite(d)) = 0;
    U(2,:) = U(2,:) + d;
  end
  f = floor(U); w = U - f;
  inside = all(f >= 1 & f <= dim(:) - 1,1);
  ii = find(inside); m = numel(ii); r = zeros(8*m,1); c = r; v = r;
  for o = 0:7
    a  = [bitget(o,1); bitget(o,2); bitget(o,3)];
    p  = f(:,ii) + a; k = o*m + (1:m);
    r(k) = ii; c(k) = sub2ind(dim,p(1,:),p(2,:),p(3,:)); v(k) = prod(a .* w(:,ii) + (1 - a) .* (1 - w(:,ii)),1);
  end
  W = sparse(r,c,v,n,prod(dim));
end
% =========================================================================
function [vdm,dc] = distortion(i,Pdbnii,Vjson,b,dwi,Vt,Yt,opts)
%distortion. Voxel displacement map of the DWI scan i by SCOPE of its b0
%  volumes and the b0 volumes of the closest DWI scan (series number) with
%  opposite phase-encoding direction along j, the same image size, and the
%  same total readout time (see scopeVDM in cat_io_dcm2bids_pp_anat).
  vdm = [];
  pe  = @(k) fieldValue(Vjson{k},'PhaseEncodingDirection','');
  trt = @(k) fieldValue(Vjson{k},'TotalReadoutTime',nan);
  snr = @(k) fieldValue(Vjson{k},'SeriesNumber',nan);
  dm  = @(k) reshape(fieldValue(Vjson{k},'Dimensions',nan(1,3)),1,[]); % (with the number of volumes)
  if ~any(strcmp(pe(i),{'j','j-'})), dc = 'none (phase-encoding axis not j)'; return; end
  p = []; di = dm(i);
  for k = dwi
    dk = dm(k);
    if k ~= i && any(strcmp(pe(k),{'j','j-'})) && ~strcmp(pe(k),pe(i)) && isequal(dk(1:3),di(1:3)) && ...
        abs(trt(k) - trt(i)) <= 0.01 * trt(i) && any(b{k} < 100)
      p(end+1) = k; %#ok<AGROW>
    end
  end
  if isempty(p), dc = 'none (no DWI scan with opposite phase-encoding direction)'; return; end
  [~,pk] = min(abs(arrayfun(snr,p) - snr(i))); p = p(pk);
  pair = [i p]; if strcmp(pe(p),'j-'), pair = [p i]; end % first image with negative polarity (cache)
  [vdm,dc] = cat_io_dcm2bids_pp_anat('scopeVDM',Pdbnii(pair),{find(b{pair(1)} < 100),find(b{pair(2)} < 100)}, ...
    Vt,Yt,pair(1) ~= i,[trt(i) trt(p)],opts);
end
% =========================================================================
function [x,res] = tensorFit(Y,b,g,s)
%tensorFit. Log-linear least squares fit of the diffusion tensor of the
%  signal Y (voxels x volumes) with the b-values b, gradient directions g
%  (3 x volumes, b0 with zero vectors), and scan s of each volume (S0 of
%  each scan, default one scan): x (6 + scans x voxels) with ln(S0) of the
%  first scan, the tensor elements Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, and ln(S0)
%  of further scans, and the RMS res of the residuals of each volume (log
%  signal).
  if ~exist('s','var'), s = ones(1,numel(b)); end
  g = g ./ max(eps,sqrt(sum(g.^2,1)));
  E = double(s(:) == (1:max(s)));                            % (volumes x scans)
  B = [E(:,1) -b(:) .* [g(1,:)'.^2 g(2,:)'.^2 g(3,:)'.^2 ...
    2*g(1,:)'.*g(2,:)' 2*g(1,:)'.*g(3,:)' 2*g(2,:)'.*g(3,:)'] E(:,2:end)];
  y   = log(max(1,double(Y)))';
  x   = B \ y;
  res = sqrt(mean((y - B*x).^2,2));
end
% =========================================================================
function [FA,MD,AD,RD,V1] = tensorMaps(x)
%tensorMaps. Fractional anisotropy (FA), mean (MD), axial (AD), and radial
%  (RD) diffusivity of the tensors x (see tensorFit) by the closed-form
%  eigenvalues of symmetric 3x3 matrices (vectorized), and the principal
%  direction V1 (3 x voxels, unit vectors in the frame of the gradient
%  directions) by the cross product of the two rows of D - l1*I with the
%  largest cross product.
  a  = x(2,:); d = x(3,:); c = x(4,:); e = x(5,:); f = x(6,:); h = x(7,:); % Dxx Dyy Dzz Dxy Dxz Dyz
  q  = (a + d + c) / 3;
  p  = sqrt(((a-q).^2 + (d-q).^2 + (c-q).^2 + 2*(e.^2 + f.^2 + h.^2)) / 6); p(p == 0) = eps;
  Ba = (a-q)./p; Bd = (d-q)./p; Bc = (c-q)./p; Be = e./p; Bf = f./p; Bh = h./p;
  r  = (Ba.*(Bd.*Bc - Bh.^2) - Be.*(Be.*Bc - Bh.*Bf) + Bf.*(Be.*Bh - Bd.*Bf)) / 2;
  ph = acos(min(1,max(-1,r))) / 3;
  l1 = q + 2*p.*cos(ph); l3 = q + 2*p.*cos(ph + 2*pi/3); l2 = 3*q - l1 - l3;
  MD = q; AD = l1; RD = (l2 + l3) / 2;
  FA = min(1,sqrt(0.5) * sqrt((l1-l2).^2 + (l2-l3).^2 + (l3-l1).^2) ./ max(eps,sqrt(l1.^2 + l2.^2 + l3.^2)));
  if nargout > 4
    n  = numel(a); r1 = [a - l1; e; f]; r2 = [e; d - l1; h]; r3 = [f; h; c - l1];
    C  = cat(3,cross(r1,r2,1),cross(r1,r3,1),cross(r2,r3,1)); % 3 x voxels x row pairs
    [m,k] = max(reshape(sum(C.^2,1),n,3),[],2);              % row pair with the largest cross product
    C  = reshape(permute(C,[2 3 1]),n*3,3);                   % (voxel, row pair) x component
    V1 = (C(sub2ind([n 3],(1:n)',k),:) ./ max(eps,sqrt(m)))';
  end
end
% =========================================================================
function P = atlasFile(name)
%atlasFile. Atlas of CAT in MNI space (templates_MNI152NLin2009cAsym).
  P = fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym',[name '.nii']);
end
% =========================================================================
function [id,abbr] = atlasRegions(name,col)
%atlasRegions. Regions (ROIid and ROIabbr) of a CAT atlas with the value 1
%  in the column col of its csv file (e.g. Vwm for white matter regions).
  C = regexp(fileread(spm_file(atlasFile(name),'ext','.csv')),'\r?\n','split');
  C = C(~cellfun('isempty',strtrim(C)));
  H = strsplit(C{1},';'); D = cellfun(@(l) strsplit(l,';'),C(2:end),'UniformOutput',false);
  D = vertcat(D{:});
  ok = strcmp(strtrim(D(:,strcmp(H,col))),'1');
  id = str2double(D(ok,strcmp(H,'ROIid'))); abbr = matlab.lang.makeValidName(strtrim(D(ok,strcmp(H,'ROIabbr'))));
end
% =========================================================================
function s = denoising(opts)
%denoising. Description of the denoising of the volumes.
  if opts.denoise, s = 'sanlm (each volume)'; else, s = 'none'; end
end
% =========================================================================
function v = fieldValue(S,fn,default)
%fieldValue. Field fn of the structure S (default if not available).
  if isstruct(S) && isfield(S,fn) && ~isempty(S.(fn)), v = S.(fn); else, v = default; end
end
% =========================================================================
function s = scanName(BIDSfile)
%scanName. BIDS name of a scan without subject and session.
  s = regexprep(spm_file(BIDSfile,'basename'),'^sub-[^_]+_(ses-[^_]+_)?','');
end
