function varargout = cat_io_dcm2bids_pp_dwi(action,varargin)
%cat_io_dcm2bids_pp_dwi. Basic preprocessing of diffusion-weighted images.
%  Diffusion tensor imaging (DTI) of each DWI scan of a session with at
%  least 12 diffusion-weighted volumes with b <= 1300 s/mm2 (DTI shell),
%  with the parameters of the preprocessing level opts.pplevel (see 
%  levelParameters in cat_io_dcm2bids_pp_anat, i.e. the resolution and settings
%  of SCOPE, and only the first scan for higher levels, whereas subsets of
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
%  derivatives (the maps only for opts.preprocessing > 1, rendered with 
%  their scan besides wFA and DECt1w, see cat_io_dcm2bids_render). The
%  segmentation is the one of the session (T1w, otherwise UNIT1, T2w, or
%  FLAIR, see anatScan in cat_io_dcm2bids_pp_anat). The gradient
%  directions are not rotated by the realignment (small effect for typical
%  motion of a few degree).
%
%  varargout = cat_io_dcm2bids_pp_dwi(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_pp_dwi('processSession',Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
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
  dti = dwi(cellfun(@(x) sum(x > 500 & x <= 1300) >= 12, b(dwi)));
  if isempty(dti), return; end
  A = cat_io_dcm2bids_pp_anat('sessionAnat',Pdbnii,datatype,suffix);
  if ~isempty(A) % BIDS file of the segmented scan (e.g. the T1w as underlay of the rendering)
    A.bids = ''; if ~isempty(BIDSfile{A.i}), A.bids = spm_file(BIDSfile{A.i},'basename'); end
  end
  p = cat_io_dcm2bids_pp_anat('levelParameters',opts.pplevel);
  for i = dti
    try
      % only the first scan for higher levels (others only if processed before)
      processScan(i,Pdbnii,Vjson,b,dwi,A,BIDSpathd{i},BIDSfile{i},opts,~p.onescan || i == dti(1));
    catch e
      cat_io_cprintf('err','%60s : DTI failed: %s (%s, line %d)\n',scanName(BIDSfile{i}), ...
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
  Pjson = fullfile(pth,['catDCM2BIDSdwi_' nam '.json']);
  maps  = {'FA','DEC','wFA','DECt1w'}; if opts.dtimaps > 1, maps = [maps {'MD','AD','RD'}]; end
  % results of older versions or without the maps of the options (e.g. the
  % wFA, or MD, AD, and RD for opts.dtimaps = 2) are processed again (also
  % without compute), where wFA needs the deformation of the segmentation
  % (and DECt1w also the non-linear coregistration)
  req     = setdiff(maps,{'DECt1w'},'stable'); if isempty(A), req = setdiff(req,{'wFA'},'stable'); end
  isnii   = @(m) exist(fullfile(pth,['catDCM2BIDSdwi_' m '_' nam '.nii']),'file') || ...
    exist(fullfile(pth,['catDCM2BIDSdwi_' m '_' nam '.nii.gz']),'file');
  missing = exist(Pjson,'file') && ~all(cellfun(isnii,req));
  compute = compute || (exist(Pjson,'file') && cat_io_dcm2bids_pp_anat('levelOutdated',Pjson)) || missing;
  if ~compute && ~exist(Pjson,'file'), return; end
  if compute && (missing || cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts))
    stime = tic;
    P  = spm_file(Pdbnii{i},'prefix','r'); if ~exist(P,'file'), P = Pdbnii{i}; end % realigned by the QC
    Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
    try
      [J,Ymap] = estimateDTI(i,Pu{1},Pdbnii,Vjson,b,dwi,A,opts);
    catch e
      cat_io_dcm2bids_helper('prepNiigz',Pu,opts); rethrow(e);
    end
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
    for mi = 1:numel(Ymap) % (maps of the scan and wFA in MNI space, DEC with the volumes red, green, and blue)
      V = Ymap(mi).V; V.fname = fullfile(pth,['catDCM2BIDSdwi_' Ymap(mi).name '_' nam '.nii']);
      V.descrip = sprintf('DTI %s (cat_io_dcm2bids_pp_dwi)',Ymap(mi).name);
      if exist(V.fname,'file'), delete(V.fname); end
      for ci = 1:size(Ymap(mi).Y,4), V.n = [ci 1]; spm_write_vol(V,Ymap(mi).Y(:,:,:,ci)); end
      cat_io_dcm2bids_helper('prepNiigz',{V.fname},opts);
    end
    J.method.duration = round(toc(stime));
    cat_io_json(Pjson,J);
  else
    J = cat_io_json(Pjson);
  end

  % copy to the BIDS derivatives (the maps for opts.preprocessing > 1)
  bnam = spm_file(BIDSfile,'basename');
  if ~exist(BIDSpathd,'dir'), mkdir(BIDSpathd); end
  copyfile(Pjson,fullfile(BIDSpathd,['catDCM2BIDSdwi_' bnam '.json']));
  if opts.preprocessing > 1
    for mi = 1:numel(maps)
      Pmap = fullfile(pth,['catDCM2BIDSdwi_' maps{mi} '_' nam '.nii']);
      if ~exist(Pmap,'file'), Pmap = fullfile(pth,['catDCM2BIDSdwi_' maps{mi} '_' nam '.nii.gz']); end
      if exist(Pmap,'file')
        cat_io_dcm2bids_helper('exportNii',Pmap,fullfile(BIDSpathd,['catDCM2BIDSdwi_' maps{mi} '_' bnam]),opts);
      end
    end
  end
  G = J.global;
  fprintf('%60s : DTI WM FA %0.3f, MD %0.2f, AD %0.2f, RD %0.2f (10^-3 mm2/s), outlier volumes %d, %s\n', ...
    scanName(BIDSfile),G.FA_WM,G.MD_WM*1000,G.AD_WM*1000,G.RD_WM*1000,J.qc.outlier_volumes, ...
    J.method.distortion_correction);
end
% =========================================================================
function [J,Ymap] = estimateDTI(i,P,Pdbnii,Vjson,b,dwi,A,opts)
%estimateDTI. DTI of the (unzipped) image P of the DWI scan i with the
%  result structure J (JSON, with the coregistration to the T1w) and the
%  maps Ymap (structure array with the name, image Y, and grid V: FA, MD,
%  AD, and RD for opts.dtimaps > 1, DEC with 4 dimensions, and wFA in MNI
%  space, i.e. not modulated, as the FA is no volume measure). The
%  principal direction of the tensor (in the frame of the gradient 
%  directions of dcm2niix, i.e. of FSL: voxel axes with a flipped first 
%  axis for images with a positive determinant) is transformed into world
%  (RAS) coordinates by the rotation of the image and into the head (MNI)
%  orientation by the rotation of the registration of the T1w (see 
%  headRotation in cat_io_dcm2bids_pp_anat), i.e. the colors of the DEC map
%  are independent of the slice orientation and of the head position. The
%  DEC map is weighted by a soft brain mask (eroded by one voxel and 
%  smoothed by 1 mm).
  Vi  = spm_vol(P); dim = Vi(1).dim; bi = b{i};
  if numel(bi) ~= numel(Vi)
    error('cat_io_dcm2bids_pp_dwi:bval','%d b-values for %d volumes.',numel(bi),numel(Vi));
  end
  g   = load([regexprep(Pdbnii{i},'\.nii(\.gz)?$','') '.bvec']);
  ib0 = find(bi < 100);
  sel = find(bi <= 1300);                                   % b0 and DTI shell
  if isempty(ib0), error('cat_io_dcm2bids_pp_dwi:b0','No b0 volume.'); end
  Yb0 = 0; for vi = ib0, Yb0 = Yb0 + spm_read_vols(Vi(vi)); end; Yb0 = Yb0 / numel(ib0);

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

  % volumes of the DTI (denoised for opts.denoise) and tensor fit of the
  % voxels with signal in all volumes (e.g. not the edge slices of the field
  % of view that are missing in some volumes after the realignment, with
  % high FA in the direction of the slices)
  Y = zeros(nnz(brain),numel(sel),'single'); bad = zeros(nnz(brain),1);
  for si = 1:numel(sel)
    Yv  = single(spm_read_vols(Vi(sel(si))));
    bad = bad + W * double(~(Yv(:) > 0));                  % (also NaN, before the denoising)
    if opts.denoise, cat_sanlm(Yv,3,1,1); end
    Y(:,si) = W * double(Yv(:));
  end
  vb = find(brain); brain(vb(bad > 0)) = false; Y = Y(bad == 0,:);
  [x,res]          = tensorFit(Y,bi(sel),g(:,sel));
  [FA,MD,AD,RD,V1] = tensorMaps(x);
  Rv = Vi(1).mat(1:3,1:3); Rv = Rv ./ sqrt(sum(Rv.^2,1)); % voxel axes in world coordinates
  if det(Rv) > 0, V1(1,:) = -V1(1,:); end                  % FSL convention of the bvec of dcm2niix
  V1  = Rh * Rv * V1; V1 = V1 ./ max(eps,sqrt(sum(V1.^2,1))); % head (MNI) orientation
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
  J.method = struct('level',opts.pplevel,'version',p.version,'DTI','log-linear least squares fit (b0 and b <= 1300 s/mm2)', ...
    'volumes',numel(sel),'bvalues',unique(50 * round(bi(sel) / 50)),'denoising',denoising(opts), ...
    'distortion_correction',dc,'brain',mask,'atlas','mori (CAT)','gradient_rotation','none', ...
    'DEC','absolute principal direction in the head (MNI) orientation times FA and soft brain mask (red x, green y, blue z)', ...
    'coregistration',Qc);
  if ~isempty(A), J.method.T1w = A.bids; end
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
function [x,res] = tensorFit(Y,b,g)
%tensorFit. Log-linear least squares fit of the diffusion tensor of the
%  signal Y (voxels x volumes) with the b-values b and gradient directions
%  g (3 x volumes, b0 with zero vectors): x (7 x voxels) with ln(S0) and
%  the tensor elements Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, and the RMS res of the
%  residuals of each volume (log signal).
  g = g ./ max(eps,sqrt(sum(g.^2,1)));
  B = [ones(numel(b),1) -b(:) .* [g(1,:)'.^2 g(2,:)'.^2 g(3,:)'.^2 ...
    2*g(1,:)'.*g(2,:)' 2*g(1,:)'.*g(3,:)' 2*g(2,:)'.*g(3,:)']];
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
