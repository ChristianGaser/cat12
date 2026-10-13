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
%        the brain (T1w segmentation of the session by the rigid
%        coregistration of the b0 image, otherwise b0 threshold) with the
%        fractional anisotropy (FA), mean (MD), axial (AD), and radial (RD)
%        diffusivity (in mm2/s), and the fit residual of each volume
%        (outlier volumes with more than 3 SD of the residuals)
%    (4) global values: medians of FA, MD, AD, and RD in the white matter
%        (WM of the T1w segmentation with FA > 0.2), the peak width of the
%        MD in the white matter (95th - 5th percentile, similar to the PSMD
%        marker of small vessel disease but without skeleton), and the
%        median MD of the brain
%    (5) regional values: median FA and MD of the white matter regions of
%        the Mori atlas of CAT (by the deformation of the T1w segmentation)
%  The results are saved in the database (catDCM2BIDSdwi_<scan>.json and
%  the maps catDCM2BIDSdwi_<FA|MD|AD|RD>_<scan>.nii in the undistorted
%  space of the scan) and copied to the BIDS derivatives (the maps only for
%  opts.preprocessing > 1). The gradient directions are not rotated by the
%  realignment (small effect for typical motion of a few degree).
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
%  exported.
  [pth,nam] = fileparts(regexprep(Pdbnii{i},'\.gz$',''));
  Pjson = fullfile(pth,['catDCM2BIDSdwi_' nam '.json']);
  maps  = {'FA','MD','AD','RD'};
  if ~compute && ~exist(Pjson,'file'), return; end
  if compute && cat_io_dcm2bids_pp_anat('levelUpdate',Pjson,opts)
    stime = tic;
    P  = spm_file(Pdbnii{i},'prefix','r'); if ~exist(P,'file'), P = Pdbnii{i}; end % realigned by the QC
    Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
    try
      [J,Ymap,Vo] = estimateDTI(i,Pu{1},Pdbnii,Vjson,b,dwi,A,opts);
    catch e
      cat_io_dcm2bids_helper('prepNiigz',Pu,opts); rethrow(e);
    end
    cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
    for mi = 1:numel(maps)
      Vo.fname   = fullfile(pth,['catDCM2BIDSdwi_' maps{mi} '_' nam '.nii']);
      Vo.descrip = sprintf('DTI %s (cat_io_dcm2bids_pp_dwi)',maps{mi});
      spm_write_vol(Vo,Ymap{mi});
      cat_io_dcm2bids_helper('prepNiigz',{Vo.fname},opts);
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
function [J,Ymap,Vo] = estimateDTI(i,P,Pdbnii,Vjson,b,dwi,A,opts)
%estimateDTI. DTI of the (unzipped) image P of the DWI scan i with the
%  result structure J (JSON), the maps Ymap (FA, MD, AD, RD), and the
%  header of the maps Vo.
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
  if ~isempty(A)
    [Ym,Vm] = readImage(A.Pm,opts);
    M     = cat_io_dcm2bids_pp_anat('coregImages',Yb0u,Vi(1),Ym,Vm);
    L     = cat_io_dcm2bids_pp_anat('sampleAnat',Vi(1),M,A,{atlasFile('mori')},opts);
    brain = L.c0 > 0.5; wmT1 = L.c0 > 2.5; atlas = L.atlas{1}; mask = 'T1w segmentation';
  else
    brain = Yb0u > 0.2 * prctile(Yb0u(:),99); wmT1 = true(dim); atlas = zeros(dim); mask = 'b0 threshold';
  end
  W = Wall(brain(:),:); clear Wall

  % volumes of the DTI (denoised for opts.denoise) and tensor fit
  Y = zeros(nnz(brain),numel(sel),'single');
  for si = 1:numel(sel)
    Yv = single(spm_read_vols(Vi(sel(si))));
    if opts.denoise, cat_sanlm(Yv,3,1,1); end
    Y(:,si) = W * double(Yv(:));
  end
  [x,res]       = tensorFit(Y,bi(sel),g(:,sel));
  [FA,MD,AD,RD] = tensorMaps(x);

  % (4) global values in the white matter and (5) Mori WM regions
  wm = wmT1(brain)' & FA > 0.2;
  G  = struct('FA_WM',median(FA(wm)),'MD_WM',median(MD(wm)),'AD_WM',median(AD(wm)),'RD_WM',median(RD(wm)), ...
    'MD_WM_peakwidth',diff(prctile(MD(wm),[5 95])),'MD_brain',median(MD),'WM_voxels',nnz(wm));
  R  = struct();
  [id,abbr] = atlasRegions('mori','Vwm');
  lab = atlas(brain)';
  for ri = 1:numel(id)
    v = lab == id(ri) & wmT1(brain)';
    if nnz(v) >= 10, R.(abbr{ri}) = struct('FA',median(FA(v)),'MD',median(MD(v)),'voxels',nnz(v)); end
  end
  mr = median(res); sr = 1.4826 * median(abs(res - mr));

  J.method = struct('level',opts.pplevel,'DTI','log-linear least squares fit (b0 and b <= 1300 s/mm2)', ...
    'volumes',numel(sel),'bvalues',unique(50 * round(bi(sel) / 50)),'denoising',denoising(opts), ...
    'distortion_correction',dc,'brain',mask,'atlas','mori (CAT)','gradient_rotation','none');
  J.global   = G;
  J.regional = R;
  J.qc       = struct('outlier_volumes',sum(res > mr + 3*sr),'residual',reshape(res,1,[]));
  Vo = struct('dim',dim,'mat',Vi(1).mat,'dt',[16 0],'pinfo',[1 0 0]');
  Ymap = {FA,MD,AD,RD};
  for mi = 1:numel(Ymap), Y = nan(dim,'single'); Y(brain) = Ymap{mi}; Ymap{mi} = Y; end
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
function [FA,MD,AD,RD] = tensorMaps(x)
%tensorMaps. Fractional anisotropy (FA), mean (MD), axial (AD), and radial
%  (RD) diffusivity of the tensors x (see tensorFit) by the closed-form
%  eigenvalues of symmetric 3x3 matrices (vectorized).
  a  = x(2,:); d = x(3,:); c = x(4,:); e = x(5,:); f = x(6,:); h = x(7,:); % Dxx Dyy Dzz Dxy Dxz Dyz
  q  = (a + d + c) / 3;
  p  = sqrt(((a-q).^2 + (d-q).^2 + (c-q).^2 + 2*(e.^2 + f.^2 + h.^2)) / 6); p(p == 0) = eps;
  Ba = (a-q)./p; Bd = (d-q)./p; Bc = (c-q)./p; Be = e./p; Bf = f./p; Bh = h./p;
  r  = (Ba.*(Bd.*Bc - Bh.^2) - Be.*(Be.*Bc - Bh.*Bf) + Bf.*(Be.*Bh - Bd.*Bf)) / 2;
  ph = acos(min(1,max(-1,r))) / 3;
  l1 = q + 2*p.*cos(ph); l3 = q + 2*p.*cos(ph + 2*pi/3); l2 = 3*q - l1 - l3;
  MD = q; AD = l1; RD = (l2 + l3) / 2;
  FA = min(1,sqrt(0.5) * sqrt((l1-l2).^2 + (l2-l3).^2 + (l3-l1).^2) ./ max(eps,sqrt(l1.^2 + l2.^2 + l3.^2)));
end
% =========================================================================
function [Y,V] = readImage(P,opts)
%readImage. Image P (.nii or .nii.gz) with its header V.
  Pu = cat_io_dcm2bids_helper('prepNii',{P},opts,0);
  V  = spm_vol(Pu{1}); V = V(1); Y = spm_read_vols(V);
  cat_io_dcm2bids_helper('prepNiigz',Pu,opts);
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
