function varargout = cat_io_dcm2bids_pp_surf(action,varargin)
%cat_io_dcm2bids_pp_surf. Quick surface preprocessing of cat_io_dcm2bids.
%  Fast surface reconstruction of the segmented anatomical scan of each
%  session (SPM label map c0, bias-corrected image m, and deformation y_
%  of the basic preprocessing, see segmentanat and sessionAnat in
%  cat_io_dcm2bids_pp_anat) for regional surface values that are
%  independent of the CAT version and for the review of the data
%  (opts.surf.run, not by default). It is a simplified version of the CAT
%  surface reconstruction (cat_surf_createCS4) without topology correction,
%  spherical mapping, and surface registration, as the regions are given
%  by the volume atlas (lotter2025 in the space of the scan):
%    (1) label map of each hemisphere with the deep structures filled as WM
%        (basal ganglia, thalamus, ventricles, and corpus callosum, but not
%        the hippocampus and amygdala) and without cerebellum and brainstem
%        (see atlasMaps in cat_io_dcm2bids_pp_anat), interpolated to the
%        resolution opts.surf.res (e.g. 0.8 mm, 0.5 mm as in CAT)
%    (2) thickness and position map by the projection-based thickness
%        (cat_vol_pbtsimpleCS4)
%    (3) central surface by the marching cubes of the position map without
%        topology correction (CAT_VolMarchingCubes -fast), reduced to
%        opts.surf.reduction times fewer faces than in CAT (2 * area in mm2
%        faces of CAT, see cat_surf_createCS4, CAT_SurfReduce)
%    (4) thickness and T1w intensity (scaled to CSF 1/3, GM 2/3, and WM 1
%        as in CAT, by the medians of the tissues of the label map) at the
%        vertices with the folding correction of CAT (Demirci et al. 2025,
%        CAT_SurfCorrectThicknessFolding) and the thickness offset of the
%        PBT resolution (see thicknessOffset), the area of the vertices, and
%        the GM region of the vertices (without the vertices of the medial
%        wall, i.e. of the filled structures or without thickness)
%    (5) absolute mean curvature (CAT_DumpCurv, in degrees, averaged over 
%        3 mm) and sulcal depth (distance to the convex hull of the 
%        hemisphere in mm, see hullDepth) of the vertices
%    (6) regional values of the cortical regions (Desikan-Killiany regions
%        of lotter2025): mean thickness, intensity, curvature, and depth,
%        area (sum), and number of vertices, and the global values of each
%        hemisphere (with the number of topology defects of its surface,
%        see topologyDefects)
%    (7) surfaces of both hemispheres reduced to opts.surf.reduce vertices
%        (e.g. 32k, CAT_SurfReduce) with the values of the nearest vertex 
%        (columns of the cdata: thickness, intensity, curvature, depth, and
%        ROIid of lotter2025, 0 for the medial wall), e.g. for the rendering
%  The results are saved in the database (catDCM2BIDSsurf_<scan>.json and 
%  catDCM2BIDSsurf_<lh|rh>_<scan>.gii) and copied to the BIDS derivatives
%  (the surfaces for opts.preprocessing > 1, otherwise only registered in
%  the derivatives index of the database). Results of another version of
%  this function (method.version) are processed again.
%
%  varargout = cat_io_dcm2bids_pp_surf(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_pp_surf('processSession',Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts)
%    J = cat_io_dcm2bids_pp_surf('processScan',A,opts,Pjson[,Pwork])
%    d = cat_io_dcm2bids_pp_surf('hullDepth',V)
%
%  See also cat_io_dcm2bids_pp_anat, cat_surf_createCS4, cat_vol_pbtsimpleCS4.

  switch action
    case {'processSession','processScan','hullDepth'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_pp_surf:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function processSession(Pdbnii,Vjson,datatype,suffix,BIDSpathd,BIDSfile,opts) %#ok<INUSL>
%processSession. Surface values of the segmented anatomical scan of a
%  session (see sessionAnat in cat_io_dcm2bids_pp_anat) with the database
%  images Pdbnii, sidecars Vjson, BIDS datatype and suffix, derivatives
%  directory BIDSpathd, and BIDS file of all scans of the session.
  A = cat_io_dcm2bids_pp_anat('sessionAnat',Pdbnii,datatype,suffix);
  if isempty(A), return; end
  [pth,nam] = fileparts(regexprep(Pdbnii{A.i},'\.gz$',''));
  Pjson = fullfile(pth,['catDCM2BIDSsurf_' regexprep(nam,'^anon_','') '.json']);
  run = ~exist(Pjson,'file') || opts.rerun;
  if ~run, J = cat_io_json(Pjson); run = ~isfield(J,'method') || ~isfield(J.method,'version') || J.method.version ~= surfVersion; end
  if run
    try
      J = processScan(A,opts,Pjson); cat_io_json(Pjson,J);
    catch e
      cat_io_cprintf('err','%60s : surface failed: %s (%s, line %d)\n',spm_file(BIDSfile{A.i},'basename'), ...
        e.message,e.stack(1).name,e.stack(1).line);
      return
    end
  end
  if ~isempty(BIDSfile{A.i}) && ~isempty(BIDSpathd{A.i})
    if ~exist(BIDSpathd{A.i},'dir'), mkdir(BIDSpathd{A.i}); end
    bnam = spm_file(BIDSfile{A.i},'basename');
    copyfile(Pjson,fullfile(BIDSpathd{A.i},['catDCM2BIDSsurf_' bnam '.json']));
    for side = {'lh','rh'}
      Pgii = surfFile(Pjson,side{1}); if ~exist(Pgii,'file'), continue; end
      Pbgii = fullfile(BIDSpathd{A.i},['catDCM2BIDSsurf_' side{1} '_' bnam '.gii']);
      cat_io_dcm2bids_helper('registerDerivative',Pgii,regexprep(Pbgii,'\.gii$',''));
      if opts.preprocessing > 1, copyfile(Pgii,Pbgii); end
    end
  end
  G = J.global;
  fprintf('%60s : surface thickness %0.2f/%0.2f mm, intensity %0.3f/%0.3f, area %0.0f/%0.0f cm2 (left/right)\n', ...
    spm_file(BIDSfile{A.i},'basename'),G.lh.thickness,G.rh.thickness,G.lh.intensity,G.rh.intensity,G.lh.area/100,G.rh.area/100);
end
% =========================================================================
function P = surfFile(Pjson,side)
%surfFile. Reduced surface of the hemisphere side (lh or rh) of the result
%  Pjson (catDCM2BIDSsurf_<side>_<scan>.gii, see processScan).
  P = strrep(regexprep(Pjson,'\.json$','.gii'),[filesep 'catDCM2BIDSsurf_'],[filesep 'catDCM2BIDSsurf_' side '_']);
end
% =========================================================================
function v = surfVersion
%surfVersion. Version of the surface values (results of other versions are
%  processed again).
  v = 3;
end
% =========================================================================
function o = thicknessOffset(res)
%thicknessOffset. Offset in mm of the projection-based thickness of the PBT
%  resolution res (linear between 0.5, 0.8, and 1.0 mm, limited to this
%  range) as the median difference of the nominal thickness of thickness
%  phantoms (mri_simulate, 1.5-3.0 mm, voxel size 0.6-1.0 mm, with the SPM 
%  segmentation of the basic preprocessing) and the measured thickness, 
%  i.e. the PBT is thinner at lower resolution (0.5 mm about as CAT). 
%  The offset of a scan also depends on its segmentation (SD about 0.17 mm
%  between the phantoms, 0.05 mm with the CAT segmentation), which is not
%  corrected.
  o = interp1([0.5 0.8 1.0],[-0.04 0.12 0.17],min(1,max(0.5,res)));
end
% =========================================================================
function J = processScan(A,opts,Pjson,Pwork)
%processScan. Surface values of the segmentation A (see sessionAnat in
%  cat_io_dcm2bids_pp_anat) as JSON structure (see the help of this file)
%  with the reduced surfaces next to the result Pjson (see surfFile), and
%  the temporary files in the directory Pwork (kept if given, e.g. for 
%  tests, otherwise a temporary directory that is removed).
  stime = tic; dt = zeros(1,4); t = tic; % (durations of the preparation, PBT, surface, and values)
  keep  = exist('Pwork','var') && ~isempty(Pwork);
  if ~keep, Pwork = tempname; end
  if ~exist(Pwork,'dir'), mkdir(Pwork); end
  P   = cat_io_dcm2bids_helper('prepNii',{A.Pc0; A.Pm; A.Py},opts,0);
  V   = spm_vol(P{1}); Yp0 = single(spm_read_vols(V)); Ym = single(spm_read_vols(spm_vol(P{2})));
  M   = cat_io_dcm2bids_pp_anat('atlasMaps',Yp0,V,P{3});
  cat_io_dcm2bids_helper('prepNiigz',P,opts);
  vx  = sqrt(sum(V.mat(1:3,1:3).^2));

  % intensity scaling (CSF 1/3, GM 2/3, WM 1) and remaining noise of the
  % WM (as in cat_surf_createCS4)
  Ymn = scaleIntensity(Ym,Yp0); clear Ym
  Ywm = cat_vol_resize(Ymn .* (Yp0 > 2.75),'reduceV',vx,2,16,'meanm');
  Ysd = cat_vol_localstat(Ywm,Ywm > 0,2,4); wmnoise = min(1/6,cat_stat_nanmean(Ysd(Ysd > 0))); clear Ywm Ysd

  % filling of the deep structures as WM (not the hippocampus and amygdala
  % and their surrounding, as in cat_surf_createCS4)
  Yfill = cat_vol_morph(M.fill,'dc',2,vx) & ~cat_vol_morph(M.hc,'dd',2,vx) & M.side > 0;
  Yp0f  = max(Yp0,3 * single(cat_vol_smooth3X(single(Yfill),1) > 0.5));
  Ywall = single(cat_vol_morph(Yfill,'dd',1,vx)); % (medial wall, not cortex)
  res   = opts.surf.res; dt(1) = toc(t);

  side = {'lh','rh'}; T = cat_io_dcm2bids_pp_anat('atlasTable','lotter2025_aparc'); G = struct(); R = struct(); nf = [0 0];
  cortex = T.id(ismember(T.lobe,{'frontal','temporal','parietal','occipital','cingulate','insula'}));
  for si = 1:2
    % (1) label map of the hemisphere at the resolution res
    t  = tic; Ys = Yp0f .* (M.side == si);
    bb = boundingBox(Ys > 0.5,round(5 ./ vx));
    [Yi,Mi] = interpCrop(max(1,Ys),V.mat,bb,res); clear Ys

    % (2) thickness and position map
    [Yth,Ypp] = cat_vol_pbtsimpleCS4(Yi,res,struct('verb',0,'gyrusrecon',1,'eidist',1,'NVBC',1, ...
      'denoise',1,'wmnoise',wmnoise)); clear Yi
    dt(2) = dt(2) + toc(t); t = tic;

    % (3) central surface (without topology correction) and its reduction
    Ppp = fullfile(Pwork,[side{si} '.ppm.nii']); Pcs = fullfile(Pwork,[side{si} '.central.gii']);
    spm_write_vol(struct('fname',Ppp,'dim',size(Ypp),'mat',Mi,'dt',[16 0],'pinfo',[1;0;0], ...
      'descrip','position map (cat_io_dcm2bids_pp_surf)'),Ypp); clear Ypp
    cat_system(sprintf('CAT_VolMarchingCubes "%s" "%s" -thresh 0.5 -fast',Ppp,Pcs),0);
    CS = loadSurf(Pcs);
    nf(si) = max(20000,round(2 * sum(faceArea(CS)) / opts.surf.reduction));
    if nf(si) < size(CS.faces,1)
      cat_system(sprintf('CAT_SurfReduce -ratio %d "%s" "%s"',nf(si),Pcs,Pcs),0); CS = loadSurf(Pcs);
    end
    nf(si) = size(CS.faces,1); dt(3) = dt(3) + toc(t); t = tic;

    % (4) values of the vertices with the folding correction of the
    % thickness and intensity
    X   = [CS.vertices'; ones(1,size(CS.vertices,1))];
    Xi  = Mi \ X; th  = spm_sample_vol(Yth,Xi(1,:),Xi(2,:),Xi(3,:),1)'; clear Yth
    Xv  = V.mat \ X; in = spm_sample_vol(Ymn,Xv(1,:),Xv(2,:),Xv(3,:),1)';
    lab = spm_sample_vol(single(M.gm),Xv(1,:),Xv(2,:),Xv(3,:),0)';
    lab(spm_sample_vol(Ywall,Xv(1,:),Xv(2,:),Xv(3,:),0)' > 0 | th < 0.1) = 0; % (medial wall)
    th  = foldingCorrection(Pcs,th,1,cat_get_defaults('extopts.thick_limit'),fullfile(Pwork,[side{si} '.thickness']));
    th(th > 0) = th(th > 0) + thicknessOffset(res);
    in  = foldingCorrection(Pcs,in,0,10,fullfile(Pwork,[side{si} '.intensity']));
    va  = vertexArea(CS);

    % (5) absolute mean curvature and sulcal depth
    Pcu = fullfile(Pwork,[side{si} '.curvature']);
    cat_system(sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pcs,Pcu),0); 
    cu  = double(cat_io_FreeSurfer('read_surf_data',Pcu)); cu = cu(:);
    sd  = hullDepth(CS.vertices);

    % (6) regional and global values of the cortical vertices
    c = ismember(lab,cortex);
    G.(side{si}) = struct('thickness',mean(th(c)),'intensity',mean(in(c)),'curvature',mean(cu(c)),'depth',mean(sd(c)), ...
      'area',sum(va(c)),'area_surface',sum(va),'vertices',nnz(c),'faces',nf(si),'defects',topologyDefects(CS));
    for ri = find(ismember(T.id,cortex) & startsWith(T.key,side{si}(1)))'
      v = lab == T.id(ri);
      if any(v)
        R.(T.key{ri}) = struct('thickness',mean(th(v)),'intensity',mean(in(v)),'curvature',mean(cu(v)), ...
          'depth',mean(sd(v)),'area',sum(va(v)),'vertices',nnz(v)); 
      end
    end

    % (7) reduced surface with the values of the nearest vertex
    Pst = fullfile(Pwork,[side{si} '.reduced.gii']);
    cat_system(sprintf('CAT_SurfReduce -ratio %d "%s" "%s"',round(2 * opts.surf.reduce),Pcs,Pst),0);
    CR  = loadSurf(Pst); ni = dsearchn(CS.vertices,CR.vertices);
    GR  = gifti(struct('vertices',single(CR.vertices),'faces',int32(CR.faces),'cdata',single([th(ni) in(ni) cu(ni) sd(ni) lab(ni)])));
    save(GR,surfFile(Pjson,side{si}),'GZipBase64Binary'); %#ok<USENS> (save method of gifti)
    dt(4) = dt(4) + toc(t);
  end
  if ~keep, rmdir(Pwork,'s'); end

  J.method = struct('version',surfVersion,'resolution',res,'reduction',opts.surf.reduction,'faces',nf, ...
    'surface','marching cubes of the PBT position map without topology correction (CAT_VolMarchingCubes -fast, CAT_SurfReduce)', ...
    'thickness','projection-based thickness (cat_vol_pbtsimpleCS4) at the central surface','thickness_offset',thicknessOffset(res), ...
    'intensity','bias-corrected T1w scaled to CSF 1/3, GM 2/3, and WM 1 at the central surface', ...
    'folding_correction','CAT_SurfCorrectThicknessFolding (Demirci et al. 2025; thickness slope 1, intensity slope 0)', ...
    'curvature','absolute mean curvature in degrees (averaged over 3 mm, CAT_DumpCurv)', ...
    'depth','sulcal depth as distance to the convex hull of the hemisphere in mm', ...
    'atlas','lotter2025_aparc (CAT, nearest GM region of the vertices, cortical regions)','wmnoise',wmnoise, ...
    'reduce',opts.surf.reduce,'surface_data',{{'thickness','intensity','curvature','depth','ROIid'}}, ...
    'duration',round(toc(stime)),'durations',struct('preparation',round(dt(1)),'pbt',round(dt(2)), ...
    'surface',round(dt(3)),'values',round(dt(4))));
  J.software = cat_io_dcm2bids_helper('software');
  J.global   = G;
  J.regional = R;
end
% =========================================================================
function Ymn = scaleIntensity(Ym,Yp0)
%scaleIntensity. Bias-corrected image Ym scaled to CSF 1/3, GM 2/3, and WM
%  1 (as the T1w intensities of CAT) by the medians of the tissues of the
%  label map Yp0 (CSF 1, GM 2, WM 3) with linear interpolation (and
%  extrapolation above the WM), by the partial volume of the tissues if
%  there are no pure tissue voxels (e.g. a failed segmentation).
  t = [cat_stat_nanmedian(Ym(Yp0 > 0.95 & Yp0 < 1.05)) cat_stat_nanmedian(Ym(Yp0 > 1.95 & Yp0 < 2.05)) ...
       cat_stat_nanmedian(Ym(Yp0 > 2.95))];
  for c = find(~isfinite(t)), t(c) = cat_stat_nanmedian(Ym(abs(Yp0 - c) < 0.5)); end
  if ~all(isfinite(t)) || any(diff([0 t]) <= 0), error('no tissue contrast of the segmentation'); end
  Ymn = single(interp1([0 t],[0 1 2 3]/3,double(Ym),'linear','extrap'));
  Ymn(Ym < 0) = 0;
end
% =========================================================================
function bb = boundingBox(Y,margin)
%bboundingBox. Bounding box of the voxels of the mask Y ([first; last]
%  voxel of each dimension) with a margin (voxels).
  bb = zeros(2,3);
  for d = 1:3
    p = any(any(permute(Y,[d setdiff(1:3,d)]),2),3);
    bb(:,d) = [max(1,find(p,1) - margin(d)); min(size(Y,d),find(p,1,'last') + margin(d))];
  end
end
% =========================================================================
function [Yi,Mi] = interpCrop(Y,M,bb,res)
%interpCrop. Cubic interpolation of the bounding box bb of the image Y
%  (voxel-to-world matrix M) to an isotropic resolution res in mm with the
%  voxel-to-world matrix Mi of the interpolated image (voxel centers at
%  (i - 0.5) / f + 0.5 of the original voxels for the factor f).
  Yc = Y(bb(1,1):bb(2,1),bb(1,2):bb(2,2),bb(1,3):bb(2,3));
  f  = sqrt(sum(M(1:3,1:3).^2)) / res;
  n  = max(1,floor(size(Yc) .* f));
  g  = arrayfun(@(d) ((1:n(d)) - 0.5) / f(d) + 0.5,1:3,'UniformOutput',false);
  F  = griddedInterpolant(Yc,'cubic','nearest');
  Yi = single(F(g));
  Mi = M * [eye(3) bb(1,:)' - 1; 0 0 0 1] * [diag(1 ./ f) (0.5 - 0.5 ./ f)'; 0 0 0 1];
end
% =========================================================================
function CS = loadSurf(P)
%loadSurf. Vertices and faces (double) of the GIfTI surface P.
  S  = gifti(P);
  CS = struct('vertices',double(S.vertices),'faces',double(S.faces));
end
% =========================================================================
function a = faceArea(CS)
%faceArea. Area of the faces of the surface CS in mm2.
  v1 = CS.vertices(CS.faces(:,1),:);
  a  = 0.5 * sqrt(sum(cross(CS.vertices(CS.faces(:,2),:) - v1,CS.vertices(CS.faces(:,3),:) - v1,2).^2,2));
end
% =========================================================================
function va = vertexArea(CS)
%vertexArea. Area of the vertices of the surface CS in mm2 (a third of the
%  area of each face).
  a  = faceArea(CS);
  va = accumarray(CS.faces(:),repmat(a / 3,3,1),[size(CS.vertices,1) 1]);
end
% =========================================================================
function g = topologyDefects(CS)
%topologyDefects. Number of topology defects of the surface CS, i.e. the 
%  handles of its components by the Euler characteristic (V - E + F = 2 
%  - 2 handles of a closed surface without defects).
  E  = unique(sort([CS.faces(:,[1 2]); CS.faces(:,[2 3]); CS.faces(:,[3 1])],2),'rows');
  cc = conncomp(graph(E(:,1),E(:,2),[],size(CS.vertices,1))); u = unique(CS.faces(:)); % (without unused vertices)
  g  = round((2 * numel(unique(cc(u))) - (numel(u) - size(E,1) + size(CS.faces,1))) / 2);
end
% =========================================================================
function d = hullDepth(V)
%hullDepth. Sulcal depth of the vertices V (n x 3) as Euclidean distance 
%  to the convex hull of the vertices in mm, i.e. the minimal distance to 
%  the planes of the facets of the hull (for points inside a convex hull).
  K = convhulln(V); p = V(K(:,1),:);
  n = cross(V(K(:,2),:) - p,V(K(:,3),:) - p,2); n = n ./ max(eps,sqrt(sum(n.^2,2)));
  n = n .* sign(sum(n .* (p - mean(V,1)),2)); b = sum(n .* p,2);  % (outward normals and offsets)
  d = zeros(size(V,1),1);
  for k = 1:10000:size(V,1)
    i = k:min(size(V,1),k + 9999); d(i) = max(0,min(b' - V(i,:) * n',[],2)); 
  end
end
% =========================================================================
function v = foldingCorrection(Pcs,v,slope,limit,Pv)
%foldingCorrection. Correction of the folding bias of the values v of the
%  vertices of the surface Pcs (Demirci et al. 2025, see
%  CAT_SurfCorrectThicknessFolding) with the thickness-based weighting
%  slope (0 - none, e.g. for intensities) and the upper limit of the
%  values, by the temporary file Pv (uncorrected values for non-finite
%  results, e.g. of degenerated vertices).
  cat_io_FreeSurfer('write_surf_data',Pv,v);
  cat_system(sprintf('CAT_SurfCorrectThicknessFolding -slope %g -max %g "%s" "%s" "%s"',slope,limit,Pcs,Pv,Pv),0);
  vc = double(cat_io_FreeSurfer('read_surf_data',Pv)); vc = vc(:);
  ok = isfinite(vc); v(ok) = vc(ok);
end
