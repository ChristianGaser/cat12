function cat_io_dcm2bids_render(PBIDS,Preport,ropts)
%cat_io_dcm2bids_render. Render slices of each scan of a BIDS directory in MNI space.
%
%  Aim: The rendering supports the visual identification (or confirmation,
%  e.g. from the report table) of outliers in large sets of similar images.
%  Therefore, the scans are rendered separately for each image type
%  (dataset/datatype/subtype, e.g. anat/t1w or fmap/epi_part-phase with the
%  suffix, the image part and inversion, and the prefix of derivatives) and
%  are comparable, as all images are sampled at the same positions in MNI 
%  space using the registration of the session (catDCM2BIDSaffine_sub-*_
%  ses-*.json in the derivatives directory of the session, see 
%  writeSessionAffines in cat_io_dcm2bids_pp_anat), i.e., independent of the 
%  matrix size, slice
%  orientation, and the position of the origin. Spatially normalized
%  derivatives (prefix w or mw, e.g. mwc1) are already in MNI space and are
%  sampled without the registration. Scans without registration are not
%  comparable and are not rendered (black tile marked by 'no MNI').
%  Besides this focus on one
%  image type, a subject-specific report shows all scans of a subject.
%  The default user should not have to change anything (the rendering is
%  relatively fast and prepares all outputs at once), whereas the (rare)
%  expert user can adapt the general organization.
%
%  Slice selection (MNI space):
%    axial z=+10       .. basal ganglia, ventricles, frontal/occipital lobe
%    coronal y=0       .. subcortical structures
%    sagittal x=0      .. between the hemispheres, shows the origin offset
%                         (x=+30 for EPI images, i.e. func, dwi, and perf with
%                         their maps and composites, as the midline shows 
%                         little of their lower resolution, see epiSlices)
%  One slice is limited to understand the whole volume, whereas the three
%  orientations already show much better what is happening in the volume.
%  Further interesting slices (e.g. for the overview pages) are
%    coronal y=-60     .. symmetric cut through the cerebellum
%    sagittal x=+/-30  .. characterize the brain (frontal, parietal,
%                         temporal, and cerebellar areas, and especially the
%                         hippocampus) of both hemispheres to avoid
%                         side-wise bias
%
%  Outputs:
%   R1) Overview pages with m x n tiles of all subjects/sessions with one
%       slice per scan (one set of pages for each slice of ropts.slices):
%         [dataset/]datatype/subtype/render_<subtype>_slice_<orientation>_<axis><pos>_p##.png
%       This gives a very brief overview of many images but is biased by
%       the selected slice, i.e., it is useful for a rough global review.
%       Each tile shows the subject, the session with the quality rating 
%       SQR, and the scan (BIDS entities, e.g. the acquisition label with 
%       the series number and protocol).
%   R2) Scan-row pages with one row per scan with a header (BIDS filename
%       and an underline over the whole page), an information panel
%       (quality ratings, image and acquisition parameters, registration),
%       and the slices axial z=+10, sagittal x=0 (x=+30 for EPI), and 
%       coronal y=0:
%         [dataset/]datatype/subtype/render_<subtype>_rows_p##.png
%       with 6 scans per page.
%   R3) Subject and session reports with the scan rows (as R2) of selected images of a
%       subject (subject perspective, see renderSubjectOrder): original 
%       anatomical scan of the segmentation (e.g. T1w), its label map (c0),
%       normalized GM, DEC map with the bias-corrected T1w as underlay, and 
%       mean EPI with the network of the rs-fMRI (similar images of several
%       scans together), otherwise all  
%       scans, as one PNG file per page of each subject with the page number
%       first (i.e. the same page of all subjects one after another, e.g. 
%       the T1w rows of the first page) and a summary page with the best 
%       image of each selected type (p00), and as one PDF per session (scans
%       of the session in the order of the subject report):
%         [dataset/]subjects/p<page>_<sub>.png
%         [dataset/]sessions/render_<sub>_<ses>.pdf
%   rows/  .. slice row of each scan without text as separate image
%             (e.g. for a later interactive HTML report)
%   cache/ .. sampled slices of each scan (see renderScanSlices)
%  The A4 portrait format of SPM is used. It allows to compare multiple
%  pages on one screen and limits the number of output pages.
%  4D data is represented by the denoised average of its first volumes of
%  the QC (see cat_io_dcm2bids_qc), otherwise by its first volume. The 
%  denoised average is not rendered as own image type, realigned 4D data 
%  (prefix r) not at all, as it equals its scan, as well as field maps
%  (fmap) and the normalized FA (wFA, for voxel-based analyses). 
%  The DTI maps of cat_io_dcm2bids_pp_dwi (FA, the direction-encoded color
%  map DEC, and MD, AD, and RD for opts.dtimaps > 1) are rendered as own 
%  image types (e.g. dwi/dwi_FA) with fixed intensity windows (FA 0-1, 
%  diffusivities 0-0.003 mm2/s, see renderWindow), the DEC map in color 
%  (principal direction in head orientation: x red, y green, z blue, with 
%  the FA as brightness).
%  Composite images are rendered as own image types (e.g. for the subject
%  report): the DEC map on the darkened bias-corrected image of the 
%  segmentation of the session (e.g. dwi_DECt1w, FA as opacity, the DEC map
%  in the space of the T1w by the non-linear coregistration if available,
%  see toT1w in cat_io_dcm2bids_pp_anat), and the mean EPI with the network of the  
%  rs-fMRI (bold_conn: schematic network as vector graphics on the pages,
%  i.e. axial slices with the 7 networks of both hemispheres, sagittal 
%  slices with the parts of the networks (hemispheres averaged, ring for 
%  the connectivity between the hemispheres), both at schematic positions
%  where the edges do not pass other nodes, and the matrix of the 
%  connectivity within and between the
%  networks instead of the coronal slice, nodes by the connectivity within
%  and edges by the connectivity between them, and the segregation of each
%  network in the legend, see networkGeometry). DWI and fMRI scans with a     
%  coregistration to the T1w of their preprocessing (see coregAnat in 
%  cat_io_dcm2bids_pp_anat) use the registration of the T1w and this 
%  coregistration (e.g. after a repositioning of the subject), i.e. the
%  images are aligned by the rigid coregistration and the affine 
%  registration of the session (the DEC composite by the non-linear 
%  coregistration, without normalization), as well as other scans that 
%  are not aligned with the anatomical scan of the session (coregistration
%  of the alignment test in their QC JSON, see scannerAlignment in 
%  cat_io_dcm2bids_pp_anat).
%
%  cat_io_dcm2bids_render(PBIDS,Preport,ropts)
%
%  PBIDS     .. main BIDS directory (e.g. outdir/study/BIDS/BIDS)
%  Preport   .. output directory (e.g. outdir/study/BIDS-report/render/BIDS)
%  ropts     .. render options (defaults in brackets, see the opts.render
%                fields of cat_io_dcm2bids_defaults)
%   .source    .. 1-raw data, 2-derivatives (catDCM2BIDS), 3-both  (3)
%   .sessions  .. 1-only the first session per subject, 0-all sessions (0)
%   .slicemode .. 'mni'      - affine registration, i.e., also scaled (default)
%                 'mnirigid' - rigid part of the registration (original size)
%   .sort      .. order of the scans in R1 and R2:
%                 'name' - by subject and session name (default)
%                 'SQR'  - subjects by their worst quality rating (worst
%                          first), the sessions of a subject stay together
%                 (the order by rating can differ between image types)
%   .slices    .. R1 slices as [orientation position] rows with orientation
%                 3-axial, 2-coronal, 1-sagittal and MNI position in mm
%                 ([3 10; 2 0; 1 0], sagittal x=0 of EPI images at x=+30),
%                 empty for no R1 (also as structure array with the fields
%                 orient and pos)
%   .tiles     .. R1 tiles per page: [3 4], [4 5], [5 7], or [6 8] ([4 5])
%   .rows      .. R2 scan-row pages (1)
%   .subjects  .. R3 subject and session reports (1)
%   internal:
%   .infoside  .. side of the information panel in R2/R3: 'left' or 'right'
%                 ('left', i.e., a row reads like a table row: identity,
%                 quality, and then the images)
%   .rowimg    .. write the slice row of each scan as image (1)
%   .cache     .. use and write the slice cache (1)
%   .fov       .. in-plane field of view in mm (200)
%   .res       .. in-plane resolution in mm (default 300/dpi, i.e., 1 mm 
%                 for 300 dpi and 0.5 mm for 600 dpi)
%   .interp    .. interpolation: 0-nearest neighbor, 1-trilinear, 2-7 B-spline
%                 degree (default: trilinear for dpi<=300 (fast), otherwise 
%                 4th degree B-spline for the finer sampling (slower))
%   .dpi       .. print resolution (300)
%   .subfilter .. only the subject and session reports (R3) of these 
%                 subjects (BIDS labels, e.g. sub-01) with the slice cache 
%                 and the slice rows for R1/R2, e.g. of a subject package of
%                 a parallel process (see cat_io_dcm2bids_parallel), whereas 
%                 R1/R2 of all subjects are rendered after all packages ({})
%   .dbdir     .. database directory with the derivatives index, i.e. the 
%                 derivatives that were not exported to the BIDS 
%                 derivatives (opts.preprocessing = 1, see registerDerivative
%                 in cat_io_dcm2bids_helper) are rendered from the database
%                 as if they were exported ('' - catDCM2BIDSdb of the output
%                 directory of PBIDS if available)
%
%  The intensities of a scan are scaled for all its slices together
%  (1-99% percentile) to keep the slices of a scan comparable. The colored
%  lines show the coordinate planes of the original (scanner) space, i.e.,
%  x=0 (light red), y=0 (green), and z=0 (blue). Labels are colored by the
%  overall quality rating SQR. Unreadable (e.g. corrupted) files give
%  black tiles ('read error').

  % defaults of the render options (see cat_io_dcm2bids_defaults)
  def = cat_io_dcm2bids_defaults;
  def = def.opts.render;
  if ~exist('ropts','var'), ropts = struct(); end
  if isfield(ropts,'orient') || isfield(ropts,'slice')
    warning('cat_io_dcm2bids_render:oldopts', ...
      'The render options "orient" and "slice" were replaced by "slices" and are ignored.');
    ropts = rmfield(ropts,intersect(fieldnames(ropts),{'orient','slice'}));
  end
  % the slice list is handled separately as cat_io_checkinopt does not
  % replace matrices by structures (and keeps the default for empty input)
  if isfield(ropts,'slices'), slices = renderSliceList(ropts.slices); ropts = rmfield(ropts,'slices'); else, slices = def.slices; end
  ropts = cat_io_checkinopt(ropts,def);
  ropts.slices = slices;
  ropts.rows   = ropts.rows > 0;  % former option with two rows (2)
  % the sampling resolution follows the print resolution (300 dpi - 1 mm), 
  % i.e., only higher resolutions use the slower B-spline interpolation 
  if isempty(ropts.res),    ropts.res    = 300 / ropts.dpi; end
  if isempty(ropts.interp), ropts.interp = 1 + 3*(ropts.dpi > 300); end

  if numel(ropts.tiles)~=2 || ~any( all( [3 4; 4 5; 5 7; 6 8] == ropts.tiles(:)' , 2 ) )
    error('cat_io_dcm2bids_render:tiles','Tiles has to be [3 4], [4 5], [5 7], or [6 8].');
  end
  if ~isempty(ropts.slices) && ( size(ropts.slices,2)~=2 || ...
     ~all(ismember(ropts.slices(:,1),1:3)) || ~all(isfinite(ropts.slices(:,2))) )
    error('cat_io_dcm2bids_render:slices',['Slices have to be [orientation position] rows ' ...
      'with orientation 3 (axial), 2 (coronal), or 1 (sagittal).']);
  end
  if ~any(strcmp(ropts.slicemode,{'mni','mnirigid'}))
    error('cat_io_dcm2bids_render:slicemode','Slicemode has to be ''mni'' or ''mnirigid''.');
  end
  if ~any(strcmpi(ropts.sort,{'name','SQR'}))
    error('cat_io_dcm2bids_render:sort','Sort has to be ''name'' or ''SQR''.');
  end
  if ~any(strcmp(ropts.infoside,{'left','right'}))
    error('cat_io_dcm2bids_render:infoside','Infoside has to be ''left'' or ''right''.');
  end
  if ~exist(PBIDS,'dir'), return; end

  % derivatives index of the database (see sourceFile)
  if isempty(ropts.dbdir), ropts.dbdir = databaseDir(PBIDS); end
  Ddb = cat_io_dcm2bids_helper('derivativeIndex',ropts.dbdir); 
  sourceFile('',Ddb);

  % all slices (K) of R1 and of the scan rows of R2/R3 to sample each scan
  % only once (the sagittal slice of EPI images at x=+30, see epiSlices)
  rowslices = [3 10; 1 0; 2 0];
  userows   = ropts.rows || ropts.subjects || ropts.rowimg;
  K = unique([ropts.slices; rowslices(1:3*userows,:)],'rows','stable');
  if isempty(K), return; end

  % common render settings
  G.pw  = 21; G.ph = 29.7;                                % A4 portrait (as SPM) in cm
  G.ovc = [1 0.25 0.25; 0.25 1 0.25; 0.35 0.55 1];        % overlay colors of the x/y/z-planes
  G.orientnam = {'sagittal','coronal','axial'}; G.axisnam = 'xyz';
  MarkColor   = cat_io_colormaps('marks+',40);            % quality colors (as command line)
  G.col2mark  = @(val) min(1,max(0,MarkColor(min(size(MarkColor,1)-3, ...
    max(1,floor( val/9.5 * size(MarkColor,1)))),:)));


  %% find and describe images
  % (with the derivatives that are only in the database, see sourceFile)
  P  = cat_vol_findfiles(PBIDS,'*sub-*.nii*');
  Pv = keys(Ddb); Pv = Pv(strncmp(Pv,[PBIDS filesep],numel(PBIDS)+1)); 
  Pv = Pv(~cellfun(@(f) exist(f,'file')>0,Pv)); 
  P  = [P(:); Pv(:)];
  P = P( ~cellfun('isempty',regexp(P,'\.nii(\.gz)?$','once')) );
  P = P( cellfun('isempty',regexp(spm_file(P,'filename'),'^catDCM2BIDS(qc|dwi_wFA|dwi_DECt1w|dwimix_wFA|dwimix_DECt1w)_','once')) ); % QC maps (shown with their scan), wFA (analyses), and DECt1w (composite)
  P = P( cellfun('isempty',regexp(spm_file(P,'filename'),'^rsub-','once')) );          % realigned 4D data (equals its scan)
  P = P( cellfun('isempty',regexp(P,'[\\/]fmap[\\/][^\\/]*$','once')) );                   % field maps
  isderiv = cat_io_contains(P,[filesep 'derivatives' filesep]);
  P = P( ~isderiv | cat_io_contains(P,[filesep 'derivatives' filesep 'catDCM2BIDS' filesep]) ); % (no derivatives of other pipelines)
  isderiv = cat_io_contains(P,[filesep 'derivatives' filesep]);
  switch ropts.source
    case 1, P = P(~isderiv); isderiv = isderiv(~isderiv);
    case 2, P = P( isderiv); isderiv = isderiv( isderiv);
  end
  if ~isempty(ropts.subfilter) % (subject packages, see cat_io_dcm2bids_parallel)
    psub = cellfun(@(f) regexp(f,'sub-[^_]+','match','once'),spm_file(P,'filename'),'UniformOutput',false); 
    P    = P(ismember(psub,cellstr(ropts.subfilter))); isderiv = isderiv(ismember(psub,cellstr(ropts.subfilter))); 
  end
  if isempty(P), return; end

  F = struct('file',P,'name','','bidsname','','prefix','','sub','','ses','','dataset','', ...
    'datatype','','subtype','','group','','json','','SQR',nan,'QR',[],'T',[],'reg',[], ...
    'readerr',0,'dim',[],'vx',[],'nvol',nan,'qcmaps',struct(),'qcused',0,'pp',[],'kind','','under','','over','','Tu',[], ...
    'seg',[],'dwi','','c0','','foreign','','t1other','');
  for fi = 1:numel(P)
    [pp,ff]  = fileparts(P{fi});
    ff       = regexprep(ff,'\.nii$','');
    si       = strfind(ff,'sub-');
    prefix   = regexprep(regexprep(ff(1:si(1)-1),'_$',''),{'^catDCM2BIDSdwimix_','^catDCM2BIDS(dwi|rsfmri)_'},{'mix',''}); % (e.g. FA, mixFA of the mixed DTI, or mean)
    bidsname = ff(si(1):end);
    F(fi).name = ff; F(fi).bidsname = bidsname; F(fi).prefix = prefix;

    % BIDS entities
    tok = regexp(bidsname,'sub-([^_]+)','tokens','once'); F(fi).sub = ['sub-' tok{1}];
    tok = regexp(bidsname,'ses-([^_]+)','tokens','once');
    if ~isempty(tok), F(fi).ses = ['ses-' tok{1}]; end
    parts = strsplit(bidsname,'_');
    F(fi).datatype = spm_file(pp,'basename');
    F(fi).subtype  = lower(parts{end});
    % other image parts (e.g. phase) and inversions (MP2RAGE) as own types
    tok = regexp(bidsname,'_(part-(?:phase|real|imag)|inv-\d+)(?=_)','tokens');
    for ti = 1:numel(tok), F(fi).subtype = [F(fi).subtype '_' tok{ti}{1}]; end
    if ~isempty(prefix), F(fi).subtype = [F(fi).subtype '_' prefix]; end

    % dataset (protocol subdirectory) relative to the main BIDS directory
    si     = strfind(P{fi},[filesep 'sub-']);
    root   = P{fi}(1:si(1)-1);
    reldir = pp(numel(root)+2:end);
    if isderiv(fi), droot = root; rroot = fileparts(fileparts(root)); else, droot = cat_io_dcm2bids_helper('derivativesDir',root); rroot = root; end
    dataset = regexprep(root(min(numel(root),numel(PBIDS)+2):end), ...
      ['(^|' regexptranslate('escape',filesep) ')derivatives.*$'],'');
    if strcmp(root,PBIDS), dataset = ''; end
    F(fi).dataset = dataset;
    F(fi).group   = fullfile(dataset,F(fi).datatype,F(fi).subtype);
    F(fi).json    = fullfile(rroot,reldir,[bidsname '.json']);  % sidecar of the raw scan

    % session registration (world to MNI), whereas spatially normalized
    % derivatives (prefix w or mw, e.g. wc1/mwc1 of SPM or wp1/mwp1 of CAT)
    % are already in MNI space (identity)
    if ~isempty(regexp(prefix,'^m?w','once'))
      F(fi).T   = eye(4);
      F(fi).reg = struct('method','template','ll',nan,'imat',spm_imatrix(eye(4)),'SOR',nan);
    else
      if isempty(F(fi).ses), sesname = F(fi).sub; else, sesname = [F(fi).sub '_' F(fi).ses]; end
      Paff = fullfile(droot,fileparts(reldir),['catDCM2BIDSaffine_' sesname '.json']);
      if exist(Paff,'file')
        try
          A = cat_io_json(Paff);
          if strcmp(ropts.slicemode,'mni'), T = A.Affine; else, T = A.Rigid; end
          if isnumeric(T) && isequal(size(T),[4 4])
            F(fi).T   = double(T);
            F(fi).reg = struct('method',A.method,'ll',A.ll,'imat',spm_imatrix(double(A.Affine)),'SOR',nan);
            if isfield(A,'qualityrating') && isfield(A.qualityrating,'SOR') && isnumeric(A.qualityrating.SOR)
              F(fi).reg.SOR = double(A.qualityrating.SOR); % session orientation rating
            end
          end
        end
      end
    end

    % denoised average of 4D data of the QC (see cat_io_dcm2bids_qc), but
    % not of phase images (e.g. part-phase or phasediff) 
    if ~isderiv(fi) && isempty(regexp(bidsname,'_part-phase_|_phase(diff|1|2)$','once'))
      Pm = niiOf(fullfile(droot,reldir,['catDCM2BIDSqc_denoisedavg_' bidsname]));
      if ~isempty(Pm), F(fi).qcmaps.denoisedavg = Pm; end
    end

    % results of the DTI and rs-fMRI preprocessing of the scan (see 
    % cat_io_dcm2bids_pp_dwi/_rsfmri): the coregistration to the T1w, i.e.
    % the scan and its maps in the space of the scan use the registration of
    % the T1w and the coregistration (e.g. after a repositioning of the 
    % subject), and the global values (information panel)
    pps = {'dwi','rsfmri'}; if startsWith(prefix,'mix'), pps = {'dwimix'}; end % (mixed DTI, see cat_io_dcm2bids_pp_dwi)
    ppcoreg = false; 
    for pp = pps
      Pj = fullfile(droot,reldir,['catDCM2BIDS' pp{1} '_' bidsname '.json']);
      if exist(Pj,'file')
        try
          F(fi).pp = cat_io_json(Pj);
          if isempty(regexp(prefix,'^m?w','once')) && isfield(F(fi).pp,'coreg') && ...
              isnumeric(F(fi).pp.coreg) && isequal(size(F(fi).pp.coreg),[4 4])
            % (the T1w of another session, i.e. a session without anatomical
            % scan, with the registration of the session of the T1w, 
            % otherwise the registration of the session by the coregistration,
            % see coregSessionAffine in cat_io_dcm2bids_pp_anat)
            so = regexp(otherT1w(F(fi)),'_(ses-[^_]+)_','tokens','once');
            if ~isempty(so)
              T2 = sessionTransform(fullfile(droot,F(fi).sub,so{1},['catDCM2BIDSaffine_' F(fi).sub '_' so{1} '.json']),ropts);
              if ~isempty(T2), F(fi).T = T2 * double(F(fi).pp.coreg); ppcoreg = true; end
            elseif ~isempty(F(fi).T)
              F(fi).T = F(fi).T * double(F(fi).pp.coreg); ppcoreg = true;
            end
          end
        end
      end
    end

    % segmentation (TIV and relative volumes, see cat_io_dcm2bids_pp_anat)
    % of the label map and the normalized GM, and the diffusion directions
    % of the b-value shells of DWI scans and their maps (bval of the scan)
    if any(strcmp(prefix,{'c0','mwc1','neuromorphometrics'}))
      Ps = fullfile(droot,reldir,['catDCM2BIDSsegus_' bidsname '.json']);
      if exist(Ps,'file'), try, F(fi).seg = cat_io_json(Ps); end, end %#ok<TRYNC>
    end
    % label map of the segmented anatomical scan for the scan and its bias-
    % corrected image (intensity window of the WM, see renderScanSlices)
    if strcmp(F(fi).datatype,'anat') && any(strcmp(prefix,{'','m','neuromorphometrics'}))
      F(fi).c0 = niiOf(fullfile(droot,reldir,['c0' bidsname]));
    end
    % atlas label map of the segmentation (see regionalVolumes in 
    % cat_io_dcm2bids_pp_anat): colored labels on the bias-corrected image 
    % (otherwise the scan) with the registration of the session
    if strcmp(F(fi).datatype,'anat') && strcmp(prefix,'neuromorphometrics')
      F(fi).kind  = 'atlas'; F(fi).Tu = F(fi).T; 
      F(fi).under = niiOf(fullfile(droot,reldir,['m' bidsname]));
      if isempty(F(fi).under), F(fi).under = niiOf(fullfile(rroot,reldir,bidsname)); end
    end
    if strcmp(F(fi).datatype,'dwi')
      Pb = regexprep(F(fi).json,'\.json$','.bval');
      if exist(Pb,'file'), try, F(fi).dwi = shellString(load(Pb)); end, end %#ok<TRYNC>
    end

    % quality ratings and the coregistration of scans that are not aligned
    % with the anatomical scan of the session (e.g. after a repositioning of
    % the subject, see scannerAlignment in cat_io_dcm2bids_pp_anat), if the
    % scan has no coregistration of its preprocessing
    Pqc = fullfile(droot,reldir,['catDCM2BIDSqc_' bidsname '.json']);
    if exist(Pqc,'file')
      try
        QC = cat_io_json(Pqc);
        if ~ppcoreg && ~isempty(F(fi).T) && isempty(regexp(prefix,'^m?w','once')) && isfield(QC,'alignment') && ...
            isfield(QC.alignment,'coreg') && isnumeric(QC.alignment.coreg) && isequal(size(QC.alignment.coreg),[4 4])
          F(fi).T = F(fi).T * double(QC.alignment.coreg);
        end
        F(fi).QR = QC.qualityrating;
        SQR = QC.qualityrating.SQR;
        if isnumeric(SQR) && isscalar(SQR) && isfinite(SQR), F(fi).SQR = double(SQR); end % e.g. empty if QC failed
      end
    end
  end
  % derivatives of the database without their raw scan (e.g. of former 
  % runs with other BIDS names)
  idb = find(ismember({F.file},Pv)); 
  F(idb(~arrayfun(@(f) exist(f.json,'file')>0,F(idb)))) = []; 
  if isempty(F), return; end


  % sessions without anatomical scan: the T1w of another session (BIDS 
  % name) for all images of the session that rely on it, i.e. the DTI and 
  % rs-fMRI results of the session (method.T1w_other_session, see otherT1w) 
  % and the scans with the registration of the session by their 
  % coregistration (see coregSessionAffine in cat_io_dcm2bids_pp_anat), 
  % shown as warning (see anatOther)
  t1o = containers.Map('KeyType','char','ValueType','char');
  for fi = 1:numel(F), t1 = otherT1w(F(fi)); if ~isempty(t1), t1o([F(fi).dataset '|' F(fi).sub '|' F(fi).ses]) = t1; end, end
  for fi = 1:numel(F)
    key = [F(fi).dataset '|' F(fi).sub '|' F(fi).ses];
    if isKey(t1o,key) && (~isempty(otherT1w(F(fi))) || (isstruct(F(fi).reg) && startsWith(F(fi).reg.method,'coreg ')))
      F(fi).t1other = t1o(key);
    end
  end

  %% composite images (own image types, e.g. for the subject report): the
  % DEC map on the bias-corrected image of the segmentation of the session
  % (see method.T1w of the DTI) as underlay with its own registration (e.g.
  % dwi_DECt1w), where the DEC map in the space of the T1w (non-linear 
  % coregistration, see toT1w in cat_io_dcm2bids_pp_anat) is used if 
  % available (otherwise the DEC map of the scan with the rigid 
  % coregistration), and the mean EPI with the network of the rs-fMRI 
  % (bold_conn)
  Fc = F([]);
  for fi = 1:numel(F)
    pp = fileparts(F(fi).file);
    switch F(fi).prefix
      case {'DEC','mixDEC'} % (DEC map of a scan or of the mixed DTI)
        if isstruct(F(fi).pp) && isfield(F(fi).pp,'method') && isfield(F(fi).pp.method,'T1w') && ~isempty(F(fi).pp.method.T1w)
          ba = regexprep(char(F(fi).pp.method.T1w),'\.nii$','');
          Fm = F(strcmp({F.prefix},'m') & strcmp({F.bidsname},ba) & strcmp({F.dataset},F(fi).dataset));
          if ~isempty(Fm) && ~isempty(Fm(1).T)
            st = regexp(ba,'_([^_]+)$','tokens','once'); % suffix of the underlay (e.g. T1w)
            C  = composite(F(fi),'dec',Fm(1).file,'',[F(fi).prefix lower(st{1})]); C.Tu = Fm(1).T;
            C.c0 = Fm(1).c0; % (window of the underlay by its WM, see renderScanSlices)
            pf = 'catDCM2BIDSdwi'; if strcmp(F(fi).prefix,'mixDEC'), pf = 'catDCM2BIDSdwimix'; end
            Pt = niiOf(fullfile(pp,[pf '_DECt1w_' F(fi).bidsname]));
            if ~isempty(Pt), C.file = Pt; C.T = Fm(1).T; end % (DEC map in the space of the T1w)
            Fc(end+1) = C; %#ok<AGROW>
          end
        end
      case 'mean'
        Pc = fullfile(pp,['catDCM2BIDSrsfmri_conn_' F(fi).bidsname '.tsv']);
        if exist(sourceFile(Pc),'file'), Fc(end+1) = composite(F(fi),'conn','',Pc,'conn'); end %#ok<AGROW>
    end
  end
  F = [F(:); Fc(:)];


  %% render each group (image type)
  % To keep the graphics load low (the interactive desktop crashed with many
  % axes), all tiles of a page are composed into one RGB image that is shown
  % in one axes, and one invisible figure is used for all pages.
  n  = round(ropts.fov / ropts.res);                         % slice size in pixel
  fh = figure('Visible','off','Color','k','InvertHardcopy','off','MenuBar','none', ... % keep the black page and white text in print (e.g. R2023b inverts them)
    'ToolBar','none','Units','centimeters','Position',[1 1 G.pw G.ph], ...
    'PaperUnits','centimeters','PaperSize',[G.pw G.ph],'PaperPosition',[0 0 G.pw G.ph]);
  fhclean = onCleanup(@() close(fh));
  Pcache  = @(F) fullfile(Preport,F.group,'cache',[F.name '.mat']);

  [groups,~,gid] = unique({F.group});
  for gi = 1:numel(groups)
    Fg = renderSortScans(F(gid==gi),ropts);
    % slices of the image type (Kg) with the slices of R1 (k1) and of the scan rows (k2)
    Kg = unique(epiSlices(Fg(1),K),'rows','stable');
    [~,k1] = ismember(unique(epiSlices(Fg(1),ropts.slices),'rows','stable'),Kg,'rows');
    [~,k2] = ismember(epiSlices(Fg(1),rowslices),Kg,'rows');

    % prepare output directory (remove old pages, not for the subject 
    % packages that only fill the cache, see ropts.subfilter)
    Pout = fullfile(Preport,groups{gi});
    if ~exist(Pout,'dir'), mkdir(Pout); end
    if isempty(ropts.subfilter)
      Pold = cat_vol_findfiles(Pout,'render_*.png',struct('depth',0));
      for oi = 1:numel(Pold), delete(Pold{oi}); end
    end

    % sample all slices of each scan once (or load them from the cache)
    Y = zeros(n,n,size(Kg,1),numel(Fg),'uint8'); Yov = Y; Yc = cell(1,numel(Fg)); nnew = 0;
    for fi = 1:numel(Fg)
      S = renderScanSlices(Fg(fi),Kg,ropts,Pcache(Fg(fi)));
      Y(:,:,:,fi) = S.Y; Yov(:,:,:,fi) = S.Yov; Yc{fi} = S.C; nnew = nnew + S.new;
      Fg(fi).readerr = S.readerr; Fg(fi).dim = S.dim; Fg(fi).vx = S.vx; Fg(fi).nvol = S.nvol; Fg(fi).qcused = S.qc;

      % slice row without text (e.g. for a later HTML report)
      if ropts.rowimg
        Prow = fullfile(Pout,'rows',[Fg(fi).name '_row.png']);
        if S.new || ~exist(Prow,'file')
          if ~exist(fileparts(Prow),'dir'), mkdir(fileparts(Prow)); end
          Yr = cell(1,3);
          for ci = 1:3, Yr{ci} = renderTileRGB(S.Y(:,:,k2(ci)),S.Yov(:,:,k2(ci)),G.ovc,colorSlice(Yc,fi,k2(ci))); end
          imwrite(cat(2,Yr{:}),Prow);
        end
      end
    end

    np1 = 0; np2 = 0;
    if ~isempty(ropts.subfilter), continue; end % (R1/R2 of all subjects after all packages)
    if ~isempty(k1)
      np1 = renderOverviewPages(fh,Fg,Y(:,:,k1,:),Yov(:,:,k1,:),colorSlices(Yc,k1),Kg(k1,:),groups{gi},Pout,ropts,G);
    end
    if ropts.rows
      np2 = renderRowPages(fh,Fg,Y(:,:,k2,:),Yov(:,:,k2,:),colorSlices(Yc,k2),Kg(k2,:),groups{gi}, ...
        fullfile(Pout,['render_' Fg(1).subtype '_rows_p%02d.png']),0,ropts,G);
    end
    cat_io_cprintf('blue',sprintf('  Rendered %4d images (%4d new) in %3d overview and %3d row page(s): %s\n', ...
      numel(Fg), nnew, np1, np2, Pout));
  end

  % remove the pages of image types of former runs that are not rendered
  % anymore (e.g. maps of former versions), i.e. other directories with a
  % slice cache (and their then empty datatype directories), only if all
  % images were rendered (raw images and derivatives, ropts.source = 3)
  Pc = {}; if ropts.source == 3 && isempty(ropts.subfilter), Pc = cat_vol_findfiles(Preport,'cache',struct('dirs',1)); end
  for ci = 1:numel(Pc)
    Pg = fileparts(Pc{ci});
    if strncmp(Pg,Preport,numel(Preport)) && ~any(endsWith(Pg,strcat(filesep,groups))) % (keeps all rendered groups)
      rmdir(Pg,'s'); Pp = fileparts(Pg); D = dir(Pp);
      if ~strcmp(Pp,Preport) && all(ismember({D.name},{'.','..'})), rmdir(Pp); end
    end
  end


  %% subject and session reports (all scans of a subject and of each of
  % its sessions, slices from the cache)
  if ropts.subjects
    [subs,~,sid] = unique(strcat({F.dataset},'|',{F.sub}));
    for si = 1:numel(subs)
      Fs = renderSubjectOrder(F(sid==si),ropts);
      Y  = zeros(n,n,3,numel(Fs),'uint8'); Yov = Y; Yc = cell(1,numel(Fs)); Ks = zeros(3,2,numel(Fs));
      for fi = 1:numel(Fs)
        Ks(:,:,fi) = epiSlices(Fs(fi),rowslices); % (slices of each scan, e.g. sagittal x=+30 for EPI)
        S = renderScanSlices(Fs(fi),Ks(:,:,fi),ropts,Pcache(Fs(fi)));
        Y(:,:,:,fi) = S.Y; Yov(:,:,:,fi) = S.Yov; Yc{fi} = S.C;
        Fs(fi).readerr = S.readerr; Fs(fi).dim = S.dim; Fs(fi).vx = S.vx; Fs(fi).nvol = S.nvol; Fs(fi).qcused = S.qc;
      end
      Pout = fullfile(Preport,Fs(1).dataset,'subjects');
      if ~exist(Pout,'dir'), mkdir(Pout); end
      Pold = [cat_vol_findfiles(Pout,['p*_' Fs(1).sub '.png'],struct('depth',0)); {fullfile(Pout,['render_' Fs(1).sub '.pdf'])}];
      for oi = 1:numel(Pold), if exist(Pold{oi},'file'), delete(Pold{oi}); end, end % (former runs and versions)
      % summary page with the best image of each selected type (page 0) and
      % the pages of the subject report as separate files with the page 
      % number first (e.g. p01_<sub>.png, i.e. the same page of all subjects
      % one after another)
      bi = summaryScans(Fs);
      if ~isempty(bi)
        renderRowPages(fh,Fs(bi),Y(:,:,:,bi),Yov(:,:,:,bi),Yc(bi),Ks(:,:,bi),[fullfile(Fs(1).dataset,Fs(1).sub) ' summary'], ...
          fullfile(Pout,['p00_' Fs(1).sub '.png']),1,ropts,G);
      end
      np = renderRowPages(fh,Fs,Y,Yov,Yc,Ks,fullfile(Fs(1).dataset,Fs(1).sub),fullfile(Pout,['p%02d_' Fs(1).sub '.png']),1,ropts,G);
      cat_io_cprintf('blue',sprintf('  Rendered %4d images in %3d subject page(s) and a summary: %s\n', ...
        numel(Fs), np, fullfile(Pout,['p*_' Fs(1).sub '.png'])));
      % session reports (scans of each session in the order of the subject report)
      [uses,~,ise] = unique({Fs.ses},'stable');
      for ei = find(~cellfun('isempty',uses(:)'))
        k    = find(ise(:)' == ei);
        % sessions without anatomical scan: the anatomical rows (scan, label
        % map, normalized GM) of the T1w of another session that the DTI
        % and rs-fMRI of this session used (with a label and a warning, see
        % anatOther and renderInfoText)
        Fk = Fs;
        if ~any(strcmp({Fs(k).datatype},'anat'))
          t1 = unique(arrayfun(@otherT1w,Fs(k),'UniformOutput',false)); t1 = t1(~cellfun('isempty',t1));
          if ~isempty(t1)
            ka = find(strcmp({Fs.datatype},'anat') & strcmp({Fs.bidsname},t1{1}) & ~strcmp({Fs.ses},uses{ei}));
            if ~isempty(ka), [Fk(ka).foreign] = deal(uses{ei}); k = [ka k]; end
          end
        end
        Pout = fullfile(Preport,Fs(1).dataset,'sessions');
        if ~exist(Pout,'dir'), mkdir(Pout); end
        Ppdf = fullfile(Pout,['render_' Fs(1).sub '_' uses{ei} '.pdf']);
        if exist(Ppdf,'file'), delete(Ppdf); end
        np = renderRowPages(fh,Fk(k),Y(:,:,:,k),Yov(:,:,:,k),Yc(k),Ks(:,:,k),fullfile(Fs(1).dataset,[Fs(1).sub '_' uses{ei}]), ...
          Ppdf,1,ropts,G);
        cat_io_cprintf('blue',sprintf('  Rendered %4d images in %3d session page(s): %s\n',numel(k),np,Ppdf));
      end
    end
  end
end
% =========================================================================
function bi = summaryScans(Fs)
%summaryScans. Scans of the summary page of a subject (see 
%  renderSubjectOrder): the best scan (lowest SQR, the first for equal
%  ratings) of each selected image type, i.e. the anatomical scan of the
%  segmentation, its label map, atlas labels, and normalized GM, the DEC map (of the 
%  mixed DTI if available), and the network of the rs-fMRI.
  st  = regexprep({Fs(~cellfun('isempty',regexp({Fs.subtype},'_c0$','once'))).subtype},'_c0$','');
  dec = {Fs(strcmp({Fs.kind},'dec')).subtype}; mix = dec(contains(dec,'mix')); if ~isempty(mix), dec = mix; end
  types = {st, strcat(st,'_c0'), strcat(st,'_neuromorphometrics'), strcat(st,'_mwc1'), dec, {'bold_conn'}}; bi = [];
  for ti = 1:numel(types)
    k = find(ismember({Fs.subtype},types{ti}));
    if isempty(k), continue; end
    q = [Fs(k).SQR]; q(isnan(q)) = inf; [~,m] = min(q); bi(end+1) = k(m); %#ok<AGROW>
  end
end
% =========================================================================
function str = shellString(b)
%shellString. Number of volumes of each b-value shell, e.g. '3x b0, 8x 
%  b50, 64x b1200' (b-values rounded to 50, values within 10% combined, and
%  shells with less than 3 volumes counted to the closest shell).
  b = 50 * round(b(:)' / 50); u = unique(b); c = arrayfun(@(v) nnz(b == v),u);
  k = 1;
  while k < numel(u) % (combine values within 10%)
    if u(k) >= 200 && u(k+1) - u(k) <= 0.1 * u(k+1)
      [~,m] = max(c(k:k+1)); v = u(k + m - 1); u(k) = v; c(k) = c(k) + c(k+1); u(k+1) = []; c(k+1) = [];
    else
      k = k + 1;
    end
  end
  while numel(u) > 1 && any(c < 3) % (small shells to the closest shell)
    [~,k] = min(c); d = abs(u - u(k)); d(k) = inf; [~,m] = min(d);
    c(m) = c(m) + c(k); u(k) = []; c(k) = [];
  end
  str = strjoin(arrayfun(@(i) sprintf('%dx b%d',c(i),u(i)),1:numel(u),'UniformOutput',false),', ');
end
% =========================================================================
function np = renderOverviewPages(fh,Fg,Y,Yov,Yc,K,group,Pout,ropts,G)
%renderOverviewPages. R1 pages with one slice per scan in m x n tiles.
%  Y and Yov are the slices (n x n x slices x scans) of the slices K, Yc 
%  the color slices of each scan (empty for gray-scale images). Each
%  tile shows the subject, the session with the quality rating SQR, and 
%  the scan (BIDS entities besides the subject, session, and suffix, e.g.
%  the acquisition label with series number and protocol) in its header.
  nx = ropts.tiles(1); ny = ropts.tiles(2); nt = nx*ny;
  fs = [8 7 6 5]; fs = fs( all( [3 4; 4 5; 5 7; 6 8] == ropts.tiles(:)' , 2 ) ); % font size
  mg = 0.002; ht = 0.02; % margin and title height (normalized)
  tw = (1 - 2*mg) / nx * G.pw; th = (1 - 2*mg - ht) / ny * G.ph; % tile size in cm
  n  = size(Y,1);
  % free header above each slice for about three text lines to reduce the
  % overlay of the labels (sub, ses with SQR, scan) with the image
  hd = round(3 * 1.2 * fs/72*2.54 * n/min(tw,th));            % header in pixel
  cw = round(max(n, (n+hd)*tw/th)); ch = round(max(n+hd, n*th/tw)); % tile size in pixel
  np = ceil(numel(Fg)/nt);

  % labels: subject, session with the quality rating, and scan (e.g. 
  % acq-007T1wMPR, i.e. several scans or protocols of a subject) in rows
  lab = cell(1,numel(Fg));
  for fi = 1:numel(Fg)
    lab{fi} = {Fg(fi).sub};
    if isempty(Fg(fi).T), lab{fi}{1} = [lab{fi}{1} ' (no MNI)']; end
    if Fg(fi).readerr,    lab{fi}{1} = [lab{fi}{1} ' (read error)']; end
    str = Fg(fi).ses;
    if ~isnan(Fg(fi).SQR), str = strtrim(sprintf('%s  SQR %0.1f',str,Fg(fi).SQR)); end
    if ~isempty(str), lab{fi}{end+1} = str; end
    scan = regexprep(regexprep(Fg(fi).bidsname,'^sub-[^_]+_(ses-[^_]+_)?',''),'_?[^_]+$','');
    if ~isempty(scan), lab{fi}{end+1} = scan; end
    % (within the tile width, e.g. long subject names)
    maxc = floor(tw / (0.6 * fs/72*2.54));
    for li = 1:numel(lab{fi})
      if numel(lab{fi}{li}) > maxc, lab{fi}{li} = [lab{fi}{li}(1:maxc-1) char(8230)]; end
    end
  end

  Gn = cell(1,numel(Fg)); % schematic networks of the rs-fMRI (vector graphics, see networkGeometry)
  for fi = find(strcmp({Fg.kind},'conn')), Gn{fi} = networkGeometry(K,Fg(fi).over,ropts,n); end
  Lv = arrayfun(@(F) overlayLines(F.T,K,ropts,n),Fg,'UniformOutput',false); % lines of the scanner planes

  % pages for each orientation and slice position
  for ki = 1:size(K,1)
    for pgi = 1:np
      clf(fh);

      % page title and legend of the overlay lines
      ax = axes('Parent',fh,'Position',[mg 1-mg-ht 1-2*mg ht],'Color','k','Visible','off');
      text(ax,0,0.5,sprintf('%s  (%s, %s, %s=%+0.0fmm, page %d/%d)', group, ...
        G.orientnam{K(ki,1)}, ropts.slicemode, G.axisnam(K(ki,1)), K(ki,2), pgi, np), ...
        'Color','w','FontSize',fs+1,'Interpreter','none','VerticalAlignment','middle');
      renderLegend(ax,1,0.5,fs,G.ovc);

      % compose all tiles of the page into one RGB image with the overlay
      ntp  = min(nt, numel(Fg) - (pgi-1)*nt);
      page = zeros(ny*ch, nx*cw, 3, 'single');
      pos  = zeros(ntp,2);
      for ti = 1:ntp
        fi = (pgi-1)*nt + ti;
        % slice below the header and centered in the remaining space
        [ix,iy]   = ind2sub([nx ny],ti);
        pos(ti,:) = [(ix-1)*cw, (iy-1)*ch]; % upper left corner of the tile
        page(pos(ti,2) + hd + floor((ch-hd-n)/2) + (1:n), pos(ti,1) + floor((cw-n)/2) + (1:n), :) = ...
          renderTileRGB(Y(:,:,ki,fi),0*Yov(:,:,ki,fi),G.ovc,colorSlice(Yc,fi,ki)) * ~(strcmp(Fg(fi).kind,'conn') && K(ki,1) == 2); % (lines as vector graphics, matrix of the networks)
      end
      ax = axes('Parent',fh,'Position',[mg mg 1-2*mg 1-2*mg-ht],'Color','k');
      image(ax,page); axis(ax,'image','off'); hold(ax,'on');
      for ti = 1:ntp % (lines of the scanner planes and networks as vector graphics on the tiles)
        fi = (pgi-1)*nt + ti; tx = pos(ti,1) + floor((cw-n)/2); ty = pos(ti,2) + hd + floor((ch-hd-n)/2);
        if ~(strcmp(Fg(fi).kind,'conn') && K(ki,1) == 2), drawOverlayLines(ax,Lv{fi},ki,tx,ty,G.ovc); end
        if ~isempty(Gn{fi}), drawNetworkVector(ax,Gn{fi},ki,tx,ty); end
      end

      % labels (in page coordinates) colored by the quality rating (white if not available)
      for ti = 1:ntp
        fi = (pgi-1)*nt + ti;
        if isnan(Fg(fi).SQR), col = [1 1 1]; else, col = G.col2mark(Fg(fi).SQR); end
        text(ax,pos(ti,1) + cw*0.02,pos(ti,2) + 1,lab{fi},'Color',col,'FontSize',fs, ...
          'Interpreter','none','VerticalAlignment','top');
      end

      % print page and give the graphics system time to finish
      Ppng = fullfile(Pout,sprintf('render_%s_slice_%s_%s%+0.0f_p%02d.png', ...
        Fg(1).subtype,G.orientnam{K(ki,1)},G.axisnam(K(ki,1)),K(ki,2),pgi));
      print(fh, Ppng, '-dpng', sprintf('-r%d',ropts.dpi));
      drawnow;
    end
  end
  np = np * size(K,1);
end
% =========================================================================
function np = renderRowPages(fh,Fg,Y,Yov,Yc,K,ttl,Ppage,showgroup,ropts,G)
%renderRowPages. Pages with one row of 3 slices per scan (R2 and R3), with
%  the color slices Yc of each scan (empty for gray-scale images).
%  Each scan starts with a header with its BIDS filename (and its image
%  type for showgroup=1) and an underline over the whole page, followed by
%  an information panel (see renderInfoText) and the slices K (common, or
%  per scan as 3 x 2 x scans, e.g. the sagittal slice x=+30 of EPI images,
%  see epiSlices). The tile size is given by the page height with 6 scans 
%  per page, whereas the
%  panel uses the remaining width. Ppage is either a PNG filename pattern
%  with the page number (%02d) or a PDF file to which all pages are
%  appended (one file, e.g. per subject).
  nps = 6;                                       % scans per page
  n   = size(Y,1); fs = 7;
  mg  = 0.2; ht = 0.6; hh = 0.55;                % margin, page title, and scan header in cm
  s   = (G.ph - 2*mg - ht)/nps - hh;             % tile size in cm
  ppc = n / s;                                   % pixel per cm
  W   = round((G.pw - 2*mg)*ppc); H = round((G.ph - 2*mg)*ppc);  % page in pixel
  pht = round(ht*ppc); phh = round(hh*ppc);      % title and header in pixel
  bh  = phh + n;                                 % height of a scan block in pixel
  if strcmp(ropts.infoside,'left'), x0 = W - 3*n; xi = 0; else, x0 = 0; xi = 3*n; end % slices and panel
  np  = ceil(numel(Fg)/nps);
  Kf  = @(fi) K(:,:,min(fi,size(K,3)));           % slices of a scan
  ism = @(fi,ci) strcmp(Fg(fi).kind,'conn') && K(ci,1,min(fi,size(K,3))) == 2; % matrix of the networks instead of the coronal slice
  ispdf = strcmpi(spm_file(Ppage,'ext'),'pdf');
  Gn = cell(1,numel(Fg)); % schematic networks of the rs-fMRI (vector graphics, see networkGeometry)
  for fi = find(strcmp({Fg.kind},'conn')), Gn{fi} = networkGeometry(Kf(fi),Fg(fi).over,ropts,n); end
  Lv = arrayfun(@(fi) overlayLines(Fg(fi).T,Kf(fi),ropts,n),1:numel(Fg),'UniformOutput',false); % lines of the scanner planes
  lh = 1.02 * fs/72*2.54 * ppc;                  % line height of the information panel (about the font size)

  for pgi = 1:np
    clf(fh);
    nsp  = min(nps, numel(Fg) - (pgi-1)*nps);
    page = zeros(H,W,3,'single');
    for si = 1:nsp
      fi = (pgi-1)*nps + si;
      y0 = pht + (si-1)*bh;                        % top of the scan block
      page(y0 + phh - (2:-1:1), :, :) = 0.6;       % underline over the whole page
      for ci = 1:3
        page(y0 + phh + (1:n), x0 + (ci-1)*n + (1:n), :) = renderTileRGB(Y(:,:,ci,fi),0*Yov(:,:,ci,fi),G.ovc,colorSlice(Yc,fi,ci)); % (lines as vector graphics)
        if ism(fi,ci), page(y0 + phh + (1:n), x0 + (ci-1)*n + (1:n), :) = 0; end % (matrix of the networks)
      end
    end
    ax = axes('Parent',fh,'Position',[mg/G.pw mg/G.ph 1-2*mg/G.pw 1-2*mg/G.ph],'Color','k');
    image(ax,page); axis(ax,'image','off'); hold(ax,'on');
    for si = 1:nsp % (lines of the scanner planes and networks as vector graphics on the tiles)
      fi = (pgi-1)*nps + si;
      for ci = 1:3
        if ~ism(fi,ci), drawOverlayLines(ax,Lv{fi},ci,x0 + (ci-1)*n,pht + (si-1)*bh + phh,G.ovc); end
        if ~isempty(Gn{fi}), drawNetworkVector(ax,Gn{fi},ci,x0 + (ci-1)*n,pht + (si-1)*bh + phh); end
      end
    end

    % page title and legend of the overlay lines
    text(ax,1,pht/2,sprintf('%s  (%s, page %d/%d)',ttl,ropts.slicemode,pgi,np), ...
      'Color','w','FontSize',fs+1,'Interpreter','none','VerticalAlignment','middle');
    renderLegend(ax,W,pht/2,fs,G.ovc);

    % scan header, slice positions, and information panel
    for si = 1:nsp
      fi = (pgi-1)*nps + si;
      y0 = pht + (si-1)*bh;
      if isnan(Fg(fi).SQR), col = [1 1 1]; else, col = G.col2mark(Fg(fi).SQR); end
      text(ax,1,y0 + phh/2,Fg(fi).name,'Color',col,'FontSize',fs,'FontWeight','bold', ...
        'Interpreter','none','VerticalAlignment','middle');
      if showgroup
        text(ax,W,y0 + phh/2,Fg(fi).group,'Color',[.6 .6 .6],'FontSize',fs, ...
          'Interpreter','none','VerticalAlignment','middle','HorizontalAlignment','right');
      end
      Ki = Kf(fi);
      for ci = 1:3
        sl = sprintf('%s=%+0.0f',G.axisnam(Ki(ci,1)),Ki(ci,2)); if ism(fi,ci), sl = 'networks (mean z)'; end
        text(ax,x0 + (ci-1)*n + 3,y0 + phh + n - 2,sl,'Color',[.6 .6 .6], ...
          'FontSize',fs-1,'Interpreter','none','VerticalAlignment','bottom');
      end
      % session without anatomical scan: warning on top of the sagittal 
      % slice of all images that rely on the T1w of another session (see 
      % anatOther)
      ao = anatOther(Fg(fi)); ci = find(Ki(:,1) == 1,1); 
      if ~isempty(ao) && ~isempty(ci)
        text(ax,x0 + (ci-1)*n + n/2,y0 + phh + 2,ao,'Color',[1 0.3 0.3],'FontSize',fs-1,'Interpreter','tex', ...
          'HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','k','Margin',1);
      end
      ti = renderInfoText(Fg(fi),G.col2mark); % (one text per line for the line height)
      for li = 1:numel(ti)
        text(ax,xi + 3,y0 + phh + 3 + (li-1)*lh,ti{li},'Color','w','FontSize',fs,'Interpreter','tex','VerticalAlignment','top');
      end
    end

    if ispdf
      exportgraphics(fh, Ppage, 'Append', true, 'ContentType', 'vector', ...
        'Resolution', ropts.dpi, 'BackgroundColor', 'k');
    else
      Pp = Ppage; if contains(Ppage,'%'), Pp = sprintf(Ppage,pgi); end % (pattern with the page number)
      print(fh, Pp, '-dpng', sprintf('-r%d',ropts.dpi));
    end
    drawnow;
  end
end
% =========================================================================
function txt = renderInfoText(F,col2mark)
%renderInfoText. TeX lines of the information panel of a scan in R2/R3:
%  subject and session (labels bold), quality ratings (the measures of the
%  image quality to the SQR, and further ratings), acquisition (scanner 
%  and timing of the sidecar of the raw scan, only for scans, their 
%  bias-corrected images, and maps of DWI and fMRI, i.e. not for the 
%  segmentation), image size, and the results of the preprocessing, with
%  the MNI registration (with the session orientation rating SOR, the 
%  coverage BBR, and the phase-encoding reversal PRR) after the quality.
  tex    = @(s) regexprep(s,'([_\\{}^])','\\$1');
  colstr = @(v) sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}',col2mark(v));
  wht    = '\color[rgb]{1 1 1}';
  red    = '\color[rgb]{1 0.3 0.3}';
  if isnan(F.SQR), c = wht; else, c = colstr(F.SQR); end

  % subject and session (labels bold)
  txt = {[c 'sub-{\bf' tex(regexprep(F.sub,'^sub-','')) '}']};
  if ~isempty(F.ses), txt{1} = [txt{1} '  ses-{\bf' tex(regexprep(F.ses,'^ses-','')) '}']; end

  % quality: the measures of the image quality to the SQR (bold) and the
  % further ratings (protocol, coverage, and phase-encoding reversal)
  rt  = @(fn) [colstr(F.QR.(fn)) tex(fn) sprintf(' %0.1f',F.QR.(fn))];
  isr = @(fn) isstruct(F.QR) && isfield(F.QR,fn) && isnumeric(F.QR.(fn)) && isscalar(F.QR.(fn)) && isfinite(F.QR.(fn));
  qm  = {'BSM','ISR','NSR','WSM','RES'}; qm = qm(cellfun(isr,qm));
  str = strjoin(cellfun(rt,qm,'UniformOutput',false),'  ');
  if ~isnan(F.SQR)
    if ~isempty(str), str = [str '  ' wht '\rightarrow '];  end
    str = [str colstr(F.SQR) sprintf('{\\bfSQR %0.1f}',F.SQR)];
  end
  if ~isempty(str), txt{end+1} = str; end
  % registration (rotation of the scanner space to MNI) before the further
  % ratings (e.g. the coverage BBR)
  if ~isempty(F.reg) && strcmp(F.reg.method,'template')
    txt{end+1} = sprintf('%sMNI: template space (spatially normalized)',wht);
  elseif ~isempty(F.reg)
    ro = round(F.reg.imat(4:6)*1800/pi)/10; ro(ro==0) = 0; % avoid "-0"
    txt{end+1} = sprintf('%sMNI: %s, rotation %+0.1f %+0.1f %+0.1f%s',wht,tex(F.reg.method),ro,char(176));
    if isfinite(F.reg.SOR), txt{end} = sprintf('%s  %sSOR %0.1f',txt{end},colstr(F.reg.SOR),F.reg.SOR); end
  end
  if isempty(F.T), txt{end+1} = [red 'no MNI registration']; end
  % further ratings: coverage and phase-encoding reversal at the end of
  % the registration line, the others (e.g. PFR) in their own line
  if isstruct(F.QR)
    fn = setdiff(fieldnames(F.QR),[qm {'BSM','ISR','NSR','WSM','RES','SQR','vx_vol'}],'stable'); fn = fn(cellfun(isr,fn));
    fg = fn(ismember(fn,{'BBR','PRR'})); fn = fn(~ismember(fn,{'BBR','PRR'}));
    if ~isempty(fg) && (~isempty(F.reg) || isempty(F.T))
      txt{end} = [txt{end} '  ' strjoin(cellfun(rt,fg(:)','UniformOutput',false),'  ')];
    elseif ~isempty(fg)
      fn = [fg; fn];
    end
    if ~isempty(fn), txt{end+1} = strjoin(cellfun(rt,fn(:)','UniformOutput',false),'  '); end
  end
  if F.readerr,    txt{end+1} = [red 'read error']; end

  % acquisition (scanner and timing of the raw scan) of scans, their 
  % bias-corrected images, and maps of DWI and fMRI (not the segmentation)
  if ~(strcmp(F.datatype,'anat') && ~any(strcmp(F.prefix,{'','m'})))
    J = struct();
    if exist(F.json,'file'), try, J = cat_io_json(F.json); end, end %#ok<TRYNC>
    str = {};
    sc  = '';
    if isfield(J,'Manufacturer') && ischar(J.Manufacturer), sc = [sc J.Manufacturer ' ']; end
    if isfield(J,'MagneticFieldStrength') && isnumeric(J.MagneticFieldStrength), sc = [sc sprintf('%gT ',J.MagneticFieldStrength)]; end
    if isfield(J,'MRAcquisitionType') && ischar(J.MRAcquisitionType), sc = [sc J.MRAcquisitionType]; end
    if ~isempty(strtrim(sc)), str{end+1} = tex(strtrim(sc)); end
    tim = {'RepetitionTime','TR'; 'EchoTime','TE'; 'InversionTime','TI'};
    for i = 1:size(tim,1)
      if isfield(J,tim{i,1}) && isnumeric(J.(tim{i,1})) && ~isempty(J.(tim{i,1}))
        str{end+1} = sprintf('%s %0.5g ms',tim{i,2},round(1000*J.(tim{i,1})(1),2)); %#ok<AGROW>
      end
    end
    if isfield(J,'FlipAngle') && isnumeric(J.FlipAngle) && ~isempty(J.FlipAngle), str{end+1} = sprintf('flip %g%s',J.FlipAngle(1),char(176)); end
    if ~isempty(str), txt{end+1} = [wht strjoin(str,', ')]; end
  end

  % image (and QC maps)
  if ~isempty(F.dim)
    txt{end+1} = sprintf('%smatrix: %dx%dx%d (%d volume(s)), voxel: %0.2fx%0.2fx%0.2f mm',wht,F.dim(1:3),F.nvol,F.vx);
  end
  if F.qcused, txt{end+1} = [wht 'denoised average of the QC']; end

  % segmentation (label map and normalized GM) and diffusion directions
  sg = @(f) fieldValue(F.seg,f,nan);
  if isstruct(F.seg) && strcmp(F.prefix,'c0')
    txt{end+1} = sprintf('%sTIV %0.0f ml, rGMV %0.1f%%, rWMV %0.1f%%, rCSFV %0.1f%%, ll %0.0f',wht,sg('TIV'), ...
      100*sg('rGMV'),100*sg('rWMV'),100*sg('rCSFV'),fieldValue(fieldValue(F.seg,'qc',struct()),'TPMll',nan));
  elseif isstruct(F.seg) && strcmp(F.prefix,'mwc1')
    txt{end+1} = sprintf('%sTIV %0.0f ml, rGMV %0.1f%%, aGMV %0.0f ml, ll %0.0f',wht,sg('TIV'),100*sg('rGMV'), ...
      sg('aGMV'),fieldValue(fieldValue(F.seg,'qc',struct()),'TPMll',nan));
  elseif strcmp(F.kind,'atlas')
    nr = numel(fieldnames(fieldValue(F.seg,'regional',struct())));
    txt{end+1} = sprintf('%satlas: %s (CAT), %d regions',wht,tex(F.prefix),nr);
  end
  if startsWith(F.prefix,'mix') && isstruct(F.pp) && isfield(F.pp,'method') && isfield(F.pp.method,'mix')
    txt{end+1} = sprintf('%smixed DTI of %d scans: %d volumes (b <= 1300)',wht, ...
      1 + numel(fieldValue(F.pp.method.mix,'further',[])),fieldValue(F.pp.method,'volumes',nan));
  elseif ~isempty(F.dwi)
    txt{end+1} = [wht 'directions: ' F.dwi];
  end

  % DTI and rs-fMRI: coregistration to the T1w (Dice of the brain masks),
  % global values, and the colors of the DEC map, segmentation, and network
  if isstruct(F.pp)
    if isfield(F.pp,'method') && isfield(F.pp.method,'coregistration') && isfield(F.pp.method.coregistration,'Dice')
      Qc = F.pp.method.coregistration; dn = nan; if isfield(Qc,'Dice_nonlinear'), dn = Qc.Dice_nonlinear; end
      cc = wht; if max([Qc.Dice Qc.Dice_scanner dn]) < 0.8, cc = red; end
      txt{end+1} = sprintf('%sDice to T1w %0.2f (scanner %0.2f)',cc,Qc.Dice,Qc.Dice_scanner);
      if isfinite(dn) && isfield(Qc,'displacement')
        txt{end} = sprintf('%s, non-linear %0.2f (%0.1f mm)',txt{end},dn,Qc.displacement(1));
        if isfield(Qc,'method') && ischar(Qc.method) && contains(Qc.method,'not used'), txt{end} = [txt{end} ' not used']; end
      end
    end
    if isfield(F.pp,'global') && isfield(F.pp.global,'FA_WM')
      g = F.pp.global;
      txt{end+1} = sprintf('%sWM: FA %0.3f, MD %0.2f, AD %0.2f, RD %0.2f (10^{-3} mm^2/s)',wht,g.FA_WM,1000*[g.MD_WM g.AD_WM g.RD_WM]);
    end
    if isfield(F.pp,'qc') && isfield(F.pp.qc,'FD_mean') && isfield(F.pp,'connectivity')
      txt{end+1} = sprintf('%sFD %0.2f mm, tSNR %0.0f, GCOR %0.3f, segregation %0.2f',wht,F.pp.qc.FD_mean, ...
        F.pp.qc.tSNR,F.pp.connectivity.GCOR,F.pp.connectivity.segregation);
    end
  end
  if any(strcmp(F.prefix,{'DEC','mixDEC'})) || strcmp(F.kind,'dec')
    st = regexp(F.under,'_([^_\\/]+)\.nii(\.gz)?$','tokens','once'); un = '';
    if strcmp(F.kind,'dec') && ~isempty(st), un = [' on ' tex(st{1})]; end % (underlay: bias-corrected T1w)
    txt{end+1} = ['\color[rgb]{1 0.3 0.3}x left-right  \color[rgb]{0.3 1 0.3}y post.-ant.  ' ...
      '\color[rgb]{0.45 0.6 1}z inf.-sup.' wht ' (FA)' un];
  end
  if strcmp(F.kind,'conn')
    % segregation of each network in its color (legend of the nodes) and the
    % description of the nodes and edges (see drawNetwork)
    N = networkNodes; nets = {'Vis','SMN','DAN','VAN','Lim','FPN','DMN'}; % (visual, somatomotor, dorsal and
    lc = 0.65 * N.col + 0.35;                        % ventral attention, limbic, frontoparietal control, default mode; lighter text)
    R  = readConn(F.over); sg = nan(1,7); if ~isempty(R) && size(R,1) == numel(N.net), sg = networkSegregation(R,N.net); end
    txt{end+1} = [wht 'seg. ' strjoin(arrayfun(@(i) sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s %s',lc(i,:),nets{i}, ...
      regexprep(sprintf('%0.2f',sg(i)),'^(-?)0\.','$1.')),1:numel(nets),'UniformOutput',false),' ')]; % (e.g. .83)
    txt{end+1} = [wht 'axial: networks per side, sagittal: parts (ring: left-right), matrix: z'];
  end
end
% =========================================================================
function Fs = renderSubjectOrder(Fs,ropts)
%renderSubjectOrder. Scans of one subject in the subject report: the 
%  selected images of the subject perspective (if available) in the order
%  (1) original anatomical scans of the type of the segmentation (e.g. 
%  T1w), (2) label map of the segmentation (c0) and the atlas labels 
%  (neuromorphometrics), (3) normalized GM (mwc1), (4) DEC map with the 
%  bias-corrected image as underlay, and (5)  
%  mean EPI with the network, i.e. similar images of several scans 
%  together (e.g. first all original T1w), ordered by session and BIDS 
%  name. Without these images (e.g. without anatomical scan), all scans by
%  session, datatype (anat, dwi, func, perf, others), and BIDS name, where
%  derivatives follow their scan (first the bias-corrected image m and the
%  segmentation c0, the DTI maps FA, DEC, MD, AD, and RD, then the others 
%  by prefix). Only the first session for ropts.sessions.
  st  = regexprep({Fs(~cellfun('isempty',regexp({Fs.subtype},'_c0$','once'))).subtype},'_c0$',''); % (e.g. t1w)
  if isempty(st), st = {'t1w'}; end
  dec = unique({Fs(strcmp({Fs.kind},'dec')).subtype});
  sel = [st strcat(st,'_c0') strcat(st,'_neuromorphometrics') strcat(st,'_mwc1') dec {'bold_conn'}];
  nr  = [ones(1,numel(st)) 2*ones(1,numel(st)) 3*ones(1,numel(st)) 4*ones(1,numel(st)) 5*ones(1,numel(dec)) 6];
  [iss,ri] = ismember({Fs.subtype},sel);
  if any(iss)
    Fs = Fs(iss); row = nr(ri(iss));
    [~,so] = sort(strcat(cellstr(num2str(row(:)))','|',{Fs.ses},'|',{Fs.bidsname})); Fs = Fs(so);
    if ropts.sessions && ~isempty(Fs), Fs = Fs( strcmp({Fs.ses},Fs(1).ses) ); end
    return
  end
  dto = {'anat','dwi','func','perf'};
  [~,dti] = ismember({Fs.datatype},dto); dti(dti==0) = numel(dto) + 1;
  pro = {'','m','c0','FA','DEC','MD','AD','RD'}; 
  [~,pri] = ismember({Fs.prefix},pro); pri(pri==0) = numel(pro) + 1;
  key = strcat({Fs.ses},'|',cellstr(num2str(dti(:)))','|',{Fs.datatype},'|',{Fs.bidsname},'|', ...
    cellstr(num2str(pri(:)))','|',{Fs.prefix});
  [~,so] = sort(key); Fs = Fs(so);
  if ropts.sessions && ~isempty(Fs)
    Fs = Fs( strcmp({Fs.ses},Fs(1).ses) );
  end
end
% =========================================================================
function S = renderScanSlices(F,K,ropts,Pcache)
%renderScanSlices. Sampled slices of a scan (with cache).
%  Each scan is read only once and all slices K ([orientation position]
%  rows) are sampled in MNI space and scaled together to uint8. The result
%  is saved in Pcache and reused if the source file, its QC maps, the 
%  registration, and the sampling parameters were not changed and the cache
%  includes all requested slices. This avoids reading all images again (e.g.
%  for large datasets) and keeps the slices for further outputs (e.g. HTML).
%  With the denoised average of the QC (F.qcmaps, see cat_io_dcm2bids_qc),
%  it is sampled instead of the first volume of 4D data. 
%  The DTI maps have fixed intensity windows (see renderWindow) and the DEC
%  map (volumes red, green, blue) and the composite images give color 
%  slices S.C (uint8 n x n x 3 x slices, otherwise empty).
%  S.Y (uint8 slices), S.Yov (overlay labels: 1-3 lines of the x/y/z-planes),
%  S.K (slices), S.readerr, header information of the scan (S.dim, S.vx, 
%  S.nvol), S.qc (1 - denoised average), and S.new (newly sampled).
  n   = round(ropts.fov / ropts.res);
  srcdate = dateOf(F.file);
  qcmaps = struct2cell(F.qcmaps); qcdate = nan(1,numel(qcmaps)); 
  for mi = 1:numel(qcmaps), D = dir(qcmaps{mi}); if ~isempty(D), qcdate(mi) = D(1).datenum; end; end
  par = struct('version',13,'srcdate',srcdate,'qc',{fieldnames(F.qcmaps)'},'qcdate',qcdate, ...
    'kind',F.kind,'part',{{F.under F.over}},'partdate',[dateOf(F.under) dateOf(F.over)], ...
    'T',F.T,'Tu',F.Tu,'fov',ropts.fov,'res',ropts.res,'interp',ropts.interp, ...
    'c0',F.c0,'c0date',dateOf(F.c0),'other',~isempty(otherT1w(F)));

  % cached slices
  if ropts.cache && exist(Pcache,'file')
    try
      C = load(Pcache,'S','par');
      [isK,ki] = ismember(K,C.S.K,'rows');
      if isequal(C.par,par) && all(isK)
        S = C.S; S.K = K; S.Y = S.Y(:,:,ki); S.Yov = S.Yov(:,:,ki); S.new = 0;
        if ~isempty(S.C), S.C = S.C(:,:,:,ki); end
        return
      end
    end
  end

  S = struct('K',K,'Y',zeros(n,n,size(K,1),'uint8'),'Yov',zeros(n,n,size(K,1),'uint8'),'C',[], ...
    'readerr',0,'dim',[],'vx',[],'nvol',nan,'qc',0,'new',1);
  try
    % header (for .nii.gz without decompressing the whole file), where 4D 
    % data is represented by its denoised average or first volume (see
    % below), i.e. only images with up to 3 volumes are read completely 
    % (spm_vol of a .nii.gz reads all volumes as double, e.g. 2.6 GB for a
    % 112x112x64x403 BOLD), e.g. the DEC map
    Pf = sourceFile(F.file);
    [~,dim,mat] = cat_io_dcm2bids_helper('niiHeaderInfo',Pf);
    S.nvol = dim(4); S.dim = dim(1:3); S.vx = sqrt(sum(mat(1:3,1:3).^2));
    if S.nvol > 3
      V = struct('dim',S.dim,'mat',mat); Vall = V; % (geometry of the volume, see renderSlice)
    else
      txt = evalc('V = spm_vol(Pf);'); %#ok<NASGU> % avoid gz-messages
      Vall = V; V = V(1);
    end
  catch
    S.readerr = 1; % e.g. corrupted files
  end
  if ~S.readerr && ~isempty(F.T)
    try
      % denoised average of 4D data of the QC (if it has finite values)
      Yvol = []; 
      if isfield(F.qcmaps,'denoisedavg')
        [Yvol,Vu] = renderAverage(F.qcmaps.denoisedavg); 
        if ~isempty(Yvol), V = Vu; S.qc = 1; end
      end
      % first volume of other 4D data (without reading the whole file)
      if isempty(Yvol) && S.nvol > 3
        try
          Yvol = cat_io_dcm2bids_helper('niiVolume',Pf,1); 
        catch % (e.g. NIfTI-2)
          txt = evalc('V = spm_vol(Pf);'); V = V(1); %#ok<NASGU>
        end
      end

      % B-spline coefficients of the volume (non-finite values set to 0) 
      % for higher interpolation degrees or volumes in memory, otherwise 
      % direct sampling
      if ropts.interp > 1 || ~isempty(Yvol)
        if isempty(Yvol), Yvol = spm_read_vols(V); end
        Yvol(~isfinite(Yvol)) = 0;
        C    = spm_bsplinc(Yvol,[ropts.interp*[1 1 1] 0 0 0]); clear Yvol
      else
        C    = [];
      end
      Ys   = zeros(n,n,size(K,1),'single');
      for ki = 1:size(K,1)
        [Ys(:,:,ki),S.Yov(:,:,ki)] = renderSlice(C,V,F.T,K(ki,1),K(ki,2),ropts);
      end
      clear C
      % common intensity scaling of all slices of the scan (fixed for DTI 
      % maps), for anatomical scans with segmentation the 98th percentile of
      % the WM (label map c0 > 2.5) as upper limit, i.e. bright fat or skull
      % do not darken the brain and the noise of the WM stays visible
      clim = renderWindow(F);
      if isempty(clim)
        Yv   = Ys(Ys~=0);
        clim = [cat_stat_prctile(Yv,1) cat_stat_prctile(Yv,99)];
        if ~isempty(F.c0) && ~any(strcmp(F.kind,{'dec','atlas'})) % (composites: the underlay, see below)
          try
            wm = sampleSlices(F.c0,1,F.T,K,ropts) > 2.5 & Ys ~= 0;
            if nnz(wm) > 100, clim(2) = prctile(Ys(wm),98); end
          end
        end
        if ~(numel(clim)==2 && all(isfinite(clim)) && diff(clim)>0), clim = [0 max(1,max(Ys(:)))]; end
      end
      S.Y  = uint8(255 * min(1,max(0,(Ys - clim(1)) / diff(clim))));

      % color slices of the DEC map of DTI (volumes red, green, and blue)
      if strcmp(F.prefix,'DEC') && numel(Vall) == 3
        S.C = zeros(n,n,3,size(K,1),'uint8');
        for ci = 1:3
          Cc = [];
          if ropts.interp > 1
            Yd = spm_read_vols(Vall(ci)); Yd(~isfinite(Yd)) = 0; Cc = spm_bsplinc(Yd,[ropts.interp*[1 1 1] 0 0 0]); 
          end
          for ki = 1:size(K,1)
            S.C(:,:,ci,ki) = uint8(255 * min(1,max(0,(renderSlice(Cc,Vall(ci),F.T,K(ki,1),K(ki,2),ropts) - clim(1)) / diff(clim))));
          end
        end
      end

      % composite images (see the help of this file)
      switch F.kind
        case 'dec'  % DEC map (direction as color, FA as opacity) on the darkened underlay (own registration,
                    % window by the WM of its label map as for anatomical scans, e.g. for bright fat)
          u  = sampleSlices(F.under,1,F.Tu,K,ropts); g = windowSlices(u);
          if ~isempty(F.c0)
            try
              wm = sampleSlices(F.c0,1,F.Tu,K,ropts) > 2.5 & u ~= 0; v = u(u ~= 0);
              if nnz(wm) > 100, g = min(1,max(0,(u - prctile(v,1)) / max(eps,prctile(u(wm),98) - prctile(v,1)))); end
            end
          end
          g  = 0.6 * g;
          d  = cat(4,sampleSlices(F.file,1,F.T,K,ropts),sampleSlices(F.file,2,F.T,K,ropts),sampleSlices(F.file,3,F.T,K,ropts));
          fa = sqrt(sum(d.^2,4)); a = 0.9 * min(1,max(0,(fa - 0.1) / 0.35)); % (no colors for low FA, e.g. GM)
          S.C = uint8(255 * permute((1-a).*g + a.*(d ./ max(eps,max(d,[],4))),[1 2 4 3]));
        case 'atlas' % labels of the atlas (colors of its table, nearest neighbor, 50% opacity) on the 
                     % underlay (window by the WM of the label map of the segmentation as for anatomical scans)
          u  = sampleSlices(F.under,1,F.Tu,K,ropts); g = windowSlices(u);
          if ~isempty(F.c0)
            try
              wm = sampleSlices(F.c0,1,F.Tu,K,ropts) > 2.5 & u ~= 0; v = u(u ~= 0);
              if nnz(wm) > 100, g = min(1,max(0,(u - prctile(v,1)) / max(eps,prctile(u(wm),98) - prctile(v,1)))); end
            end
          end
          L  = round(sampleSlices(F.file,1,F.T,K,ropts,0));
          [lc,la] = atlasColors(F.prefix,max(L(:)));
          c  = zeros([size(L) 3],'single'); 
          for d = 1:3, cd = lc(:,d); c(:,:,:,d) = reshape(cd(max(1,L(:))),size(L)); end
          a  = 0.5 * single(L > 0) .* reshape(la(max(1,L(:))),size(L)); % (white matter uncolored)
          S.C = uint8(255 * permute((1-a).*g + a.*c,[1 2 4 3]));
        case 'conn' % mean EPI (the network is drawn as vector graphics on the pages, see networkGeometry)
          S.C = [];
      end
    catch
      S.readerr = 1; % e.g. unsupported datatype
    end
  end
  if ropts.cache
    try
      if ~exist(fileparts(Pcache),'dir'), mkdir(fileparts(Pcache)); end
      save(Pcache,'S','par');
    end
  end
end
% =========================================================================
function [Y,V] = renderAverage(P)
%renderAverage. Denoised average P of 4D data of the QC (see 
%  cat_io_dcm2bids_qc) as image Y with the header V (same space as the 
%  first volume of the scan, empty without finite values, e.g. of a failed
%  QC).
  txt = evalc('V = spm_vol(P);'); %#ok<NASGU>
  Y   = single(spm_read_vols(V)); 
  if ~any(isfinite(Y(:)) & Y(:)~=0), Y = []; V = []; end
end
% =========================================================================
function Fg = renderSortScans(Fg,ropts)
%renderSortScans. Order of the scans of one image type.
%  By subject and session name, optionally only the first session of each
%  subject. For sort='SQR', the subjects are ordered by their worst
%  rating (unrated subjects at the end), whereas the sessions of a subject
%  stay together in their original order.
  [~,so] = sort( strcat({Fg.sub},'_',{Fg.ses},'_',{Fg.file}) ); Fg = Fg(so);
  if ropts.sessions
    [~,ui] = unique({Fg.sub},'first'); Fg = Fg(sort(ui));
  end
  if strcmpi(ropts.sort,'SQR')
    [subs,~,si] = unique({Fg.sub},'stable');
    wSQR = -inf(1,numel(subs));
    for i = 1:numel(subs)
      v = [Fg(si==i).SQR]; v = v(~isnan(v));
      if ~isempty(v), wSQR(i) = max(v); end
    end
    [~,so]   = sort(wSQR,'descend');   % stable for equal ratings
    rk(so)   = 1:numel(so);
    [~,fo]   = sort(rk(si));            % stable, i.e., sessions keep their order
    Fg = Fg(fo);
  end
end
% =========================================================================
function sl = renderSliceList(s)
%renderSliceList. Slices as [orientation position] rows given by a matrix,
%  a (GUI) structure array or cell with the fields orient and pos, or a
%  cell of [orientation position] pairs.
  if isempty(s)
    sl = zeros(0,2);
  elseif isnumeric(s)
    sl = double(s);
  elseif isstruct(s)
    sl = [reshape([s.orient],[],1) reshape([s.pos],[],1)];
  elseif iscell(s)
    sl = cell2mat(cellfun(@renderSliceList,s(:),'UniformOutput',false));
  else
    error('cat_io_dcm2bids_render:slices','Unknown slice definition.');
  end
end
% =========================================================================
function K = epiSlices(F,K)
%epiSlices. Slices K ([orientation position] rows) of the image F, where 
%  the sagittal slice x=0 of EPI images (datatypes func, dwi, and perf, 
%  also their maps and composites) is replaced by x=+30, as the midline 
%  shows little of their lower resolution (e.g. the cortex, insula, and
%  hippocampus at x=+30).
  if any(strcmp(F.datatype,{'func','dwi','perf'})), K(K(:,1) == 1 & K(:,2) == 0,2) = 30; end
end
% =========================================================================
function t1 = otherT1w(F)
%otherT1w. BIDS name of the T1w of another session that the DTI or rs-fMRI
%  of the image F used (session without anatomical scan, see 
%  method.T1w_other_session of cat_io_dcm2bids_pp_dwi/_rsfmri), otherwise ''.
  t1 = '';
  if isstruct(F.pp) && isfield(F.pp,'method') && isfield(F.pp.method,'T1w_other_session') && ...
      isequal(F.pp.method.T1w_other_session,true) && isfield(F.pp.method,'T1w')
    t1 = regexprep(char(F.pp.method.T1w),'\.nii(\.gz)?$','');
  end
end
% =========================================================================
function s = anatOther(F)
%anatOther. Warning of a session without anatomical scan in the image F
%  (TeX, e.g. 'no ses-w0001 T1w -> use ses-w0000 T1w', otherwise ''): the 
%  anatomical rows of the other session in the session report (F.foreign)
%  and all images of the session that rely on the T1w of the other session
%  (F.t1other, i.e. the DTI and rs-fMRI results and the scans with the
%  registration by their coregistration).
  s = ''; t = {}; 
  if ~isempty(F.foreign)
    t = [{F.foreign F.ses} regexp(F.bidsname,'_([^_]+)$','tokens','once')]; 
  elseif ~isempty(F.t1other)
    t = [{F.ses} regexp(F.t1other,'_(ses-[^_]+)_','tokens','once') regexp(F.t1other,'_([^_]+)$','tokens','once')]; 
  end
  if numel(t) == 3, s = sprintf('no %s %s \\rightarrow use %s %s',t{1},t{3},t{2},t{3}); end
end
% =========================================================================
function T = sessionTransform(Paff,ropts)
%sessionTransform. Transformation from world to MNI space of the session
%  registration Paff (catDCM2BIDSaffine_*.json, affine or rigid part by
%  ropts.slicemode), empty if not available.
  T = [];
  if ~exist(Paff,'file'), return; end
  try
    A = cat_io_json(Paff);
    if strcmp(ropts.slicemode,'mni'), T = A.Affine; else, T = A.Rigid; end
    if isnumeric(T) && isequal(size(T),[4 4]), T = double(T); else, T = []; end
  catch
    T = [];
  end
end
% =========================================================================
function Yrgb = renderTileRGB(Y,Yov,ovc,Yc)
%renderTileRGB. Gray-scale slice (uint8) with colored overlay lines (Yov 
%  1-3), or the color slice Yc (uint8 n x n x 3, e.g. of the DEC map of 
%  DTI) with the overlay lines.
  if exist('Yc','var') && ~isempty(Yc)
    Yrgb = single(Yc)/255;
  else
    Yg   = single(Y)/255; Yrgb = cat(3,Yg,Yg,Yg);
  end
  for ci = 1:3
    Yc = Yrgb(:,:,ci);
    for k = 1:3, Yc(Yov==k) = ovc(k,ci); end
    Yrgb(:,:,ci) = Yc;
  end
end
% =========================================================================
function Yc = colorSlice(C,fi,k)
%colorSlice. Color slice k of scan fi of the color slices C (cell of
%  n x n x 3 x slices or empty for gray-scale images).
  if isempty(C{fi}), Yc = []; else, Yc = C{fi}(:,:,:,k); end
end
% =========================================================================
function C = colorSlices(C,k)
%colorSlices. Color slices k of all scans (see colorSlice).
  for fi = 1:numel(C), if ~isempty(C{fi}), C{fi} = C{fi}(:,:,:,k); end, end
end
% =========================================================================
function C = composite(F,kind,under,over,tag)
%composite. Composite image of the image F (see the composite images in 
%  cat_io_dcm2bids_render): kind (dec, conn) with the underlay and 
%  overlay files (or the connectivity matrix) as own image type with the 
%  tag (e.g. dwi_DEC to dwi_DECt1w).
  C = F; C.kind = kind; C.under = under; C.over = over; C.prefix = tag; C.name = [F.name '_' tag];
  C.subtype = regexprep(F.subtype,['_' F.prefix '$'],['_' tag]);
  C.group   = fullfile(F.dataset,F.datatype,C.subtype);
end
% =========================================================================
function P = niiOf(P)
%niiOf. Existing NIfTI of the name P without extension (.nii.gz or .nii,
%  also in the database, see sourceFile), otherwise ''.
  if exist(sourceFile([P '.nii.gz']),'file'), P = [P '.nii.gz']; elseif exist(sourceFile([P '.nii']),'file'), P = [P '.nii']; else, P = ''; end
end
% =========================================================================
function P = sourceFile(P,D)
%sourceFile. File of the BIDS derivative P that is read, i.e. its database
%  file for a derivative that was not exported (opts.preprocessing = 1, 
%  see registerDerivative in cat_io_dcm2bids_helper), otherwise P. All 
%  names and outputs of the rendering use the BIDS name P, only the reading
%  uses the database file. The derivatives index D (map of the BIDS 
%  derivatives to their database files, see derivativeIndex) is set once 
%  per call of cat_io_dcm2bids_render by sourceFile('',D).
  persistent index
  if nargin > 1, index = D; return; end
  if ~isempty(index) && ~isempty(P) && isKey(index,P) && ~exist(P,'file'), P = index(P); end
end
% =========================================================================
function Pdb = databaseDir(PBIDS)
%databaseDir. Database directory catDCM2BIDSdb of the output directory of
%  the BIDS directory PBIDS (outdir/BIDSdir/BIDSsubdir), otherwise ''.
  Pdb = ''; P = PBIDS;
  for i = 1:3
    P = fileparts(P);
    if exist(fullfile(P,'catDCM2BIDSdb'),'dir'), Pdb = fullfile(P,'catDCM2BIDSdb'); return; end
  end
end
% =========================================================================
function d = dateOf(P)
%dateOf. Modification date of the file P (also in the database, see 
%  sourceFile, NaN if not available).
  d = nan; if ~isempty(P), D = dir(sourceFile(P)); if ~isempty(D), d = D(1).datenum; end, end
end
% =========================================================================
function Ys = sampleSlices(P,vol,T,K,ropts,hold)
%sampleSlices. Slices K of the volume vol of the image P (linear 
%  interpolation, or hold 0 for nearest neighbor, e.g. of labels, 0 
%  outside) by the transformation T to MNI space (see renderSlice).
  if ~exist('hold','var'), hold = 1; end
  n = round(ropts.fov / ropts.res); Ys = zeros(n,n,size(K,1),'single');
  P = sourceFile(P); 
  txt = evalc('V = spm_vol(P);'); V = V(min(vol,numel(V))); o = ropts; o.interp = hold; %#ok<NASGU>
  for ki = 1:size(K,1), Ys(:,:,ki) = renderSlice([],V,T,K(ki,1),K(ki,2),o); end
end
% =========================================================================
function [lc,la] = atlasColors(name,n)
%atlasColors. RGB colors lc (0-1) of the labels 1 to n of the atlas name 
%  of CAT (ROIcolor of its table templates_MNI152NLin2009cAsym/<name>.csv,
%  gray for labels without color) and their opacity la (0 for the white 
%  matter, i.e. White Matter in the name, otherwise 1).
  persistent cache
  if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
  if ~isKey(cache,name)
    C  = {}; Pc = fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym',[name '.csv']);
    if exist(Pc,'file')
      C = regexp(fileread(Pc),'\r?\n','split'); C = C(~cellfun('isempty',strtrim(C)));
      C = cellfun(@(l) strsplit(l,';'),C,'UniformOutput',false);
    end
    id = []; rgb = zeros(0,3); wm = [];
    if numel(C) > 1
      hdr = C{1}; ic = find(strcmp(hdr,'ROIcolor')); ii = find(strcmp(hdr,'ROIid')); in = find(strcmp(hdr,'ROIname'));
      if ~isempty(ic) && ~isempty(ii)
        for k = 2:numel(C)
          v = str2double(strsplit(strtrim(C{k}{ic}),' '));
          if numel(v) == 3 && all(isfinite(v)), id(end+1) = str2double(C{k}{ii}); rgb(end+1,:) = v / 255; end %#ok<AGROW>
          if ~isempty(in) && ~isempty(regexpi(C{k}{in},'white.?matter','once')), wm(end+1) = str2double(C{k}{ii}); end %#ok<AGROW>
        end
      end
    end
    cache(name) = struct('id',id,'rgb',rgb,'wm',wm);
  end
  A  = cache(name); lc = 0.6 * ones(max(1,n),3); la = ones(max(1,n),1);
  ok = A.id >= 1 & A.id <= n; lc(A.id(ok),:) = A.rgb(ok,:);
  w  = A.wm(A.wm >= 1 & A.wm <= n); la(w) = 0;
end
% =========================================================================
function g = windowSlices(Y)
%windowSlices. Slices Y scaled to 0-1 by the 1st and 99th percentile.
  v = Y(Y ~= 0); clim = [cat_stat_prctile(v,1) cat_stat_prctile(v,99)];
  if ~(numel(clim) == 2 && all(isfinite(clim)) && diff(clim) > 0), clim = [0 max(1,max(Y(:)))]; end
  g = min(1,max(0,(Y - clim(1)) / diff(clim)));
end
% =========================================================================
function G = networkGeometry(K,Pconn,ropts,n)
%networkGeometry. Schematic network of the connectivity matrix Pconn (Fisher
%  z of the 100 cortical parcels of the Schaefer2018 atlas, see 
%  cat_io_dcm2bids_pp_rsfmri) for the slices K of n x n pixel at schematic
%  positions (see networkNodes): axial slices with the 14 nodes of the 7 
%  networks of both hemispheres, sagittal slices with the parts of the 
%  networks (both hemispheres averaged) and a ring for the connectivity 
%  between the hemispheres (homotopic), and coronal slices with the matrix
%  of the mean connectivity within and between the networks (both 
%  hemispheres averaged). The node diameter (8
%  to 16 mm) and opacity (0.5 to 1) are given by the internal connectivity 
%  (mean z of the parcels within a hemisphere, 0.2 to 0.8), the edge width
%  (0 to 8 mm) and opacity (0 to 1) by the mean z between the nodes (0.1
%  to 0.7, within the hemispheres for the parts). G.v{slice} describes 
%  each slice (type 'nodes' with xy, d, a, ring, col, edges, ew, ea, or 
%  'matrix' with M). Empty without matrix.
  G = []; N = networkNodes; R = readConn(Pconn);
  if isempty(R) || size(R,1) ~= numel(N.node), return; end
  d0 = 8 / ropts.res; Fz = isfinite(R); Rz = R; Rz(~Fz) = 0;
  bm = @(A,B) (A' * Rz * B) ./ max(1,A' * double(Fz) * B);            % mean z between the parcels of groups (columns)
  % axial: nodes of the networks of each hemisphere
  A  = double(N.node(:) == (1:14)); Z = bm(A,A);
  ax = nodeStruct(N.pos,diag(Z),Z,N.col,d0,true(14,1),zeros(14,1));
  % sagittal: parts of the networks (left and right hemisphere averaged)
  np = max(N.part); L = double(N.part(:) == (1:np) & ~N.right(:)); Rh = double(N.part(:) == (1:np) & N.right(:));
  Zs = (bm(L,L) + bm(Rh,Rh)) / 2;                                        % within the hemispheres
  sa = nodeStruct(N.ppos,diag(Zs),Zs,N.col(N.pnet,:),d0,true(np,1),min(1,max(0,(diag(bm(L,Rh)) - 0.1) / 0.6)));
  % matrix: mean z within (diagonal) and between the networks (hemispheres averaged)
  Ln = double(N.net(:) == (1:7) & ~N.right(:)); Rn = double(N.net(:) == (1:7) & N.right(:));
  M  = (bm(Ln,Ln) + bm(Rn,Rn)) / 2;
  G  = struct('n',n,'col',N.col(1:7,:),'v',{cell(1,size(K,1))});
  c0 = [0; -18; 18];                                                     % center of the slices (see renderSlice)
  for ki = 1:size(K,1)
    an = K(ki,1); ap = setdiff(1:3,an); c = c0; c(an) = K(ki,2);
    if an == 2, G.v{ki} = struct('type','matrix','M',M); continue; end
    if an == 1, V = sa; else, V = ax; end
    V.xy = [(V.xy(:,1) - c(ap(1))) / ropts.res + (n+1)/2, n + 1 - ((V.xy(:,2) - c(ap(2))) / ropts.res + (n+1)/2)];
    G.v{ki} = V;
  end
end
% =========================================================================
function V = nodeStruct(P,zin,Z,col,d0,show,ring)
%nodeStruct. Nodes at the positions P (mm) with the internal connectivity
%  zin, the mean z between the nodes Z, the colors col, the minimum 
%  diameter d0 (pixel), and the ring weight (0 to 1, e.g. of the
%  homotopic connectivity), see networkGeometry.
  w = min(1,max(0,(zin(:) - 0.2) / 0.6)); m = size(P,1);
  [ea,eb] = find(triu(true(m),1)); we = min(1,max(0,(Z(sub2ind([m m],ea,eb)) - 0.1) / 0.6));
  [we,o] = sort(we); k = we > 0.02; o = o(k); we = we(k);
  V = struct('type','nodes','xy',P,'d',d0 * (1 + w),'a',0.5 + 0.5 * w,'ring',ring(:),'col',col,'show',show, ...
    'edges',[ea(o) eb(o)],'ew',d0 * we,'ea',we);
end
% =========================================================================
function drawNetworkVector(ax,G,ki,x0,y0)
%drawNetworkVector. Network G (see networkGeometry) of the slice ki as 
%  vector graphics in the axes ax of a page image (data units in pixel of
%  the page, tile with the upper left corner x0,y0): edges as bands with
%  the colors of their nodes (interpolated), nodes as disks with a ring 
%  (homotopic connectivity), or the matrix of the networks.
  if isempty(G) || isempty(G.v{ki}), return; end
  V = G.v{ki};
  if strcmp(V.type,'matrix'), drawNetworkMatrix(ax,V.M,G.col,G.n,x0,y0); return; end
  xy = V.xy + [x0 y0]; s = V.show;
  for e = 1:size(V.edges,1)
    i = V.edges(e,1); j = V.edges(e,2);
    if ~(s(i) && s(j)), continue; end
    v = xy(j,:) - xy(i,:); L = norm(v); if L < eps, continue; end
    q = [-v(2) v(1)] / L * V.ew(e) / 2;               % half width perpendicular to the edge
    patch('Parent',ax,'XData',[xy(i,1)+q(1) xy(j,1)+q(1) xy(j,1)-q(1) xy(i,1)-q(1)], ...
      'YData',[xy(i,2)+q(2) xy(j,2)+q(2) xy(j,2)-q(2) xy(i,2)-q(2)],'FaceVertexCData',V.col([i j j i],:), ...
      'FaceColor','interp','EdgeColor','none','FaceAlpha',V.ea(e));
  end
  t = linspace(0,2*pi,41);
  for k = find(s(:)')
    r = V.d(k)/2;
    if V.ring(k) > 0.02 % (ring: annulus outside the node)
      rw = r + 0.4 * min(V.d) * V.ring(k);             % (up to 40% of the smallest node diameter)
      patch('Parent',ax,'XData',xy(k,1) + [rw*cos(t) r*cos(fliplr(t))],'YData',xy(k,2) + [rw*sin(t) r*sin(fliplr(t))], ...
        'FaceColor',0.5 * V.col(k,:) + 0.5,'FaceAlpha',0.3 + 0.7 * V.ring(k),'EdgeColor','none');
    end
    patch('Parent',ax,'XData',xy(k,1) + r*cos(t),'YData',xy(k,2) + r*sin(t),'FaceColor',V.col(k,:), ...
      'FaceAlpha',V.a(k),'EdgeColor','none');
  end
end
% =========================================================================
function drawNetworkMatrix(ax,M,col,n,x0,y0)
%drawNetworkMatrix. Matrix M of the mean Fisher z within (diagonal) and
%  between the 7 networks (colors col) in the tile (n pixel, upper left 
%  corner x0,y0) with the names of the networks and the values in the 
%  cells (black for z <= 0, red to yellow for z of 0 to 0.8, blue for 
%  negative z). The font size is given by the cell size on the page (about
%  6 pt for the scan rows).
  nets = {'Vis','SMN','DAN','VAN','Lim','FPN','DMN'}; cs = 0.12 * n; gx = x0 + 0.14 * n; gy = y0 + 0.09 * n;
  u = get(ax,'Units'); set(ax,'Units','centimeters'); p = get(ax,'Position'); set(ax,'Units',u);
  fs = 0.40 * cs * min(p(3) / diff(get(ax,'XLim')),p(4) / diff(get(ax,'YLim'))) * 72/2.54; % (cell in cm of the page)
  for i = 1:7
    for j = 1:7
      z = M(i,j); t = min(1,max(0,z / 0.8)); tn = min(1,max(0,-z / 0.2));
      c = [min(1,2*t) max(0,2*t - 1) * 0.9 0]; if z < 0, c = [0.1 0.3 1] * tn; end
      if ~isfinite(z), c = [0.2 0.2 0.2]; end
      patch('Parent',ax,'XData',gx + (j-1)*cs + [0 cs cs 0],'YData',gy + (i-1)*cs + [0 0 cs cs],'FaceColor',c, ...
        'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.3);
      if isfinite(z)
        tc = [1 1 1]; if t > 0.55, tc = [0 0 0]; end
        text(ax,gx + (j-0.5)*cs,gy + (i-0.5)*cs,regexprep(sprintf('%0.2f',round(z,2) + 0),'^(-?)0\.','$1.'),'Color',tc,'FontSize',fs, ... % (+0: no -0.00)
          'HorizontalAlignment','center','VerticalAlignment','middle');
      end
    end
    lc = 0.65 * col(i,:) + 0.35;
    text(ax,gx - 0.02*n,gy + (i-0.5)*cs,nets{i},'Color',lc,'FontSize',fs,'HorizontalAlignment','right','VerticalAlignment','middle');
    text(ax,gx + (i-0.5)*cs,gy - 0.02*n,nets{i},'Color',lc,'FontSize',fs,'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
end
% =========================================================================
function s = networkSegregation(R,net)
%networkSegregation. System segregation of each of the 7 networks of the
%  connectivity matrix R (Fisher z of the parcels with the network net),
%  i.e. (within - between) / within with the mean z of the edges within 
%  the network and of the edges between the network and all other networks
%  (Chan et al. 2014, as the segregation of the rs-fMRI of all networks).
  s = nan(1,7);
  for ni = 1:7
    in = net(:) == ni; W = R(in,in); W = W(triu(true(nnz(in)),1)); B = R(in,~in);
    w  = mean(W(isfinite(W))); bt = mean(B(isfinite(B))); s(ni) = (w - bt) / w;
  end
end
% =========================================================================
function N = networkNodes
%networkNodes. Schematic network of the rs-fMRI (persistent) of the 100
%  parcels of the Schaefer2018 atlas of CAT (see networks in 
%  cat_io_dcm2bids_pp_rsfmri) with the network net (7 networks of Yeo et 
%  al. 2011) and hemisphere (right) of each parcel: 
%   - the 14 nodes of the networks of both hemispheres (node of each 
%     parcel, nodes 1 to 7 left and 8 to 14 right, colors col) with their
%     axial positions pos (MNI x,y in mm)
%   - the parts of the networks in the sagittal projection (part of each
%     parcel, network pnet of each part, both hemispheres): 2 parts of the
%     visual, somatomotor, dorsal and ventral attention, and frontoparietal
%     network (e.g. anterior and posterior), one of the limbic, and 3 of the
%     default mode network (frontal, posterior, and temporal), by k-means 
%     of the parcel centers (y,z, deterministic initialization), with their
%     sagittal positions ppos (MNI y,z in mm)
%  The positions are schematic (designed once): close to the main spot of
%  the network or part within the cortex (and within the brain of the axial
%  slice z=+10 and of the sagittal slice x=+30 of EPI images), so that the
%  edges of typical rs-fMRI data (pilot data of several sites) do not pass
%  other nodes, i.e. the frontal and the superior parietal regions are used
%  as well (axial nodes symmetric).
  persistent S
  if isempty(S)
    parc = 'Schaefer2018_100Parcels_17Networks_order';
    P = cat_io_dcm2bids_pp_rsfmri('atlasFile',parc); V = spm_vol(P); Y = round(spm_read_vols(V));
    [~,net] = cat_io_dcm2bids_pp_rsfmri('networks',parc); net = net(:); np = numel(net);
    C = regexp(fileread(spm_file(P,'ext','.csv')),'\r?\n','split'); C = C(~cellfun('isempty',strtrim(C)));
    D = cellfun(@(l) strsplit(l,';'),C(2:end),'UniformOutput',false); D = vertcat(D{:});
    [~,so] = sort(str2double(D(:,1))); D = D(so,:);
    right = contains(D(:,3),'_RH_'); node = net + 7 * right; node(net == 0) = 0;
    X = nan(np,3);
    for p = 1:np
      [i,j,k] = ind2sub(size(Y),find(Y == p));
      if ~isempty(i), x = V.mat * [mean(i); mean(j); mean(k); 1]; X(p,:) = x(1:3)'; end
    end
    % parts of the networks in the sagittal projection
    kn = [2 2 2 2 1 2 3]; part = zeros(np,1); pnet = zeros(0,1);
    for ni = 1:7
      p = find(net == ni & isfinite(X(:,1))); id = kmeansDet(X(p,2:3),kn(ni));
      for c = 1:kn(ni), pnet(end+1,1) = ni; part(p(id == c)) = numel(pnet); end %#ok<AGROW>
    end
    % schematic positions (MNI, mm): axial nodes of the left hemisphere (x,y,
    % right mirrored) and sagittal parts (y,z) in the order of the parts, i.e.
    % superior and inferior visual, superior and lateral somatomotor, posterior
    % and anterior dorsal attention, anterior and posterior ventral attention,
    % limbic, anterior and posterior frontoparietal, and posterior, anterior,
    % and temporal default mode network
    pos  = [-26 -94; -62 -23; -48 -75; -55 12; -12 24; -38 51; -16 66];   % Vis SMN DAN VAN Lim FPN DMN
    ppos = [-86 11; -54 -21; -10 68; 1 -4; -82 36; -40 65; 22 21; 11 54; 12 -37; 51 32; -67 56; -55 19; 47 -2; -15 -36];
    col  = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] / 255;
    S = struct('net',net,'right',right,'node',node,'col',[col; col],'pos',[pos; -pos(:,1) pos(:,2)], ...
      'part',part,'pnet',pnet,'ppos',ppos);
  end
  N = S;
end
% =========================================================================
function id = kmeansDet(Z,k)
%kmeansDet. k-means clusters id of the points Z (rows) with a deterministic
%  initialization (the point farthest from the mean, then the points 
%  farthest from the chosen centers).
  id = ones(size(Z,1),1); if k < 2, return; end
  dist = @(A,B) sqrt(max(0,sum(A.^2,2) + sum(B.^2,2)' - 2 * A * B'));
  [~,m] = max(dist(Z,mean(Z,1))); C = Z(m,:);
  for c = 2:k, [~,m] = max(min(dist(Z,C),[],2)); C(c,:) = Z(m,:); end
  for it = 1:100
    [~,id] = min(dist(Z,C),[],2); C0 = C;
    for c = 1:k, if any(id == c), C(c,:) = mean(Z(id == c,:),1); end, end
    if isequal(C,C0), break; end
  end
end
% =========================================================================
function v = fieldValue(S,fn,default)
%fieldValue. Field fn of the structure S (default if not available or empty).
  if isstruct(S) && isfield(S,fn) && ~isempty(S.(fn)), v = S.(fn); else, v = default; end
end
% =========================================================================
function R = readConn(P)
%readConn. Connectivity matrix of the table P (parcel names in the first
%  row and column, see cat_io_dcm2bids_pp_rsfmri).
  R = [];
  try
    C = cat_io_csv(sourceFile(P),'','',struct('delimiter','\t','convert2double',1)); C = C(2:end,2:end);
    R = nan(size(C));
    for i = 1:numel(C), if isnumeric(C{i}) && ~isempty(C{i}), R(i) = double(C{i}(1)); end, end
  end
end
% =========================================================================
function clim = renderWindow(F)
%renderWindow. Fixed intensity window of the DTI maps (comparable scans): 
%  FA 0-1, DEC 0-0.8 (FA weighted directions), MD/AD/RD 0-0.003 mm2/s, 
%  otherwise empty (percentiles of the scan).
  switch F.prefix
    case 'FA',              clim = [0 1];
    case 'DEC',             clim = [0 0.8];
    case {'MD','AD','RD'},  clim = [0 0.003];
    otherwise,              clim = [];
  end
end
% =========================================================================
function renderLegend(ax,x,y,fs,ovc)
%renderLegend. Legend of the overlay lines (right aligned at x,y).
  str = sprintf(['original (scanner) planes:  \\color[rgb]{%g %g %g}x=0  ' ...
    '\\color[rgb]{%g %g %g}y=0  \\color[rgb]{%g %g %g}z=0'], ovc'); 
  text(ax,x,y,str,'Color','w','FontSize',fs,'Interpreter','tex', ...
    'HorizontalAlignment','right','VerticalAlignment','middle');
end
% =========================================================================
function L = overlayLines(T,K,ropts,n)
%overlayLines. Lines where the slices K (n x n pixel, see renderSlice)
%  cross the coordinate planes x=0, y=0, and z=0 of the original (scanner)
%  space with the transformation T from world to MNI space, as rows
%  [plane x1 y1 x2 y2] in the pixel coordinates of the displayed tile 
%  (column and row) for each slice (cell, planes that are nearly parallel
%  to the slice are ignored), e.g. for vector graphics on the pages.
  L = cell(1,size(K,1)); if isempty(T), return; end
  c0 = [0; -18; 18];
  for ki = 1:size(K,1)
    an = K(ki,1); ap = setdiff(1:3,an); c = c0; c(an) = K(ki,2);
    Mw = zeros(4); Mw(4,4) = 1;
    Mw(ap(1),1) = ropts.res; Mw(ap(1),4) = c(ap(1)) - (n+1)/2*ropts.res;
    Mw(ap(2),2) = ropts.res; Mw(ap(2),4) = c(ap(2)) - (n+1)/2*ropts.res;
    Mw(an,3)    = 1;         Mw(an,4)    = c(an) - 1;
    A  = T \ Mw; L{ki} = zeros(0,5);
    for k = 1:3 % plane a*i + b*j + d = 0 of the slice voxel (i,j), displayed at x = i, y = n + 1 - j (rot90)
      a = A(k,1); b = A(k,2); d = A(k,3) + A(k,4);
      if norm([a b]) < 0.3*ropts.res, continue; end
      e = [0.5 n+0.5]; P = zeros(0,2);
      if abs(b) > eps, for x = e, j = -(a*x + d)/b; y = n + 1 - j; if y >= e(1) && y <= e(2), P(end+1,:) = [x y]; end, end, end %#ok<AGROW>
      if abs(a) > eps, for y = e, j = n + 1 - y; x = -(b*j + d)/a; if x >= e(1) && x <= e(2), P(end+1,:) = [x y]; end, end, end %#ok<AGROW>
      if size(P,1) >= 2
        [~,i1] = min(P(:,1) + P(:,2)); [~,i2] = max(P(:,1) + P(:,2)); if i1 == i2, i2 = find((1:size(P,1))' ~= i1,1); end
        L{ki}(end+1,:) = [k P(i1,:) P(i2,:)];
      end
    end
  end
end
% =========================================================================
function drawOverlayLines(ax,L,ki,x0,y0,ovc)
%drawOverlayLines. Lines L (see overlayLines) of the slice ki as vector
%  graphics in the axes ax of a page image (tile with the upper left corner
%  x0,y0 in page pixels) in the colors ovc of the planes.
  if isempty(L) || isempty(L{ki}), return; end
  for li = 1:size(L{ki},1)
    l = L{ki}(li,:);
    plot(ax,x0 + l([2 4]),y0 + l([3 5]),'-','Color',ovc(l(1),:),'LineWidth',0.5);
  end
end
% =========================================================================
function [Ys,Yov] = renderSlice(C,V,T,an,pos,ropts)
%renderSlice. Slice of the (first) volume V sampled in MNI space.
%  C are the B-spline coefficients of V (spm_bsplinc) of degree ropts.interp
%  (>1), otherwise (empty C) spm_slice_vol samples V directly.
%  T is the transformation from world (scanner) to MNI space (identity for
%  spatially normalized images, i.e., world coordinates in MNI space). The slice
%  normal is the MNI axis an (1-x sagittal, 2-y coronal, 3-z axial) at the
%  position pos (mm), the in-plane axes are the other two axes centered
%  at the typical brain bounding box. Yov labels the lines where the slice
%  crosses the coordinate planes x=0 (1), y=0 (2), and z=0 (3) of the
%  original (scanner) space.
  n  = round(ropts.fov / ropts.res);
  ap = setdiff(1:3,an);              % in-plane axes (horizontal, vertical)
  c0 = [0; -18; 18];                 % center of the typical brain bounding box

  % slice voxel (i,j,1) to MNI (mm) and then to image voxel coordinates
  c  = c0; c(an) = pos;
  Mw = zeros(4); Mw(4,4) = 1;
  Mw(ap(1),1) = ropts.res; Mw(ap(1),4) = c(ap(1)) - (n+1)/2*ropts.res;
  Mw(ap(2),2) = ropts.res; Mw(ap(2),4) = c(ap(2)) - (n+1)/2*ropts.res;
  Mw(an,3)    = 1;         Mw(an,4)    = c(an) - 1;
  M  = (T * V.mat) \ Mw;
  [ii,jj] = ndgrid(1:n,1:n);
  if isempty(C)
    Ys = spm_slice_vol(V, M, [n n], [ropts.interp NaN]);
    Ys(isnan(Ys)) = 0;
  else
    x  = M(1,1)*ii + M(1,2)*jj + M(1,3) + M(1,4);
    y  = M(2,1)*ii + M(2,2)*jj + M(2,3) + M(2,4);
    z  = M(3,1)*ii + M(3,2)*jj + M(3,3) + M(3,4);
    Ys = spm_bsplins(C,x,y,z,[ropts.interp*[1 1 1] 0 0 0]);
    Ys( ~isfinite(Ys) | x<0.99 | y<0.99 | z<0.99 | x>V.dim(1)+0.01 | y>V.dim(2)+0.01 | z>V.dim(3)+0.01 ) = 0;
  end

  % original coordinate axes: pixels where the original world coordinate k
  % is (about) zero, i.e., lines of about one pixel but at least 1 mm width 
  % (planes that are nearly parallel to the slice are ignored)
  A   = T \ Mw;                      % slice voxel to world (scanner) coordinates
  Yov = zeros(n,'uint8');
  for k = 1:3
    g = A(k,1:2);
    if norm(g) < 0.3*ropts.res, continue; end
    Yov( abs(A(k,1)*ii + A(k,2)*jj + A(k,3) + A(k,4)) <= max(abs(g))/2 * max(1,1/ropts.res) ) = k;
  end

  % first in-plane axis left to right and second upwards, i.e., axial with
  % anterior on top (neurological view), coronal/sagittal with superior on
  % top (sagittal with anterior on the right)
  Ys  = rot90(Ys);
  Yov = rot90(Yov);
end
