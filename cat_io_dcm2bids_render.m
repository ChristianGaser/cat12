function cat_io_dcm2bids_render(PBIDS,Preport,ropts)
%cat_io_dcm2bids_render. Render slices of each scan of a BIDS directory in MNI space.
%
%  Aim: The rendering supports the visual identification (or confirmation,
%  e.g. from the report table) of outliers in large sets of similar images.
%  Therefore, the scans are rendered separately for each image type
%  (dataset/datatype/subtype, e.g. anat/t1w) and are comparable, as all
%  images are sampled at the same positions in MNI space using the
%  registration of the session (catDCM2BIDSaffine_sub-*_ses-*.json in the
%  derivatives directory of the session, see writeSessionAffines in
%  cat_io_dcm2bids_pp), i.e., independent of the matrix size, slice
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
%         [dataset/]datatype/subtype/render_<subtype>_<orientation>_<axis><pos>_p##.png
%       This gives a very brief overview of many images but is biased by
%       the selected slice, i.e., it is useful for a rough global review.
%   R2) Scan-row pages with one row per scan with a header (BIDS filename
%       and an underline over the whole page), an information panel
%       (quality ratings, image and acquisition parameters, registration),
%       and the slices axial z=+10, sagittal x=0, and coronal y=0:
%         [dataset/]datatype/subtype/render_<subtype>_rows_p##.png
%       with 6 scans per page.
%   R3) Subject reports with the scan rows (as R2) of all scans (image
%       types) of a subject ordered by session, datatype (anat, dwi, func,
%       fmap, perf, others), and name (i.e., derivatives follow their scan)
%       as one multi-page PDF per subject:
%         [dataset/]subjects/render_<sub>.pdf
%   rows/  .. slice row of each scan without text as separate image
%             (e.g. for a later interactive HTML report)
%   cache/ .. sampled slices of each scan (see renderScanSlices)
%  The A4 portrait format of SPM is used. It allows to compare multiple
%  pages on one screen and limits the number of output pages.
%  4D data is represented by its first volume. Further session reports
%  with a compressed representation of 4D data (e.g. FA maps for DWI)
%  should be prepared in the QC.
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
%                 ([3 10; 2 0; 1 0]), empty for no R1 (also as structure
%                 array with the fields orient and pos)
%   .tiles     .. R1 tiles per page: [3 4], [4 5], [5 7], or [6 8] ([4 5])
%   .rows      .. R2 scan-row pages (1)
%   .subjects  .. R3 subject reports (1)
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
%
%  The intensities of a scan are scaled for all its slices together
%  (1-99% percentile) to keep the slices of a scan comparable. The colored
%  lines show the coordinate planes of the original (scanner) space, i.e.,
%  x=0 (red), y=0 (green), and z=0 (blue). Labels are colored by the
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

  % all slices (K) of R1 (k1) and of the scan rows of R2/R3 (k2) to sample
  % each scan only once
  rowslices = [3 10; 1 0; 2 0];
  userows   = ropts.rows || ropts.subjects || ropts.rowimg;
  K = unique([ropts.slices; rowslices(1:3*userows,:)],'rows','stable');
  if isempty(K), return; end
  [~,k1] = ismember(ropts.slices,K,'rows');
  [~,k2] = ismember(rowslices,K,'rows');

  % common render settings
  G.pw  = 21; G.ph = 29.7;                                % A4 portrait (as SPM) in cm
  G.ovc = [1 0.25 0.25; 0.25 1 0.25; 0.35 0.55 1];        % overlay colors of the x/y/z-planes
  G.orientnam = {'sagittal','coronal','axial'}; G.axisnam = 'xyz';
  MarkColor   = cat_io_colormaps('marks+',40);            % quality colors (as command line)
  G.col2mark  = @(val) min(1,max(0,MarkColor(min(size(MarkColor,1)-3, ...
    max(1,floor( val/9.5 * size(MarkColor,1)))),:)));


  %% find and describe images
  P = cat_vol_findfiles(PBIDS,'*sub-*.nii*');
  P = P( ~cellfun('isempty',regexp(P,'\.nii(\.gz)?$','once')) );
  isderiv = cat_io_contains(P,[filesep 'derivatives' filesep]);
  switch ropts.source
    case 1, P = P(~isderiv); isderiv = isderiv(~isderiv);
    case 2, P = P( isderiv); isderiv = isderiv( isderiv);
  end
  if isempty(P), return; end

  F = struct('file',P,'name','','bidsname','','prefix','','sub','','ses','','dataset','', ...
    'datatype','','subtype','','group','','json','','SQR',nan,'QR',[],'T',[],'reg',[], ...
    'readerr',0,'dim',[],'vx',[],'nvol',nan);
  for fi = 1:numel(P)
    [pp,ff]  = fileparts(P{fi});
    ff       = regexprep(ff,'\.nii$','');
    si       = strfind(ff,'sub-');
    prefix   = regexprep(ff(1:si(1)-1),'_$','');
    bidsname = ff(si(1):end);
    F(fi).name = ff; F(fi).bidsname = bidsname; F(fi).prefix = prefix;

    % BIDS entities
    tok = regexp(bidsname,'sub-([^_]+)','tokens','once'); F(fi).sub = ['sub-' tok{1}];
    tok = regexp(bidsname,'ses-([^_]+)','tokens','once');
    if ~isempty(tok), F(fi).ses = ['ses-' tok{1}]; end
    parts = strsplit(bidsname,'_');
    F(fi).datatype = spm_file(pp,'basename');
    F(fi).subtype  = lower(parts{end});
    if ~isempty(prefix), F(fi).subtype = [F(fi).subtype '_' prefix]; end

    % dataset (protocol subdirectory) relative to the main BIDS directory
    si     = strfind(P{fi},[filesep 'sub-']);
    root   = P{fi}(1:si(1)-1);
    reldir = pp(numel(root)+2:end);
    if isderiv(fi), droot = root; rroot = fileparts(root); else, droot = fullfile(root,'derivatives'); rroot = root; end
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

    % quality ratings
    Pqc = fullfile(droot,reldir,['catDCM2BIDSqc_' bidsname '.json']);
    if exist(Pqc,'file')
      try
        QC = cat_io_json(Pqc);
        F(fi).QR = QC.qualityrating;
        SQR = QC.qualityrating.SQR;
        if isnumeric(SQR) && isscalar(SQR) && isfinite(SQR), F(fi).SQR = double(SQR); end % e.g. empty if QC failed
      end
    end
  end


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

    % prepare output directory (remove old pages)
    Pout = fullfile(Preport,groups{gi});
    if ~exist(Pout,'dir'), mkdir(Pout); end
    Pold = cat_vol_findfiles(Pout,'render_*.png',struct('depth',0));
    for oi = 1:numel(Pold), delete(Pold{oi}); end

    % sample all slices of each scan once (or load them from the cache)
    Y = zeros(n,n,size(K,1),numel(Fg),'uint8'); Yov = Y; nnew = 0;
    for fi = 1:numel(Fg)
      S = renderScanSlices(Fg(fi),K,ropts,Pcache(Fg(fi)));
      Y(:,:,:,fi) = S.Y; Yov(:,:,:,fi) = S.Yov; nnew = nnew + S.new;
      Fg(fi).readerr = S.readerr; Fg(fi).dim = S.dim; Fg(fi).vx = S.vx; Fg(fi).nvol = S.nvol;

      % slice row without text (e.g. for a later HTML report)
      if ropts.rowimg
        Prow = fullfile(Pout,'rows',[Fg(fi).name '_row.png']);
        if S.new || ~exist(Prow,'file')
          if ~exist(fileparts(Prow),'dir'), mkdir(fileparts(Prow)); end
          Yr = cell(1,3);
          for ci = 1:3, Yr{ci} = renderTileRGB(S.Y(:,:,k2(ci)),S.Yov(:,:,k2(ci)),G.ovc); end
          imwrite(cat(2,Yr{:}),Prow);
        end
      end
    end

    np1 = 0; np2 = 0;
    if ~isempty(k1)
      np1 = renderOverviewPages(fh,Fg,Y(:,:,k1,:),Yov(:,:,k1,:),K(k1,:),groups{gi},Pout,ropts,G);
    end
    if ropts.rows
      np2 = renderRowPages(fh,Fg,Y(:,:,k2,:),Yov(:,:,k2,:),K(k2,:),groups{gi}, ...
        fullfile(Pout,['render_' Fg(1).subtype '_rows_p%02d.png']),0,ropts,G);
    end
    cat_io_cprintf('blue',sprintf('  Rendered %4d images (%4d new) in %3d overview and %3d row page(s): %s\n', ...
      numel(Fg), nnew, np1, np2, Pout));
  end


  %% subject reports (all scans of a subject, slices from the cache)
  if ropts.subjects
    [subs,~,sid] = unique(strcat({F.dataset},'|',{F.sub}));
    for si = 1:numel(subs)
      Fs = renderSubjectOrder(F(sid==si),ropts);
      Y  = zeros(n,n,3,numel(Fs),'uint8'); Yov = Y;
      for fi = 1:numel(Fs)
        S = renderScanSlices(Fs(fi),K(k2,:),ropts,Pcache(Fs(fi)));
        Y(:,:,:,fi) = S.Y; Yov(:,:,:,fi) = S.Yov;
        Fs(fi).readerr = S.readerr; Fs(fi).dim = S.dim; Fs(fi).vx = S.vx; Fs(fi).nvol = S.nvol;
      end
      Pout = fullfile(Preport,Fs(1).dataset,'subjects');
      if ~exist(Pout,'dir'), mkdir(Pout); end
      Ppdf = fullfile(Pout,['render_' Fs(1).sub '.pdf']);
      if exist(Ppdf,'file'), delete(Ppdf); end
      np = renderRowPages(fh,Fs,Y,Yov,K(k2,:),fullfile(Fs(1).dataset,Fs(1).sub),Ppdf,1,ropts,G);
      cat_io_cprintf('blue',sprintf('  Rendered %4d images in %3d subject page(s): %s\n', ...
        numel(Fs), np, Ppdf));
    end
  end
end
% =========================================================================
function np = renderOverviewPages(fh,Fg,Y,Yov,K,group,Pout,ropts,G)
%renderOverviewPages. R1 pages with one slice per scan in m x n tiles.
%  Y and Yov are the slices (n x n x slices x scans) of the slices K. Each
%  tile shows the subject, session, and quality rating SQR in its header.
  nx = ropts.tiles(1); ny = ropts.tiles(2); nt = nx*ny;
  fs = [8 7 6 5]; fs = fs( all( [3 4; 4 5; 5 7; 6 8] == ropts.tiles(:)' , 2 ) ); % font size
  mg = 0.002; ht = 0.02; % margin and title height (normalized)
  tw = (1 - 2*mg) / nx * G.pw; th = (1 - 2*mg - ht) / ny * G.ph; % tile size in cm
  n  = size(Y,1);
  % free header above each slice for about two text lines to reduce the
  % overlay of the labels (sub, ses, SQR) with the image
  hd = round(2 * 1.2 * fs/72*2.54 * n/min(tw,th));            % header in pixel
  cw = round(max(n, (n+hd)*tw/th)); ch = round(max(n+hd, n*th/tw)); % tile size in pixel
  np = ceil(numel(Fg)/nt);

  % labels: subject, session and quality rating in separate rows
  lab = cell(1,numel(Fg));
  for fi = 1:numel(Fg)
    lab{fi} = {Fg(fi).sub};
    if isempty(Fg(fi).T), lab{fi}{1} = [lab{fi}{1} ' (no MNI)']; end
    if Fg(fi).readerr,    lab{fi}{1} = [lab{fi}{1} ' (read error)']; end
    if ~isempty(Fg(fi).ses), lab{fi}{end+1} = Fg(fi).ses; end
    if ~isnan(Fg(fi).SQR),   lab{fi}{end+1} = sprintf('SQR: %0.1f',Fg(fi).SQR); end
  end

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
          renderTileRGB(Y(:,:,ki,fi),Yov(:,:,ki,fi),G.ovc);
      end
      ax = axes('Parent',fh,'Position',[mg mg 1-2*mg 1-2*mg-ht],'Color','k');
      image(ax,page); axis(ax,'image','off');

      % labels (in page coordinates) colored by the quality rating (white if not available)
      for ti = 1:ntp
        fi = (pgi-1)*nt + ti;
        if isnan(Fg(fi).SQR), col = [1 1 1]; else, col = G.col2mark(Fg(fi).SQR); end
        text(ax,pos(ti,1) + cw*0.02,pos(ti,2) + 1,lab{fi},'Color',col,'FontSize',fs, ...
          'Interpreter','none','VerticalAlignment','top');
      end

      % print page and give the graphics system time to finish
      Ppng = fullfile(Pout,sprintf('render_%s_%s_%s%+0.0f_p%02d.png', ...
        Fg(1).subtype,G.orientnam{K(ki,1)},G.axisnam(K(ki,1)),K(ki,2),pgi));
      print(fh, Ppng, '-dpng', sprintf('-r%d',ropts.dpi));
      drawnow;
    end
  end
  np = np * size(K,1);
end
% =========================================================================
function np = renderRowPages(fh,Fg,Y,Yov,K,ttl,Ppage,showgroup,ropts,G)
%renderRowPages. Pages with one row of 3 slices per scan (R2 and R3).
%  Each scan starts with a header with its BIDS filename (and its image
%  type for showgroup=1) and an underline over the whole page, followed by
%  an information panel (see renderInfoText) and the slices K. The tile
%  size is given by the page height with 6 scans per page, whereas the
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
  slab = arrayfun(@(ki) sprintf('%s=%+0.0f',G.axisnam(K(ki,1)),K(ki,2)),1:size(K,1),'UniformOutput',false);
  ispdf = strcmpi(spm_file(Ppage,'ext'),'pdf');

  for pgi = 1:np
    clf(fh);
    nsp  = min(nps, numel(Fg) - (pgi-1)*nps);
    page = zeros(H,W,3,'single');
    for si = 1:nsp
      fi = (pgi-1)*nps + si;
      y0 = pht + (si-1)*bh;                        % top of the scan block
      page(y0 + phh - (2:-1:1), :, :) = 0.6;       % underline over the whole page
      for ci = 1:3
        page(y0 + phh + (1:n), x0 + (ci-1)*n + (1:n), :) = renderTileRGB(Y(:,:,ci,fi),Yov(:,:,ci,fi),G.ovc);
      end
    end
    ax = axes('Parent',fh,'Position',[mg/G.pw mg/G.ph 1-2*mg/G.pw 1-2*mg/G.ph],'Color','k');
    image(ax,page); axis(ax,'image','off');

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
      for ci = 1:3
        text(ax,x0 + (ci-1)*n + 3,y0 + phh + n - 2,slab{ci},'Color',[.6 .6 .6], ...
          'FontSize',fs-1,'Interpreter','none','VerticalAlignment','bottom');
      end
      text(ax,xi + 3,y0 + phh + 3,renderInfoText(Fg(fi),G.col2mark),'Color','w', ...
        'FontSize',fs,'Interpreter','tex','VerticalAlignment','top');
    end

    if ispdf
      exportgraphics(fh, Ppage, 'Append', true, 'ContentType', 'vector', ...
        'Resolution', ropts.dpi, 'BackgroundColor', 'k');
    else
      print(fh, sprintf(Ppage,pgi), '-dpng', sprintf('-r%d',ropts.dpi));
    end
    drawnow;
  end
end
% =========================================================================
function txt = renderInfoText(F,col2mark)
%renderInfoText. TeX lines of the information panel of a scan in R2/R3:
%  subject, session, quality ratings, image size, acquisition parameters
%  (from the sidecar of the raw scan), and MNI registration (with the session
%  orientation rating SOR).
  tex    = @(s) regexprep(s,'([_\\{}^])','\\$1');
  colstr = @(v) sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}',col2mark(v));
  wht    = '\color[rgb]{1 1 1}';
  red    = '\color[rgb]{1 0.3 0.3}';
  if isnan(F.SQR), c = wht; else, c = colstr(F.SQR); end

  % subject and quality
  txt = {[c tex(F.sub)]};
  if ~isempty(F.ses), txt{end+1} = [c tex(F.ses)]; end
  if ~isnan(F.SQR),   txt{end+1} = sprintf('%sSQR: %0.1f',c,F.SQR); end
  if isstruct(F.QR) % further ratings
    fn = setdiff(fieldnames(F.QR),{'SQR','vx_vol'},'stable'); str = '';
    for i = 1:numel(fn)
      v = F.QR.(fn{i});
      if isnumeric(v) && isscalar(v) && isfinite(v), str = [str sprintf('%s%s %0.1f  ',colstr(v),tex(fn{i}),v)]; end %#ok<AGROW>
    end
    if ~isempty(str), txt{end+1} = str; end
  end
  if isempty(F.T), txt{end+1} = [red 'no MNI registration']; end
  if F.readerr,    txt{end+1} = [red 'read error']; end

  % image
  if ~isempty(F.dim)
    txt{end+1} = sprintf('%smatrix: %dx%dx%d, %d volume(s)',wht,F.dim(1:3),F.nvol);
    txt{end+1} = sprintf('%svoxel: %0.2fx%0.2fx%0.2f mm',wht,F.vx);
  end

  % acquisition
  J = struct();
  if exist(F.json,'file')
    try, J = cat_io_json(F.json); end
  end
  str = '';
  if isfield(J,'Manufacturer') && ischar(J.Manufacturer), str = [str J.Manufacturer ' ']; end
  if isfield(J,'MagneticFieldStrength') && isnumeric(J.MagneticFieldStrength), str = [str sprintf('%gT ',J.MagneticFieldStrength)]; end
  if isfield(J,'MRAcquisitionType') && ischar(J.MRAcquisitionType), str = [str J.MRAcquisitionType]; end
  if ~isempty(str), txt{end+1} = [wht tex(strtrim(str))]; end
  tim = {'RepetitionTime','TR'; 'EchoTime','TE'; 'InversionTime','TI'}; str = '';
  for i = 1:size(tim,1)
    if isfield(J,tim{i,1}) && isnumeric(J.(tim{i,1})) && ~isempty(J.(tim{i,1}))
      str = [str sprintf('%s %0.4g ms  ',tim{i,2},1000*J.(tim{i,1})(1))]; %#ok<AGROW>
    end
  end
  if isfield(J,'FlipAngle') && isnumeric(J.FlipAngle) && ~isempty(J.FlipAngle), str = [str sprintf('FA %g%s',J.FlipAngle(1),char(176))]; end
  if ~isempty(str), txt{end+1} = [wht str]; end

  % registration (rotation of the scanner space to MNI)
  if ~isempty(F.reg) && strcmp(F.reg.method,'template')
    txt{end+1} = sprintf('%sMNI: template space (spatially normalized)',wht);
  elseif ~isempty(F.reg)
    rt = round(F.reg.imat(4:6)*1800/pi)/10; rt(rt==0) = 0; % avoid "-0"
    txt{end+1} = sprintf('%sMNI: %s, rotation %+0.1f %+0.1f %+0.1f%s',wht,tex(F.reg.method),rt,char(176));
    if isfinite(F.reg.SOR), txt{end} = sprintf('%s  %sSOR %0.1f',txt{end},colstr(F.reg.SOR),F.reg.SOR); end
  end
end
% =========================================================================
function Fs = renderSubjectOrder(Fs,ropts)
%renderSubjectOrder. Order of the scans of one subject in the subject
%  report: by session, datatype (anat, dwi, func, fmap, perf, others), and
%  BIDS name, where derivatives follow their scan. Only the first session
%  for ropts.sessions.
  dto = {'anat','dwi','func','fmap','perf'};
  [~,dti] = ismember({Fs.datatype},dto); dti(dti==0) = numel(dto) + 1;
  key = strcat({Fs.ses},'|',cellstr(num2str(dti(:)))','|',{Fs.datatype},'|',{Fs.bidsname},'|',{Fs.prefix});
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
%  is saved in Pcache and reused if the source file, the registration, and
%  the sampling parameters were not changed and the cache includes all
%  requested slices. This avoids reading all images again (e.g. for large
%  datasets) and keeps the slices for further outputs (e.g. HTML).
%  S.Y (uint8 slices), S.Yov (overlay labels), S.K (slices), S.readerr,
%  header information (S.dim, S.vx, S.nvol), and S.new (newly sampled).
  n   = round(ropts.fov / ropts.res);
  D   = dir(F.file);
  if isempty(D), srcdate = nan; else, srcdate = D(1).datenum; end
  par = struct('version',2,'srcdate',srcdate,'T',F.T,'fov',ropts.fov,'res',ropts.res,'interp',ropts.interp);

  % cached slices
  if ropts.cache && exist(Pcache,'file')
    try
      C = load(Pcache,'S','par');
      [isK,ki] = ismember(K,C.S.K,'rows');
      if isequal(C.par,par) && all(isK)
        S = C.S; S.K = K; S.Y = S.Y(:,:,ki); S.Yov = S.Yov(:,:,ki); S.new = 0;
        return
      end
    end
  end

  S = struct('K',K,'Y',zeros(n,n,size(K,1),'uint8'),'Yov',zeros(n,n,size(K,1),'uint8'), ...
    'readerr',0,'dim',[],'vx',[],'nvol',nan,'new',1);
  try
    txt = evalc('V = spm_vol(F.file);'); %#ok<NASGU> % avoid gz-messages
    S.nvol = numel(V); V = V(1);
    S.dim  = V.dim; S.vx = sqrt(sum(V.mat(1:3,1:3).^2));
  catch
    S.readerr = 1; % e.g. corrupted files
  end
  if ~S.readerr && ~isempty(F.T)
    try
      % B-spline coefficients of the volume (non-finite values set to 0) 
      % for higher interpolation degrees, otherwise direct sampling
      if ropts.interp > 1
        Yvol = spm_read_vols(V); Yvol(~isfinite(Yvol)) = 0;
        C    = spm_bsplinc(Yvol,[ropts.interp*[1 1 1] 0 0 0]); clear Yvol
      else
        C    = [];
      end
      Ys   = zeros(n,n,size(K,1),'single');
      for ki = 1:size(K,1)
        [Ys(:,:,ki),S.Yov(:,:,ki)] = renderSlice(C,V,F.T,K(ki,1),K(ki,2),ropts);
      end
      clear C
      % common intensity scaling of all slices of the scan
      Yv   = Ys(Ys~=0);
      clim = [prctile(Yv,1) prctile(Yv,99)];
      if ~(numel(clim)==2 && all(isfinite(clim)) && diff(clim)>0), clim = [0 max(1,max(Ys(:)))]; end
      S.Y  = uint8(255 * min(1,max(0,(Ys - clim(1)) / diff(clim))));
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
function Yrgb = renderTileRGB(Y,Yov,ovc)
%renderTileRGB. Gray-scale slice (uint8) with colored overlay lines.
  Yrgb = repmat(single(Y)/255,1,1,3);
  for ci = 1:3
    Yc = Yrgb(:,:,ci);
    for k = 1:3, Yc(Yov==k) = ovc(k,ci); end
    Yrgb(:,:,ci) = Yc;
  end
end
% =========================================================================
function renderLegend(ax,x,y,fs,ovc)
%renderLegend. Legend of the overlay lines (right aligned at x,y).
  text(ax,x,y,sprintf(['original (scanner) planes:  \\color[rgb]{%g %g %g}x=0  ' ...
    '\\color[rgb]{%g %g %g}y=0  \\color[rgb]{%g %g %g}z=0'], ovc'), ...
    'Color','w','FontSize',fs,'Interpreter','tex', ...
    'HorizontalAlignment','right','VerticalAlignment','middle');
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
