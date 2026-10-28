function varargout = cat_io_dcm2bids_reportfig(action,varargin)
%cat_io_dcm2bids_reportfig. Report figure of cat_io_dcm2bids.
%  Diagrams of the report tables (see cat_io_dcm2bids_report) as A4 pages
%  for the project manager (data per site and time), the principal
%  investigator (sample and brain markers), and the scanner operator
%  (ratings, outliers, and systematic scan issues): an overview page of the
%  BIDS directory and one page per site (for 2 to opts.reportfig.maxsites 
%  sites, i.e. a large number of sites, e.g. of a traveling-subject study, 
%  is only shown together, without site columns, site labels, and site 
%  pages) with the header and the panels (rows)
%    header        logos of the design, project, BIDS directory (site and
%                  scanner), numbers of subjects, sessions, and scans, and
%                  the median ratings of the sessions with icons (protocol
%                  fit rating PFR, protocol conformity percentage PCP,
%                  brain coverage rating BBR, session orientation rating
%                  SOR, phase-encoding reversal rating PRR, scan quality
%                  rating SQR)
%    data          sessions per month (by site, from the date of the session
%                  label, i.e. not for coded session labels of anonymize > 1)
%                  and the site table (sites as columns, numbers as blue
%                  boxes, ratings in the colors of the marks)
%    sample        age/sex table with pictograms (numbers as blue boxes, the
%                  mean age, its SD, and the age range as yellow-orange
%                  boxes for older subjects and more variable ages) and the
%                  sex pie, age pyramid (age at the first session of each
%                  subject, males left, females right, other or unknown sex
%                  half on each side in the center, see opts.reportfig.
%                  agerange and agestep, by default the age range of the
%                  sample, see ageRange), and the pyramid of the subjects
%                  per site with the mean age and its SD
%    markers       relative tissue volumes (best T1w scan), white matter FA
%                  (best DWI scan), and system segregation (best BOLD scan)
%                  by age with the robust age/sex model of the sample (see
%                  cat_io_dcm2bids_stats, points by site, with grid)
%    homogeneity   sample homogeneity (quartic mean Z-score, corrected for
%                  age and sex for at least 10 sessions) of the mwc1 maps,
%                  the normalized FA maps, and the connectivity matrices by
%                  the quality rating of their scans (points by site)
%    QC            session ratings (in the order of the header, points by
%                  site; for one site the ratings of all scans with the
%                  protocol status of each scan as PCP, points by datatype)
%                  and image quality SQR of the best scan of each scan type
%                  as box plots (marks 0.5 to 6.5 with their colors in the
%                  background, worse ratings are clipped as in the CAT
%                  report), and the completeness and protocol status of the
%                  sessions
%  The values are given in or next to the diagrams. The pages are saved for
%  each design of opts.reportfig.design (one or several, by default all) as
%    <outdir>/<subdir>/<BIDSdir>-report/report_<BIDSsubdir>_<design>.pdf
%  and the overview of the first design is shown in the SPM Graphics window
%  in interactive sessions. The statistics are columns of the session table
%  (see cat_io_dcm2bids_stats).
%  A design (DCM2BIDS/colormaps/<design>.json, e.g. CAT_bright, CAT_dark,
%  or DZPG_bright) defines the logos (one or two image files of the CAT
%  directory, top left and top right of the header), the icon of
%  the ratings, the colors of the text and background (foreground full and
%  light, background light and full), and the colormaps of discrete values
%  (sites without color of the site dictionary, see the color column of the
%  project table) and of continuous values (name of cat_io_colormaps or of
%  MATLAB, or RGB matrix), whereas the ratings use the colors of the marks of the
%  QC (cat_io_colormaps('marks+'), as the command line output). The font
%  size opts.reportfig.fontsize (8 pt) is scaled by opts.reportfig.fontscale
%  for header-like elements (0.96, titles of the diagrams) and for axes,
%  table entries, and legends (0.8), whereas the header of the page has
%  its own sizes.
%
%  varargout = cat_io_dcm2bids_reportfig(action,varargin)
%
%  Actions (see the help of the local functions):
%    Ppdf = cat_io_dcm2bids_reportfig('writeReportFigure',Poutdir,job)
%
%  See also cat_io_dcm2bids_report, cat_io_dcm2bids_stats, cat_plot_boxplot.

  switch action
    case {'writeReportFigure'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_reportfig:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Ppdf = writeReportFigure(Poutdir,job)
%writeReportFigure. Report figures of the report tables of the BIDS
%  directory of the job for each design of opts.reportfig.design (see the
%  help of this file), returns the PDF files (cell).
  Prep = fullfile(Poutdir,[job.BIDSdir '-report']);
  Ptab = fullfile(Prep,'tables',job.BIDSsubdir);
  fmt  = job.opts.tableformat;
  S    = readTable(fullfile(Ptab,['scans.' fmt]));
  E    = readTable(fullfile(Ptab,['sessions.' fmt]));
  Q    = readTable(fullfile(Ptab,['project.' fmt]));
  designs = cellstr(job.opts.reportfig.design);
  Ppdf = strcat(fullfile(Prep,['report_' job.BIDSsubdir '_']),designs(:)','.pdf');
  if isempty(E), Ppdf = {}; return; end % (the PDFs of former runs are kept without tables)
  Pold = [{fullfile(Prep,['report_' job.BIDSsubdir '.pdf'])}, ...   % (all designs of former runs)
    strcat(fullfile(Prep,['report_' job.BIDSsubdir '_']),designNames(),'.pdf')];
  for oi = 1:numel(Pold), if exist(Pold{oi},'file'), delete(Pold{oi}); end, end
  ws = warning('off','MATLAB:print:ContentTypeImageSuggested'); % vector output of many objects
  for di = 1:numel(designs)
    job.opts.reportfig.design = designs{di};
    F = setupFigure(job,E,Q);

    % overview page and one page per site (for several sites, but not for 
    % more than opts.reportfig.maxsites sites, see manySites)
    pages = {''};
    if numel(F.sites) > 1 && ~manySites(F), pages = [pages F.sites]; end %#ok<AGROW>
    fh = figure('Visible','off','Color',F.bg,'Units','centimeters','Position',[1 1 F.pw F.ph], ...
      'PaperUnits','centimeters','PaperSize',[F.pw F.ph],'PaperPosition',[0 0 F.pw F.ph], ...
      'MenuBar','none','ToolBar','none');
    fhclean = onCleanup(@() close(fh));
    for pgi = 1:numel(pages)
      clf(fh); set(fh,'Color',F.bg);
      if isempty(pages{pgi})
        si = true(numel(S),1); ei = true(numel(E),1);
      else
        si = strcmp(textCol(S,'site'),pages{pgi}); ei = strcmp(textCol(E,'site'),pages{pgi});
      end
      drawPage(fh,F,S(si),E(ei),pages{pgi},pgi,numel(pages));
      try
        exportgraphics(fh,Ppdf{di},'ContentType','vector','BackgroundColor',F.bg,'Padding','figure','Append',pgi > 1);
      catch % (no Padding in older MATLAB versions, e.g. R2023b)
        exportgraphics(fh,Ppdf{di},'ContentType','vector','BackgroundColor',F.bg,'Append',pgi > 1);
      end
      if pgi == 1 && di == 1, showSPM(fh,F); end
    end
    clear fhclean
    fprintf('  Report figure: %s\n',Ppdf{di});
  end
  warning(ws);
end
% =========================================================================
function names = designNames()
%designNames. Names of the designs of the report figure (JSON files of
%  DCM2BIDS/colormaps, see loadDesign).
  D = dir(fullfile(fileparts(mfilename('fullpath')),'DCM2BIDS','colormaps','*.json'));
  names = regexprep({D.name},'\.json$','');
end
% =========================================================================
function F = setupFigure(job,E,Q)
%setupFigure. Common settings of all pages: page size, design (colors,
%  logos, and colormaps, see loadDesign), fonts, scan types, sites, and the
%  statistics of the whole sample (models and homogeneity thresholds),
%  i.e., the site pages show the same models.
  F.pw = 21; F.ph = 29.7; F.mx = 0.3;               % A4 portrait in cm and side margin
  F.project = job.subdir; F.dataset = job.BIDSsubdir;
  F.ropts   = job.opts.reportfig;
  [F.ropts.agerange,F.ropts.agestep] = ageRange(numCol(E,'age'),F.ropts); % (the same for all pages)
  F.fs      = F.ropts.fontsize;                     % text
  F.fsh     = F.fs * F.ropts.fontscale(1);          % header-like elements (titles of the diagrams)
  F.fss     = F.fs * F.ropts.fontscale(2);          % axes, table entries, and legends
  F.fst     = F.fs * 1.5;                           % title of the page (see drawHeader)
  D = loadDesign(F.ropts.design);
  F.design  = D.Name; F.logos = D.Logos; F.icon = D.Icons;
  F.fg  = D.Colors(1,:); F.fgl = D.Colors(2,:);    % foreground full and light (text)
  F.bgl = D.Colors(3,:); F.bg  = D.Colors(4,:);    % background light (panels) and full (page)
  F.mix  = @(c,a) a .* c + (1 - a) .* F.bg;        % color c with the opacity a on the page
  F.mixp = @(c,a) a .* c + (1 - a) .* F.bgl;       % color c with the opacity a on a panel
  F.marks   = [0.5 6.5];                            % rating range (worse ratings are clipped)
  MarkColor = cat_io_colormaps('marks+',40);        % quality colors (as command line)
  F.col2mark = @(val) min(1,max(0,MarkColor(min(size(MarkColor,1)-3, ...
    max(1,floor( val/9.5 * size(MarkColor,1)))),:)));
  [F.types,F.tnam] = scanTypes(job.opts);
  F.sites   = unique(textCol(E,'site'))'; F.sites(cellfun('isempty',F.sites)) = [];
  F.sitecol = colormapOf(D.Discrete,max(1,numel(F.sites)));
  qs = textCol(Q,'site'); qc = textCol(Q,'color');  % colors of the site dictionary (project table)
  for si = 1:numel(F.sites)
    h = qc(strcmp(qs,F.sites{si}));
    if ~isempty(h) && ~isempty(regexp(h{1},'^#[0-9a-fA-F]{6}$','once')), F.sitecol(si,:) = sscanf(h{1}(2:end),'%2x')' / 255; end
  end
  F.dtnam   = {'anat','func','dwi','fmap','perf'};  % datatypes of the scan ratings of one site (see drawRatings)
  F.dtcol   = colormapOf(D.Discrete,numel(F.dtnam) + 1); % (last: other datatypes)
  F.ccmap   = colormapOf(D.Continuous,64);          % continuous values (not used yet)
  F.cntcol  = [0.15 0.42 0.82];                     % numbers of the tables (more blue for larger values)
  F.agecol  = [0.96 0.62 0.12];                     % ages of the tables (more yellow-orange for older subjects, see ageColor)
  F.markbg  = 0.2;                                  % opacity of the mark colors in the background of the rating diagrams
  F.red     = [0.85 0.10 0.10];                     % outliers and thresholds (lighter for dark designs)
  if mean(F.bg) < 0.5, F.red = [1 0.35 0.35]; end
  F.sexnam  = {'M','F','O',''}; F.sexlab = {'male','female','other','n/a'};
  F.sexcol  = [0 0.443 0.737; 0.745 0.086 0.133; 0.55 0.55 0.55; 0.78 0.78 0.78]; % #0071bc, #be1622
  F.statnam = {'conform','semiconform','nonconform','untested','missing'};
  F.statcol = [F.col2mark(1); F.col2mark(3); F.col2mark(5.5); 0.62 0.62 0.62; F.mix([0.6 0.6 0.6],0.3)];
  F.tisnam  = {'GM','WM','CSF'};
  F.tiscol  = [0.20 0.60 0.30; 0.90 0.55 0.15; 0.25 0.45 0.85];
  F.Q       = Q;

  % models of the volumes and markers and thresholds of the homogeneity of
  % the whole sample (see cat_io_dcm2bids_stats)
  age = numCol(E,'age'); sex = textCol(E,'sex');
  Y   = 100 * [numCol(E,'T1w_rGMV') numCol(E,'T1w_rWMV') numCol(E,'T1w_rCSFV')];
  [~,F.M]  = cat_io_dcm2bids_stats('volumeModel',Y,age,sex);
  [~,F.Mm] = cat_io_dcm2bids_stats('volumeModel',[numCol(E,'dwi_FA_WM') numCol(E,'bold_segregation')],age,sex);
  F.qthr   = nan(1,3); fq = {'T1w_QMZ','dwi_QMZ','bold_QMZ'};
  for qi = 1:3
    q = numCol(E,fq{qi}); q = q(isfinite(q));
    if numel(q) >= 3, F.qthr(qi) = mean(q) + 2*std(q); end
  end
  % SQR of light and severe artifacts of the whole sample (see sqrLimits)
  F.yellow = [0.85 0.60 0.00]; if mean(F.bg) < 0.5, F.yellow = [1 0.80 0.30]; end
  F.sqrlim = [sqrLimits(numCol(E,'T1w_SQR')); sqrLimits(numCol(E,'dwi_SQR')); sqrLimits(numCol(E,'bold_SQR'))];
end
% =========================================================================
function L = sqrLimits(q)
%sqrLimits. Typical protocol quality (the 25th percentile of the SQR, i.e.
%  the border of the best 25% of the sample) and the borders of light and
%  severe artifacts [typical light severe]: one and two interquartile 
%  ranges (IQR) above the typical quality, at least 0.25 and 0.75 marks 
%  (e.g. motion artifacts of T1w data appear slightly about 0.5 and 
%  clearly about 1.0 marks below the typical quality), i.e. the IQR (robust
%  to the artifacts, two IQR above the 25th percentile are about two SD 
%  above the mean of normally distributed data) widens the borders of a 
%  heterogeneous sample (e.g. several sites or protocols), whereas the 
%  minimum keeps them relevant for a very homogeneous sample. NaN for less
%  than 5 values.
  q = q(isfinite(q)); L = nan(1,3);
  if numel(q) < 5, return; end
  q25 = prctile(q,25); iqr = prctile(q,75) - q25;
  L = [q25 q25 + max(0.25,iqr) q25 + max(0.75,2*iqr)];
end
% =========================================================================
function E = sessionPCP(E,S)
%sessionPCP. Protocol conformity percentage of the scans of each session
%  (with protocol test) as mark (100% as 1, 0% as 6, see the header), NaN 
%  without tested scans.
  if isempty(E), return; end
  [E.PCP] = deal(nan);
  st = textCol(S,'protocol_status'); ok = ismember(st,{'conform','semiconform','nonconform'});
  sk = strcat(textCol(S,'participant_id'),'|',textCol(S,'session_id'));
  ek = strcat(textCol(E,'participant_id'),'|',textCol(E,'session_id'));
  for ei = 1:numel(E)
    m = ok & strcmp(sk,ek{ei});
    if any(m), E(ei).PCP = 1 + 5 * (1 - mean(strcmp(st(m),'conform'))); end
  end
end
% =========================================================================
function D = loadDesign(name)
%loadDesign. Design of the report figure (DCM2BIDS/colormaps/<name>.json,
%  see the help of this file) with the logos and the icon as files of the
%  CAT directory (only existing ones), and the bright CAT design for
%  unknown designs.
  Pcat = fileparts(mfilename('fullpath'));
  Pd   = fullfile(Pcat,'DCM2BIDS','colormaps');
  D    = struct('Name','CAT_bright','Logos',{{}},'Icons','', ...
    'Colors',[0.15 0.15 0.17; 0.45 0.45 0.50; 0.93 0.94 0.96; 1 1 1],'Discrete','nejm','Continuous','cool');
  P    = fullfile(Pd,[char(name) '.json']);
  if ~exist(P,'file')
    cat_io_cprintf('warn','  Report figure: unknown design "%s" (CAT_bright is used).\n',char(name));
    P = fullfile(Pd,'CAT_bright.json');
  end
  if exist(P,'file')
    J = jsondecode(fileread(P));
    for fn = fieldnames(D)'
      if isfield(J,fn{1}), D.(fn{1}) = J.(fn{1}); end
    end
  end
  D.Logos = cellstr(D.Logos);
  D.Logos = fullfile(Pcat,D.Logos(~cellfun('isempty',D.Logos)));
  D.Logos = D.Logos(cellfun(@(f) exist(f,'file') > 0,D.Logos));
  if ~isempty(D.Icons), D.Icons = fullfile(Pcat,char(D.Icons)); end
  if ~exist(D.Icons,'file'), D.Icons = ''; end
  D.Colors = double(D.Colors);
  if ~isequal(size(D.Colors),[4 3]), error('cat_io_dcm2bids_reportfig:design','The colors of the design need 4 RGB rows.'); end
end
% =========================================================================
function C = colormapOf(c,n)
%colormapOf. Colormap c (name of cat_io_colormaps or of a MATLAB colormap,
%  or RGB rows that are repeated if needed) with n colors.
  if ischar(c) || isstring(c)
    try
      C = cat_io_colormaps(char(c),n);
    catch
      try, C = feval(char(c),n); catch, C = lines(n); end
    end
  else
    C = double(c); C = C(mod(0:n-1,size(C,1)) + 1,:);
  end
  C = min(1,max(0,C(1:n,1:3)));
end
% =========================================================================
function drawPage(fh,F,S,E,site,pgi,npg)
%drawPage. One page of the report figure with the scans S and sessions E
%  of all sites (site = '') or of one site in the rows of the help.
  w = F.pw - 2*F.mx;
  drawHeader(fh,F,S,E,site,[F.mx F.ph-F.mx-2.1 w 2.1]);

  % data: sessions per month and sites (panel)
  % (the site table gets more space for more site columns, about 1.3 cm 
  % each, down to a third of the panel for the timeline)
  p = [F.mx 22.5 w 4.6]; drawPanel(fh,F,p);
  nc = ~isempty(site) + numel(F.sites) * (isempty(site) && ~manySites(F)) + (numel(F.sites) > 1); % (see drawSiteTable)
  ws = min(p(3) * 2/3, max(p(3) - 12.2, 1.8 + 1.3 * nc));
  drawTimeline(fh,F,E,[p(1)+0.15 p(2)+0.1 p(3)-ws-0.6 p(4)-0.2]);
  drawSiteTable(fh,F,site,[p(1)+p(3)-ws-0.15 p(2)+0.15 ws p(4)-0.3]);

  % sample: age/sex table with the sex pie, age pyramid, and subjects per
  % site (panel)
  p = [F.mx 17.65 w 4.6]; drawPanel(fh,F,p);
  drawSexTable(fh,F,E,[p(1)+0.15 p(2)+0.1 6.6 p(4)-0.2]);
  drawAge(fh,F,E,[p(1)+6.95 p(2)+0.1 6.1 p(4)-0.2]);
  drawSiteSample(fh,F,E,[p(1)+13.25 p(2)+0.1 p(3)-13.4 p(4)-0.2]);

  % brain markers by age (panel)
  c = panelCells(fh,F,[F.mx 12.0 w 5.4],3);
  drawVolumes(fh,F,E,c(1,:));
  drawMarker(fh,F,E,c(2,:),'dwi_FA_WM',F.Mm(1),'White matter FA (best DWI scan)','FA',3);
  drawMarker(fh,F,E,c(3,:),'bold_segregation',F.Mm(2),'System segregation (best BOLD scan)','segregation',2);

  % sample homogeneity (corrected for age and sex, panel)
  c = panelCells(fh,F,[F.mx 6.35 w 5.4],3);
  drawHomogeneity(fh,F,E,c(1,:),'T1w_QMZ','T1w_SQR',F.qthr(1),'Homogeneity of the GM maps (mwc1)','T1w');
  drawHomogeneity(fh,F,E,c(2,:),'dwi_QMZ','dwi_SQR',F.qthr(2),'Homogeneity of the FA maps (wFA)','DWI');
  drawHomogeneity(fh,F,E,c(3,:),'bold_QMZ','bold_SQR',F.qthr(3),'Homogeneity of the connectivity','BOLD');

  % quality control (panel)
  c = panelCells(fh,F,[F.mx 0.6 w 5.5],3);
  rn = {'PFR','PCP','BBR','SOR','PRR','SQR'}; % (the order of the header)
  if isscalar(unique(textCol(E,'site'))) % one site: the ratings of all scans by datatype
    drawRatings(fh,F,scanPCP(S),rn,rn,c(1,:),'Scan ratings','datatype');
  else
    drawRatings(fh,F,sessionPCP(E,S),rn,rn,c(1,:),'Session ratings');
  end
  drawRatings(fh,F,E,strcat(F.tnam,'_SQR'),F.tnam,c(2,:),'Image quality SQR (best scans)');
  drawCompleteness(fh,F,S,E,c(3,:));

  % footer
  ax = textAxes(fh,F,[F.mx 0.15 w 0.35]);
  text(ax,0,0.5,sprintf('cat_io_dcm2bids report figure of %s/%s, %s (design %s)',F.project,F.dataset, ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm')),F.design),'FontSize',F.fss,'Color',F.fgl,'Interpreter','none');
  text(ax,1,0.5,sprintf('page %d/%d',pgi,npg),'FontSize',F.fss,'Color',F.fgl,'HorizontalAlignment','right');
end
% =========================================================================
function c = cells(p,n)
%cells. n equal columns [x y w h] (rows of c, in cm) of the row p.
  w = p(3) / n; c = [p(1) + (0:n-1)'*w repmat(p(2),n,1) repmat(w,n,1) repmat(p(4),n,1)];
end
% =========================================================================
function [ax,tax] = cellAxes(fh,F,c,ttl,side)
%cellAxes. Axes of a diagram in the cell c (in cm) with the title ttl and
%  margins for the tick labels, and invisible axes of the width side (cm,
%  0 for none) right of the diagram for statistics and legends.
  th  = 0.42;                                       % title height
  text(textAxes(fh,F,[c(1) c(2)+c(4)-th c(3) th]),0.02,0.5,ttl,'FontSize',F.fsh,'FontWeight','bold', ...
    'Color',F.fg,'Interpreter','none','VerticalAlignment','middle');
  ml = 1.0; mb = 0.85; mr = 0.15 + side;           % left, bottom, and right margin
  ax = newAxes(fh,F,[c(1)+ml c(2)+mb c(3)-ml-mr c(4)-th-mb-0.1]);
  tax = [];
  if side > 0, tax = textAxes(fh,F,[c(1)+c(3)-side c(2)+mb side c(4)-th-mb-0.1]); end
end
% =========================================================================
function drawPanel(fh,F,p)
%drawPanel. Background of a panel (light background of the design).
  ax = textAxes(fh,F,p);
  patch(ax,[0 1 1 0],[0 0 1 1],F.bgl,'EdgeColor','none');
end
% =========================================================================
function c = panelCells(fh,F,p,n)
%panelCells. Panel p (see drawPanel) with n equal cells (see cells) inside.
  drawPanel(fh,F,p); c = cells([p(1)+0.15 p(2)+0.1 p(3)-0.3 p(4)-0.2],n);
end
% =========================================================================
function w = drawImage(fh,F,P,pos,align,square)
%drawImage. Image file P (e.g. PNG with transparency) of the height pos(4)
%  (cm) at the position pos (x is the right edge for align = 'right') with
%  its aspect ratio, or for square fitted into the square of this height
%  (at the top and at the side of align, e.g. the logos), returns the width
%  w (0 without image).
  w = 0;
  if isempty(P) || ~exist(P,'file'), return; end
  if ~exist('align','var') || isempty(align), align = 'left'; end
  if ~exist('square','var'), square = false; end
  try
    [I,map,A] = imread(P);
    if ~isempty(map), I = ind2rgb(I,map); end
    ar = size(I,2) / size(I,1); h = pos(4); wi = h * ar; w = wi;   % height and width of the image
    if square, w = pos(4); if ar > 1, wi = w; h = w / ar; end, end
    x = pos(1); if strcmp(align,'right'), x = pos(1) - wi; pos(1) = pos(1) - w; end
    ax = axes('Parent',fh,'Units','normalized','Position',[x/F.pw (pos(2) + pos(4) - h)/F.ph wi/F.pw h/F.ph], ...
      'Visible','off','HitTest','off');
    if isempty(A), image(ax,I); else, image(ax,I,'AlphaData',double(A) / double(intmax(class(A)))); end
    axis(ax,'image','off'); ax.Toolbar = []; disableDefaultInteractivity(ax);
  catch e
    cat_io_cprintf('warn','  Report figure: image "%s" not shown: %s\n',P,e.message); w = 0;
  end
end
% =========================================================================
function drawHeader(fh,F,S,E,site,p)
%drawHeader. Logos of the design (left and right), title, numbers of
%  subjects, sessions, and scans, the scanner of site pages, and the table
%  of the median ratings of the sessions with icons (fonts of 85% of the
%  page).
  if isempty(site) && isscalar(F.sites), site = F.sites{1}; end % project of one site
  hs = 0.85;                                        % font scale of the header
  lh = 1.6; x0 = p(1); x1 = p(1) + p(3); yl = p(2) + p(4) - lh - 0.1; % (square space of the logos)
  if ~isempty(F.logos), x0 = x0 + drawImage(fh,F,F.logos{1},[x0 yl 0 lh],'left',true) + 0.35; end
  if numel(F.logos) > 1, x1 = x1 - drawImage(fh,F,F.logos{2},[x1 yl 0 lh],'right',true) - 0.35; end

  % title, BIDS directory, numbers, and scanner
  tw = 6.4;                                         % width of the rating table
  ax = textAxes(fh,F,[x0 p(2) x1 - x0 - tw - 0.2 p(4)]);
  text(ax,0,0.86,sprintf('DCM2BIDS report: %s',F.project),'FontSize',F.fst*hs,'FontWeight','bold', ...
    'Interpreter','none','Color',F.fg);
  if isempty(site)
    str = sprintf('BIDS directory %s: overview of %d sites',F.dataset,numel(F.sites));
  else
    str = sprintf('BIDS directory %s: site %s',F.dataset,site);
  end
  text(ax,0,0.62,str,'FontSize',F.fs*1.2*hs,'Interpreter','none','Color',F.fg);
  nsub = numel(unique(textCol(E,'participant_id')));
  text(ax,0,0.41,sprintf('%d subjects, %d sessions, %d scans',nsub,numel(E),numel(S)), ...
    'FontSize',F.fs*hs,'Interpreter','none','Color',F.fg);
  if ~isempty(site) % scanner of the site
    sc = strtrim(strcat(textCol(S,'manufacturer'),{' '},textCol(S,'model')));
    fs = numCol(S,'field_strength'); fs = fs(isfinite(fs));
    if ~isempty(fs), sc = strcat(sc,sprintf(' %gT',median(fs))); end
    text(ax,0,0.22,sprintf('%s, software %s',mostFrequent(sc),valueCounts(textCol(S,'software'))), ...
      'FontSize',F.fss*hs,'Interpreter','none','Color',F.fgl);
  end

  % median ratings of the sessions and protocol conformity of the scans
  % (fully written names with the icon of the design)
  st = textCol(S,'protocol_status'); st = st(~cellfun('isempty',st));
  if isempty(st), pcp = nan; else, pcp = 100 * mean(strcmp(st,'conform')); end
  R  = {'Protocol fit rating (PFR)',             median(numCol(E,'PFR'),'omitnan'), 0
        'Protocol conformity percentage (PCP)',  pcp,                               1
        'Brain coverage rating (BBR)',           median(numCol(E,'BBR'),'omitnan'), 0
        'Session orientation rating (SOR)',      median(numCol(E,'SOR'),'omitnan'), 0
        'Phase-encoding reversal rating (PRR)',  median(numCol(E,'PRR'),'omitnan'), 0
        'Scan quality rating (SQR)',             median(numCol(E,'SQR'),'omitnan'), 0};
  tx = x1 - tw; rh = (p(4) - 0.1) / size(R,1); y0 = p(2) + p(4) - 0.1 - rh;
  for ri = 1:size(R,1)
    y = y0 - (ri-1) * rh;
    wi = 0; if ~isempty(F.icon), wi = drawImage(fh,F,F.icon,[tx y+0.06 0 rh-0.12]); end
    ax = textAxes(fh,F,[tx + wi + 0.12 y tw - wi - 0.12 rh]);
    v  = R{ri,2};
    if R{ri,3}, mark = 1 + 5 * (1 - v/100); vstr = sprintf('%0.0f%%',v); % 100% as 1, 0% as 6
    else,       mark = v; vstr = sprintf('%0.2f',v);
    end
    if isfinite(v), col = F.mix(F.col2mark(mark),0.6); else, col = F.mix(F.fgl,0.15); vstr = 'n/a'; end
    text(ax,0,0.5,R{ri,1},'FontSize',F.fs*hs,'Color',F.fg,'VerticalAlignment','middle');
    bx = 1 - 1.0 / (tw - wi - 0.12);               % value box of 1 cm
    patch(ax,[bx 1 1 bx],[0.1 0.1 0.9 0.9],col,'EdgeColor','none');
    text(ax,(bx+1)/2,0.5,vstr,'FontSize',F.fs*hs,'FontWeight','bold','HorizontalAlignment','center', ...
      'VerticalAlignment','middle','Color',F.fg);
  end
  ax = textAxes(fh,F,[p(1) p(2)-0.02 p(3) 0.04]);
  line(ax,[0 1],[0.5 0.5],'Color',F.fgl,'LineWidth',0.5);
end
% =========================================================================
function drawTimeline(fh,F,E,c)
%drawTimeline. Sessions per month (stacked by site) of the session date of
%  the session labels (ses-YYYYMMDD) with the number of each month.
  sid = textCol(E,'session_id');
  tok = regexp(sid,'^ses-(\d{4})(\d{2})\d{2}','tokens','once');
  ym  = nan(numel(sid),2);
  for i = 1:numel(tok), if ~isempty(tok{i}), ym(i,:) = str2double(tok{i}); end, end
  ok  = ym(:,1) >= 1980 & ym(:,1) <= 2100 & ym(:,2) >= 1 & ym(:,2) <= 12;
  str = 'Sessions per month';
  if any(ok)
    mi  = (ym(:,1) - 1) * 12 + ym(:,2);              % month index
    m0  = min(mi(ok)); nm = max(mi(ok)) - m0 + 1;
    lab = arrayfun(@(m) sprintf('%04d-%02d',floor((m-1)/12)+1,mod(m-1,12)+1),m0:m0+nm-1,'UniformOutput',false);
    str = sprintf('%s (%d sessions, %s to %s)',str,sum(ok),lab{1},lab{end});
    if any(~ok), str = sprintf('%s, %d without date',str,sum(~ok)); end
    if c(3) < 10, str = sprintf('Sessions (%d, %s to %s)',sum(ok),lab{1},lab{end}); end % (narrow cell for many sites)
  end
  ax = cellAxes(fh,F,c,str,0);
  if ~any(ok)
    axis(ax,'off');
    text(ax,0.5,0.5,'no session dates (e.g. coded session labels of anonymize > 1)','FontSize',F.fs, ...
      'HorizontalAlignment','center','Color',F.fgl);
    return
  end
  site = textCol(E,'site');
  sites = F.sites; if isempty(sites), sites = {''}; end
  cnt = @(sel) accumarray(reshape(mi(sel),[],1) - m0 + 1,1,[nm 1]); % sessions per month
  C = zeros(nm,numel(sites));
  for si = 1:numel(sites), C(:,si) = cnt(ok & strcmp(site,sites{si})); end
  C(:,end+1) = cnt(ok & ~ismember(site,sites));   % sessions without site
  hb = bar(ax,1:nm+1,[C; zeros(1,size(C,2))],1,'stacked','EdgeColor',F.bg,'LineWidth',0.2); % (empty row to stack one month)
  for si = 1:numel(hb)
    if si <= numel(sites), set(hb(si),'FaceColor',siteColor(F,sites{si})); else, set(hb(si),'FaceColor',[.6 .6 .6]); end
  end
  hold(ax,'on');
  tot = sum(C,2);
  if nm <= 48
    for i = find(tot' > 0)
      text(ax,i,tot(i),num2str(tot(i)),'FontSize',F.fss,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','Color',F.fg);
    end
  end
  step = max(1,ceil(nm/max(2,round(c(3)/1.6))));
  set(ax,'XLim',[0.4 nm+0.6],'YLim',[0 max(tot)*1.3 + 0.5],'XTick',1:step:nm,'XTickLabel',lab(1:step:nm));
  ylabel(ax,'sessions','FontSize',F.fss,'Color',F.fg);
end
% =========================================================================
function drawSiteTable(fh,F,site,c)
%drawSiteTable. Site table of the project table with the sites as columns
%  (overview page: all sites and all; site page: the site and all) and the
%  numbers of subjects, sessions, and scans (blue boxes, more blue for more
%  data of a site), the median ratings (colors of the marks), and the
%  protocol conformity percentage as rows, i.e. also the legend of the site
%  colors.
  Q  = F.Q; qs = textCol(Q,'site');
  if isempty(site), cols = F.sites; else, cols = {site}; end
  if isempty(site) && manySites(F), cols = {}; end % (only all sites together)
  if any(strcmp(qs,'all')) && ~isscalar(F.sites), cols = [cols {'all'}]; end
  if ~isfield(Q,'PCP') && isfield(Q,'conform'), [Q.PCP] = Q.conform; end % (tables of former versions)
  rows = {'subjects','subjects'; 'sessions','sessions'; 'scans','scans'; ...
          'SQR (median)','SQR'; 'PFR (median)','PFR'; 'SOR (median)','SOR'; 'PCP','PCP'};
  ax = textAxes(fh,F,c);
  ttl = 'Sites'; if isempty(site) && manySites(F), ttl = sprintf('Sites (%d, all together)',numel(F.sites)); end
  text(ax,0.0,1.0,ttl,'FontSize',F.fsh,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
  nc = numel(cols); x0 = 1.8 / c(3); cw = (1 - x0) / max(1,nc); rh = 0.88 / (size(rows,1) + 1.3); % (row names 1.8 cm)
  V  = nan(size(rows,1),nc);
  for ci = 1:nc
    qi = find(strcmp(qs,cols{ci}),1);
    for ri = 1:size(rows,1)
      if ~isempty(qi) && isfield(Q,rows{ri,2}) && isnumeric(Q(qi).(rows{ri,2})) && ~isempty(Q(qi).(rows{ri,2}))
        V(ri,ci) = double(Q(qi).(rows{ri,2}));
      end
    end
  end
  ns = nc - any(strcmp(cols,'all'));               % site columns (scaling of the numbers)
  y  = 0.88 - rh*0.65;                             % header (1.3 rows) with the site colors
  hh = gobjects(1,nc);
  for ci = 1:nc
    xc = x0 + (ci - 0.5) * cw;
    if ~strcmp(cols{ci},'all')
      patch(ax,xc + cw*[-0.42 0.42 0.42 -0.42],y + rh*[0.47 0.47 0.62 0.62],siteColor(F,cols{ci}),'EdgeColor','none');
    end
    hh(ci) = text(ax,xc,y - rh*0.05,cols{ci},'FontSize',F.fss,'FontWeight','bold','Color',F.fg, ...
      'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none');
  end
  fitTexts(hh,cw*0.94,F.fss,true);                 % (two lines for long site names)
  hv = gobjects(0);
  for ri = 1:size(rows,1)
    y = 0.88 - (ri + 0.8) * rh;
    text(ax,0,y,rows{ri,1},'FontSize',F.fss,'Color',F.fg,'VerticalAlignment','middle');
    vmax = max([V(ri,1:ns) eps]);
    for ci = 1:nc
      v = V(ri,ci); xc = x0 + (ci - 0.5) * cw; col = [];
      if ri <= 3 && isfinite(v)                    % numbers (blue boxes)
        str = sprintf('%d',round(v));
        if ~strcmp(cols{ci},'all'), col = F.mixp(F.cntcol,0.6 * min(1,v / vmax)); end
      elseif ri <= 6 && isfinite(v)                % ratings
        str = sprintf('%0.2f',v); col = F.mixp(F.col2mark(v),0.6);
      elseif isfinite(v)                           % protocol conformity percentage
        str = sprintf('%0.0f%%',v); col = F.mixp(F.col2mark(1 + 5*(1 - v/100)),0.6);
      else
        str = 'n/a';
      end
      if ~isempty(col), patch(ax,xc + cw*[-0.45 0.45 0.45 -0.45],y + rh*[-0.42 -0.42 0.42 0.42],col,'EdgeColor','none'); end
      fw = 'normal'; if strcmp(cols{ci},'all'), fw = 'bold'; end
      hv(end+1) = text(ax,xc,y,str,'FontSize',F.fss,'Color',F.fg,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontWeight',fw); %#ok<AGROW>
    end
  end
  fitTexts(hv,cw*0.86,F.fss,false);
end
% =========================================================================
function fitTexts(h,w,fs,split)
%fitTexts. Common font size of the texts h (one line, default size fs) for
%  the width w (data units of their axes): first two lines for split
%  (at the change of case, letters and digits, or separators next to the
%  middle, with 90% of the font size), then smaller fonts down to 60% of
%  fs, and finally shortened texts.
  if isempty(h), return; end
  ext = @() max(arrayfun(@(hi) hi.Extent(3),h));
  if ext() <= w, return; end
  fs1 = fs;
  if split
    for hi = 1:numel(h), set(h(hi),'String',splitLabel(get(h(hi),'String'))); end
    fs1 = 0.9 * fs; set(h,'FontSize',fs1);
  end
  e = ext();
  if e > w, set(h,'FontSize',max(0.6 * fs,fs1 * w / e)); end
  for hi = 1:numel(h) % (shortened by the ratio of the widths, i.e. few steps)
    s = cellstr(get(h(hi),'String')); e = h(hi).Extent(3);
    while e > w && max(cellfun('length',s)) > 3
      [n,li] = max(cellfun('length',s)); s0 = s{li};
      s{li} = shortLabel(s0,max(3,min(n - 1,floor(n * w / e))));
      if strcmp(s{li},s0), break; end % (no shorter label, e.g. 'x...' for 3 characters)
      set(h(hi),'String',s); e = h(hi).Extent(3);
    end
  end
end
% =========================================================================
function many = manySites(F)
%manySites. True for more than opts.reportfig.maxsites sites (default 20), 
%  which are only shown together (without site columns, site labels, and 
%  site pages), as the diagrams are not readable for so many sites.
  maxs = 20; if isfield(F.ropts,'maxsites'), maxs = F.ropts.maxsites; end
  many = numel(F.sites) > maxs;
end
% =========================================================================
function s = splitLabel(s)
%splitLabel. Label s in two lines (cell) at the change of case (not after
%  the first capital letter of a word), letters and digits, or at a
%  separator next to its middle (unchanged for labels below 5 characters).
  s = char(s); n = numel(s);
  if n < 5, return; end
  sep = find(ismember(s,'_- .'));
  if ~isempty(sep)
    [~,i] = min(abs(sep - (n+1)/2)); b = sep(i);
    s = {s(1:b-1), s(b+1:end)}; return
  end
  cls = 1*isstrprop(s,'lower') + 2*isstrprop(s,'upper') + 3*isstrprop(s,'digit');
  k   = find(cls(2:end) ~= cls(1:end-1) & ~(cls(1:end-1) == 2 & cls(2:end) == 1)) + 1;
  if isempty(k), k = ceil(n/2) + 1; end
  [~,i] = min(abs(k - (n/2 + 1))); b = k(i);
  s = {s(1:b-1), s(b:end)};
end
% =========================================================================
function [age,sex] = subjects(E)
%subjects. Age and sex of the subjects (first session of each subject) of
%  the sessions E (sex O for other and '' for unknown).
  [~,ui] = unique(textCol(E,'participant_id'),'stable');
  age = numCol(E(ui),'age'); sex = textCol(E(ui),'sex');
  sex(~ismember(sex,{'M','F','O'})) = {''};
end
% =========================================================================
function drawSexTable(fh,F,E,c)
%drawSexTable. Table of the subjects of each sex (pictograms) in the row
%  height of the site table with the number and percentage (blue boxes, 
%  more blue for more subjects), and the mean age, its SD, and the age 
%  range (yellow-orange boxes for older subjects and larger SD and range, 
%  see ageColor), and the pie of the sex of the subjects.
  [age,sex] = subjects(E); n = numel(sex);
  ax = textAxes(fh,F,c);
  text(ax,0,1,sprintf('Subjects (%d)',n),'FontSize',F.fsh,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
  grp = {'all',true(size(sex)),0; 'male',strcmp(sex,'M'),1; 'female',strcmp(sex,'F'),2; ...
         'other/n/a',~ismember(sex,{'M','F'}),3};
  % columns (left edge and width in cm) and rows (row 0 is the header)
  X   = [1.32 1.15; 2.55 0.62; 3.25 0.72; 4.05 0.72]; hdr = {'n (%)','age','SD','range'};
  rh  = 0.45; yc = @(r) 1 - (0.46 + (r + 0.5) * rh) / c(4);
  box = @(ci,r,col) patch(ax,(X(ci,1) + [0 X(ci,2) X(ci,2) 0]) / c(3),yc(r) + 0.42 * rh / c(4) * [-1 -1 1 1], ...
    col,'EdgeColor','none');
  val = @(ci,r,str,fw) text(ax,(X(ci,1) + X(ci,2)/2) / c(3),yc(r),str,'FontSize',F.fss,'FontWeight',fw, ...
    'Color',F.fg,'HorizontalAlignment','center','VerticalAlignment','middle');
  for ci = 1:4, val(ci,0,hdr{ci},'bold'); end
  nmax = max([cellfun(@sum,grp(2:end,2)); 1]);
  for gi = 1:size(grp,1)
    g = grp{gi,2}; a = age(g & isfinite(age)); y = yc(gi); fw = 'normal';
    if grp{gi,3} > 0
      drawPerson(ax,c,0.15 / c(3),y,0.34,grp{gi,3},F.sexcol(grp{gi,3},:));
      col = F.sexcol(grp{gi,3},:); if grp{gi,3} == 3, col = F.fgl; end
      box(1,gi,F.mixp(F.cntcol,0.6 * sum(g) / nmax));
    else
      col = F.fg; fw = 'bold';                       % (all subjects without box as in the site table)
    end
    text(ax,0.35 / c(3),y,grp{gi,1},'FontSize',F.fss,'Color',col,'VerticalAlignment','middle','FontWeight',fw);
    val(1,gi,sprintf('%d (%0.0f%%)',sum(g),100 * sum(g) / max(1,n)),fw);
    if ~isempty(a), box(2,gi,ageColor(F,mean(a),'mean')); val(2,gi,sprintf('%0.1f',mean(a)),fw); end
    if numel(a) > 1
      box(3,gi,ageColor(F,std(a),'sd'));             val(3,gi,sprintf('\\pm%0.1f',std(a)),fw);
      box(4,gi,ageColor(F,max(a) - min(a),'range')); val(4,gi,sprintf('%0.0f-%0.0f',min(a),max(a)),fw);
    end
  end
  na = sum(~isfinite(age));
  if na > 0, text(ax,0,0.02,sprintf('%d subject(s) without age',na),'FontSize',F.fss,'Color',F.fgl,'VerticalAlignment','bottom'); end
  % pie right of the table (centered on its rows)
  cnt = cellfun(@(s) sum(strcmp(sex,s)),F.sexnam);
  d   = min(1.5,c(3) - sum(X(end,:)) - 0.2);       % diameter in cm
  pax = newAxes(fh,F,[c(1) + c(3) - d - 0.05, c(2) + c(4) - 0.46 - 2.5 * rh - d/2, d, d]);
  if sum(cnt) > 0
    nz = find(cnt > 0); h = pie(pax,cnt(nz));
    for i = 1:numel(nz) % percentage in the slice
      set(h(2*i-1),'FaceColor',F.sexcol(nz(i),:),'EdgeColor',F.bgl);
      if isscalar(nz), set(h(2*i-1),'EdgeColor','none'); end % (no edge of a full pie)
      pc = 100 * cnt(nz(i)) / sum(cnt); str = '';
      if pc >= 8, str = sprintf('%0.0f%%',pc); end
      set(h(2*i),'String',str,'Position',0.55 * get(h(2*i),'Position'),'FontSize',F.fss, ...
        'FontWeight','bold','Color','w','HorizontalAlignment','center');
    end
  end
  axis(pax,'equal','off');
end
% =========================================================================
function col = ageColor(F,v,kind)
%ageColor. Box color (on a panel) of the age value v in years of the kind 
%  'mean' (mean age), 'sd' (its SD), or 'range' (range of the ages) with 
%  more yellow-orange for older subjects (relative to 0 to 100 years, i.e.
%  not the age range of the diagrams of the sample, see ageRange), a larger
%  SD (up to 25 years), and a larger range (up to 100 years), i.e. 
%  transparent for young subjects with similar ages.
  ar = [0 100]; r = diff(ar);
  switch kind
    case 'mean', a = (v - ar(1)) / r;
    case 'sd',   a = v / (r / 4);
    otherwise,   a = v / r;
  end
  col = F.mixp(F.agecol,0.6 * min(1,max(0,a)));
end
% =========================================================================
function drawPerson(ax,c,x,y,s,kind,col)
%drawPerson. Pictogram of a person (kind 1 male, 2 female, 3 other) of the
%  height s (cm) centered at x,y (0-1) of the invisible axes ax of the cell
%  c (cm, for the aspect ratio).
  sx = s / c(3); sy = s / c(4);
  t  = linspace(0,2*pi,24);
  patch(ax,x + 0.12*sx*cos(t),y + (0.37 + 0.12*sin(t))*sy,col,'EdgeColor','none');   % head
  switch kind
    case 1, bx = [-0.17 -0.17 0.17 0.17]; by = [-0.5 0.2 0.2 -0.5];                   % trousers/shirt
    case 2, bx = [-0.25 -0.09 0.09 0.25]; by = [-0.5 0.2 0.2 -0.5];                   % dress
    otherwise, bx = [-0.21 -0.13 0.13 0.21]; by = [-0.5 0.2 0.2 -0.5];
  end
  patch(ax,x + bx*sx,y + by*sy,col,'EdgeColor','none');
end
% =========================================================================
function drawAge(fh,F,E,c)
%drawAge. Age pyramid of the subjects (age at the first session, males
%  left, females right, other or unknown sex half on each side in the 
%  center) with the counts of each bin.
  [age,sex] = subjects(E);
  edges = F.ropts.agerange(1):F.ropts.agestep:F.ropts.agerange(2);
  cen   = edges(1:end-1) + diff(edges)/2;
  agec  = min(edges(end),max(edges(1),age));      % ages outside the range in the first/last bin
  cM = histcounts(agec(strcmp(sex,'M')),edges);
  cF = histcounts(agec(strcmp(sex,'F')),edges);
  cU = histcounts(agec(~ismember(sex,{'M','F'})),edges);
  ax = cellAxes(fh,F,c,sprintf('Age at the first session (%g-year bins)',F.ropts.agestep),0);
  hold(ax,'on');
  xm = max([cU/2 + cM, cU/2 + cF, 1]) * 1.3;
  plot(ax,[0 0],edges([1 end]),'Color',F.fgl,'LineWidth',0.5); % (below the bars and numbers)
  pyramidBars(ax,F,[cen cen(end) + F.ropts.agestep],[cM; cF; cU]',1,1,xm,c(3) - 1.15); % (empty row to stack one bin)
  xt = 0:niceStep(max([cM cF 1])/3):max([cM cF 1]);
  set(ax,'XLim',[-xm xm],'YLim',edges([1 end]),'YTick',ageTicks(F), ...
    'XTick',unique([-xt xt]),'XTickLabel',arrayfun(@(x) num2str(abs(x)),unique([-xt xt]),'UniformOutput',false));
  xlabel(ax,'subjects','FontSize',F.fss,'Color',F.fg); ylabel(ax,'age in years','FontSize',F.fss,'Color',F.fg);
  text(ax,-xm*0.95,edges(end) - diff(edges([1 end]))*0.04,F.sexlab{1},'FontSize',F.fss,'Color',F.sexcol(1,:), ...
    'FontWeight','bold','VerticalAlignment','top');
  text(ax, xm*0.95,edges(end) - diff(edges([1 end]))*0.04,F.sexlab{2},'FontSize',F.fss,'Color',F.sexcol(2,:), ...
    'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','right');
  if any(cU), text(ax,0,edges(end) - diff(edges([1 end]))*0.04,'other/n/a','FontSize',F.fss,'Color',F.sexcol(3,:), ...
    'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','center','BackgroundColor',F.bg,'Margin',0.5); end
end
% =========================================================================
function pyramidBars(ax,F,y,C,bw,labels,xm,aw)
%pyramidBars. Bars of a pyramid at the positions y with the counts C (rows
%  of males, females, and other or unknown sex) of the width bw (1 for 
%  touching bars): males left, females right, and other or unknown sex half
%  on each side in the center (as they are neither male nor female), with 
%  the numbers at the bars for labels (and in the center bar if it is wide
%  enough for the x-limit xm and the axes width aw in cm).
  n  = size(C,1); C(end+1:numel(y),:) = 0; h = C(:,3) / 2; 
  hl = barh(ax,y,-[h C(:,1)],bw,'stacked','EdgeColor',F.bg,'LineWidth',0.3);
  hr = barh(ax,y, [h C(:,2)],bw,'stacked','EdgeColor',F.bg,'LineWidth',0.3);
  set([hl(1) hr(1)],'FaceColor',F.sexcol(3,:)); set(hl(2),'FaceColor',F.sexcol(1,:)); set(hr(2),'FaceColor',F.sexcol(2,:));
  if ~labels, return; end
  cmpu = aw / (2 * xm);                            % cm per subject
  for i = 1:n
    if C(i,1) > 0, text(ax,-h(i) - C(i,1) - xm*0.02,y(i),num2str(C(i,1)),'FontSize',F.fss, ...
        'HorizontalAlignment','right','VerticalAlignment','middle','Color',F.fg); end
    if C(i,2) > 0, text(ax, h(i) + C(i,2) + xm*0.02,y(i),num2str(C(i,2)),'FontSize',F.fss, ...
        'HorizontalAlignment','left','VerticalAlignment','middle','Color',F.fg); end
    if C(i,3) > 0 && C(i,3) * cmpu > 0.14 * numel(num2str(C(i,3))) + 0.08
      text(ax,0,y(i),num2str(C(i,3)),'FontSize',F.fss,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Color','w');
    end
  end
end
% =========================================================================
function drawSiteSample(fh,F,E,c)
%drawSiteSample. Pyramid of the subjects per site (males left, females
%  right, gray other or unknown sex half on each side in the center, see
%  pyramidBars) with the mean age and its SD of each site in boxes next to
%  the pyramid (colors as the age/sex table, see ageColor).
  sites = F.sites; site = textCol(E,'site');
  if isempty(sites), sites = {''}; end
  sites = sites(cellfun(@(s) any(strcmp(site,s)),sites));
  if isempty(sites), sites = {''}; end
  ns = numel(sites); C = zeros(ns,3); A = nan(ns,2);
  many = manySites(F) && ns > 1;                   % (only the bars, without site labels, numbers, and ages)
  for si = 1:ns
    [a,s] = subjects(E(strcmp(site,sites{si})));
    C(si,:) = [sum(strcmp(s,'M')) sum(strcmp(s,'F')) sum(~ismember(s,{'M','F'}))];
    a = a(isfinite(a)); if ~isempty(a), A(si,1) = mean(a); end
    if numel(a) > 1, A(si,2) = std(a); end
  end
  sw = 1.5; [ax,tax] = cellAxes(fh,F,c,'Subjects per site',sw);
  xm = max([C(:,3)/2 + C(:,1); C(:,3)/2 + C(:,2); 1]) * 1.5; % (room for the numbers)
  % bars of the thickness of the age pyramid bars (see drawAge) with at
  % least 0.3 cm per site (for the age boxes), i.e. thinner bars only for 
  % many sites
  th  = c(4) - 1.45;                               % height of the axes in cm (see cellAxes)
  nb  = max(1,round(diff(F.ropts.agerange) / F.ropts.agestep)); 
  sp  = min(th / (ns + 0.2),max(th / nb / 0.85,0.3)); % space per site in cm
  bw  = min(0.85,th / nb / sp);                    % bar width (of the space per site)
  pad = max(0,(th / sp - ns - 0.2) / 2);           % (empty space above and below the sites)
  plot(ax,[0 0],[0.4 - pad ns + 0.6 + pad],'Color',F.fgl,'LineWidth',0.5); % (below the bars and numbers)
  pyramidBars(ax,F,1:ns+1,C,bw,~many,xm,c(3) - 1.15 - sw); % (empty row to stack one site)
  xt = 0:niceStep(xm/1.5/2):xm/1.5; xt = unique([-xt xt]);
  set(ax,'YDir','reverse','YLim',[0.4 - pad ns + 0.6 + pad],'YTick',1:ns,'XLim',[-xm xm],'TickLabelInterpreter','tex', ...
    'YTickLabel',strcat({'\bf '},regexprep(sites,'([_^{}\\])','\\$1')), ...   % (bold site codes)
    'XTick',xt,'XTickLabel',arrayfun(@(x) num2str(abs(x)),xt,'UniformOutput',false));
  xlabel(ax,'subjects','FontSize',F.fss,'Color',F.fg);
  if many
    set(ax,'YTick',[]); ylabel(ax,sprintf('%d sites',ns),'FontSize',F.fss,'Color',F.fg); 
    return
  end
  % mean age and SD of each site (boxes aligned with the bars)
  X  = [0.04 0.62; 0.72 0.72] / sw;               % boxes (left edge and width, normalized)
  bh = min(0.36,0.85 * th / (ns + 0.2 + 2*pad)) / th;
  kind = {'mean','sd'}; hdr = {'age','SD'};
  for ci = 1:2
    text(tax,X(ci,1) + X(ci,2)/2,1.0,hdr{ci},'FontSize',F.fss,'FontWeight','bold','Color',F.fg, ...
      'HorizontalAlignment','center','VerticalAlignment','bottom');
  end
  for si = 1:ns
    y = 1 - (si - 0.4 + pad) / (ns + 0.2 + 2*pad);
    for ci = 1:2
      if ~isfinite(A(si,ci)), continue; end
      patch(tax,X(ci,1) + [0 X(ci,2) X(ci,2) 0],y + bh/2 * [-1 -1 1 1],ageColor(F,A(si,ci),kind{ci}),'EdgeColor','none');
      fmt = '%0.1f'; if ci == 2, fmt = '\\pm%0.1f'; end
      text(tax,X(ci,1) + X(ci,2)/2,y,sprintf(fmt,A(si,ci)),'FontSize',F.fss,'Color',F.fg,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
    end
  end
end
% =========================================================================
function drawVolumes(fh,F,E,c)
%drawVolumes. Relative tissue volumes of the best T1w scan of each session
%  by age (marker by sex) with the robust age/sex model of the sample and
%  the volume statistics next to the diagram (z-scores of the session
%  table, see cat_io_dcm2bids_stats).
  [ax,tax] = cellAxes(fh,F,c,'Tissue volumes (best T1w scan)',1.7);
  age = numCol(E,'age'); sex = textCol(E,'sex');
  fn  = {'T1w_rGMV','T1w_rWMV','T1w_rCSFV'};
  Y   = 100 * [numCol(E,fn{1}) numCol(E,fn{2}) numCol(E,fn{3})];
  Z   = [numCol(E,[fn{1} '_z']) numCol(E,[fn{2} '_z']) numCol(E,[fn{3} '_z'])];
  ok  = isfinite(age) & any(isfinite(Y),2);
  if ~any(ok), noData(ax,F,{'no T1w segmentation','(preprocessing)'}); return; end
  hold(ax,'on');
  mk = {'s','o','^','^'}; xr = [min(age(ok)) max(age(ok))] + [-1 1];
  for ti = 1:3
    subjectLines(ax,F,E,age,Y(:,ti),ok & isfinite(Y(:,ti)),F.tiscol(ti,:));
    for gi = 1:4
      g = ok & isfinite(Y(:,ti)) & strcmp(sex,F.sexnam{gi});
      if gi == 4, g = ok & isfinite(Y(:,ti)) & ~ismember(sex,F.sexnam(1:3)); end
      points(ax,age(g),Y(g,ti),mk{gi},F.tiscol(ti,:));
    end
    o = ok & abs(Z(:,ti)) > 3; % outliers of the model
    plot(ax,age(o),Y(o,ti),'o','MarkerSize',5,'MarkerEdgeColor',F.red,'LineWidth',0.5);
    modelLine(ax,F,F.M(ti),xr,F.tiscol(ti,:));
  end
  set(ax,'XLim',F.ropts.agerange,'XTick',ageTicks(F),'YLim',[0 max(60,ceil(max(Y(:))/10)*10)]);
  gridOn(ax,F);
  xlabel(ax,'age in years','FontSize',F.fss,'Color',F.fg); ylabel(ax,'volume in % of TIV','FontSize',F.fss,'Color',F.fg);
  % statistics next to the diagram
  tiv = numCol(E,'T1w_TIV'); y = 0.98; dy = 0.095; ds = 0.062; % (lines and lines of a measure)
  for ti = 1:3
    v = Y(isfinite(Y(:,ti)),ti); M = F.M(ti);
    text(tax,0.08,y,sprintf('%s %0.1f\\pm%0.1f%%',F.tisnam{ti},mean(v),std(v)),'FontSize',F.fss, ...
      'Color',F.tiscol(ti,:),'FontWeight','bold');
    if startsWith(M.method,'robust linear')
      y = y - ds; text(tax,0.08,y,sprintf(' %+0.2f%%/10 y',10*M.b(2)),'FontSize',F.fss,'Color',F.fg);
    end
    y = y - dy;
  end
  t = tiv(isfinite(tiv));
  if ~isempty(t), text(tax,0.08,y,sprintf('TIV %0.0f\\pm%0.0f ml',mean(t),std(t)),'FontSize',F.fss,'Color',F.fg); y = y - dy; end
  text(tax,0.08,y,sprintf('n = %d',sum(ok)),'FontSize',F.fss,'Color',F.fg); y = y - dy;
  sexLegend(tax,F,y,sum(any(abs(Z(ok,:)) > 3,2)));
end
% =========================================================================
function drawMarker(fh,F,E,c,fn,M,ttl,ylab,nd)
%drawMarker. Brain marker fn of the sessions E (e.g. the white matter FA
%  of the best DWI scan) by age (points by site, marker by sex) with the
%  robust age/sex model M of the sample, the outliers of the model (|z| >
%  3), and the statistics next to the diagram (nd decimals).
  [ax,tax] = cellAxes(fh,F,c,ttl,1.7);
  age = numCol(E,'age'); sex = textCol(E,'sex'); site = textCol(E,'site');
  v   = numCol(E,fn); z = numCol(E,[fn '_z']);
  ok  = isfinite(age) & isfinite(v);
  if ~any(ok), noData(ax,F,{sprintf('no %s',ylab),'(preprocessing)'}); return; end
  hold(ax,'on');
  mk = {'s','o','^','^'}; xr = [min(age(ok)) max(age(ok))] + [-1 1];
  subjectLines(ax,F,E,age,v,ok,site);
  [us,~,ui] = unique(site);
  for si = 1:numel(us)
    col = siteColor(F,us{si});
    for gi = 1:4
      g = ok & ui == si & strcmp(sex,F.sexnam{gi});
      if gi == 4, g = ok & ui == si & ~ismember(sex,F.sexnam(1:3)); end
      points(ax,age(g),v(g),mk{gi},col);
    end
  end
  o = ok & abs(z) > 3;
  plot(ax,age(o),v(o),'o','MarkerSize',5.5,'MarkerEdgeColor',F.red,'LineWidth',0.5);
  modelLine(ax,F,M,xr,F.fgl);
  yr = [min(v(ok)) max(v(ok))]; yr = yr + [-0.15 0.15] * max(diff(yr),abs(mean(yr))*0.05 + eps);
  set(ax,'XLim',F.ropts.agerange,'XTick',ageTicks(F),'YLim',yr);
  gridOn(ax,F);
  xlabel(ax,'age in years','FontSize',F.fss,'Color',F.fg); ylabel(ax,ylab,'FontSize',F.fss,'Color',F.fg);
  % statistics next to the diagram
  y = 0.98; dy = 0.095; ds = 0.062; fmt = sprintf('%%0.%df',nd); % (lines and lines of the measure)
  text(tax,0.08,y,sprintf(['mean ' fmt '\\pm' fmt],mean(v(ok)),std(v(ok))),'FontSize',F.fss,'Color',F.fg,'FontWeight','bold');
  if startsWith(M.method,'robust linear')
    y = y - ds; text(tax,0.08,y,sprintf([' %+' fmt(2:end) '/10 y'],10*M.b(2)),'FontSize',F.fss,'Color',F.fg);
    if contains(M.method,'sex'), y = y - ds; text(tax,0.08,y,sprintf([' %+' fmt(2:end) ' male'],M.b(3)),'FontSize',F.fss,'Color',F.fg); end
  end
  y = y - dy; text(tax,0.08,y,sprintf('n = %d',sum(ok)),'FontSize',F.fss,'Color',F.fg); y = y - dy;
  sexLegend(tax,F,y,sum(o));
end
% =========================================================================
function t = ageTicks(F)
%ageTicks. Ticks of the age range of the age pyramid and the brain markers
%  (e.g. 20 years for 0 to 100 years).
  ar = F.ropts.agerange; t = ar(1):niceStep(diff(ar)/5):ar(2);
end
% =========================================================================
function [ar,st] = ageRange(age,ropts)
%ageRange. Age range ar of the diagrams and bin size st of the age pyramid
%  in years (opts.reportfig.agerange and agestep), by default the range of
%  the ages of the sample rounded to 10 years for a range of at least 20
%  years, to 5 years for 10 to 20 years, and to 1 year below (0 to 100 
%  years without ages), and a bin size of 1, 2, 5, or 10 years for at 
%  most 20 bins.
  ar = []; st = [];
  if isfield(ropts,'agerange'), ar = double(ropts.agerange); end
  if isfield(ropts,'agestep'),  st = double(ropts.agestep);  end
  if numel(ar) ~= 2 || ~(ar(2) > ar(1))
    a = age(isfinite(age)); if isempty(a), a = [0 100]; end
    r = max(a) - min(a);
    if r >= 20, u = 10; elseif r >= 10, u = 5; else, u = 1; end
    ar = floor(min(a) / u) * u; ar(2) = max(ar(1) + u,ceil(max(a) / u) * u);
  end
  if isempty(st) || ~(st > 0)
    m = [1 2 5 10]; st = m(find(diff(ar) ./ m <= 20,1));
    if isempty(st), st = niceStep(diff(ar) / 20); end
  end
end
% =========================================================================
function gridOn(ax,F)
%gridOn. Light grid of the diagram ax below the data (brain markers).
  set(ax,'XGrid','on','YGrid','on','GridColor',F.fgl,'GridAlpha',0.25,'Layer','bottom');
end
% =========================================================================
function subjectLines(ax,F,E,x,y,ok,col)
%subjectLines. Thin lines between the sessions (points x,y of the valid 
%  sessions ok) of each subject in the order of the sessions (drawn before
%  the points, i.e. behind them) with the color col (RGB, or the sites of 
%  the sessions for the color of the site).
  sub = textCol(E,'participant_id'); ses = textCol(E,'session_id'); idx = find(ok(:));
  if numel(idx) < 2, return; end
  [us,~,ui] = unique(sub(idx));
  for si = 1:numel(us)
    k = idx(ui == si); if numel(k) < 2, continue; end
    [~,so] = sort(ses(k)); k = k(so);
    if iscell(col), c = siteColor(F,col{k(1)}); else, c = col; end
    plot(ax,x(k),y(k),'-','Color',F.mix(c,0.45),'LineWidth',0.3);
  end
end
% =========================================================================
function points(ax,x,y,mk,col)
%points. Data points of the brain markers (marker mk and color col) with
%  30% opacity (overlapping points are darker).
  if isempty(x), return; end
  scatter(ax,x,y,12,col,mk,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeColor','none');
end
% =========================================================================
function modelLine(ax,F,M,xr,col)
%modelLine. Line of the robust age model M at the mean sex (median for the
%  robust z-score) in the age range xr.
  if M.n == 0, return; end
  xa = [0 0]; if isfinite(M.mage), xa = xr - M.mage; end
  plot(ax,xr,M.b(1) + M.b(2) * xa,'-','Color',col,'LineWidth',1);
end
% =========================================================================
function sexLegend(tax,F,y,nout)
%sexLegend. Legend of the markers of the sex and of the outliers of the
%  model (red rings) with their number nout.
  plot(tax,0.12,y,'s','MarkerSize',2.8,'MarkerFaceColor',F.fgl,'MarkerEdgeColor','none');
  text(tax,0.22,y,'male','FontSize',F.fss,'Color',F.fg);
  plot(tax,0.62,y,'o','MarkerSize',2.8,'MarkerFaceColor',F.fgl,'MarkerEdgeColor','none');
  text(tax,0.72,y,'female','FontSize',F.fss,'Color',F.fg);
  plot(tax,0.12,y - 0.095,'o','MarkerSize',5,'MarkerEdgeColor',F.red,'LineWidth',0.5);
  text(tax,0.22,y - 0.095,sprintf('|z| > 3: %d',nout),'FontSize',F.fss,'Color',F.red);
end
% =========================================================================
function drawHomogeneity(fh,F,E,c,fq,fs,thr,ttl,typ)
%drawHomogeneity. Sample homogeneity fq (quartic mean Z-score QMZ, e.g. of
%  the mwc1 maps of the best T1w scan) by the quality rating fs of the
%  scans (points by site) with the threshold thr of atypical maps (mean + 2
%  SD of the whole sample, as cat_stat_homogeneity), the SQR borders of 
%  light and severe artifacts (yellow and red, see sqrLimits), and the 
%  labels (without sub- and ses-) of the (at most 3) sessions above the 
%  threshold and of the (at most 3) worst sessions with severe artifacts,
%  and thin lines between the sessions of a subject (see subjectLines).
  [ax,tax] = cellAxes(fh,F,c,ttl,1.7);
  qmz = numCol(E,fq); sqr = numCol(E,fs); site = textCol(E,'site');
  ok  = isfinite(qmz) & isfinite(sqr);
  if ~any(ok), noData(ax,F,{'no sample homogeneity',sprintf('(at least 3 %s scans with preprocessing)',typ)}); return; end
  hold(ax,'on');
  yr = [min(qmz(ok)) max([qmz(ok); thr])]; yr = yr + [-0.1 0.2] * max(diff(yr),0.05);
  for mi = 1:6 % mark colors of the quality
    patch(ax,mi + [-0.5 0.5 0.5 -0.5],yr([1 1 2 2]),F.mix(F.col2mark(mi),F.markbg),'EdgeColor','none');
  end
  if isfinite(thr), plot(ax,F.marks,[thr thr],'--','Color',F.red,'LineWidth',0.6); end
  L = F.sqrlim(find(strcmp({'T1w','DWI','BOLD'},typ),1),:); % (typical quality, light and severe artifacts)
  if all(isfinite(L))
    plot(ax,L([2 2]),yr,':','Color',F.yellow,'LineWidth',0.8);
    plot(ax,L([3 3]),yr,'--','Color',F.red,'LineWidth',0.6);
  end
  x = min(F.marks(2),max(F.marks(1),sqr));
  subjectLines(ax,F,E,x,qmz,ok,site);
  [us,~,ui] = unique(site);
  for si = 1:numel(us)
    g = ok & ui == si; col = siteColor(F,us{si});
    if any(g), plot(ax,x(g),qmz(g),'o','MarkerSize',3,'MarkerFaceColor',col,'MarkerEdgeColor',F.mix(col,0.5),'LineWidth',0.2); end
  end
  set(ax,'XLim',F.marks,'YLim',yr,'XTick',1:6,'Layer','top'); % (before the labels for their extent)
  hi  = find(ok & qmz > thr);
  [~,so] = sort(qmz(hi),'descend'); hl = hi(so(1:min(end,3)));
  hs  = find(ok & sqr > L(3));                      % (severe artifacts by the SQR)
  [~,so] = sort(sqr(hs),'descend'); hs = hs(so(1:min(end,3)));
  hl  = union(hl,hs); [~,so] = sort(qmz(hl),'descend'); hl = hl(so);
  sub = regexprep(textCol(E,'participant_id'),'^sub-',''); ses = regexprep(textCol(E,'session_id'),'^ses-','');
  yi = inf;
  for i = 1:numel(hl) % right of the point, or left if it does not fit (below the previous label)
    xi = x(hl(i)); yi = min(qmz(hl(i)),yi - 0.065*diff(yr));
    h  = text(ax,xi + 0.15,yi,shortLabel(strtrim([sub{hl(i)} ' ' ses{hl(i)}]),31),'FontSize',F.fss*0.85, ...
      'Color',F.fg,'Interpreter','none','HorizontalAlignment','left','Clipping','on');
    if sum(h.Extent([1 3])) > F.marks(2)
      set(h,'Position',[xi - 0.15 yi 0],'HorizontalAlignment','right');
      if h.Extent(1) < F.marks(1), set(h,'Position',[F.marks(1) + 0.05 yi 0],'HorizontalAlignment','left'); end
    end
  end
  xlabel(ax,sprintf('%s quality SQR (mark)',typ),'FontSize',F.fss,'Color',F.fg); ylabel(ax,'QMZ','FontSize',F.fss,'Color',F.fg);
  % statistics next to the diagram
  y = 0.98; dy = 0.095;
  text(tax,0.08,y,sprintf('QMZ %0.3f\\pm%0.3f',mean(qmz(ok)),std(qmz(ok))),'FontSize',F.fss,'Color',F.fg,'FontWeight','bold'); y = y - dy;
  text(tax,0.08,y,sprintf('n = %d',sum(isfinite(qmz))),'FontSize',F.fss,'Color',F.fg); y = y - dy;
  if isfinite(thr)
    text(tax,0.08,y,sprintf('mean + 2 SD %0.3f',thr),'FontSize',F.fss,'Color',F.red); y = y - dy;
    text(tax,0.08,y,sprintf('above: %d',numel(hi)),'FontSize',F.fss,'Color',F.fg); y = y - dy;
  end
  if all(isfinite(L))
    text(tax,0.08,y,sprintf('SQR typical %0.2f',L(1)),'FontSize',F.fss,'Color',F.fg); y = y - dy;
    text(tax,0.08,y,sprintf('light > %0.2f: %d',L(2),sum(ok & sqr > L(2) & sqr <= L(3))),'FontSize',F.fss,'Color',F.yellow); y = y - dy;
    text(tax,0.08,y,sprintf('severe > %0.2f: %d',L(3),sum(ok & sqr > L(3))),'FontSize',F.fss,'Color',F.red); y = y - dy;
  end
  if sum(isfinite(qmz)) >= 10, str = {'corrected for','age and sex'}; else, str = {'not corrected','(< 10 sessions)'}; end
  text(tax,0.08,y,str{1},'FontSize',F.fss,'Color',F.fgl); text(tax,0.08,y - dy,str{2},'FontSize',F.fss,'Color',F.fgl);
end
% =========================================================================
function drawCompleteness(fh,F,S,E,c)
%drawCompleteness. Protocol status of the best scan of each typical scan
%  type and session (conform, semiconform, nonconform, without protocol
%  test, or missing scan) in percent of the sessions.
  ax  = cellAxes(fh,F,[c(1) c(2)+0.6 c(3) c(4)-0.6],'Completeness and protocol status',0.75);
  sk  = strcat(textCol(S,'participant_id'),'|',textCol(S,'session_id'));
  ek  = strcat(textCol(E,'participant_id'),'|',textCol(E,'session_id'));
  st  = textCol(S,'protocol_status'); ty = textCol(S,'type');
  rnk = 4 * ones(numel(S),1);                    % without protocol test
  for ri = 1:3, rnk(strcmp(st,F.statnam{ri})) = ri; end
  nt  = numel(F.types); P = zeros(nt,5); nses = numel(E);
  for ti = 1:nt
    for ei = 1:nses
      r = rnk(strcmp(sk,ek{ei}) & strcmp(ty,F.types{ti}));
      if isempty(r), P(ti,5) = P(ti,5) + 1; else, P(ti,min(r)) = P(ti,min(r)) + 1; end
    end
  end
  P  = 100 * P / max(1,nses);
  hb = barh(ax,1:nt+1,[P; zeros(1,5)],0.7,'stacked','EdgeColor',F.bg,'LineWidth',0.3); % (empty row to stack one type)
  for ci = 1:5, set(hb(ci),'FaceColor',F.statcol(ci,:)); end
  hold(ax,'on');
  for ti = 1:nt
    x0 = 0;
    for ci = 1:5
      if P(ti,ci) >= 12
        text(ax,x0 + P(ti,ci)/2,ti,sprintf('%0.0f%%',P(ti,ci)),'FontSize',F.fss*0.9, ...
          'HorizontalAlignment','center','Color',F.fg);
      end
      x0 = x0 + P(ti,ci);
    end
    text(ax,103,ti,sprintf('%0.0f/%d',nses*(1 - P(ti,5)/100),nses),'FontSize',F.fss,'Color',F.fg);
  end
  set(ax,'YDir','reverse','YLim',[0.4 nt+0.6],'YTick',1:nt,'YTickLabel',F.tnam,'XLim',[0 100], ...
    'TickLabelInterpreter','none');
  xlabel(ax,'sessions in %','FontSize',F.fss,'Color',F.fg);
  % legend below the axes (two rows)
  lax = textAxes(fh,F,[c(1)+1.0 c(2) c(3)-1.2 0.6]);
  xl  = [0 0.33 0.66 0 0.33]; yl = [0.75 0.75 0.75 0.25 0.25];
  for ci = 1:5
    patch(lax,xl(ci) + [0 0.04 0.04 0],yl(ci) + [-0.18 -0.18 0.18 0.18],F.statcol(ci,:),'EdgeColor',F.fgl,'LineWidth',0.2);
    text(lax,xl(ci) + 0.06,yl(ci),F.statnam{ci},'FontSize',F.fss,'Color',F.fg,'VerticalAlignment','middle');
  end
end
% =========================================================================
function drawRatings(fh,F,E,fields,names,c,ttl,colorby)
%drawRatings. Box plots of the ratings of the fields of the records E
%  (sessions, or the scans of one site) (marks 0.5 to 6.5 with the colors
%  of the marks, worse ratings are clipped) with the data points (colored
%  by site, or by datatype for colorby 'datatype' with a legend), the
%  median, and the number of records of each rating.
  if ~exist('colorby','var'), colorby = 'site'; end
  ax = cellAxes(fh,F,c,ttl,0);
  n  = numel(fields); data = cell(1,n); site = textCol(E,colorby);
  for fi = 1:n, data{fi} = numCol(E,fields{fi}); end
  nv = cellfun(@(d) sum(isfinite(d)),data);
  if ~any(nv), noData(ax,F,{'no ratings'}); return; end
  datac = cellfun(@(d) min(F.marks(2),max(F.marks(1),d(isfinite(d)))),data,'UniformOutput',false);
  datac(nv == 0) = {nan};
  set(0,'CurrentFigure',fh); set(fh,'CurrentAxes',ax);
  cat_plot_boxplot(datac,struct('names',{names},'ylim',F.marks,'ygrid',0,'fontsize',F.fss, ...
    'groupcolor',repmat([.55 .55 .62],n,1),'showdata',0,'outliers',0,'median',2,'trans',0.4,'sat',0.4));
  hold(ax,'on');
  % mark colors in the background (as the QC ratings)
  hp = gobjects(1,6);
  for mi = 1:6
    hp(mi) = patch(ax,[0.4 n+1 n+1 0.4],mi + [-0.5 -0.5 0.5 0.5],F.mix(F.col2mark(mi),F.markbg),'EdgeColor','none');
  end
  uistack(hp,'bottom');
  % data points (by site, behind the boxes, clipped ratings as triangles),
  % median, and number of each rating
  hd = gobjects(0);
  for fi = find(nv > 0)
    d = data{fi}; ok = isfinite(d); d = d(ok); s = site(ok);
    [~,so] = sort(d); x = zeros(size(d));
    x(so) = fi + 0.32 * (mod((0:numel(d)-1)' * 0.618,1) - 0.5);
    dc = min(F.marks(2),max(F.marks(1),d));
    mk = repmat('o',size(d)); mk(d > F.marks(2)) = '^'; mk(d < F.marks(1)) = 'v';
    [us,~,ui] = unique(s);
    for si = 1:numel(us)
      if strcmp(colorby,'datatype'), col = dtColor(F,us{si}); else, col = siteColor(F,us{si}); end
      for m = 'o^v'
        g = ui == si & mk == m;
        if any(g)
          hd(end+1) = plot(ax,x(g),dc(g),m,'MarkerSize',2.5,'MarkerFaceColor',col, ...
            'MarkerEdgeColor',F.mix(col,0.6),'LineWidth',0.2); %#ok<AGROW>
        end
      end
    end
    md = median(d);
    text(ax,fi + 0.40,min(F.marks(2)-0.2,max(F.marks(1)+0.2,md)),sprintf('%0.2f',md),'FontSize',F.fss, ...
      'Color',F.fg,'FontWeight','bold');
  end
  if ~isempty(hd), uistack(hd,'bottom'); end
  uistack(hp,'bottom');
  set(ax,'XLim',[0.4 n+1],'YLim',F.marks,'YTick',1:6,'YDir','normal','Layer','top','Color',F.bg, ...
    'XColor',F.fgl,'YColor',F.fgl,'FontSize',F.fss, ...
    'XTick',1:n,'TickLabelInterpreter','tex','XTickLabel',cellfun(@(nm,k) ...
    sprintf('%s\\newline{\\fontsize{%g}n=%d}',strrep(nm,'_','\_'),F.fss*0.85,k),names,num2cell(nv),'UniformOutput',false));
  ylabel(ax,'mark','FontSize',F.fss,'Color',F.fg);
  if strcmp(colorby,'datatype') % legend of the datatypes (top left)
    dt = unique(site(~cellfun('isempty',site))); [~,o] = ismember(dt,F.dtnam); o(o==0) = numel(F.dtnam) + 1; [~,o] = sort(o); dt = dt(o);
    str = strjoin(cellfun(@(d) sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s',dtColor(F,d),d),dt(:)','UniformOutput',false),'  ');
    text(ax,0.55,F.marks(2) - 0.15,str,'FontSize',F.fss,'Interpreter','tex','VerticalAlignment','top');
  end
end
% =========================================================================
function c = dtColor(F,dt)
%dtColor. Color of the datatype (scan ratings of one site, see drawRatings).
  di = find(strcmp(F.dtnam,dt),1); if isempty(di), di = numel(F.dtnam) + 1; end
  c  = F.dtcol(di,:);
end
% =========================================================================
function S = scanPCP(S)
%scanPCP. Protocol status of each scan as mark of the protocol conformity
%  (conform 1, semiconform 3.5, nonconform 6, NaN without protocol test),
%  e.g. for the scan ratings of one site.
  if isempty(S), return; end
  st = textCol(S,'protocol_status'); m = nan(numel(S),1);
  m(strcmp(st,'conform')) = 1; m(strcmp(st,'semiconform')) = 3.5; m(strcmp(st,'nonconform')) = 6;
  for si = 1:numel(S), S(si).PCP = m(si); end
end
% =========================================================================
function noData(ax,F,str)
%noData. Message in the empty axes ax (e.g. without preprocessing).
  axis(ax,'off');
  text(ax,0.5,0.5,str,'FontSize',F.fss,'HorizontalAlignment','center','Color',F.fgl);
end
% =========================================================================
function showSPM(fh,F)
%showSPM. Copy of the page in the SPM Graphics window in interactive
%  sessions (fonts scaled to the window height).
  if ~usejava('desktop') || ~feature('ShowFigureWindows'), return; end
  try
    fg = spm_figure('GetWin','Graphics'); spm_figure('Clear',fg); set(fg,'Color',F.bg);
    copyobj(findobj(allchild(fh),'flat','Type','axes'),fg); % all diagrams and texts are axes
    u0 = get(fg,'Units'); set(fg,'Units','centimeters'); p = get(fg,'Position'); set(fg,'Units',u0);
    sc = p(4) / F.ph;
    h  = findall(fg,'-property','FontSize');
    for hi = 1:numel(h), set(h(hi),'FontSize',get(h(hi),'FontSize') * sc); end
  catch e
    cat_io_cprintf('warn','  Report figure not shown in the SPM Graphics window: %s\n',e.message);
  end
end
% =========================================================================
function ax = newAxes(fh,F,pos)
%newAxes. Axes of a diagram at the position pos ([x y w h] in cm) in the
%  style of the design, also for high-level plots (NextPlot add, i.e. 
%  without reset of the axes, e.g. by bar or pie).
  ax = axes('Parent',fh,'Units','normalized','Position',[pos(1)/F.pw pos(2)/F.ph pos(3)/F.pw pos(4)/F.ph], ...
    'FontSize',F.fss,'Box','on','Layer','top','Color',F.bg,'XColor',F.fgl,'YColor',F.fgl, ...
    'TickLength',[0.008 0.008],'LineWidth',0.5,'NextPlot','add');
  ax.Toolbar = []; disableDefaultInteractivity(ax); % no interactive elements (faster)
end
% =========================================================================
function ax = textAxes(fh,F,pos)
%textAxes. Invisible axes (0 to 1) for text at the position pos (in cm).
  ax = axes('Parent',fh,'Units','normalized','Position',[pos(1)/F.pw pos(2)/F.ph pos(3)/F.pw pos(4)/F.ph], ...
    'XLim',[0 1],'YLim',[0 1],'Visible','off','HitTest','off');
  ax.Toolbar = []; disableDefaultInteractivity(ax);
  hold(ax,'on');
end
% =========================================================================
function st = niceStep(v)
%niceStep. Step of 1, 2, or 5 times a power of 10 (at least 1) of the value v.
  e  = 10^floor(log10(max(v,1))); f = max(v,1) / e;
  m  = [1 2 5 10]; st = e * m(find(m >= f - eps,1));
end
% =========================================================================
function c = siteColor(F,site)
%siteColor. Color of the site (gray without site).
  si = find(strcmp(F.sites,site),1);
  if isempty(si), c = [.5 .5 .5]; else, c = F.sitecol(si,:); end
end
% =========================================================================
function s = shortLabel(s,n)
%shortLabel. Label of at most n characters (shortened in the middle).
  if numel(s) > n, s = [s(1:ceil(n/2)-1) '...' s(end-floor(n/2)+3:end)]; end
end
% =========================================================================
function [types,names] = scanTypes(opts)
%scanTypes. Typical scan types (datatype/suffix, see opts.scantypes) and
%  their column names (as in cat_io_dcm2bids_report).
  types = cellstr(opts.scantypes);
  names = regexprep(types,'^[^/]*/','');
  dup   = cellfun(@(n) sum(strcmp(names,n)) > 1, names);
  names(dup) = strrep(types(dup),'/','_');
  names = regexprep(names,'[^a-zA-Z0-9_]','_');
end
% =========================================================================
function R = readTable(P)
%readTable. Report table (csv or tsv, see writeTable of
%  cat_io_dcm2bids_helper) as structure array (one element per row) with
%  numeric columns (n/a as NaN) and text columns (n/a as '').
  R = struct([]);
  if ~exist(P,'file'), return; end
  if strcmpi(spm_file(P,'ext'),'tsv'), dl = sprintf('\t'); else, dl = ','; end
  C = readcell(P,'FileType','text','Delimiter',dl,'DatetimeType','text','DurationType','text');
  if size(C,1) < 2, return; end
  hdr = matlab.lang.makeValidName(cellfun(@char,C(1,:),'UniformOutput',false));
  V   = C(2:end,:); n = size(V,1);
  isna = cellfun(@(x) isa(x,'missing') || (ischar(x) && strcmp(x,'n/a')),V);
  isnm = cellfun(@(x) isnumeric(x) || islogical(x),V);
  vals = cell(1,numel(hdr));
  for ci = 1:numel(hdr)
    if all(isna(:,ci) | isnm(:,ci))
      x = nan(n,1); x(isnm(:,ci)) = cellfun(@double,V(isnm(:,ci),ci));
      vals{ci} = num2cell(x);
    else
      v = V(:,ci); v(isna(:,ci)) = {''};
      v(isnm(:,ci)) = cellfun(@(x) sprintf('%.15g',x),v(isnm(:,ci)),'UniformOutput',false);
      vals{ci} = v;
    end
  end
  R = cell2struct([vals{:}],hdr,2);
end
% =========================================================================
function x = numCol(R,fn)
%numCol. Numeric column fn of the records R as column vector (NaN if not
%  available or not numeric).
  x = nan(numel(R),1);
  if isempty(R) || ~isfield(R,fn), return; end
  for i = 1:numel(R)
    v = R(i).(fn);
    if (isnumeric(v) || islogical(v)) && ~isempty(v), x(i) = double(v(1)); end
  end
end
% =========================================================================
function c = textCol(R,fn)
%textCol. Text column fn of the records R as column cell ('' if not
%  available), also of numeric values (e.g. a site named by the serial
%  number of the scanner).
  c = repmat({''},numel(R),1);
  if isempty(R) || ~isfield(R,fn), return; end
  for i = 1:numel(R)
    v = R(i).(fn);
    if ischar(v), c{i} = v;
    elseif isnumeric(v) && isscalar(v) && isfinite(v), c{i} = sprintf('%.15g',v);
    end
  end
end
% =========================================================================
function v = mostFrequent(c)
%mostFrequent. Most frequent (non-empty) text of the cell c.
  c = c(~cellfun('isempty',c));
  if isempty(c), v = ''; return; end
  [u,~,ui] = unique(c);
  [~,mi] = max(accumarray(ui(:),1));
  v = u{mi};
end
% =========================================================================
function v = valueCounts(c)
%valueCounts. Different (non-empty) texts of the cell c with their counts
%  (most frequent first).
  c = c(~cellfun('isempty',c));
  if isempty(c), v = 'n/a'; return; end
  [u,~,ui] = unique(c);
  n = accumarray(ui(:),1);
  [n,so] = sort(n,'descend'); u = u(so);
  v = strjoin(arrayfun(@(i) sprintf('%s (%d)',u{i},n(i)),1:numel(u),'UniformOutput',false),', ');
end
