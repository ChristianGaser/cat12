function varargout = cat_io_dcm2bids_participantreport(action,varargin)
%cat_io_dcm2bids_participantreport. Participant report of cat_io_dcm2bids.
%  Illustrated A4 report for the scanned person (non-expert readers, not
%  the researcher or physician) of each session with a segmented T1w,
%  saved only in the private directory of the BIDS directory
%    <outdir>/<subdir>/<BIDSdir>-private/<BIDSsubdir>/participantreports/
%      <sub>_<ses>_participantreport_<lang>.pdf
%  (opts.participantreport.run, English and German by default) with
%    header   title, study logo, and identity of the participant (ID, name,
%             birth date, sex, weight, scan date, and age of the original
%             DICOM header in the database, i.e. to confirm the person)
%    intro    brain (GM, WM, CSF), MRI scanner (with illustration), and
%             images and measures, with the abbreviations
%    slices   bias-corrected T1w in the head orientation (rigid MNI
%             registration of the session) through the right hippocampus
%             (sagittal, coronal, axial) with the DEC map of the DWI in the
%             WM (in the space of the T1w, see cat_io_dcm2bids_pp_dwi), the
%             cortical lobes lightly in the GM (lotter2025, discrete colors
%             of the design), the hippocampus, and the removed face (side
%             of the plane of spm_deface, hatched)
%    surface  pial surface of CAT (derivatives/CAT<version>, otherwise the
%             quick surfaces of cat_io_dcm2bids_pp_surf, otherwise the GM
%             boundary of the segmentation without thickness) from the
%             right-front-top with the cortical thickness (dark gray-pink
%             2 mm to bright 4 mm), the sulci darkened by the sulcal depth,
%             the cerebellum and brainstem (atlas), and the transparent head
%             (view-dependent transparency of cat_surf_render2, i.e. with
%             the removed face)
%    values   GMV, WMV, and HCV of the segmentation of the session and CT
%             of the surface in normative band plots by age and sex (5-95%
%             and 25-75%) of the reference models opts.participantreport.
%             norms (DCM2BIDS/participantreport/norms.json, placeholders),
%             and the overall quality rating of the T1w (SQR of the QC,
%             e.g. reduced by motion, as grade and percentage)
%    notes    image quality (QC), individuality of the brain (diagnosis
%             only by trained physicians), and privacy (defacing)
%    footer   disclaimer, study, organization, website, institute, and the
%             CAT logo
%  The layout opts.participantreport.variant ('A' magazine rows, 'B'
%  sidebar, 'C' infographic, 'D' landscape) is drawn as a language-free
%  base (images, numbers, abbreviations) with the texts of the language
%  files DCM2BIDS/participantreport/<lang>.json (opts.participantreport.
%  lang, one PDF each, a <lang>.json of the study directory replaces it,
%  missing fields are English), where scripts without reliable PDF fonts
%  (Arabic, Chinese, Japanese: Render "image") are drawn as Java2D text
%  images. The design is the study design of the report figure (colormap.
%  json of the study directory, see cat_io_dcm2bids_reportfig), the
%  organization and website are the Name and the first link of the study
%  JSON, and the institute is given by the DICOM header.
%
%  varargout = cat_io_dcm2bids_participantreport(action,varargin)
%
%  Actions (see the help of the local functions):
%    Ppdf = cat_io_dcm2bids_participantreport('writeReports',Poutdir,job[,Psub])
%    Ppdf = cat_io_dcm2bids_participantreport('writeReport',R,Pout,popts,G)
%    Ppdf = cat_io_dcm2bids_participantreport('demo',Pout[,variants,langs,Pstudy])
%    L    = cat_io_dcm2bids_participantreport('language',lang[,Pstudydir])
%    R    = cat_io_dcm2bids_participantreport('sessionData',Poutdir,job[,Psub])
%    I    = cat_io_dcm2bids_participantreport('sessionImages',R,G)
%    G    = cat_io_dcm2bids_participantreport('designOf',job)
%
%  See also cat_io_dcm2bids_reportfig, cat_io_dcm2bids_render, cat_surf_render2.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  switch action
    case {'writeReports','writeReport','demo','language','sessionData','designOf'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    case 'sessionImages'
      varargout{1} = sessionImages(varargin{:},atlasIds);
    otherwise
      error('cat_io_dcm2bids_participantreport:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Ppdf = writeReports(Poutdir,job,Psub)
%writeReports. Participant reports of the sessions of the BIDS directory
%  of the job (see sessionData) in its private directory, returns the PDFs.
%  Psub limits the reports to these subject directories of the database 
%  (e.g. the subject package of a parallel process, see processPackage in
%  cat_io_dcm2bids_parallel, without the title line).
  if ~exist('Psub','var'), Psub = {}; end
  if isempty(Psub), fprintf('\nDCM2BIDS - Participant reports.\n'); end
  Ppdf  = {};
  popts = reportOptions(job);
  R     = sessionData(Poutdir,job,Psub);
  if isempty(R)
    if isempty(Psub), fprintf('  No sessions with segmented T1w.\n'); end
    return
  end
  Pout  = fullfile(Poutdir,[job.BIDSdir '-private'],job.BIDSsubdir,'participantreports');
  G     = designOf(job);
  for ri = 1:numel(R)
    try
      Ppdf = [Ppdf writeReport(R(ri),Pout,popts,G)]; %#ok<AGROW>
    catch e
      cat_io_cprintf('err','  Participant report of %s failed: %s (%s, line %d)\n',R(ri).name,e.message,e.stack(1).name,e.stack(1).line);
    end
  end
end
% =========================================================================
function Ppdf = writeReport(R,Pout,popts,G)
%writeReport. Report PDFs of the session R (see sessionData) for each
%  variant and language of popts (see reportOptions) with the design G
%  (see designOf) in the directory Pout as <name>_participantreport_<lang>
%  .pdf (with _<variant> for several variants), returns the PDF files.
  if ~exist(Pout,'dir'), mkdir(Pout); end
  A = atlasIds;
  I = sessionImages(R,G,A);   % (language-free images of all pages)
  M = measureValues(R,I,normModels(popts.norms));
  I.qc = qcRating(R);
  variants = cellstr(popts.variant); Ppdf = {};
  for vi = 1:numel(variants)
    for li = 1:numel(popts.lang)
      L = language(popts.lang{li},popts.studydir);
      F = newPage(G,L,R,M,I,variants{vi});
      feval(['variant' upper(variants{vi})],F);
      P = fullfile(Pout,sprintf('%s_participantreport_%s',R.name,L.Language));
      if numel(variants) > 1, P = sprintf('%s_%s',P,upper(variants{vi})); end
      cat_io_dcm2bids_helper('printPDF',F.fh,[P '.pdf']); % (A4 page, see printPDF)
      close(F.fh); Ppdf{end+1} = [P '.pdf']; %#ok<AGROW>
      fprintf('  Participant report: %s\n',Ppdf{end});
    end
  end
end
% =========================================================================
function Ppdf = demo(Pout,variants,langs,Pstudy)
%demo. Participant reports of the CAT templates with placeholder identity
%  and values (see demoSession) for the selection of the layout and the
%  test of the languages: variants ('A' to 'D', default all), languages
%  (default {'en'}), and study JSON of the design (default DZPG), in Pout.
  Pcat = fileparts(mfilename('fullpath'));
  if nargin < 2 || isempty(variants), variants = {'A','B','C','D'}; end
  if nargin < 3 || isempty(langs), langs = {'en'}; end
  if nargin < 4 || isempty(Pstudy), Pstudy = fullfile(Pcat,'DCM2BIDS','studies','DZPG','study.json'); end
  job   = struct('subdir','DEMO','dataset',struct('study',Pstudy),'opts',struct());
  popts = reportOptions(job); popts.variant = cellstr(variants); popts.lang = cellstr(langs);
  Ppdf  = writeReport(demoSession,Pout,popts,designOf(job));
end
% =========================================================================
function popts = reportOptions(job)
%reportOptions. Options of the report (job.opts.participantreport, see
%  cat_io_dcm2bids_defaults): variant, languages, reference models (''
%  for DCM2BIDS/participantreport/norms.json), and the study directory
%  (directory of the study JSON, for the language files).
  popts = struct('variant','B','lang',{{'en','de'}},'norms','','studydir','');
  if isfield(job,'opts') && isfield(job.opts,'participantreport')
    for fn = fieldnames(job.opts.participantreport)', popts.(fn{1}) = job.opts.participantreport.(fn{1}); end
  end
  if ischar(popts.lang), popts.lang = regexp(strtrim(popts.lang),'[\s,;]+','split'); end % (GUI: "en de")
  popts.lang = cellstr(popts.lang);
  if isfield(job,'dataset') && isfield(job.dataset,'study') && ~isempty(job.dataset.study)
    popts.studydir = fileparts(char(job.dataset.study));
  end
end
% =========================================================================
function G = designOf(job)
%designOf. Design of the report: colors of the study design of the report
%  figure (see reportDesigns in cat_io_dcm2bids_reportfig, the bright CAT
%  colors for dark designs, as the report is printed), the discrete colors
%  (accent, cortical lobes, hippocampus), the study logo (last logo of the
%  design that is not the CAT logo), the CAT logo, the thickness colormap,
%  and the study texts (project, organization and website of the study
%  JSON).
  Pcat = fileparts(mfilename('fullpath'));
  j = job; j.opts.reportfig.design = {'study'};
  P = cat_io_dcm2bids_reportfig('reportDesigns',j);
  D = cat_io_dcm2bids_reportfig('loadDesign',P{1});
  if mean(D.Colors(4,:)) < 0.5
    B = cat_io_dcm2bids_reportfig('loadDesign',fullfile(Pcat,'DCM2BIDS','colormaps','CAT_bright.json')); D.Colors = B.Colors;
  end
  G.fg  = D.Colors(1,:); G.fgl = D.Colors(2,:); G.bgl = D.Colors(3,:); G.bg = [1 1 1];
  G.disc    = cat_io_dcm2bids_reportfig('colormapOf',D.Discrete,8);
  G.acc     = G.disc(1,:);
  G.lobecol = G.disc([2 5 4 8 6 7],:);      % frontal, temporal, parietal, occipital, cingulate, insula
  G.hccol   = G.disc(3,:);
  G.headcol = 0.55 * G.fgl + 0.45;
  G.Pcat    = fullfile(Pcat,'doc','images','CAT.png');
  Ps        = D.Logos(~strcmpi(cellfun(@(p) spm_file(p,'filename'),D.Logos,'UniformOutput',false),'CAT.png'));
  G.Plogo   = ''; if ~isempty(Ps), G.Plogo = Ps{end}; end
  G.cm      = interp1([0 0.5 1],[0.36 0.29 0.32; 0.70 0.58 0.62; 0.97 0.89 0.91],linspace(0,1,256)); % thickness 2 to 4 mm
  G.study   = struct('Study','','Organization','','Website','');
  if isfield(job,'subdir'), G.study.Study = job.subdir; end
  if isfield(job,'dataset') && isfield(job.dataset,'study') && exist(char(job.dataset.study),'file')
    try
      S = jsondecode(fileread(char(job.dataset.study)));
      if isfield(S,'Name'), G.study.Organization = char(S.Name); end
      if isfield(S,'ReferencesAndLinks') && ~isempty(S.ReferencesAndLinks) % (empty list: [])
        w = cellstr(S.ReferencesAndLinks); w = w(~cellfun('isempty',regexp(w,'^(https?://|www\.)','once')));
        if ~isempty(w), G.study.Website = regexprep(w{1},'^https?://|/$',''); end
      end
    catch e
      cat_io_cprintf('warn','  Participant report: study JSON not read (%s).\n',e.message);
    end
  end
end
% =========================================================================
function L = language(lang,Pstudydir)
%language. Texts of the language lang (<lang>.json of the study directory
%  Pstudydir, otherwise of DCM2BIDS/participantreport) with the English
%  texts for missing fields (and unknown languages).
  Pd = fullfile(fileparts(mfilename('fullpath')),'DCM2BIDS','participantreport');
  L  = jsondecode(fileread(fullfile(Pd,'en.json')));
  P  = {fullfile(Pd,[lang '.json'])};
  if nargin > 1 && ~isempty(Pstudydir), P = [{fullfile(Pstudydir,[lang '.json'])} P]; end
  P  = P(cellfun(@(f) exist(f,'file') == 2,P));
  if isempty(P)
    if ~strcmp(lang,'en'), cat_io_cprintf('warn','  Participant report: no language file "%s" (English is used).\n',lang); end
    return
  end
  L = mergeStruct(L,jsondecode(fileread(P{1})));
end
% =========================================================================
function a = mergeStruct(a,b)
%mergeStruct. Fields of the structure b in a (recursively).
  for fn = fieldnames(b)'
    if isfield(a,fn{1}) && isstruct(a.(fn{1})) && isstruct(b.(fn{1})), a.(fn{1}) = mergeStruct(a.(fn{1}),b.(fn{1}));
    else, a.(fn{1}) = b.(fn{1}); end
  end
end
% =========================================================================
function N = normModels(P)
%normModels. Reference models of the measures (JSON P, '' for DCM2BIDS/
%  participantreport/norms.json): mean (quadratic in age - AgeCenter) and
%  SD of each measure and sex.
  if isempty(P), P = fullfile(fileparts(mfilename('fullpath')),'DCM2BIDS','participantreport','norms.json'); end
  N = jsondecode(fileread(P));
  if iscell(N.Measures), N.Measures = [N.Measures{:}]; end
end
% =========================================================================
% data of the sessions
% =========================================================================
function R = sessionData(Poutdir,job,Psub)
%sessionData. Sessions of the BIDS directory of the job with a bias-
%  corrected T1w of the basic preprocessing (index of the derivatives of
%  the database, catDCM2BIDSderivatives.json, i.e. the BIDS derivative or
%  its database source for opts.preprocessing = 1), one T1w per session
%  (the one with CAT surfaces, otherwise the first): T1w (msub*), label
%  map (c0), atlas (lotter2025), segmentation values (catDCM2BIDSsegus),
%  registration (catDCM2BIDSaffine, rigid part for the head orientation
%  and the face plane of spm_deface by the affine part), DEC map in the
%  space of the T1w (catDCM2BIDSdwi_DECt1w), surfaces (CAT, quick
%  surfaces), and the identity, institute, and field strength of the
%  original DICOM header (sidecar of the database). Psub limits the
%  sessions to these subject directories of the database (default all).
%  Sessions without BIDS session directory (e.g. quarantined sessions, see
%  quarantineSession in cat_io_dcm2bids_process) are skipped.
  if ~exist('Psub','var') || isempty(Psub), Psub = {fullfile(Poutdir,'catDCM2BIDSdb')}; dp = 3; else, Psub = cellstr(Psub); dp = 2; end
  R    = struct([]);
  Pbr  = [strrep(fullfile(job.BIDSdir,job.BIDSsubdir),'\','/') '/'];
  Pcat = cat_io_dcm2bids_pp_cat('derivativesName',0);
  Pidx = {}; % (<db>/<sub>/<ses>/)
  for si = 1:numel(Psub), Pidx = [Pidx; cat_vol_findfiles(Psub{si},'catDCM2BIDSderivatives.json',struct('depth',dp))]; end %#ok<AGROW>
  for ii = 1:numel(Pidx)
    try
      D = jsondecode(fileread(Pidx{ii})); E = D.derivatives;
      if iscell(E), E = cellfun(@(e) struct('derivative',e.derivative,'source',e.source),E); end
    catch
      continue
    end
    if isempty(E), continue; end
    der = strrep({E.derivative},'\','/'); src = fullfile(fileparts(Pidx{ii}),{E.source});
    ok  = startsWith(der,Pbr) & contains(der,'/derivatives/catDCM2BIDS/'); der = der(ok); src = src(ok);
    fil = cellfun(@(d,s) existing({fullfile(Poutdir,d),s}),der,src,'UniformOutput',false);
    nam = regexprep(der,'^.*/','');
    for ti = find(~cellfun('isempty',regexp(nam,'^msub-.*_T1w\.nii(\.gz)?$','once')))
      bn  = regexprep(nam{ti},'^m|\.nii(\.gz)?$','');
      tok = regexp(bn,'^(sub-[^_]+)_?(ses-[^_]+)?','tokens','once'); sub = tok{1}; ses = tok{2};
      if ~exist(fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,sub,ses),'dir'), continue; end % (not exported, e.g. quarantine)
      Pd  = fileparts(der{ti});                         % BIDS derivatives directory of the T1w
      byName = @(pre,ext) fil(strcmp(nam,[pre bn ext]) | strcmp(nam,[pre bn ext '.gz']));
      name = sub; if ~isempty(ses), name = [sub '_' ses]; end
      r = struct('sub',sub,'ses',ses,'name',name,'bn',bn, ...
        'T1w',fil{ti},'c0',first(byName('c0','.nii')),'atlas',first(byName('lotter2025_','.nii')),'DEC','', ...
        'segus','','qc','','sidecar','','Rigid',eye(4),'face',[0 -1.1 0.98 100],'surf',struct('kind','','P',{{}}), ...
        'id',struct(),'inst',struct('Institution','','Address',''),'field',NaN,'age',NaN);
      % segmentation values and original sidecar (database source of the T1w)
      base = regexprep(spm_file(src{ti},'filename'),'^m(anon_)?|\.nii(\.gz)?$','');
      r.segus   = existing({fullfile(Poutdir,Pd,['catDCM2BIDSsegus_' bn '.json']), fullfile(fileparts(src{ti}),['catDCM2BIDSsegus_' base '.json'])});
      r.sidecar = existing({fullfile(fileparts(src{ti}),[base '.json'])});
      r.qc      = existing({fullfile(Poutdir,Pd,['catDCM2BIDSqc_' bn '.json']), ...
        fullfile(fileparts(src{ti}),['catDCM2BIDSqc_anon_' base '.json']), fullfile(fileparts(src{ti}),['catDCM2BIDSqc_' base '.json'])});
      % registration of the session (rigid head orientation and the face plane by the affine part)
      Psd = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,'derivatives','catDCM2BIDS',sub,ses);
      Paf = fullfile(Psd,['catDCM2BIDSaffine_' r.name '.json']);
      if exist(Paf,'file')
        try
          T = jsondecode(fileread(Paf)); r.Rigid = double(T.Rigid); r.face = [0 -1.1 0.98 100] * double(T.Affine);
        catch e
          cat_io_cprintf('warn','  Participant report: registration of %s not read (%s).\n',r.name,e.message);
        end
      end
      % DEC map in the space of the T1w (DWI of the session)
      dec = find(~cellfun('isempty',regexp(der,['/' sub '/' ses '/?dwi/catDCM2BIDSdwi_DECt1w_'],'once')) & ~cellfun('isempty',fil),1);
      if ~isempty(dec), r.DEC = fil{dec}; end
      % surfaces: CAT, otherwise the quick surfaces
      Pc = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,'derivatives',Pcat,sub,ses,'anat');
      Ph = fullfile(Pc,strcat({'lh','rh'},'.central.',bn,'.gii'));
      if all(cellfun(@(f) exist(f,'file') == 2,Ph))
        r.surf = struct('kind','CAT','P',{Ph});
      else
        Pq = [byName('catDCM2BIDSsurf_lh_','.gii') byName('catDCM2BIDSsurf_rh_','.gii')];
        if numel(Pq) == 2 && all(~cellfun('isempty',Pq)), r.surf = struct('kind','quick','P',{Pq}); end
      end
      r = sidecarData(r);
      % one T1w per session (the one with CAT surfaces)
      si = []; if ~isempty(R), si = find(strcmp({R.name},r.name),1); end
      if isempty(si), R = [R r]; %#ok<AGROW>
      elseif strcmp(r.surf.kind,'CAT') && ~strcmp(R(si).surf.kind,'CAT'), R(si) = r;
      end
    end
  end
end
% =========================================================================
function r = sidecarData(r)
%sidecarData. Identity (ID, name, birth date, sex, weight, scan date, age),
%  institute, and field strength of the original DICOM sidecar of the
%  session r (database).
  r.id = struct('PID','','Name','','Birth','','Sex','','Weight',NaN,'Scan','');
  if isempty(r.sidecar), return; end
  try J = jsondecode(fileread(r.sidecar)); catch, return; end
  f = @(fn) fieldOf(J,fn,'');
  nm = strsplit(char(string(f('PatientName'))),'^');                 % (DICOM PN: family^given^middle)
  if numel(nm) > 1, nm = strtrim(strjoin([nm(2:min(3,end)) nm(1)],' ')); else, nm = nm{1}; end
  r.id = struct('PID',char(string(f('PatientID'))),'Name',nm,'Birth',char(string(f('PatientBirthDate'))), ...
    'Sex',char(string(f('PatientSex'))),'Weight',double(fieldOf(J,'PatientWeight',NaN)),'Scan','');
  t = char(string(f('AcquisitionDateTime'))); if isempty(t), t = char(string(f('AcquisitionDate'))); end
  if numel(t) >= 10, r.id.Scan = t(1:10); end
  r.age = ageOf(r.id.Birth,r.id.Scan);
  if ~isfinite(r.age), r.age = double(fieldOf(J,'PatientAge',NaN)); end
  r.inst  = struct('Institution',char(string(f('InstitutionName'))), ...
    'Address',regexprep(char(string(f('InstitutionAddress'))),'\s*,\s*',', '));
  r.field = double(fieldOf(J,'MagneticFieldStrength',NaN));
end
% =========================================================================
function R = demoSession
%demoSession. Session of the CAT templates (MNI152NLin2009cAsym) for the
%  demo: T1w with simulated defacing (plane of spm_deface outside of the
%  smoothed brain mask), label map of the shooting template (CSF 1, GM 2,
%  WM 3), lotter2025, DEC map of the hcp1065 tracts by their main
%  direction (FA 0.6), template surfaces with thickness and the sulcal
%  depth by the square root of the sulc, and placeholder identity, values,
%  and institute.
  C = fileparts(mfilename('fullpath')); T = fullfile(C,'templates_MNI152NLin2009cAsym');
  Vt = spm_vol(fullfile(T,'Template_T1.nii')); Yt = spm_read_vols(Vt);
  Vb = spm_vol(fullfile(T,'brainmask.nii'));   Yb = spm_read_vols(Vb); Bs = zeros(size(Yb)); spm_smooth(double(Yb > 0.5),Bs,[10 10 10]);
  [i,j,k] = ndgrid(1:size(Yt,1),1:size(Yt,2),1:size(Yt,3)); Xw = Vt.mat * [i(:) j(:) k(:) ones(numel(i),1)]';
  b  = sampleAt(struct('Y',Bs,'mat',Vb.mat),Xw(1:3,:),1);
  Yt(([0 -1.1 0.98 100] * Xw)' < 0 & b < 0.02) = 0;
  Vg = spm_vol(fullfile(T,'Template_4_GS.nii')); Yg = spm_read_vols(Vg(1)); Yw = spm_read_vols(Vg(2));
  [i,j,k] = ndgrid(1:size(Yg,1),1:size(Yg,2),1:size(Yg,3)); Xg = Vg(1).mat * [i(:) j(:) k(:) ones(numel(i),1)]';
  bg = reshape(sampleAt(struct('Y',double(Yb),'mat',Vb.mat),Xg(1:3,:),1),size(Yg));
  c0 = (1 + Yg + 2 * Yw) .* (bg > 0.5);
  % DEC placeholder: tracts colored by their main direction
  Vh = spm_vol(fullfile(T,'hcp1065.nii')); Yh = round(spm_read_vols(Vh));
  H  = regexp(fileread(fullfile(T,'hcp1065.csv')),'\r?\n','split'); H = H(2:end); H = H(~cellfun('isempty',strtrim(H)));
  H  = cellfun(@(l) strsplit(l,';'),H,'UniformOutput',false); hid = cellfun(@(c) str2double(c{1}),H); hab = cellfun(@(c) c{2},H,'UniformOutput',false);
  hd = repmat([0.15 0.80 0.25],numel(hid),1);                                         % association: A-P
  hd(~cellfun('isempty',regexp(hab,'^bCC','once')),:) = repmat([0.90 0.15 0.15],nnz(~cellfun('isempty',regexp(hab,'^bCC','once'))),1);
  iv = ~cellfun('isempty',regexp(hab,'^[lr](CBT|CPT|CST|CS[ASP]|FAT|PAT|RST|TRS|VOF)$','once'));
  hd(iv,:) = repmat([0.20 0.35 0.95],nnz(iv),1);                                      % projection: S-I
  [ok,hi] = ismember(Yh,hid); Yd = zeros([size(Yh) 3]);
  for c = 1:3, d = zeros(size(Yh)); d(ok) = 0.6 * hd(hi(ok),c); Yd(:,:,:,c) = d; end
  % surfaces of the template (pial-like by the half thickness, depth to the convex hull)
  S = struct('vertices',{},'faces',{},'thickness',{},'depth',{});
  for h = {'lh','rh'}
    g = gifti(fullfile(C,'templates_surfaces_32k',[h{1} '.central.Template_T1.gii']));
    t = cat_io_FreeSurfer('read_surf_data',fullfile(C,'templates_surfaces_32k',[h{1} '.thickness.Template_T1']));
    v = double(g.vertices); f = double(g.faces); v = v + vertexNormals(v,f) .* t/2;
    S(end+1) = struct('vertices',v,'faces',f,'thickness',t,'depth',cat_io_dcm2bids_pp_surf('hullDepth',v)); %#ok<AGROW>
  end
  R = struct('sub','sub-DEMO','ses','ses-1','name','sub-DEMO_ses-1','bn','', ...
    'T1w',struct('Y',Yt,'mat',Vt.mat),'c0',struct('Y',c0,'mat',Vg(1).mat), ...
    'atlas',fullfile(T,'lotter2025_aparc.nii'),'DEC',struct('Y',Yd,'mat',Vh.mat), ...
    'segus','','qc','','sidecar','','Rigid',eye(4),'face',[0 -1.1 0.98 100],'surf',struct('kind','data','P',{{}}), ...
    'id',struct('PID','P0012','Name','Max Mustermann','Birth','1980-01-01','Sex','M','Weight',78, ...
    'Scan',char(datetime('now','Format','yyyy-MM-dd'))),'inst',struct('Institution','Example Institute of Radiology', ...
    'Address','Example Street 1, 00000 City (placeholder)'),'field',3,'age',NaN);
  R.age = ageOf(R.id.Birth,R.id.Scan); R.surfdata = S; R.values = struct('GMV',612,'WMV',548,'CT',2.31,'HCV',7.9,'SQR',2.1);
end
% =========================================================================
function M = measureValues(R,I,N)
%measureValues. Measures of the reference models N with the values of the
%  session R (GMV, WMV, and HCV of the segmentation, CT of the surface)
%  and their percentile P by the model of the sex (mean of both models for
%  other or unknown sex) at the age of the session (n/a outside of the age
%  range of the models).
  v = struct('GMV',NaN,'WMV',NaN,'CT',I.ct,'HCV',NaN);
  if isfield(R,'values')
    for fn = fieldnames(R.values)', v.(fn{1}) = R.values.(fn{1}); end
  elseif ~isempty(R.segus)
    try
      J = jsondecode(fileread(R.segus)); v.GMV = fieldOf(J,'aGMV',NaN); v.WMV = fieldOf(J,'aWMV',NaN);
      if isfield(J,'regional') && all(isfield(J.regional,{'lHIP','rHIP'}))
        v.HCV = fieldOf(J.regional.lHIP,'Vgm',NaN) + fieldOf(J.regional.rHIP,'Vgm',NaN);
      end
    catch e
      cat_io_cprintf('warn','  Participant report: segmentation values of %s not read (%s).\n',R.name,e.message);
    end
  end
  ac = N.AgeCenter; ar = N.AgeRange(:)';
  M  = struct('name',{},'unit',{},'fmt',{},'val',{},'mu',{},'sd',{},'P',{},'age',{},'range',{});
  for i = 1:numel(N.Measures)
    Ni = N.Measures(i);
    switch upper(R.id.Sex)
      case 'M', b = Ni.Male.Mean(:)'; sd = Ni.Male.SD;
      case 'F', b = Ni.Female.Mean(:)'; sd = Ni.Female.SD;
      otherwise, b = (Ni.Male.Mean(:)' + Ni.Female.Mean(:)')/2; sd = (Ni.Male.SD + Ni.Female.SD)/2;
    end
    mu  = @(x) b(1) + b(2)*(x - ac) + b(3)*(x - ac).^2;
    val = NaN; if isfield(v,Ni.Name), val = double(v.(Ni.Name)); end
    P   = NaN; if isfinite(val) && R.age >= ar(1) && R.age <= ar(2), P = round(100 * 0.5*erfc(-((val - mu(R.age))/sd)/sqrt(2))); end
    M(end+1) = struct('name',Ni.Name,'unit',Ni.Unit,'fmt',Ni.Format,'val',val,'mu',mu,'sd',sd,'P',P,'age',R.age,'range',ar); %#ok<AGROW>
  end
end
% -------------------------------------------------------------------------
function Q = qcRating(R)
%qcRating. Overall quality rating of the T1w (SQR of the QC of
%  cat_io_dcm2bids_qc, in marks, e.g. reduced by motion artifacts) as
%  grade (A+ to F) and percentage rating points of CAT (see cat_stat_marks)
%  with the color of the marks (as the CAT reports), empty grade without.
  Q = struct('mark',NaN,'grade','','rps',NaN,'col',[0.5 0.5 0.5]);
  if isfield(R,'values') && isfield(R.values,'SQR'), Q.mark = R.values.SQR;
  elseif ~isempty(R.qc)
    try
      J = jsondecode(fileread(R.qc));
      if isfield(J,'qualityrating'), Q.mark = double(fieldOf(J.qualityrating,'SQR',NaN)); end
    catch e
      cat_io_cprintf('warn','  Participant report: QC of %s not read (%s).\n',R.name,e.message);
    end
  end
  if ~isscalar(Q.mark) || ~isfinite(Q.mark), Q.mark = NaN; return; end
  grades = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  Q.grade = grades{max(1,min(numel(grades),round((Q.mark + 2/3)*3 - 3)))};
  Q.rps   = min(100,max(0,105 - Q.mark*10));
  cm      = min(1,max(0,cat_io_colormaps('marks+',40)));
  Q.col   = cm(min(size(cm,1) - 3,max(1,floor(Q.mark/9.5*size(cm,1)))),:);
end
% =========================================================================
% language-free images
% =========================================================================
function I = sessionImages(R,G,A)
%sessionImages. Language-free images of the session R: slices (see
%  sliceImages) and surface rendering with the mean cortical thickness ct
%  (see surfaceImage), DEC map and thickness available (hasdec, hasct).
  V = loadVolumes(R);
  [I.sl,wmI] = sliceImages(R,V,G,A);
  [I.surf,I.ct,I.co] = surfaceImage(R,V,G,A,wmI);
  I.hasdec = ~isempty(V.DEC.Y); I.hasct = isfinite(I.ct);
end
% -------------------------------------------------------------------------
function V = loadVolumes(R)
%loadVolumes. Volumes of the session R (T1w, label map c0, atlas, DEC with
%  3 channels) with the data Y and the voxel-to-world matrix mat (empty Y
%  if not available, given volumes of the demo are kept).
  for fn = {'T1w','c0','atlas','DEC'}
    P = R.(fn{1});
    if isstruct(P), V.(fn{1}) = P; continue; end
    V.(fn{1}) = struct('Y',[],'mat',eye(4));
    if isempty(P) || ~exist(P,'file'), continue; end
    evalc('Vi = spm_vol(P);'); % (avoid messages of gz files)
    V.(fn{1}) = struct('Y',spm_read_vols(Vi),'mat',Vi(1).mat);
  end
end
% -------------------------------------------------------------------------
function v = sampleAt(S,Xw,hold)
%sampleAt. Values of the volume S (see loadVolumes) at the world
%  coordinates Xw (3 x n) with the interpolation hold (channels of 4D
%  volumes as columns), empty without volume.
  v = []; if isempty(S.Y), return; end
  X = S.mat \ [Xw; ones(1,size(Xw,2))]; nc = size(S.Y,4); v = zeros(size(Xw,2),nc);
  for c = 1:nc, v(:,c) = spm_sample_vol(double(S.Y(:,:,:,c)),X(1,:),X(2,:),X(3,:),hold); end
end
% -------------------------------------------------------------------------
function [S,wmI] = sliceImages(R,V,G,A)
%sliceImages. RGB images (0.3 mm, head orientation) of the sagittal,
%  coronal, and axial slice through the centroid of the right hippocampus
%  with the T1w (WM 0.87), the cortical lobes lightly in the GM, the
%  hippocampus (filled and outlined), the DEC map in the WM (by the FA),
%  the removed face (plane of spm_deface, black and hatched up to 30 mm
%  from the head), and the positions of the labels (hippocampus, GM, WM,
%  and lateral ventricles), returns also the WM intensity of the T1w.
  r = 0.3; xr = -95:r:95; yr = -130:r:95; zr = -75:r:100; c = [28 -20 -14];
  if ~isempty(V.atlas.Y)
    [i,j,k] = ind2sub(size(V.atlas.Y),find(ismember(round(V.atlas.Y),A.rhc)));
    if ~isempty(i), h = R.Rigid * V.atlas.mat * [mean(i) mean(j) mean(k) 1]'; c = h(1:3)'; end
  end
  o = {'sag','cor','axi'}; C = cell(1,3);
  for oi = 1:3
    switch o{oi}
      case 'sag', [Z,Y] = ndgrid(fliplr(zr),yr); X = c(1)*ones(size(Y)); U = Y; W = Z; xd = yr([1 end]); yd = zr([end 1]);
      case 'cor', [Z,X] = ndgrid(fliplr(zr),xr); Y = c(2)*ones(size(X)); U = X; W = Z; xd = xr([1 end]); yd = zr([end 1]);
      case 'axi', [Y,X] = ndgrid(fliplr(yr),xr); Z = c(3)*ones(size(X)); U = X; W = Y; xd = xr([1 end]); yd = yr([end 1]);
    end
    Xw = R.Rigid \ [X(:) Y(:) Z(:) ones(numel(X),1)]'; sz = size(X);
    s  = struct('I',zeros(sz),'P0',zeros(sz),'L',zeros(sz),'D',[], ...
      'pv',reshape(R.face * Xw,sz),'X',X,'U',U,'W',W,'xd',xd,'yd',yd);
    if ~isempty(V.T1w.Y),   s.I  = reshape(sampleAt(V.T1w,Xw(1:3,:),1),sz); end
    if ~isempty(V.c0.Y),    s.P0 = reshape(sampleAt(V.c0,Xw(1:3,:),1),sz); end
    if ~isempty(V.atlas.Y), s.L  = round(reshape(sampleAt(V.atlas,Xw(1:3,:),0),sz)); end
    if ~isempty(V.DEC.Y),   s.D  = abs(reshape(sampleAt(V.DEC,Xw(1:3,:),1),[sz 3])); end
    C{oi} = s;
  end
  I0  = cell2mat(cellfun(@(s) s.I(:),C,'UniformOutput',false)'); P0 = cell2mat(cellfun(@(s) s.P0(:),C,'UniformOutput',false)');
  wmI = median(I0(P0 > 2.75)); if ~isfinite(wmI) || wmI <= 0, wmI = cat_stat_prctile(I0(I0 > 0),95); end
  clip = @(x) min(1,max(0,x));
  for oi = 1:3
    s = C{oi}; I = clip(s.I / (1.15 * wmI)); rgb = repmat(I,[1 1 3]);
    gm = clip(1 - abs(s.P0 - 2)); wm = clip(s.P0 - 2);
    % cortical lobes (light colors of the design)
    lobe = A.lut(min(numel(A.lut),max(0,s.L) + 1)); ok = lobe > 0; a = 0.45 * gm .* ok;
    for ci = 1:3
      lc = zeros(size(I)); lc(ok) = G.lobecol(lobe(ok),ci); lc = 0.6 * lc + 0.4;
      rgb(:,:,ci) = rgb(:,:,ci) .* (1 - a) + a .* lc .* (0.35 + 0.75 * I);
    end
    % DEC map in the WM (by the FA, i.e. the norm of the FA-weighted directions)
    if ~isempty(s.D)
      nd = sqrt(sum(s.D.^2,3)); a = 0.85 * clip((nd - 0.08) / 0.3) .* wm; dc = s.D ./ max(eps,max(s.D,[],3));
      for ci = 1:3, rgb(:,:,ci) = rgb(:,:,ci) .* (1 - a) + a .* dc(:,:,ci) .* (0.5 + 0.6 * I); end
    end
    % hippocampus
    hc = ismember(s.L,A.hc); he = hc & ~imerode(hc,ones(2*round(0.5/r) + 1)); a = 0.35 * hc + 0.5 * he; % (outline 0.5 mm)
    for ci = 1:3, rgb(:,:,ci) = rgb(:,:,ci) .* (1 - a) + a * G.hccol(ci); end
    % removed face (black with hatch up to 30 mm from the head)
    hd = s.I > 0.08 * wmI; f = s.pv < 0 & ~hd & bwdist(hd) * r < 30;
    [ii,jj] = ndgrid(1:size(f,1),1:size(f,2)); hatch = f & mod(ii + jj,round(7/r)) < round(1.5/r); % (lines of 1.5 mm every 7 mm)
    rgb = rgb .* repmat(~f,[1 1 3]) + 0.30 * repmat(hatch,[1 1 3]);
    % positions of the labels (most interior pixel of the structure)
    rhc = hc; if ~strcmp(o{oi},'sag'), rhc = hc & s.X > 0; end
    up  = s.W > c(3) + 10; lh = s.X < 0; if strcmp(o{oi},'sag'), lh = true(size(lh)); end
    S.(o{oi}) = struct('rgb',min(1,max(0,rgb)),'xd',s.xd,'yd',s.yd,'hc',centroid(rhc,s.U,s.W), ...
      'gm',interior(gm > 0.6 & ok & lh & up,s.U,s.W),'wm',interior(wm > 0.8 & lh & up,s.U,s.W), ...
      'csf',interior(ismember(s.L,A.lv),s.U,s.W));
  end
end
% -------------------------------------------------------------------------
function p = centroid(m,U,W)
%centroid. Centroid of the mask m in the slice coordinates U, W (empty).
  p = []; if any(m(:)), p = [mean(U(m)) mean(W(m))]; end
end
% -------------------------------------------------------------------------
function p = interior(m,U,W)
%interior. Most interior pixel of the mask m (maximum distance to its
%  border) in the slice coordinates U, W (empty for an empty mask).
  p = []; if ~any(m(:)), return; end
  d = bwdist(~m); [~,i] = max(d(:)); p = [U(i) W(i)];
end
% -------------------------------------------------------------------------
function [img,ct,co] = surfaceImage(R,V,G,A,wmI)
%surfaceImage. Rendering (RGB, cropped) of the brain surfaces (see
%  brainSurfaces) in the transparent head (see headSurface, view-dependent
%  transparency of cat_surf_render2, i.e. transparent where the head faces
%  the camera) from the right-front-top, with the mean cortical thickness
%  and the callouts of thin and thick cortex (see thicknessCallouts).
  az = 118; el = 20;
  [B,ct,S] = brainSurfaces(R,V,G,A); H = headSurface(R,V,wmI);
  co = thicknessCallouts(R,S,[vertcat(B.vertices); H.vertices],az,el);
  fh = figure('Visible','off','Color','w','Units','inches','Position',[0 0 6 6],'InvertHardcopy','off', ...
    'PaperUnits','inches','PaperSize',[6 6],'PaperPosition',[0 0 6 6]);
  fhclean = onCleanup(@() close(fh));
  ax = axes(fh,'Position',[0 0 1 1]); hold(ax,'on');
  for bi = 1:numel(B)
    patch(ax,'Vertices',B(bi).vertices,'Faces',B(bi).faces,'FaceVertexCData',B(bi).col,'FaceColor','interp', ...
      'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.40,'DiffuseStrength',0.70, ...
      'SpecularStrength',0.05,'SpecularExponent',10);
  end
  hp = [];
  if ~isempty(H.vertices)
    hp = patch(ax,'Vertices',H.vertices,'Faces',H.faces,'FaceColor',G.headcol,'EdgeColor','none', ...
      'FaceLighting','gouraud','AmbientStrength',0.55,'DiffuseStrength',0.45,'SpecularStrength',0.1);
  end
  axis(ax,'equal','off'); view(ax,az,el); camlight(ax,'headlight'); camlight(ax,30,40);
  set(ax,'AmbientLightColor',[1 1 1]); lighting(ax,'gouraud'); material(ax,'dull');
  if ~isempty(hp), cat_surf_render2('ViewTransparency',hp,[0.02 0.5 2.5]); end
  img = print(fh,'-RGBImage','-r600','-image'); clear fhclean   % (3600 x 3600 pixel, exportgraphics halves it on HiDPI systems)
  m = any(img < 245,3);
  if any(m(:)), img = img(find(any(m,2),1):find(any(m,2),1,'last'),find(any(m,1),1):find(any(m,1),1,'last'),:); end
end
% -------------------------------------------------------------------------
function co = thicknessCallouts(R,S,Vall,az,el)
%thicknessCallouts. Callouts {point, label position, label} in fractions
%  of the cropped rendering (orthographic projection of the view az, el,
%  the crop is the projected bounding box of all vertices Vall) of thin
%  and thick cortex (3rd and 97th percentile of the thickness of visible
%  gyral crowns of the right hemisphere, i.e. facing the camera with a
%  sulcal depth below 1.5 mm, labeled with their thickness).
  co = {}; if numel(S) < 2 || isempty(Vall), return; end
  T  = viewmtx(az,el); vd = [sind(az)*cosd(el) -cosd(az)*cosd(el) sind(el)];
  pa = T(1:2,:) * [Vall'; ones(1,size(Vall,1))]; bx = [min(pa(1,:)) max(pa(1,:))]; by = [min(pa(2,:)) max(pa(2,:))];
  v  = toHead(R.Rigid,S(2).vertices); t = S(2).thickness; d = S(2).depth;
  ok = vertexNormals(v,S(2).faces) * vd' > 0.6 & isfinite(t) & t > 0.5; if ~isempty(d), ok = ok & d < 1.5; end
  if nnz(ok) < 20, return; end
  p  = T(1:2,:) * [v(ok,:)'; ones(1,nnz(ok))]; f = [(p(1,:) - bx(1))/diff(bx); (by(2) - p(2,:))/diff(by)]'; t = t(ok);
  [~,i1] = min(abs(t - cat_stat_prctile(t,3))); [~,i2] = min(abs(t - cat_stat_prctile(t,97)));
  lx = @(p) 0.08 + 0.84 * (p(1) > 0.5);   % (label at the side of the point, thin at the top, thick at the bottom)
  co = {f(i1,:),[lx(f(i1,:)) 0.06],sprintf('%0.1f mm',t(i1)); f(i2,:),[lx(f(i2,:)) 0.94],sprintf('%0.1f mm',t(i2))};
end
% -------------------------------------------------------------------------
function [B,ct,S] = brainSurfaces(R,V,G,A)
%brainSurfaces. Brain surfaces of the rendering in head coordinates with
%  their colors: pial surfaces with the thickness (gray-pink colormap of
%  the design, 2 to 4 mm) darkened in the sulci by the sulcal depth (see
%  thicknessColors), otherwise the GM boundary of the label map (without
%  thickness), and the cerebellum and brainstem of the atlas (constant
%  color), with the mean cortical thickness ct (NaN without) and the pial
%  surfaces S (see loadSurfaces).
  B = struct('vertices',{},'faces',{},'col',{}); ct = NaN;
  S = loadSurfaces(R);
  if ~isempty(S)
    t = [];
    for si = 1:numel(S)
      B(end+1) = struct('vertices',toHead(R.Rigid,S(si).vertices),'faces',S(si).faces, ...
        'col',thicknessColors(S(si).thickness,G.cm,S(si).depth)); %#ok<AGROW>
      t = [t; S(si).thickness(isfinite(S(si).thickness) & S(si).thickness > 0.5)]; %#ok<AGROW>
    end
    ct = mean(t);
  elseif ~isempty(V.c0.Y)
    L = zeros(size(V.c0.Y));
    if ~isempty(V.atlas.Y) && isequal(size(V.atlas.Y),size(V.c0.Y)), L = round(V.atlas.Y); end
    fv = maskSurface(V.c0.Y > 1.5 & ~ismember(L,A.cb),V.c0.mat,R.Rigid,3,120000);
    B(end+1) = struct('vertices',fv.vertices,'faces',fv.faces,'col',repmat(G.cm(128,:),size(fv.vertices,1),1));
  end
  if ~isempty(V.atlas.Y)
    fv = maskSurface(ismember(round(V.atlas.Y),A.cb),V.atlas.mat,R.Rigid,5,40000);
    if ~isempty(fv.vertices), B(end+1) = struct('vertices',fv.vertices,'faces',fv.faces,'col',repmat(G.cm(110,:),size(fv.vertices,1),1)); end
  end
end
% -------------------------------------------------------------------------
function S = loadSurfaces(R)
%loadSurfaces. Pial surfaces (world coordinates) of both hemispheres with
%  thickness and sulcal depth: CAT (pial surface, otherwise central
%  surface shifted by the half thickness, depth of surfextract), quick
%  surfaces (central surface shifted by the half thickness, columns
%  thickness and depth of the cdata, medial wall without thickness), or
%  the surfaces of the demo, otherwise empty.
  S = struct('vertices',{},'faces',{},'thickness',{},'depth',{});
  if isfield(R,'surfdata'), S = R.surfdata; return; end
  for hi = 1:numel(R.surf.P)
    try
      g = gifti(R.surf.P{hi}); v = double(g.vertices); f = double(g.faces);
      switch R.surf.kind
        case 'CAT'
          Pt = regexprep(R.surf.P{hi},'\.central\.(.*)\.gii$','.thickness.$1');
          Pd = regexprep(R.surf.P{hi},'\.central\.(.*)\.gii$','.depth.$1');
          Pp = strrep(R.surf.P{hi},'.central.','.pial.');
          t  = cat_io_FreeSurfer('read_surf_data',Pt);
          if exist(Pp,'file'), p = gifti(Pp); v = double(p.vertices); else, v = v + vertexNormals(v,f) .* t/2; end
          if exist(Pd,'file'), d = cat_io_FreeSurfer('read_surf_data',Pd); else, d = cat_io_dcm2bids_pp_surf('hullDepth',v); end
        case 'quick'
          D = double(g.cdata); t = D(:,1); d = D(:,4); t(D(:,5) == 0) = NaN;
          v = v + vertexNormals(v,f) .* max(0,D(:,1))/2;
      end
      S(end+1) = struct('vertices',v,'faces',f,'thickness',t(:),'depth',d(:)); %#ok<AGROW>
    catch e
      cat_io_cprintf('warn','  Participant report: surface %s not read (%s).\n',R.surf.P{hi},e.message);
      S = S([]); return
    end
  end
end
% -------------------------------------------------------------------------
function col = thicknessColors(t,cm,d)
%thicknessColors. Colors of the thickness t (2 to 4 mm by the colormap
%  cm, dark without thickness, e.g. the medial wall) with the sulcal depth
%  shadow (darkened by up to 50% at 15 mm depth d).
  k   = round(min(1,max(0,(t - 2)/2)) * (size(cm,1) - 1)) + 1; k(~isfinite(k)) = 1;
  col = cm(k,:); col(~isfinite(t),:) = repmat(0.8 * cm(1,:),nnz(~isfinite(t)),1);
  if ~isempty(d), d(~isfinite(d)) = 0; col = col .* (1 - 0.5 * min(1,max(0,d/15)).^0.8); end
end
% -------------------------------------------------------------------------
function H = headSurface(R,V,wmI)
%headSurface. Outer surface of the head of the T1w in head coordinates:
%  filled largest component of the smoothed T1w (4 mm) above 15% of the WM
%  intensity at about 2 mm, i.e. with the removed face of the defacing,
%  without the neck (below 78 mm under the AC in the head orientation).
  H = struct('vertices',zeros(0,3),'faces',zeros(0,3));
  if isempty(V.T1w.Y) || ~isfinite(wmI), return; end
  vx = sqrt(sum(V.T1w.mat(1:3,1:3).^2)); st = max(1,round(2 ./ vx));
  Ys = zeros(size(V.T1w.Y)); spm_smooth(double(V.T1w.Y),Ys,4 ./ vx);
  Ys = Ys(1:st(1):end,1:st(2):end,1:st(3):end) > 0.15 * wmI;
  cc = bwconncomp(Ys,6); if cc.NumObjects == 0, return; end
  [~,li] = max(cellfun(@numel,cc.PixelIdxList)); Ys(:) = false; Ys(cc.PixelIdxList{li}) = true;
  Ys = imfill(Ys,'holes');
  mat = V.T1w.mat * [diag(st) (1 - st(:)); 0 0 0 1];   % (voxels of the subsampled grid)
  [i,j,k] = ndgrid(1:size(Ys,1),1:size(Ys,2),1:size(Ys,3)); z = R.Rigid(3,:) * mat * [i(:) j(:) k(:) ones(numel(i),1)]';
  Ys(z < -78) = false;
  H  = maskSurface(Ys,mat,R.Rigid,5,40000);
end
% -------------------------------------------------------------------------
function fv = maskSurface(M,mat,Rig,sm,nf)
%maskSurface. Isosurface of the mask M (voxel-to-world matrix mat)
%  smoothed by sm voxels (odd) in head coordinates (Rig), reduced to at
%  most nf faces.
  fv = struct('vertices',zeros(0,3),'faces',zeros(0,3));
  if ~any(M(:)), return; end
  [i,j,k] = ind2sub(size(M),find(M)); b = [min(i) min(j) min(k)] - 3; e = [max(i) max(j) max(k)] + 3;
  b = max(1,b); e = min(size(M),e);
  Ms = smooth3(double(M(b(1):e(1),b(2):e(2),b(3):e(3))),'gaussian',sm);
  f  = isosurface(Ms,0.5); if isempty(f.vertices), return; end
  if size(f.faces,1) > nf, f = reducepatch(f,nf); end
  vox = f.vertices(:,[2 1 3]) + b - 1;
  fv  = struct('vertices',toHead(Rig * mat,vox),'faces',f.faces);
end
% -------------------------------------------------------------------------
function v = toHead(T,v)
%toHead. Points v (n x 3) transformed by the matrix T.
  v = (T * [v'; ones(1,size(v,1))])'; v = v(:,1:3);
end
% -------------------------------------------------------------------------
function N = vertexNormals(V,F)
%vertexNormals. Outward unit normals of the vertices V of the faces F
%  (area-weighted face normals, oriented away from the center).
  fn = cross(V(F(:,2),:) - V(F(:,1),:),V(F(:,3),:) - V(F(:,1),:),2);
  N  = zeros(size(V));
  for d = 1:3, N(:,d) = accumarray(F(:),repmat(fn(:,d),3,1),[size(V,1) 1]); end
  N  = N ./ max(eps,sqrt(sum(N.^2,2)));
  if mean(sum((V - mean(V,1)) .* N,2)) < 0, N = -N; end
end
% -------------------------------------------------------------------------
function A = atlasIds
%atlasIds. ROIids of lotter2025 (see atlasTable in cat_io_dcm2bids_pp_anat):
%  right and both hippocampi, lateral ventricles, cerebellum (GM and WM)
%  and brainstem, and the lobe index of the cortical regions (lut of
%  ROIid + 1: frontal, temporal, parietal, occipital, cingulate, insula).
  T = cat_io_dcm2bids_pp_anat('atlasTable','lotter2025_aparc'); id = @(k) T.id(ismember(T.key,k));
  A = struct('rhc',id({'rHIP'}),'hc',id({'lHIP','rHIP'}),'lv',id({'lLV','rLV'}), ...
    'cb',id({'lCB','rCB','lCBWM','rCBWM','bBS'}),'lut',zeros(max(T.id) + 1,1));
  lobes = {'frontal','temporal','parietal','occipital','cingulate','insula'};
  for li = 1:numel(lobes), A.lut(T.id(strcmp(T.lobe,lobes{li})) + 1) = li; end
end
% =========================================================================
% page
% =========================================================================
function F = newPage(G,L,R,M,I,variant)
%newPage. Figure of the A4 page (landscape for variant D) with a page
%  axes in cm (origin top-left) and the texts of the identity.
  if strcmpi(variant,'D'), pw = 29.7; ph = 21; else, pw = 21; ph = 29.7; end
  fh = figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 pw ph],'InvertHardcopy','off', ...
    'PaperUnits','centimeters','PaperSize',[pw ph],'PaperPosition',[0 0 pw ph]);
  ax = axes(fh,'Units','normalized','Position',[0 0 1 1],'XLim',[0 pw],'YLim',[0 ph],'YDir','reverse','Visible','off');
  hold(ax,'on'); rectangle(ax,'Position',[0 0 pw ph],'FaceColor','w','EdgeColor','none'); % (full page in the export)
  F = G; F.fh = fh; F.ax = ax; F.pw = pw; F.ph = ph; F.m = 1.2; F.L = L; F.M = M; F.I = I; F.R = R; F.variant = upper(variant);
  F.font = L.Font; F.rtl = strcmp(L.Direction,'rtl'); F.img = strcmp(L.Render,'image');
  F.jfont = javaFont(L.ImageFont); F.dpi = 400; F.wcache = containers.Map();
  sx = ''; if isfield(L.Sex,upper(R.id.Sex)), sx = L.Sex.(upper(R.id.Sex)); end
  F.id = struct('PID',R.id.PID,'Name',R.id.Name,'Birth',R.id.Birth,'Sex',sx,'Weight','','Scan',R.id.Scan,'Age','');
  if isnumeric(R.id.Weight) && isscalar(R.id.Weight) && R.id.Weight > 0, F.id.Weight = sprintf('%g kg',round(R.id.Weight)); end
  if isfinite(R.age), F.id.Age = sprintf('%0.1f',R.age); end
  for fn = fieldnames(F.id)', if isempty(F.id.(fn{1})), F.id.(fn{1}) = 'n/a'; end, end
  F.study = G.study; F.study.Institution = R.inst.Institution; F.study.Address = R.inst.Address;
  if isempty(F.study.Institution) && isempty(F.study.Address), F.study.Institution = 'n/a'; end
end
% =========================================================================
% variants
% =========================================================================
function variantA(F)
%variantA. Magazine rows (3 columns): brain and scanner text, scanner
%  illustration and images text, and the notes (the abbreviations are
%  defined in the texts), the large surface with its caption and the mean
%  thickness, the slices in a row with the QC rating, and the four band
%  plots in a row.
%  The rows below the surface are placed from the bottom upwards, i.e. the
%  surface gets the remaining height (longer languages).
  cw = (F.pw - 2*F.m - 2*0.6)/3; cx = F.m + (0:2)*(cw + 0.6); wf = F.pw - 2*F.m;
  yb = bottom(F); y = header(F,'plain');
  y1 = intro(F,cx(1),y,cw,{'brain','scanner'});
  scannerIcon(F,[cx(2) y cw 2.3]);
  y2 = intro(F,cx(2),y + 2.6,cw,{'images'});
  y3 = notes(F,cx(3),y,cw,1);
  y  = max([y1 y2 y3]) + 0.3;
  % band plots and slices from the bottom upwards
  yn = yb - para(F,0,0,wf,F.L.Captions.Norms,8,'gap',0,'draw',false) - 0.1;
  para(F,F.m,yn,wf,F.L.Captions.Norms,8);
  bw = (wf - 3*0.7)/4; ya = yn - 0.95 - 1.9;
  for i = 1:4, bandPlot(F,[F.m + 0.3 + (i-1)*(bw + 0.7) ya bw - 0.3 1.9],i,true); end
  ys = ya - 0.6 - para(F,0,0,wf,sliceCaption(F),8,'gap',0,'draw',false);
  para(F,F.m,ys,wf,sliceCaption(F),8);
  hl = min(4.3,max(3.2,(ys - y)/2.4)); ys = ys - 0.1 - hl - 0.5;
  slices(F,[F.m ys wf hl],'row',true);
  % surface with caption, colorbar, and mean thickness (if there is space)
  hs = min(7,max(2.5,ys - 0.35 - y));
  surfPanel(F,[cx(1) y 2*cw + 0.6 hs]);
  yc = para(F,cx(3),y + 0.4,cw,surfCaption(F),8);
  if F.I.hasct, colorbarCT(F,[cx(3) yc + 0.3 cw - 1.6 0.28]); yc = yc + 1.1; end
  if yc + 1.6 < ys, yc = bigValue(F,cx(3),yc + 0.2,3); end
  fits(F,max(yc,y + hs + 0.35),ys);
end
% -------------------------------------------------------------------------
function variantB(F)
%variantB. Sidebar (4 columns): the whole text in the first column (with
%  the scanner and the abbreviations), and the surface, slices, band
%  plots, the table of the values (with the QC rating and the caption of
%  the band plots), and the notes in the other three.
  sw = 5.0; x0 = F.m + sw + 0.6; w0 = F.pw - F.m - x0; cw = (w0 - 2*0.5)/3;
  yb = bottom(F); y = header(F,'side');
  rectangle(F.ax,'Position',[F.m - 0.2 y - 0.2 sw + 0.4 yb - y],'Curvature',[0.04 0.02],'FaceColor',F.bgl,'EdgeColor','none');
  yi = intro(F,F.m,y,sw,{'brain','scannerimg','scanner','images'},F.bgl);
  fits(F,glossary(F,F.m,yi + 0.1,sw,glossKeys(F,1),7,F.bgl),yb - 0.3);
  yn = yb - notes(F,0,0,w0,3,F.bg,false); notes(F,x0,yn,w0,3);
  sr = surfPanel(F,[x0 y 0.68*w0 6.6]); xc = sr(1) + sr(3) + 0.3; wc = F.pw - F.m - xc;
  yc = para(F,xc,y + 0.5,wc,surfCaption(F),8);
  if F.I.hasct, colorbarCT(F,[xc yc + 0.3 wc - 1.4 0.28]); yc = yc + 1.1; end
  y  = max(y + sr(4),yc) + 0.4;
  h  = slices(F,[x0 y w0 4.2],'row',false);
  y  = para(F,x0,y + h + 0.15,w0,sliceCaption(F),8) + 0.3;
  bw = (2*cw + 0.5 - 1.4)/2;
  for i = 1:4, bandPlot(F,[x0 + 0.7 + mod(i-1,2)*(bw + 0.7) y + 0.4 + floor((i-1)/2)*3.1 bw 2.2],i,i > 2); end
  yt = valuesTable(F,x0 + 2*(cw + 0.5),y,cw);
  yt = para(F,x0 + 2*(cw + 0.5),yt + 0.2,cw,F.L.Captions.Norms,7.5);
  fits(F,max(yt,y + 6.5),yn - 0.3);
end
% -------------------------------------------------------------------------
function variantC(F)
%variantC. Infographic: colored header band and numbered steps that
%  alternate between text and images, i.e. (1) brain text with the slices
%  (large sagittal) and the QC rating, (2) scanner text with the large
%  surface and callouts, and (3) measures with tiles of large numbers and
%  small band plots, with the notes in three columns at the bottom.
  cw = (F.pw - 2*F.m - 2*0.6)/3; cx = F.m + (0:2)*(cw + 0.6); w2 = 2*cw + 0.6; wf = F.pw - 2*F.m;
  yb = bottom(F); y = header(F,'band');
  yb = yb - notes(F,0,0,wf,3,F.bg,false); notes(F,F.m,yb,wf,3);
  % (1) brain text + slices
  badge(F,cx(1),y,1); y1 = intro(F,cx(1) + 0.75,y,cw - 0.75,{'brain'});
  h  = slices(F,[cx(2) y w2 4.4],'feature',true);
  y2 = para(F,cx(2),y + h + 0.1,w2,sliceCaption(F),8);
  y  = max(y1,y2) + 0.35;
  % (2) scanner text + surface with callouts
  sr = surfPanel(F,[cx(1) y 0.58*w2 5.2],F.I.co); xs = sr(1) + sr(3) + 0.3; ws = cx(1) + w2 - xs;
  ys = para(F,xs,y + 0.3,ws,surfCaption(F),8);
  if F.I.hasct, colorbarCT(F,[xs ys + 0.3 ws - 1.4 0.28]); ys = ys + 1.1; end
  ys = max(ys,sr(2) + sr(4));
  badge(F,cx(3),y,2); y2 = intro(F,cx(3) + 0.75,y,cw - 0.75,{'scanner'});
  scannerIcon(F,[cx(3) + 0.75 y2 cw - 0.75 1.8]);
  y  = max(ys,y2 + 2.0) + 0.3;
  % (3) measures text + values
  badge(F,cx(3),y,3); y3 = intro(F,cx(3) + 0.75,y,cw - 0.75,{'images'});
  para(F,cx(1),y,w2,F.L.Labels.Values,10,'fw','bold','col',F.acc);
  tw = (w2 - 0.4)/2;
  for i = 1:4, tile(F,[cx(1) + mod(i-1,2)*(tw + 0.4) y + 0.6 + floor((i-1)/2)*2.5 tw 2.3],i); end
  fits(F,max(y3,para(F,cx(1),y + 5.65,w2,F.L.Captions.Norms,8)),yb - 0.3);
end
% -------------------------------------------------------------------------
function variantD(F)
%variantD. Landscape (4 columns): one column per element, i.e. text,
%  slices with the scanner, surface with the table of the values (with
%  the QC rating), and band plots, with the notes in three columns at the
%  bottom.
  cw = (F.pw - 2*F.m - 3*0.6)/4; cx = F.m + (0:3)*(cw + 0.6); wf = F.pw - 2*F.m;
  yb = bottom(F); y = header(F,'wide');
  yb = yb - notes(F,0,0,wf,3,F.bg,false); notes(F,F.m,yb,wf,3); yb = yb - 0.3;
  fits(F,intro(F,cx(1),y,cw,{'brain','scanner','images'}),yb);
  h  = slices(F,[cx(2) y cw 9.0],'stack',false);
  yc = para(F,cx(2),y + h + 0.15,cw,sliceCaption(F),8);
  hi = min(0.42*cw,yb - yc - 0.3); if hi > 1.4, scannerIcon(F,[cx(2) yc + 0.2 cw hi]); end
  sr = surfPanel(F,[cx(3) y cw 5.6]); yc = sr(2) + sr(4) + 0.25;
  if F.I.hasct, colorbarCT(F,[cx(3) + 0.3 yc cw - 1.4 0.28]); yc = yc + 0.85; end
  yc = para(F,cx(3),yc,cw,surfCaption(F),8);
  fits(F,valuesTable(F,cx(3),yc + 0.3,cw),yb);
  hn = para(F,0,0,cw - 0.7,F.L.Captions.Norms,8,'gap',0,'draw',false);
  bh = min(2.2,(yb - y - 0.4 - 3*0.75 - 0.95 - hn)/4);
  for i = 1:4, bandPlot(F,[cx(4) + 0.7 y + 0.4 + (i-1)*(bh + 0.75) cw - 0.7 bh],i,i == 4); end
  fits(F,para(F,cx(4) + 0.7,y + 0.4 + 3*(bh + 0.75) + bh + 0.95,cw - 0.7,F.L.Captions.Norms,8),yb);
end
% -------------------------------------------------------------------------
function fits(F,y,yb)
%fits. Warning if the content (bottom y, with the gap below the last
%  paragraph) overlaps the next element (top yb, e.g. the footer).
  if y > yb + 0.3
    cat_io_cprintf('warn','  Participant report (%s, %s): content overlaps by %0.1f cm.\n',F.variant,F.L.Language,y - yb);
  end
end
% =========================================================================
% components
% =========================================================================
function y = header(F,style)
%header. Title, study logo, and identity of the participant, returns the
%  y below the header: 'plain' (logo above the identity card), 'side'
%  (logo beside the card), 'band' (title and logo in a colored band, the
%  identity in a row below), or 'wide' (landscape, identity in one row).
  L = F.L; m = F.m; hasl = ~isempty(F.Plogo);
  switch style
    case 'band'
      hb = 2.6; lw = 0;
      rectangle(F.ax,'Position',[0 0 F.pw hb],'FaceColor',F.fg,'EdgeColor','none');
      if hasl
        rectangle(F.ax,'Position',[F.pw - m - 2.0 0.45 2.0 1.7],'Curvature',[0.2 0.25],'FaceColor','w','EdgeColor','none');
        logo(F,F.Plogo,[F.pw - m - 1.85 0.55 1.7 1.5],'center',[1 1 1]); lw = 2.4;
      end
      ltext(F,m,1.05,L.Title,20,'fw','bold','col',[1 1 1],'bg',F.fg,'va','middle','maxw',F.pw - 2*m - lw);
      ltext(F,m,1.95,L.Subtitle,9.5,'col',mix(F,[1 1 1],0.8,F.fg),'bg',F.fg,'va','middle','maxw',F.pw - 2*m - lw);
      y = idrow(F,[m hb + 0.3 F.pw - 2*m 1.0],[]) + 0.45;
    case 'wide'
      lw = 0; if hasl, lw = 1.5; logo(F,F.Plogo,[F.pw - m - lw 0.6 lw lw],'right',F.bg); end
      r = [F.pw - m - lw - 0.3*hasl - 15.8 0.85 15.8 1.1];
      idrow(F,r,F.bgl);
      ltext(F,m,1.05,L.Title,19,'fw','bold','col',F.fg,'va','middle','maxw',r(1) - m - 0.4);
      ltext(F,m,1.85,L.Subtitle,9.5,'col',F.fgl,'va','middle','maxw',r(1) - m - 0.4);
      y = 3.0;
    case 'side'
      lw = 0; if hasl, lw = 1.6; logo(F,F.Plogo,[F.pw - m - lw 0.75 lw lw],'right',F.bg); end
      r = [F.pw - m - lw - 0.3*hasl - 8.6 0.75 8.6 1.6];
      idcard(F,r);
      ltext(F,m,1.25,L.Title,19,'fw','bold','col',F.fg,'va','middle','maxw',r(1) - m - 0.4);
      ltext(F,m,2.1,L.Subtitle,9.5,'col',F.fgl,'va','middle','maxw',r(1) - m - 0.4);
      y = r(2) + r(4) + 0.75;
    otherwise % plain
      if hasl, logo(F,F.Plogo,[F.pw - m - 1.2 0.45 1.2 1.2],'right',F.bg); end
      r = [F.pw - m - 8.6 1.85 8.6 1.6];
      idcard(F,r);
      ltext(F,m,1.6,L.Title,19,'fw','bold','col',F.fg,'va','middle','maxw',r(1) - m - 0.4);
      ltext(F,m,2.5,L.Subtitle,9.5,'col',F.fgl,'va','middle','maxw',r(1) - m - 0.4);
      y = r(2) + r(4) + 0.75;
  end
  plot(F.ax,[m F.pw - m],[y y] - 0.3,'-','Color',F.acc,'LineWidth',1.2);
end
% -------------------------------------------------------------------------
function idcard(F,r)
%idcard. Identity of the participant in a light card with two rows (ID,
%  name, sex | birth date, scan date, age, weight).
  rectangle(F.ax,'Position',r,'Curvature',[0.06 0.2],'FaceColor',F.bgl,'EdgeColor','none');
  fn = {'PID','Name','Sex','Birth','Scan','Age','Weight'};
  cx = [0 0.30 0.76 0 0.30 0.55 0.76]; cy = [0 0 0 1 1 1 1]; nx = [0.30 0.76 1 0.30 0.55 0.76 1];
  rh = (r(4) - 0.2)/2; w = r(3) - 0.4;
  for i = 1:numel(fn)
    x = r(1) + 0.25 + cx(i)*w; y = r(2) + 0.12 + cy(i)*rh; mw = (nx(i) - cx(i))*w - 0.15;
    ltext(F,x,y,F.L.Identity.(fn{i}),6.5,'col',F.fgl,'bg',F.bgl,'maxw',mw);
    shrink(text(F.ax,x,y + 0.3,F.id.(fn{i}),'FontName','Helvetica','FontSize',9.5,'FontWeight','bold','Color',F.fg, ...
      'VerticalAlignment','top','Interpreter','none'),mw);
  end
end
% -------------------------------------------------------------------------
function y = idrow(F,r,bg)
%idrow. Identity of the participant in one row (ID, name, birth date,
%  sex, weight, scan date, age) on the background bg (none for empty),
%  returns the bottom of the row.
  if ~isempty(bg), rectangle(F.ax,'Position',r,'Curvature',[0.04 0.2],'FaceColor',bg,'EdgeColor','none'); else, bg = F.bg; end
  fn = {'PID','Name','Birth','Sex','Weight','Scan','Age'}; rw = [1.0 2.0 1.15 0.85 0.8 1.15 0.75];
  x0 = r(1) + 0.25 + [0 cumsum(rw(1:end-1))]/sum(rw)*(r(3) - 0.4); dw = rw/sum(rw)*(r(3) - 0.4);
  for i = 1:numel(fn)
    ltext(F,x0(i),r(2) + 0.1,F.L.Identity.(fn{i}),6.5,'col',F.fgl,'bg',bg,'maxw',dw(i) - 0.2);
    shrink(text(F.ax,x0(i),r(2) + 0.42,F.id.(fn{i}),'FontName','Helvetica','FontSize',9.5,'FontWeight','bold','Color',F.fg, ...
      'VerticalAlignment','top','Interpreter','none'),dw(i) - 0.2);
    if i > 1, plot(F.ax,[x0(i) x0(i)] - 0.12,r(2) + [0.12 r(4) - 0.12],':','Color',F.fgl,'LineWidth',0.6); end
  end
  y = r(2) + r(4);
end
% -------------------------------------------------------------------------
function y0 = bottom(F)
%bottom. Footer of all variants: disclaimer (bold), study, organization,
%  website, and institute, and the CAT logo (bottom right), returns the y
%  above the footer.
  s = F.L.Footer;
  for fn = fieldnames(F.study)', s = strrep(s,['{' fn{1} '}'],F.study.(fn{1})); end
  s  = regexprep(regexprep(regexprep(s,'\(\s*,\s*','('),'\s*,\s*([\)\.])','$1'),'\s*\(\s*\)',''); % (empty fields)
  s  = sprintf('%s %s, %s.',s,cat_version,char(datetime('now','Format','yyyy-MM-dd')));
  lw = 1.75; w = F.pw - 2*F.m - lw - 0.4;
  hd = para(F,0,0,w,F.L.Intro.Disclaimer,6.5,'fw','bold','gap',0.04,'draw',false);
  hf = para(F,0,0,w,s,6,'gap',0,'draw',false);
  yf = F.ph - 0.5 - max(0.95,hd + hf) - 0.12;
  plot(F.ax,[F.m F.pw - F.m],[yf yf],'-','Color',F.fgl,'LineWidth',0.4);
  y  = para(F,F.m,yf + 0.12,w,F.L.Intro.Disclaimer,6.5,'fw','bold','col',F.fg,'gap',0.04);
  para(F,F.m,y,w,s,6,'col',F.fgl);
  logo(F,F.Pcat,[F.pw - F.m - lw yf + 0.14 lw 0.85],'right',F.bg);
  y0 = yf - 0.3;
end
% -------------------------------------------------------------------------
function y = notes(F,x,y,w,ncol,bg,draw)
%notes. "Good to know": image quality (QC rating), individuality of the
%  brain (diagnosis by physicians), and privacy (defacing) in ncol columns
%  on the background bg, returns the y below (without drawing for draw =
%  false, e.g. for the height).
  if nargin < 6, bg = F.bg; end
  if nargin < 7, draw = true; end
  N  = F.L.Notes; k = {'Quality','Individual','Privacy'}; g = 0.5; cw = (w - (ncol - 1)*g)/ncol;
  y  = para(F,x,y,w,N.Head,10,'fw','bold','col',F.acc,'gap',0.06,'bg',bg,'draw',draw); yr = y;
  for i = 1:numel(k)
    c = mod(i - 1,ncol); if c == 0, y = yr; end
    yi = para(F,x + c*(cw + g),y,cw,N.([k{i} 'Head']),8,'fw','bold','gap',0.03,'bg',bg,'draw',draw);
    yr = max(yr,para(F,x + c*(cw + g),yi,cw,N.(k{i}),8,'gap',0.15,'bg',bg,'draw',draw));
  end
  y = yr;
end
% -------------------------------------------------------------------------
function qcBadge(F,x,y)
%qcBadge. Overall quality rating of the T1w (grade in a circle with the
%  color of the marks, QC and percentage rating points) at x,y (top-left).
  Q = F.I.qc; if isempty(Q.grade), return; end
  tc = [1 1 1]; if Q.col * [0.299; 0.587; 0.114] > 0.6, tc = F.fg; end   % (dark text on light colors)
  rectangle(F.ax,'Position',[x y 0.5 0.5],'Curvature',[1 1],'FaceColor',Q.col,'EdgeColor','none');
  text(F.ax,x + 0.25,y + 0.26,Q.grade,'FontName','Helvetica','FontSize',7,'FontWeight','bold','Color',tc, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
  text(F.ax,x + 0.62,y + 0.25,sprintf('QC %0.0f%%',Q.rps),'FontName','Helvetica','FontSize',7,'FontWeight','bold', ...
    'Color',F.fg,'VerticalAlignment','middle');
end
% -------------------------------------------------------------------------
function y = intro(F,x,y,w,parts,bg)
%intro. Parts of the introduction text (flow layout on the background bg),
%  returns the y below.
  if nargin < 6, bg = F.bg; end
  I = F.L.Intro;
  for k = 1:numel(parts)
    switch parts{k}
      case 'brain',      y = para(F,x,y,w,I.BrainHead,10,'fw','bold','col',F.acc,'gap',0.06,'bg',bg); y = para(F,x,y,w,I.Brain,9,'gap',0.3,'bg',bg);
      case 'scanner',    y = para(F,x,y,w,I.ScannerHead,10,'fw','bold','col',F.acc,'gap',0.06,'bg',bg); y = para(F,x,y,w,I.Scanner,9,'gap',0.3,'bg',bg);
      case 'scannerimg', scannerIcon(F,[x y w 0.42*w]); y = y + 0.42*w + 0.15;
      case 'images',     y = para(F,x,y,w,I.ImagesHead,10,'fw','bold','col',F.acc,'gap',0.06,'bg',bg); y = para(F,x,y,w,I.Images,9,'gap',0.3,'bg',bg);
    end
  end
end
% -------------------------------------------------------------------------
function k = glossKeys(F,long)
%glossKeys. Abbreviations of the glossary (also of the measures for long),
%  DEC only with DEC map.
  k = {'MRI','T1w','DEC','GM','WM','CSF','HC','CT','QC'};
  if long, k = [k {'GMV','WMV','HCV','P'}]; end
  if isempty(F.I.qc.grade), k = setdiff(k,{'QC'},'stable'); end
  if ~F.I.hasdec, k = setdiff(k,{'DEC'},'stable'); end
end
% -------------------------------------------------------------------------
function y = glossary(F,x,y,w,keys,sz,bg)
%glossary. Abbreviations of the base layer with the definitions of the
%  language file (on the background bg), returns the y below.
  if nargin < 7, bg = F.bg; end
  y = para(F,x,y,w,F.L.Labels.Abbreviations,sz + 1,'fw','bold','col',F.acc,'gap',0.05,'bg',bg);
  kw = 0.95;
  for ki = 1:numel(keys)
    text(F.ax,x,y,keys{ki},'FontName','Helvetica','FontSize',sz,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
    y = para(F,x + kw,y,w - kw,F.L.Abbreviations.(keys{ki}),sz,'gap',0.06,'col',F.fg,'bg',bg);
  end
end
% -------------------------------------------------------------------------
function s = sliceCaption(F)
%sliceCaption. Caption of the slices (without DEC map).
  s = F.L.Captions.Slices; if ~F.I.hasdec, s = F.L.Captions.SlicesNoDEC; end
end
% -------------------------------------------------------------------------
function s = surfCaption(F)
%surfCaption. Caption of the surface (without thickness).
  s = F.L.Captions.Surface; if ~F.I.hasct, s = F.L.Captions.SurfaceNoCT; end
end
% -------------------------------------------------------------------------
function y = valuesTable(F,x,y,w)
%valuesTable. Table of the values with the typical range (25th to 75th
%  percentile) and the percentile, returns the y below.
  y = para(F,x,y,w,F.L.Labels.Values,10,'fw','bold','col',F.acc,'gap',0.15);
  c = [x x + 0.25*w x + 0.56*w x + w]; hdr = {F.L.Labels.Measure,F.L.Labels.Value,F.L.Labels.Typical};
  hh = 0.3;
  for j = 1:3, [~,h] = ltext(F,c(j),y,hdr{j},7,'col',F.fgl,'maxw',c(j+1) - c(j) - 0.1); hh = max(hh,h); end
  text(F.ax,c(4),y,'P','FontName','Helvetica','FontSize',7,'Color',F.fgl,'VerticalAlignment','top','HorizontalAlignment','right');
  y = y + hh + 0.1; plot(F.ax,[x x + w],[y y] - 0.05,'-','Color',F.fgl,'LineWidth',0.4);
  for i = 1:numel(F.M)
    M = F.M(i); mu = M.mu(M.age);
    t = {M.name,[fmtValue(M.fmt,M.val) ' ' M.unit],'',fmtP(M.P,'%d')};
    if isfinite(mu), t{3} = [sprintf(M.fmt,mu - 0.674*M.sd) '-' sprintf(M.fmt,mu + 0.674*M.sd)]; end
    for j = 1:4
      text(F.ax,c(j),y + 0.05,t{j},'FontName','Helvetica','FontSize',7.5,'FontWeight',pick(j==1 || j==4,'bold','normal'), ...
        'Color',F.fg,'VerticalAlignment','top','HorizontalAlignment',pick(j==4,'right','left'));
    end
    y = y + 0.48;
  end
  Q = F.I.qc;
  if ~isempty(Q.grade) % (overall quality rating of the T1w)
    text(F.ax,c(1),y + 0.05,'QC','FontName','Helvetica','FontSize',7.5,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
    plot(F.ax,c(2) + 0.12,y + 0.22,'o','MarkerSize',6,'MarkerFaceColor',Q.col,'MarkerEdgeColor','none');
    text(F.ax,c(2) + 0.35,y + 0.05,sprintf('%s (%0.0f%%)',Q.grade,Q.rps),'FontName','Helvetica','FontSize',7.5, ...
      'Color',F.fg,'VerticalAlignment','top');
    y = y + 0.48;
  end
  y = para(F,x,y + 0.1,w,F.L.Labels.Percentile,7,'col',F.fgl);
end
% -------------------------------------------------------------------------
function badge(F,x,y,n)
%badge. Numbered circle of a step (variant C).
  rectangle(F.ax,'Position',[x y 0.55 0.55],'Curvature',[1 1],'FaceColor',F.acc,'EdgeColor','none');
  text(F.ax,x + 0.275,y + 0.29,sprintf('%d',n),'FontName','Helvetica','FontSize',10,'FontWeight','bold','Color','w', ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
end
% -------------------------------------------------------------------------
function scannerIcon(F,r)
%scannerIcon. Illustration of an MRI scanner with a person (front view of
%  the gantry with the table and the person in side view) with the field
%  strength.
  ax = rax(F,r); axis(ax,'off'); set(ax,'XLim',[0 9],'YLim',[0 4.2]); daspect(ax,[1 1 1]);
  rectangle(ax,'Position',[0.2 0.1 4.0 4.0],'Curvature',[0.35 0.35],'FaceColor',F.bgl,'EdgeColor',F.fgl,'LineWidth',0.8);
  rectangle(ax,'Position',[0.95 0.85 2.5 2.5],'Curvature',[1 1],'FaceColor','w','EdgeColor',F.fgl,'LineWidth',0.8);
  rectangle(ax,'Position',[1.25 1.15 1.9 1.9],'Curvature',[1 1],'FaceColor','none','EdgeColor',mix(F,F.disc(2,:),0.5),'LineWidth',0.6,'LineStyle',':');
  rectangle(ax,'Position',[1.7 1.45 7.0 0.32],'Curvature',[0.3 1],'FaceColor',F.fgl,'EdgeColor','none');   % table
  patch(ax,[7.4 8.2 8.0 7.6],[1.45 1.45 0.1 0.1],F.fgl,'EdgeColor','none');                               % stand
  rectangle(ax,'Position',[1.85 1.78 0.62 0.62],'Curvature',[1 1],'FaceColor',F.acc,'EdgeColor','none');    % head
  rectangle(ax,'Position',[2.5 1.8 3.6 0.62],'Curvature',[0.5 1],'FaceColor',F.acc,'EdgeColor','none');     % body
  rectangle(ax,'Position',[6.0 1.8 1.7 0.42],'Curvature',[0.6 1],'FaceColor',F.acc,'EdgeColor','none');     % legs
  fs = 3; if isfinite(F.R.field), fs = F.R.field; end
  text(ax,2.2,3.75,sprintf('%g T',round(fs*10)/10),'FontName','Helvetica','FontSize',7.5,'FontWeight','bold','Color',F.fg,'HorizontalAlignment','center');
  text(ax,6.2,3.3,'MRI','FontName','Helvetica','FontSize',8,'FontWeight','bold','Color',F.fgl,'HorizontalAlignment','center');
end
% -------------------------------------------------------------------------
function r = surfPanel(F,r,co)
%surfPanel. Surface rendering fitted into the rectangle r (returns the
%  used rectangle) with optional callouts co = {[px py],[lx ly],label; ...}
%  in image fractions.
  im = F.I.surf; ar = size(im,2)/size(im,1); w = min(r(3),r(4)*ar); h = w/ar;
  r  = [r(1) + (r(3) - w)/2 r(2) w h];
  ax = rax(F,r); image(ax,im); axis(ax,'off');
  set(ax,'YDir','reverse','XLim',[0.5 size(im,2) + 0.5],'YLim',[0.5 size(im,1) + 0.5]);
  if nargin > 2
    sz = [size(im,2) size(im,1)];
    for ci = 1:size(co,1)
      p = co{ci,1}.*sz; q = co{ci,2}.*sz;
      plot(ax,[p(1) q(1)],[p(2) q(2)],'-','Color',F.fg,'LineWidth',0.7);
      plot(ax,p(1),p(2),'o','MarkerSize',3.5,'MarkerFaceColor',F.fg,'MarkerEdgeColor','w');
      text(ax,q(1),q(2),co{ci,3},'FontName','Helvetica','FontSize',7.5,'FontWeight','bold','Color',F.fg, ...
        'BackgroundColor','w','EdgeColor',F.fgl,'Margin',1.5,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
  end
end
% -------------------------------------------------------------------------
function colorbarCT(F,r)
%colorbarCT. Colorbar of the cortical thickness (2 to 4 mm).
  ax = rax(F,r); image(ax,[2 4],[0 1],reshape(F.cm,[1 size(F.cm)])); set(ax,'YTick',[],'XTick',2:0.5:4, ...
    'XTickLabel',{'2','','3','','4'},'FontName','Helvetica','FontSize',7,'TickLength',[0 0],'XColor',F.fgl,'Box','off','YDir','normal', ...
    'XLim',[2 4],'YLim',[0 1]);
  text(ax,4.15,0.5,'CT (mm)','FontName','Helvetica','FontSize',7,'FontWeight','bold','Color',F.fg,'VerticalAlignment','middle');
end
% -------------------------------------------------------------------------
function h = slices(F,r,arr,qc)
%slices. Slices (see sliceImages) in the rectangle r: 'row' (sag, cor,
%  axi), 'stack' (sag over cor and axi), or 'feature' (large sagittal
%  with labels, cor and axi stacked right), with the QC rating below for
%  qc (see qcBadge), returns the used height.
  ar = [225 190 190]./[175 175 225]; g = 0.25; o = {'sag','cor','axi'};
  switch arr
    case 'row'
      hh = min(r(4),(r(3) - 2*g)/sum(ar)); x = r(1) + (r(3) - hh*sum(ar) - 2*g)/2;
      for i = 1:3, sliceAx(F,[x r(2) hh*ar(i) hh],o{i},i == 2); x = x + hh*ar(i) + g; end
      decLegend(F,[r(1) + r(3) - 3.0 r(2) + hh + 0.05 3.0 0.45]); h = hh + 0.5;
      if qc, qcBadge(F,x - sum(ar)*hh - 3*g,r(2) + hh + 0.03); end
    case 'stack'
      w1 = r(3); h1 = w1/ar(1); h2 = (r(3) - g)/(ar(2) + ar(3));
      sliceAx(F,[r(1) r(2) w1 h1],'sag',false);
      sliceAx(F,[r(1) r(2) + h1 + g h2*ar(2) h2],'cor',true);
      sliceAx(F,[r(1) + h2*ar(2) + g r(2) + h1 + g h2*ar(3) h2],'axi',false);
      decLegend(F,[r(1) + r(3) - 3.0 r(2) + h1 + h2 + g + 0.05 3.0 0.45]); h = h1 + h2 + g + 0.5;
    case 'feature'
      h1 = r(4); w1 = h1*ar(1); h2 = (h1 - g)/2;
      if w1 + g + h2*ar(2) > r(3), h1 = r(3)/(ar(1) + ar(2)/2); w1 = h1*ar(1); h2 = (h1 - g)/2; end
      sliceAx(F,[r(1) r(2) w1 h1],'sag',true);
      sliceAx(F,[r(1) + w1 + g r(2) h2*ar(2) h2],'cor',false);
      sliceAx(F,[r(1) + w1 + g r(2) + h2 + g h2*ar(3) h2],'axi',false);
      decLegend(F,[r(1) + w1 + g r(2) + h1 + 0.05 3.0 0.45]); h = h1 + 0.5;
      if qc, qcBadge(F,r(1),r(2) + h1 + 0.03); end
  end
end
% -------------------------------------------------------------------------
function sliceAx(F,r,ori,labels)
%sliceAx. One slice with the orientation letters, the hippocampus, and
%  optional tissue labels (GM, WM, CSF).
  S  = F.I.sl.(ori); ax = rax(F,r); xd = S.xd; yd = S.yd; W = diff(sort(xd)); H = diff(sort(yd));
  image(ax,'XData',xd,'YData',yd,'CData',S.rgb); set(ax,'YDir','normal'); axis(ax,'off');
  set(ax,'XLim',sort(xd),'YLim',sort(yd));
  switch ori, case 'sag', lab = {'P','A','S'}; case 'cor', lab = {'L','R','S'}; otherwise, lab = {'L','R','A'}; end
  t = @(x,y,s,ha,va) text(ax,x,y,s,'FontName','Helvetica','FontSize',7,'Color',[0.85 0.85 0.85],'HorizontalAlignment',ha,'VerticalAlignment',va);
  t(min(xd) + 4,mean(yd),lab{1},'left','middle'); t(max(xd) - 4,mean(yd),lab{2},'right','middle');
  t(mean(xd),max(yd) - 3,lab{3},'center','top');
  if ~isempty(S.hc), lbl(F,ax,S.hc,[S.hc(1) + 0.22*W S.hc(2) - 0.12*H],'HC',F.hccol); end
  if labels
    x0 = min(xd); z0 = min(yd);
    if ~isempty(S.gm),  lbl(F,ax,S.gm,[x0 + 0.08*W z0 + 0.20*H],'GM',[1 1 1]); end
    if ~isempty(S.wm),  lbl(F,ax,S.wm,[x0 + 0.24*W z0 + 0.88*H],'WM',[1 1 1]); end
    if ~isempty(S.csf), lbl(F,ax,S.csf,[x0 + 0.63*W z0 + 0.88*H],'CSF',[1 1 1]); end
  end
end
% -------------------------------------------------------------------------
function lbl(~,ax,p,q,s,col)
%lbl. Label s at q with a line to p (axes data).
  plot(ax,[p(1) q(1)],[p(2) q(2)],'-','Color',col,'LineWidth',0.5);
  text(ax,q(1),q(2),s,'FontName','Helvetica','FontSize',7,'FontWeight','bold','Color',col, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0 0 0],'Margin',0.5);
end
% -------------------------------------------------------------------------
function decLegend(F,r)
%decLegend. Fibre directions of the DEC colors (L-R red, A-P green, S-I
%  blue, only with DEC map) and the hatch of the removed face.
  ax = rax(F,r); axis(ax,'off'); set(ax,'XLim',[0 3.3],'YLim',[0 1]);
  if F.I.hasdec
    c = [0.90 0.15 0.15; 0.15 0.80 0.25; 0.20 0.35 0.95]; s = {'LR','AP','SI'};
    for i = 1:3
      plot(ax,(i-1)*0.8 + [0.05 0.35],[0.5 0.5],'-','Color',c(i,:),'LineWidth',2.5);
      text(ax,(i-1)*0.8 + 0.4,0.5,s{i},'FontName','Helvetica','FontSize',6.5,'Color',F.fg,'VerticalAlignment','middle');
    end
  end
  patch(ax,[2.5 2.85 2.85 2.5],[0.25 0.25 0.75 0.75],[0.1 0.1 0.1],'EdgeColor','none');
  for k = 0:3, plot(ax,2.5 + [k*0.12 min(0.35,k*0.12 + 0.35)],[0.25 0.25 + min(0.5,(0.35 - k*0.12)/0.35*0.5)],'-','Color',[0.6 0.6 0.6],'LineWidth',0.5); end
end
% -------------------------------------------------------------------------
function bandPlot(F,r,i,showx,bg)
%bandPlot. Normative band (5-95% and 25-75%, median) of the measure i by
%  age with the value of the participant and its percentile (on the
%  background bg, e.g. of the tiles).
  if nargin < 5, bg = F.bg; end
  M  = F.M(i); a = linspace(min(M.range(1),M.age - 5),max(M.range(2),M.age + 5),150); a = a(isfinite(a));
  if isempty(a), a = linspace(M.range(1),M.range(2),150); end
  mu = M.mu(a); z = [-1.645 -0.674 0.674 1.645];
  ax = rax(F,r);
  fill(ax,[a fliplr(a)],[mu + z(1)*M.sd fliplr(mu + z(4)*M.sd)],mix(F,F.acc,0.16),'EdgeColor','none');
  fill(ax,[a fliplr(a)],[mu + z(2)*M.sd fliplr(mu + z(3)*M.sd)],mix(F,F.acc,0.34),'EdgeColor','none');
  plot(ax,a,mu,'-','Color',mix(F,F.acc,0.8),'LineWidth',0.8);
  if isfinite(M.val) && isfinite(M.age)
    plot(ax,M.age,M.val,'o','MarkerSize',5,'MarkerFaceColor',F.fg,'MarkerEdgeColor','w','LineWidth',0.8);
    text(ax,M.age + 0.04*diff(a([1 end])),M.val,fmtP(M.P,'P%d'),'FontName','Helvetica','FontSize',7,'FontWeight','bold', ...
      'Color',F.fg,'VerticalAlignment','middle');
  end
  yl = [min(mu + z(1)*M.sd) max(mu + z(4)*M.sd)]; if isfinite(M.val), yl = [min(yl(1),M.val) max(yl(2),M.val)]; end
  set(ax,'XLim',a([1 end]),'XTick',20:20:80,'YLim',yl + [-0.05 0.05]*diff(yl),'FontName','Helvetica','FontSize',7, ...
    'TickDir','out','TickLength',[0.02 0.02],'XColor',F.fgl,'YColor',F.fgl,'Box','off','Color','none','LineWidth',0.4);
  if ~showx, set(ax,'XTickLabel',{}); end
  text(ax,0,1.04,sprintf('%s (%s)',M.name,M.unit),'Units','normalized','FontName','Helvetica','FontSize',7.5,'FontWeight','bold', ...
    'Color',F.fg,'VerticalAlignment','bottom');
  if showx, ltext(F,r(1) + r(3)/2,r(2) + r(4) + 0.38,F.L.Labels.Age,7,'col',F.fgl,'bg',bg,'ha','center'); end
end
% -------------------------------------------------------------------------
function tile(F,r,i)
%tile. Large number of a measure with its percentile and a small band
%  plot (variant C).
  M = F.M(i);
  rectangle(F.ax,'Position',r,'Curvature',[0.08 0.12],'FaceColor',F.bgl,'EdgeColor','none');
  text(F.ax,r(1) + 0.25,r(2) + 0.2,M.name,'FontName','Helvetica','FontSize',8,'FontWeight','bold','Color',F.acc,'VerticalAlignment','top');
  text(F.ax,r(1) + 0.25,r(2) + 0.65,fmtValue(M.fmt,M.val),'FontName','Helvetica','FontSize',18,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
  text(F.ax,r(1) + 0.25,r(2) + 1.5,sprintf('%s | %s',M.unit,fmtP(M.P,'P%d')),'FontName','Helvetica','FontSize',7.5,'Color',F.fgl,'VerticalAlignment','top');
  bandPlot(F,[r(1) + 0.45*r(3) r(2) + 0.4 0.52*r(3) r(4) - 1.05],i,true,F.bgl);
end
% -------------------------------------------------------------------------
function y = bigValue(F,x,y,i)
%bigValue. Large number of the measure i with its unit and percentile.
  M = F.M(i);
  text(F.ax,x,y,M.name,'FontName','Helvetica','FontSize',8,'FontWeight','bold','Color',F.acc,'VerticalAlignment','top');
  text(F.ax,x,y + 0.45,fmtValue(M.fmt,M.val),'FontName','Helvetica','FontSize',20,'FontWeight','bold','Color',F.fg,'VerticalAlignment','top');
  text(F.ax,x + 2.3,y + 0.75,sprintf('%s | %s',M.unit,fmtP(M.P,'P%d')),'FontName','Helvetica','FontSize',8,'Color',F.fgl,'VerticalAlignment','top');
  y = y + 1.4;
end
% -------------------------------------------------------------------------
function logo(F,P,r,align,bg)
%logo. Image file P fitted into the rectangle r (aligned left, center, or
%  right) on the background color bg (blended transparency, as the PDF of
%  the painters renderer has no image transparency).
  if isempty(P) || ~exist(P,'file'), return; end
  [im,map,a] = imread(P); if ~isempty(map), im = ind2rgb(im,map); end
  if isinteger(im), im = double(im) / double(intmax(class(im))); end
  if isempty(a), a = ones(size(im,1),size(im,2)); elseif isinteger(a), a = double(a) / double(intmax(class(a))); end
  if size(im,3) == 1, im = repmat(im,[1 1 3]); end
  im = im .* a + (1 - a) .* reshape(bg,1,1,3);
  ar = size(im,2)/size(im,1); h = r(4); w = h*ar; if w > r(3), w = r(3); h = w/ar; end
  switch align, case 'right', x = r(1) + r(3) - w; case 'center', x = r(1) + (r(3) - w)/2; otherwise, x = r(1); end
  ax = rax(F,[x r(2) + (r(4) - h)/2 w h]); image(ax,min(1,max(0,im))); axis(ax,'off');
  set(ax,'YDir','reverse','XLim',[0.5 size(im,2) + 0.5],'YLim',[0.5 size(im,1) + 0.5]);
end
% =========================================================================
% text of the language layer
% =========================================================================
function y = para(F,x,y,w,s,sz,varargin)
%para. Paragraph s wrapped to the width w (cm) at x,y (top) with the font
%  size sz (options fw, col, gap, bg, draw), returns the y below the
%  paragraph: MATLAB text (wrapped by measured words, CJK characters), or
%  Java2D text image (Render "image", see javaText).
  o = struct('fw','normal','col',F.fg,'gap',0.15,'bg',F.bg,'draw',true);
  for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
  if isempty(s), return; end
  if F.img
    [im,~,h] = javaText(F,s,sz,w,o.col,o.bg,strcmp(o.fw,'bold'));
    if o.draw, placeImage(F,[x y w h],im); end
    y = y + h + o.gap; return
  end
  lines = wrap(F,s,w,sz,o.fw);
  if F.rtl, xx = x + w; ha = 'right'; else, xx = x; ha = 'left'; end
  t = text(F.ax,xx,y,lines,'FontName',F.font,'FontSize',sz,'FontWeight',o.fw,'Color',o.col, ...
    'VerticalAlignment','top','HorizontalAlignment',ha,'Interpreter','none');
  e = get(t,'Extent'); y = y + e(4) + o.gap;
  if ~o.draw, delete(t); end
end
% -------------------------------------------------------------------------
function [w,h] = ltext(F,x,y,s,sz,varargin)
%ltext. One line of the language layer at x,y (options fw, col, bg, ha,
%  va, maxw for a smaller font, at least 5 pt, or a smaller image), returns
%  its width and height (cm).
  o = struct('fw','normal','col',F.fg,'bg',F.bg,'ha','left','va','top','maxw',Inf);
  for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
  w = 0; h = 0; if isempty(s), return; end
  if F.img
    [im,w,h] = javaText(F,s,sz,Inf,o.col,o.bg,strcmp(o.fw,'bold'));
    if w > o.maxw, h = h*o.maxw/w; w = o.maxw; end
    switch o.ha, case 'right', x = x - w; case 'center', x = x - w/2; end
    switch o.va, case 'middle', y = y - h/2; case 'bottom', y = y - h; end
    placeImage(F,[x y w h],im); return
  end
  t = text(F.ax,x,y,s,'FontName',F.font,'FontSize',sz,'FontWeight',o.fw,'Color',o.col, ...
    'HorizontalAlignment',o.ha,'VerticalAlignment',o.va,'Interpreter','none');
  e = shrink(t,o.maxw); w = e(3); h = e(4);
end
% -------------------------------------------------------------------------
function e = shrink(t,maxw)
%shrink. Smaller font of the text t (at least 5 pt) for the width maxw
%  (cm of the page axes), returns its extent.
  e = get(t,'Extent'); sz = get(t,'FontSize');
  while e(3) > maxw && sz > 5, sz = sz - 0.5; set(t,'FontSize',sz); e = get(t,'Extent'); end
end
% -------------------------------------------------------------------------
function lines = wrap(F,s,w,sz,fw)
%wrap. Greedy line breaks of the text s for the width w (cm) by the
%  measured widths of the words (CJK: characters).
  if any(s > 11903 & s < 40960)
    tok = regexp(s,'[A-Za-z0-9.,:;()%\-]+|\S','match'); sep = '';
  else
    tok = strsplit(strtrim(s)); sep = ' ';
  end
  sw = 0; if ~isempty(sep), sw = twidth(F,'x x',sz,fw) - twidth(F,'xx',sz,fw); end
  lines = {}; cur = ''; cw = 0;
  for i = 1:numel(tok)
    tw = twidth(F,tok{i},sz,fw);
    if isempty(cur), cur = tok{i}; cw = tw;
    elseif cw + sw + tw <= w, cur = [cur sep tok{i}]; cw = cw + sw + tw; %#ok<AGROW>
    else, lines{end+1} = cur; cur = tok{i}; cw = tw; %#ok<AGROW>
    end
  end
  if ~isempty(cur), lines{end+1} = cur; end
end
% -------------------------------------------------------------------------
function w = twidth(F,s,sz,fw)
%twidth. Width of the text s in cm (cached).
  key = sprintf('%s|%g|%s|%s',F.font,sz,fw,s);
  if isKey(F.wcache,key), w = F.wcache(key); return; end
  t = text(F.ax,0,0,s,'FontName',F.font,'FontSize',sz,'FontWeight',fw,'Interpreter','none');
  e = get(t,'Extent'); delete(t); w = e(3); F.wcache(key) = w; % (handle)
end
% -------------------------------------------------------------------------
function [im,w,h] = javaText(F,s,sz,wcm,col,bg,bold)
%javaText. Text s rendered by Java2D (TextLayout with bidi and shaping,
%  line breaks of the LineBreakMeasurer for the width wcm in cm, one line
%  for Inf) with F.dpi in the color col on the background bg (RGB image,
%  width w and height h in cm).
  st = java.awt.Font.PLAIN; if bold, st = java.awt.Font.BOLD; end
  f  = java.awt.Font(F.jfont,st,round(sz/72*F.dpi));
  b0 = java.awt.image.BufferedImage(1,1,java.awt.image.BufferedImage.TYPE_INT_ARGB); g0 = b0.createGraphics(); hints(g0);
  as = java.text.AttributedString(s); as.addAttribute(java.awt.font.TextAttribute.FONT,f);
  if F.rtl, as.addAttribute(java.awt.font.TextAttribute.RUN_DIRECTION,java.awt.font.TextAttribute.RUN_DIRECTION_RTL); end
  it = as.getIterator(); frc = g0.getFontRenderContext(); tl = {};
  if isinf(wcm)
    tl{1} = java.awt.font.TextLayout(it,frc); wpx = ceil(tl{1}.getAdvance()) + 4;
  else
    wpx = round(wcm/2.54*F.dpi); lbm = java.awt.font.LineBreakMeasurer(it,frc);
    while lbm.getPosition() < it.getEndIndex(), tl{end+1} = lbm.nextLayout(single(wpx)); end %#ok<AGROW>
  end
  g0.dispose();
  hpx = ceil(sum(cellfun(@(l) l.getAscent() + l.getDescent() + l.getLeading(),tl))) + 4;
  bi = java.awt.image.BufferedImage(wpx,hpx,java.awt.image.BufferedImage.TYPE_INT_ARGB); g = bi.createGraphics(); hints(g);
  c  = int32(round(255*col)); g.setColor(java.awt.Color(c(1),c(2),c(3)));
  yy = 0;
  for li = 1:numel(tl)
    yy = yy + tl{li}.getAscent(); xx = 0; if F.rtl, xx = wpx - tl{li}.getAdvance(); end
    tl{li}.draw(g,single(xx),single(yy)); yy = yy + tl{li}.getDescent() + tl{li}.getLeading();
  end
  g.dispose();
  px = bi.getRGB(0,0,wpx,hpx,[],0,wpx); px = reshape(typecast(int32(px),'uint32'),wpx,hpx)';
  a  = double(bitshift(px,-24))/255; im = zeros(hpx,wpx,3);
  for ci = 1:3, im(:,:,ci) = a*col(ci) + (1 - a)*bg(ci); end
  w  = wpx/F.dpi*2.54; h = hpx/F.dpi*2.54;
end
% -------------------------------------------------------------------------
function hints(g)
%hints. Antialiasing and fractional metrics of the Java2D text.
  g.setRenderingHint(java.awt.RenderingHints.KEY_TEXT_ANTIALIASING,java.awt.RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
  g.setRenderingHint(java.awt.RenderingHints.KEY_FRACTIONALMETRICS,java.awt.RenderingHints.VALUE_FRACTIONALMETRICS_ON);
end
% -------------------------------------------------------------------------
function f = javaFont(name)
%javaFont. Java font family name (the logical font SansSerif with the
%  fallback fonts of the system if not installed).
  persistent fam
  f = 'SansSerif'; if isempty(name), return; end
  if isempty(fam), fam = cell(java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getAvailableFontFamilyNames()); end
  if any(strcmpi(fam,name)), f = name; end
end
% -------------------------------------------------------------------------
function placeImage(F,r,im)
%placeImage. RGB image im in the rectangle r (cm, top-left).
  ax = rax(F,r); image(ax,im); axis(ax,'off');
  set(ax,'YDir','reverse','XLim',[0.5 size(im,2) + 0.5],'YLim',[0.5 size(im,1) + 0.5]);
end
% =========================================================================
% helpers
% =========================================================================
function ax = rax(F,r)
%rax. Axes at the rectangle r = [x y w h] in cm from the top-left corner.
  ax = axes(F.fh,'Units','centimeters','Position',[r(1) F.ph - r(2) - r(4) r(3) r(4)]); hold(ax,'on');
end
% -------------------------------------------------------------------------
function c = mix(~,c,a,b)
%mix. Color c with the opacity a on the background b (default white).
  if nargin < 4, b = [1 1 1]; end
  c = a*c + (1 - a)*b;
end
% -------------------------------------------------------------------------
function v = pick(c,a,b)
%pick. a if c else b.
  if c, v = a; else, v = b; end
end
% -------------------------------------------------------------------------
function s = fmtValue(fmt,v)
%fmtValue. Value v by the format fmt, n/a for NaN.
  if isfinite(v), s = sprintf(fmt,v); else, s = 'n/a'; end
end
% -------------------------------------------------------------------------
function s = fmtP(P,fmt)
%fmtP. Percentile P by the format fmt, P n/a for NaN.
  if isfinite(P), s = sprintf(fmt,P); else, s = 'P n/a'; end
end
% -------------------------------------------------------------------------
function a = ageOf(birth,scan)
%ageOf. Age in years between the dates birth and scan (yyyy-mm-dd, NaN).
  a = NaN;
  try
    a = days(datetime(scan(1:10),'InputFormat','yyyy-MM-dd') - datetime(birth(1:10),'InputFormat','yyyy-MM-dd')) / 365.2425;
  catch
  end
end
% -------------------------------------------------------------------------
function v = fieldOf(S,fn,default)
%fieldOf. Field fn of the structure S, otherwise default (also if empty).
  v = default; if isfield(S,fn) && ~isempty(S.(fn)), v = S.(fn); end
end
% -------------------------------------------------------------------------
function P = existing(P)
%existing. First existing file of the cell P, otherwise ''.
  P = P(cellfun(@(f) exist(f,'file') == 2,P)); if isempty(P), P = ''; else, P = P{1}; end
end
% -------------------------------------------------------------------------
function P = first(P)
%first. First non-empty file of the cell P, otherwise ''.
  P = P(~cellfun('isempty',P)); if isempty(P), P = ''; else, P = P{1}; end
end
