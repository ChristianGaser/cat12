function [ok,res] = cat_DCM2BIDS_selftest(Pdata,Pout,tests,levels)
%cat_DCM2BIDS_selftest. Self-test of cat_io_dcm2bids with real DICOM data.
%  Regression tests of the DICOM import and BIDS export (with a fake
%  dcm2niix that only logs its calls to detect unnecessary conversions)
%  and the evaluation of the preprocessing levels (opts.pplevel) by the
%  processing time and the accuracy of the measures in relation to level 1.
%  The tests use the DZPG pilot data, i.e. the data directory Pdata with
%  DCM/DCM_pilots/<site> and protocols (DZPG3T and sites.json). The test
%  data are selected by their structure (and not by names):
%    site_jena_JE7                 .. the first two subjects with sessions
%                                     <subject>_MR<n>_<date> (both sessions
%                                     of the first subject, the first
%                                     session of the second subject)
%    site_munich_humanPhantom2023  .. flat DICOM folder (DICOMDIR)
%    site_tuebingen_pilot202402    .. DWI series split into two directories,
%                                     copied field-map series, T1w series
%                                     under a path with special characters,
%                                     and the derived ADC series
%    site_berlin_pilot202511       .. session of the preprocessing levels
%                                     and the DICOM export (classic DICOM)
%  All outputs are written to Pout, i.e. the data directory is not changed.
%  The fake dcm2niix is a shell script (macOS/Linux).
%
%  [ok,res] = cat_DCM2BIDS_selftest(Pdata,Pout,tests,levels)
%
%  Pdata  .. data directory (default: selection dialog)
%  Pout   .. output directory (default: <tempdir>/cat_DCM2BIDS_selftest)
%  tests  .. test groups (default: all)
%            'import'    - DICOM import, export of the database, repeated
%                          import, JSON-only output followed by the images
%            'selection' - export of selected subjects, sessions, scans and
%                          DICOM directories, rerun, and expected errors
%            'special'   - stored settings, double copies, symbolic links,
%                          relative paths, and dcm2niix options
%            'inputs'    - split and copied series, flat DICOM folder,
%                          protocol tables, output 0, and UTF-8 paths
%            'phenotypic'- phenotypic data tables (matching to the sessions,
%                          sessions files, conflicts, key columns by number,
%                          no names in BIDS, persistence, corrected tables)
%            'export'    - anonymized DICOM export (zip files, tables,
%                          reimport, repeated run, protocol scans of a
%                          JSON-only import, cleaning, mosaic defacing)
%                          and the anonymization level 2 (BIDS labels,
%                          no personal data, DICOM export and reimport)
%            'pplevels'  - preprocessing levels (time and accuracy, coregistration, DTI mask)
%            'parallel'  - subject packages by parallel processes (nproc 2)
%                          with the same outputs as one process (with QC,
%                          preprocessing, and rendering, about 30 min)
%  levels .. preprocessing levels of 'pplevels' (default: 1:0.5:4), where
%            level 1 is the reference (its results are kept in Pout and
%            only processed once)
%  ok     .. true if all tests passed
%  res    .. structure array of the tests (group, name, ok, info)
%
%  The groups 'selection' and 'special' use the import of 'import' (that
%  runs before them if it is not part of tests and its output is missing).
%
%  See also cat_io_dcm2bids, cat_tst_unittests.

  if ~exist('Pdata','var') || isempty(Pdata)
    Pdata = spm_select(1,'dir','Select the DZPG pilot data directory (with DCM and protocols)');
  end
  if ~exist('Pout','var')   || isempty(Pout),   Pout   = fullfile(tempdir,'cat_DCM2BIDS_selftest'); end
  if ~exist('tests','var')  || isempty(tests),  tests  = {'import','selection','special','inputs','phenotypic','export','pplevels','parallel'}; end
  if ~exist('levels','var') || isempty(levels), levels = 1:0.5:4; end
  tests = cellstr(tests);

  E   = setup(char(Pdata),char(Pout));
  res = struct('group',{},'name',{},'ok',{},'info',{});
  if any(strcmp(tests,'import')) || (any(ismember(tests,{'selection','special'})) && ~exist(E.PL0,'file'))
    res = [res runGroup('import',@testImport,E)];
  end
  for g = {'selection','special','inputs','phenotypic','export'}
    if any(strcmp(tests,g{1}))
      res = [res runGroup(g{1},str2func(['test' upper(g{1}(1)) g{1}(2:end)]),E)]; %#ok<AGROW>
    end
  end
  if any(strcmp(tests,'pplevels'))
    res = [res runGroup('pplevels',@(E) testLevels(E,levels),E)];
  end
  if any(strcmp(tests,'parallel'))
    res = [res runGroup('parallel',@testParallel,E)];
  end

  % summary
  fprintf('\ncat_DCM2BIDS_selftest (%s)\n%s\n',Pout,repmat('-',1,100));
  for ri = 1:numel(res)
    if res(ri).ok, st = 'ok'; else, st = 'FAILED'; end
    fprintf('%-10s %-26s %-7s %s\n',res(ri).group,res(ri).name,st,res(ri).info);
  end
  ok = all([res.ok]);
  fprintf('%s\n%d of %d tests passed\n',repmat('-',1,100),nnz([res.ok]),numel(res));
end
% =========================================================================
function R = runGroup(name,fh,E)
%runGroup. Test group with the error of the group as failed test.
  fprintf('cat_DCM2BIDS_selftest: %s (%s)\n',name,datestr(now,'HH:MM:SS')); %#ok<TNOW1,DATST>
  try
    R = fh(E);
  catch e
    R = struct('name','error','ok',false,'info',sprintf('%s (%s, line %d)',e.message,e.stack(1).name,e.stack(1).line));
  end
  [R.group] = deal(name);
  R = orderfields(R,{'group','name','ok','info'});
end
% =========================================================================
function E = setup(Pdata,Pout)
%setup. Test data, output directories, fake dcm2niix, and job template.
  E.Pdata = Pdata; E.Pout = Pout;
  P = fullfile(Pdata,'DCM','DCM_pilots');
  if ~exist(P,'dir'), error('cat_DCM2BIDS_selftest:noData','Cannot find "%s".',P); end
  if ~exist(Pout,'dir'), mkdir(Pout); end

  % JE7: the first two subjects with sessions <subject>_MR<n>_<date>
  D = dir(fullfile(P,'site_jena_JE7')); D = D([D.isdir]);
  tok = regexp({D.name},'^(.+)_MR(\d+)_(\d{8})$','tokens','once');
  D = D(~cellfun('isempty',tok)); tok = vertcat(tok{:});
  [~,si] = sort({D.name}); D = D(si); tok = tok(si,:);
  subs = unique(tok(:,1),'stable');
  if numel(subs) < 2, error('cat_DCM2BIDS_selftest:noData','Two subjects of site_jena_JE7 are required.'); end
  A = find(strcmp(tok(:,1),subs{1})); B = find(strcmp(tok(:,1),subs{2}));
  if numel(A) < 2, error('cat_DCM2BIDS_selftest:noData','Two sessions of the first subject of site_jena_JE7 are required.'); end
  E.dcm = fullfile(P,'site_jena_JE7',{D(A(1)).name; D(A(2)).name; D(B(1)).name});
  Ds = dir(E.dcm{1}); Ds = Ds([Ds.isdir] & ~startsWith({Ds.name},'.'));
  E.dcmses    = fullfile(E.dcm{1},Ds(1).name);              % DICOM session folder
  Dt = dir(fullfile(E.dcmses,'*T1w*'));
  E.dcmseries = fullfile(E.dcmses,Dt(1).name);              % DICOM T1w series

  % Munich (flat DICOM folder), Tuebingen (series), Berlin (preprocessing)
  E.MU  = fullfile(P,'site_munich_humanPhantom2023');
  TU    = dir(fullfile(P,'site_tuebingen_pilot202402','*','*','25_diff_standard_*_AP'));
  E.TU  = TU(1).folder;
  Fa    = dir(fullfile(E.TU,'*_ADC')); E.ADC = fullfile(E.TU,Fa(1).name);
  B     = dir(fullfile(P,'site_berlin_pilot202511')); B = B([B.isdir] & ~startsWith({B.name},'.'));
  E.BER = fullfile(P,'site_berlin_pilot202511',B(1).name);
  E.dicts = struct('Pprotocoldirs',{{fullfile(Pdata,'protocols','DZPG3T')}}, ...
    'Pcenterdict',{{fullfile(Pdata,'protocols','sites.json')}});

  % runs of the import tests
  E.T   = fullfile(Pout,'run');
  E.DB  = fullfile(E.T,'study','catDCM2BIDSdb');
  E.B   = fullfile(E.T,'study','BIDS-pseudonymized'); % BIDS directory of the anonymization level 1
  E.PL0 = fullfile(Pout,'L0.mat');                           % BIDS listing of the import

  % fake dcm2niix (logs its calls)
  E.flog = fullfile(Pout,'dcm2niix_calls.log');
  E.fake = fullfile(Pout,'fake_dcm2niix.sh');
  fid = fopen(E.fake,'w'); fprintf(fid,'#!/bin/sh\necho "$@" >> "%s"\n',E.flog); fclose(fid);
  system(sprintf('chmod +x "%s"',E.fake));
  % dcm2niix of the import (the tested release, see getDCM2NIIX in 
  % cat_io_dcm2bids_importDCM), e.g. for the reference conversion of the 
  % flat folder with the same version
  try
    E.dcm2niix = cat_io_dcm2bids_importDCM('getDCM2NIIX'); 
  catch
    E.dcm2niix = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix';
    if ~ismac, [~,E.dcm2niix] = system('which dcm2niix'); E.dcm2niix = strtrim(E.dcm2niix); end
  end
  [~,ver] = system(sprintf('"%s" --version',E.dcm2niix)); 
  E.dcm2niixver = regexp(ver,'v[0-9]+\.[0-9]+\.[0-9]+','match','once'); 

  E.job0 = struct('outdir',{{E.T}},'subdir','study');
  E.job0.opts = struct('runqc',0,'preprocessing',0,'render',struct('run',0),'reportfig',struct('run',0));
end
% =========================================================================
function R = testImport(E)
%testImport. DICOM import (baseline listing L0 of the BIDS export), export
%  of the whole database, repeated import, and JSON-only output (1)
%  followed by the images (2).
  R = struct('name',{},'ok',{},'info',{});
  fresh(E.T); mkdir(E.T); cleanLog(E);

  % (1) DICOM import
  job = E.job0; job.data = E.dcm; tic; evalc('cat_io_dcm2bids(job);'); t = toc;
  L0 = listing(E.B); save(E.PL0,'L0');
  snr   = cat_vol_findfiles(E.DB,'snr-*',struct('dirs',1,'depth',3));
  nside = cellfun(@(d) numel(dir(fullfile(d,'*=*=*=*.json'))) - numel(dir(fullfile(d,'catDCM2BIDS*.json'))), snr);
  nimp  = cellfun(@(d) numel(dir(fullfile(d,'catDCM2BIDSimport_*.json'))), snr);
  nimg  = cellfun(@(d) numel(dir(fullfile(d,'*=*=*=*.nii.gz'))) - numel(dir(fullfile(d,'anon_*.nii.gz'))), snr);
  I = cat_io_json(cat_vol_findfiles(snr{1},'catDCM2BIDSimport_*.json')); I = I(1);
  R(end+1) = struct('name','import','ok',~isempty(L0) && all(nside==1) && all(nimp==1) && all(nimg==1) && ...
    ~exist(fullfile(E.DB,'catDCM2BIDSimportpath'),'dir') && isempty(dir(fullfile(E.DB,'catDCM2BIDStmp*'))) && ...
    exist(I.DICOMdirectory,'dir') > 0, 'info',sprintf('%0.0f s, %d BIDS files, %d scans',t,numel(L0),numel(snr)));

  % report figure of the report tables of the import for all designs (see
  % cat_io_dcm2bids_reportfig, the default): one PDF per design with the
  % overview and one page per site
  Prep = cat_vol_findfiles(fullfile(E.T,'study'),'*-report',struct('dirs',1,'depth',0)); ok = false; info = 'no report tables';
  if ~isempty(Prep)
    Ptab = cat_vol_findfiles(fullfile(Prep{1},'tables'),'*',struct('dirs',1,'depth',0));
    def  = cat_io_dcm2bids_defaults; job = struct('opts',def.opts,'subdir','study', ...
      'BIDSdir',regexprep(spm_file(Prep{1},'filename'),'-report$',''),'BIDSsubdir',spm_file(Ptab{1},'filename'));
    T = readcell(fullfile(Ptab{1},['sessions.' def.opts.tableformat]),'FileType','text');
    s = T(2:end,strcmp(T(1,:),'site')); s = s(cellfun(@(x) (ischar(x) && ~isempty(x)) || isnumeric(x),s));
    nsite = numel(unique(cellfun(@num2str,s,'UniformOutput',false))); np = 1 + nsite*(nsite > 1);
    designs = cellstr(def.opts.reportfig.design); npg = zeros(1,numel(designs)); tic;
    evalc('Ppdf = cat_io_dcm2bids_reportfig(''writeReportFigure'',fullfile(E.T,''study''),job);');
    for di = 1:min(numel(designs),numel(Ppdf))
      if exist(Ppdf{di},'file'), npg(di) = numel(regexp(fileread(Ppdf{di}),'/Type\s*/Page[^s]')); end
    end
    ok   = all(npg == np);
    info = sprintf('%0.0f s, pages %s of %d (%s)',toc,mat2str(npg),np,strjoin(designs,', '));
  end
  R(end+1) = struct('name','report-figure','ok',ok,'info',info);

  % colors of the sites in the project table with the site dictionary (see
  % siteColors in cat_io_dcm2bids_protocols): report tables of the import
  % again with the dictionary, i.e. the colors of the cities of the scanners
  % of the site registry (the import uses no dictionary), with a color for
  % all sites ("*" entry of the dictionary), different colors (brightness
  % of sites with the same color), and a color of the dictionary
  ok = false; info = 'no report tables';
  if ~isempty(Prep)
    jobc = job; jobc.dicts = E.dicts; evalc('cat_io_dcm2bids_report(''writeReportTables'',fullfile(E.T,''study''),jobc);'); %#ok<NASGU>
    Q  = readcell(fullfile(Ptab{1},['project.' def.opts.tableformat]),'FileType','text');
    qs = cellfun(@num2str,Q(2:end,strcmp(Q(1,:),'site')),'UniformOutput',false); st = ~strcmp(qs,'all');
    qc = cellfun(@num2str,Q(2:end,strcmp(Q(1,:),'color')),'UniformOutput',false);
    J  = cat_io_json(E.dicts.Pcenterdict{1}); if isstruct(J), J = num2cell(J); end
    dc = cellfun(@(e) upper(char(string(e.Color))),J(cellfun(@(e) isfield(e,'Color'),J)),'UniformOutput',false);
    hex  = ~cellfun('isempty',regexp(qc,'^#[0-9A-F]{6}$','once'));
    ok   = any(st) && all(hex(st)) && numel(unique(qc(st))) == nnz(st) && any(ismember(qc(st),dc));
    info = sprintf('%d of %d sites with color (%d of the dictionary), %d different colors',nnz(hex & st),nnz(st), ...
      nnz(ismember(qc(st),dc)),numel(unique(qc(st))));
  end
  R(end+1) = struct('name','site-colors','ok',ok,'info',info);

  % (2) export of the whole database (no conversion)
  resetBIDS(E.T); job = E.job0; job.data = {E.DB}; job.opts.Pdcm2nii = E.fake; tic; evalc('cat_io_dcm2bids(job);'); t = toc;
  L1 = listing(E.B);
  R(end+1) = struct('name','export-all','ok',isequal(L0,L1) && ~exist(E.flog,'file'), ...
    'info',sprintf('%0.0f s, identical: %d, dcm2niix calls: %d',t,isequal(L0,L1),exist(E.flog,'file')>0));

  % (3) the same DICOM directories again (imported before, no conversion)
  resetBIDS(E.T); job = E.job0; job.data = E.dcm; job.opts.Pdcm2nii = E.fake; tic; evalc('cat_io_dcm2bids(job);'); t = toc;
  L2 = listing(E.B);
  R(end+1) = struct('name','dicom-again','ok',isequal(L0,L2) && ~exist(E.flog,'file'), ...
    'info',sprintf('%0.0f s, identical: %d, dcm2niix calls: %d',t,isequal(L0,L2),exist(E.flog,'file')>0));

  % (4) JSON-only output (1) in a new output directory, then the images (2)
  T2 = fullfile(E.Pout,'run1'); fresh(T2);
  DB2 = fullfile(T2,'study','catDCM2BIDSdb'); B2 = fullfile(T2,'study','BIDS-pseudonymized');
  job = E.job0; job.outdir = {T2}; job.data = E.dcm; job.opts.output = 1; tic; evalc('cat_io_dcm2bids(job);'); t1 = toc;
  nimg1 = numel(cat_vol_findfiles(DB2,'*.nii*')); Lj = listing(B2);
  job.opts.output = 2; tic; evalc('cat_io_dcm2bids(job);'); t2 = toc; L3 = listing(B2);
  R(end+1) = struct('name','output1-then-2','ok',nimg1==0 && ~any(contains(Lj,'.nii')) && isequal(names(L3),names(L0)), ...
    'info',sprintf('output 1: %0.0f s, %d images in the database; output 2: %0.0f s, equal names: %d', ...
    t1,nimg1,t2,isequal(names(L3),names(L0))));

  % (5) acquisition time of the sessions (acq_time of the sessions files) of
  % the anonymization levels 0 (date and time) and 1 (date and full hour),
  % the level 2 is tested by the export (anonymized-bids)
  T3 = fullfile(E.Pout,'acqtime'); fresh(T3); acq = {'',''};
  pat = {'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$','^\d{4}-\d{2}-\d{2}T\d{2}:00:00$'}; 
  for level = 0:1
    job = E.job0; job.outdir = {T3}; job.data = E.dcm(1); job.opts.output = 1; job.opts.anonymize = level;
    evalc('cat_io_dcm2bids(job);');
    Ps = cat_vol_findfiles(fullfile(T3,'study'),'sub-*_sessions.tsv'); 
    Ps = Ps(contains(Ps,[filesep 'BIDS' repmat('-pseudonymized',1,level) filesep])); 
    if ~isempty(Ps)
      C = cat_io_csv(Ps{1},'','',struct('delimiter','\t','convert2double',0)); c = find(strcmp(C(1,:),'acq_time'),1); 
      if ~isempty(c) && size(C,1) > 1, acq{level + 1} = char(string(C{2,c})); end
    end
  end
  ok = all(cellfun(@(a,p) ~isempty(regexp(a,p,'once')),acq,pat)) && strcmp(acq{1}(1:min(13,end)),acq{2}(1:min(13,end))); 
  R(end+1) = struct('name','acq-time','ok',ok,'info',sprintf(['level 0 date and time: %d, level 1 date and full ' ...
    'hour: %d, same date and hour: %d'],~isempty(regexp(acq{1},pat{1},'once')),~isempty(regexp(acq{2},pat{2},'once')), ...
    numel(acq{1}) >= 13 && numel(acq{2}) >= 13 && strcmp(acq{1}(1:min(13,end)),acq{2}(1:min(13,end))))); 
end
% =========================================================================
function R = testSelection(E)
%testSelection. Export of selected subjects, sessions, scans, and DICOM
%  directories (subsets of the import listing), rerun, and expected errors.
  R = struct('name',{},'ok',{},'info',{}); cleanLog(E);
  load(E.PL0,'L0');
  [sA,sB,sesA] = dbSubjects(E); bA = bidsSub(sA,E.DB); bB = bidsSub(sB,E.DB);
  t1w = dir(fullfile(E.DB,sA,sesA{1},'snr-*T1w*')); t1w = fullfile(t1w(1).folder,t1w(1).name);
  pat = @(s) [bA '/' bidsSes(s)]; anat = [pat(sesA{1}) '/(anat/.*_T1w|catDCM2BIDSaffine)'];
  tests = {
    'subject',        {fullfile(E.DB,sA)},                                  bA
    'session',        {fullfile(E.DB,sA,sesA{2})},                          pat(sesA{2})
    'scan',           {t1w},                                                anat
    'dicom-session',  {E.dcmses},                                           pat(sesA{1})
    'dicom-series',   {E.dcmseries},                                        anat
    'mixed',          {E.dcm{3}; fullfile(E.DB,sB)},                        bB
    'selected-twice', {fullfile(E.DB,sA); fullfile(E.DB,sA,sesA{2})},       bA
    };
  for ti = 1:size(tests,1)
    resetBIDS(E.T); job = E.job0; job.data = tests{ti,2}; job.opts.Pdcm2nii = E.fake;
    txt = evalc('cat_io_dcm2bids(job);');
    Lt = scans(listing(E.B)); Le = sel(scans(L0),tests{ti,3});
    R(end+1) = struct('name',tests{ti,1},'ok',isequal(Lt,Le) && ~isempty(Le), ...
      'info',sprintf('%d/%d files, repeated: %d, cannot find: %d',numel(Lt),numel(Le), ...
      numel(strfind(txt,'ignore repeated scan')),numel(strfind(txt,'Cannot find')))); %#ok<AGROW>
  end

  % rerun with database input (no conversion, but processing again)
  resetBIDS(E.T); job = E.job0; job.data = {E.DB}; job.opts.Pdcm2nii = E.fake; job.opts.rerun = 1;
  evalc('cat_io_dcm2bids(job);'); Lt = scans(listing(E.B));
  R(end+1) = struct('name','rerun-all','ok',isequal(names(Lt),names(scans(L0))) && ~exist(E.flog,'file'), ...
    'info',sprintf('%d files, dcm2niix calls: %d',numel(Lt),exist(E.flog,'file')>0));

  % database of another output directory: error before any output
  T2 = fullfile(E.Pout,'other'); fresh(T2);
  job = E.job0; job.outdir = {T2}; job.data = {fullfile(E.DB,sA)};
  try cat_io_dcm2bids(job); id = 'none'; catch e, id = e.identifier; end
  R(end+1) = struct('name','foreign-db','ok',strcmp(id,'cat_io_dcm2bids:foreignDatabase') && ~exist(T2,'dir'),'info',id);

  % no database directory
  mkdir(fullfile(E.DB,'foo')); job = E.job0; job.data = {fullfile(E.DB,'foo')};
  try cat_io_dcm2bids(job); id = 'none'; catch e, id = e.identifier; end
  rmdir(fullfile(E.DB,'foo'));
  R(end+1) = struct('name','no-db-dir','ok',strcmp(id,'cat_io_dcm2bids:noDatabaseDirectory'),'info',id);

  % scan without image and without DICOM directory (output 2)
  Pd = dir(fullfile(E.DB,sA,sesA{2},'snr-*T2w*')); Pd = fullfile(Pd(1).folder,Pd(1).name);
  Pi = cat_vol_findfiles(Pd,'catDCM2BIDSimport_*.json'); I0 = cat_io_json(Pi{1}); I = I0;
  I.DICOMdirectory = '/nonexistent/dicom'; cat_io_json(Pi{1},I);
  Pim = cat_vol_findfiles(Pd,'*=*=*=*.nii.gz'); Pim = Pim(~contains(Pim,'anon_')); movefile(Pim{1},[Pim{1} '.bak']);
  resetBIDS(E.T); job = E.job0; job.data = {Pd}; job.opts.Pdcm2nii = E.fake;
  try txt = evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; txt = ''; end
  movefile([Pim{1} '.bak'],Pim{1}); cat_io_json(Pi{1},I0);
  R(end+1) = struct('name','missing-dicom','ok',strcmp(id,'none') && contains(txt,'DICOM directory not available') && ...
    isempty(scans(listing(E.B))),'info',id);

  % the same directory for the input and the output, and an input within the
  % output directory (not the database): error before any output
  T3 = fullfile(E.Pout,'inout'); fresh(T3); mkdir(T3);
  job = E.job0; job.outdir = {T3}; job.data = {T3};
  try cat_io_dcm2bids(job); id = 'none'; catch e, id = e.identifier; end
  job = E.job0; job.data = {E.B};
  try cat_io_dcm2bids(job); id2 = 'none'; catch e, id2 = e.identifier; end
  R(end+1) = struct('name','input-output','ok',strcmp(id,'cat_io_dcm2bids:inputOutput') && ...
    strcmp(id2,'cat_io_dcm2bids:inputOutput') && ~exist(fullfile(T3,'study'),'dir'),'info',sprintf('%s, %s',id,id2));
  R(end+1) = struct('name','no-conversion','ok',~exist(E.flog,'file'),'info','fake dcm2niix not called');
end
% =========================================================================
function R = testSpecial(E)
%testSpecial. Stored settings of a previous version, two copies of a series
%  in the database, symbolic links, relative paths, dcm2niix options, site
%  and session labels (synthetic databases), site dictionaries with mixed
%  fields, and the minimum age of the participants file.
  R = struct('name',{},'ok',{},'info',{}); cleanLog(E);
  load(E.PL0,'L0'); Lb = names(scans(L0));
  [sA,~,sesA] = dbSubjects(E); bA = bidsSub(sA,E.DB);
  t1w = dir(fullfile(E.DB,sA,sesA{1},'snr-*T1w*')); t1w = fullfile(t1w(1).folder,t1w(1).name);
  jsonsDB = @(d) numel(dir(fullfile(d,'*=*=*=*.json'))) - numel(dir(fullfile(d,'catDCM2BIDS*.json')));
  sesLb = sel(Lb,[bA '/' bidsSes(sesA{1})]);

  % stored settings of a previous version (with longDBnames)
  resetBIDS(E.T); job = E.job0; job.data = {E.DB}; job.opts.Pdcm2nii = E.fake; evalc('cat_io_dcm2bids(job);');
  Popts = fullfile(E.B,'BIDS','catDCM2BIDSdb.mat');
  load(Popts,'opts','dicts'); opts.longDBnames = 1; save(Popts,'opts','dicts');
  try txt = evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; txt = ''; end
  R(end+1) = struct('name','old-settings','ok',strcmp(id,'none') && ~contains(txt,'does not fit') && ...
    isequal(names(scans(listing(E.B))),Lb),'info',id);

  % database scan directory with two copies of the series (older imports)
  Pr = dir(fullfile(t1w,'*=*=*=*.json')); Pr = Pr(~startsWith({Pr.name},'catDCM2BIDS'));
  Pr = fullfile(t1w,Pr(1).name); [~,nam] = fileparts(Pr);
  Pc = fullfile(t1w,['selftestcopy' nam(find(nam=='=',1):end)]);
  sn = regexp(nam,'=(\d+)$','tokens','once'); sn = sprintf('acq-%03d',str2double(sn{1}));
  copyfile(Pr,[Pc '.json']); copyfile(strrep(Pr,'.json','.nii.gz'),[Pc '.nii.gz']);
  resetBIDS(E.T); job = E.job0; job.data = {t1w}; job.opts.Pdcm2nii = E.fake; txt = evalc('cat_io_dcm2bids(job);');
  Lt = names(scans(listing(E.B)));
  delete([Pc '.json']); delete([Pc '.nii.gz']);
  R(end+1) = struct('name','db-two-copies','ok',isequal(Lt,sel(Lb,[bA '/' bidsSes(sesA{1}) '/(anat/.*_T1w|catDCM2BIDSaffine)'])) && ...
    any(contains(Lt,sn)),'info',sprintf('%d files, repeated: %d',numel(Lt),numel(strfind(txt,'ignore repeated'))));

  % DICOM session selected by another path (symbolic link): read again, but no new scans
  Ld = fullfile(E.Pout,'dicomlink'); if ~exist(Ld,'dir'), system(sprintf('ln -s "%s" "%s"',E.dcm{1},Ld)); end
  resetBIDS(E.T); job = E.job0; job.data = {Ld}; txt = evalc('cat_io_dcm2bids(job);'); Lt = names(scans(listing(E.B)));
  R(end+1) = struct('name','dicom-other-path','ok',jsonsDB(t1w)==1 && contains(txt,'0 new scans') && isequal(Lt,sesLb), ...
    'info',sprintf('sidecars of the T1w: %d, %d files',jsonsDB(t1w),numel(Lt)));

  % DICOM rerun: headers read again, images kept, no second copies
  resetBIDS(E.T); job = E.job0; job.data = E.dcm(1); job.opts.rerun = 1; evalc('cat_io_dcm2bids(job);');
  Lt = names(scans(listing(E.B))); D1 = dir(fullfile(E.DB,sA,sesA{1},'*','*=*=*=*a.json'));
  R(end+1) = struct('name','dicom-rerun','ok',jsonsDB(t1w)==1 && isempty(D1) && ...
    isempty(dir(fullfile(E.DB,'catDCM2BIDStmp*'))) && isequal(Lt,sesLb), ...
    'info',sprintf('sidecars of the T1w: %d, a-files: %d, %d files',jsonsDB(t1w),numel(D1),numel(Lt)));

  % ignoreScouts = 0: dcm2niix without -i y (JSON-only import into a new output)
  T3 = fullfile(E.Pout,'run3'); fresh(T3);
  job = E.job0; job.outdir = {T3}; job.data = {E.dcmseries}; job.opts.ignoreScouts = 0; job.opts.output = 1;
  evalc('cat_io_dcm2bids(job);');
  Pi = cat_vol_findfiles(fullfile(T3,'study','catDCM2BIDSdb'),'catDCM2BIDSimport_*.json'); I = cat_io_json(Pi{1});
  R(end+1) = struct('name','no-ignoreScouts','ok',isscalar(Pi) && contains(I.dcm2niixOptions,'-i n'),'info',I.dcm2niixOptions);

  % output path with "scout" (symbolic link) and a database selection by this path
  Ls = fullfile(E.Pout,'scoutlink'); if ~exist(Ls,'dir'), system(sprintf('ln -s "%s" "%s"',E.T,Ls)); end
  resetBIDS(E.T); job = E.job0; job.outdir = {Ls}; job.opts.Pdcm2nii = E.fake;
  job.data = {fullfile(Ls,'study','catDCM2BIDSdb',sA,sesA{1})};
  try evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; end
  Lt = names(scans(listing(E.B)));
  R(end+1) = struct('name','scout-in-outdir','ok',strcmp(id,'none') && isequal(Lt,sesLb),'info',sprintf('%s, %d files',id,numel(Lt)));

  % database selected by a symbolic link, real output directory
  resetBIDS(E.T); job = E.job0; job.opts.Pdcm2nii = E.fake; job.data = {fullfile(Ls,'study','catDCM2BIDSdb',sA)};
  try evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; end
  Lt = names(scans(listing(E.B)));
  R(end+1) = struct('name','symlink-selection','ok',strcmp(id,'none') && isequal(Lt,sel(Lb,bA)),'info',sprintf('%s, %d files',id,numel(Lt)));

  % relative output directory, absolute selection
  resetBIDS(E.T); old = cd(E.Pout); job = E.job0; job.outdir = {'run'}; job.opts.Pdcm2nii = E.fake; job.data = {fullfile(E.DB,sA)};
  try evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; end
  cd(old); Lt = names(scans(listing(E.B)));
  R(end+1) = struct('name','relative-outdir','ok',strcmp(id,'none') && isequal(Lt,sel(Lb,bA)),'info',sprintf('%s, %d files',id,numel(Lt)));
  R(end+1) = struct('name','no-conversion','ok',~exist(E.flog,'file'),'info','fake dcm2niix not called');

  % site registry without site dictionary (synthetic database): the same 
  % abbreviation (Heidelberg and Heidenheim, two Siemens scanners of Jena),
  % another manufacturer and field strength, no address, and in a second 
  % run no manufacturer and a new scanner of a registered abbreviation 
  % (Heilbronn), with unchanged abbreviations of the first run
  Ds = fullfile(E.Pout,'sites','catDCM2BIDSdb'); fresh(fileparts(Ds));
  sc = {'1001','Str,Heidelberg,BW,DE',3,'SIEMENS'; '1002','Str,Heidenheim,BW,DE',3,'Siemens Healthineers'; ...
        '1003','Str,Jena,TH,DE',3,'SIEMENS'; '1004','Str,Jena,TH,DE',3,'SIEMENS'; ...
        '1005','Str,Jena,TH,DE',3,'Philips Medical Systems'; '1006','Str,Jena,TH,DE',1.5,'SIEMENS'; ...
        '1007','',3,'GE MEDICAL SYSTEMS'; '1008','Str,Magdeburg,ST,DE',3,''; '1009','Str,Heilbronn,BW,DE',3,'SIEMENS'};
  for k = 1:7, fakeScan(Ds,sc(k,:)); end
  o = struct('sitelength',3,'sitevendor',2); S1 = {}; S2 = {}; %#ok<NASGU> (evalc)
  evalc('S1 = cat_io_dcm2bids_protocols(''siteRegistry'',Ds,o);');
  for k = 8:9, fakeScan(Ds,sc(k,:)); end
  evalc('S2 = cat_io_dcm2bids_protocols(''siteRegistry'',Ds,o);');
  e1 = {'HeASI30T','HeBSI30T','JeASI30T','JeBSI30T','JenPH30T','JenSI15T','UnkGE30T'}; e2 = [e1 {'MagXX30T','HeCSI30T'}];
  R(end+1) = struct('name','site-registry','ok',isequal(S1(:,2)',e1) && isequal(S2(:,2)',e2) && isequal(S2(:,1)',sc(:,1)'), ...
    'info',sprintf('run 1: %s; run 2: %s',strjoin(S1(:,2)',' '),strjoin(S2(:,2)',' ')));

  % site dictionary with different fields of the entries (cell array of 
  % jsondecode), a siteID (anonymization level 2), a numeric serial number,
  % an entry without serial number, and a color entry of cities (skipped, 
  % see siteColors)
  Pd = fullfile(fileparts(Ds),'sites_mixed.json'); Sd = {}; %#ok<NASGU> (evalc)
  cat_io_dcm2bids_helper('writeText',Pd,['[{"DeviceSerialNumber":"1001","InstitutionAbbreviation":"Abc","siteID":"07","InstitutionName":"X"},' ...
    '{"InstitutionAbbreviation":"Def","DeviceSerialNumber":2002},{"InstitutionName":"no serial"},' ...
    '{"Name":"site","Cities":["Jena","Halle"],"Color":"#A82117"}]']);
  evalc('Sd = cat_io_dcm2bids_protocols(''getSites'',{Pd});');
  R(end+1) = struct('name','site-dictionary-mixed','ok',isequal(Sd,{'1001','Abc','07';'2002','Def',''}), ...
    'info',sprintf('%d entries: %s',size(Sd,1),strjoin(strcat(Sd(:,1),'=',Sd(:,2))',' ')));

  % minimum age of the participants file over runs (e.g. a later run with 
  % only a follow-up session)
  Pp = fullfile(E.Pout,'participants','participants.tsv'); fresh(fileparts(Pp));
  ages = [45.3 40.1 50.2 NaN]; a = nan(size(ages));
  for k = 1:numel(ages)
    cat_io_dcm2bids_bids('writeParticipantTSV',struct('sub','sub-01','V',struct('PatientSex','F','PatientAge',ages(k)), ...
      'extra',[]),Pp,1);
    T = cat_io_csv(Pp,'','',struct('delimiter','\t','convert2double',-2)); a(k) = T{2,3};
  end
  R(end+1) = struct('name','participants-min-age','ok',isequal(a,[45.3 40.1 40.1 40.1]),'info',sprintf('ages %s',mat2str(a)));

  % session labels (synthetic database): two studies of a day, a study two
  % days later, and a study over midnight (two directories) by the 
  % anonymization level (0-2, weeks), and the labels of the import of the 
  % anonymized DICOM export (StudyID with the session label, the dates of 
  % the labels, and the times)
  Dl = fullfile(E.Pout,'sessions'); fresh(Dl);
  S = {'ses-1-20250101','20250101120000'; 'ses-2-20250101','20250101150000'; 'ses-3-20250103','20250103090000'; ...
       'ses-4-20250105','20250105235000'; 'ses-4-20250106','20250105235000'};
  T = {'ses-w0000-19000101','19000101120000'; 'ses-w0000b-19000101','19000101150000'; 'ses-w0000c-19000101','19000101090000'; ...
       'ses-w0001-19000108','19000108235000'; 'ses-w0001-19000109','19000108235000'};
  for k = 1:5, fakeSession(Dl,'sub-A-1',S(k,:)); fakeSession(Dl,'sub-B-1',T(k,:)); end
  off = days(datetime(1900,1,1) - datetime(2025,1,1)); o = struct('BIDSsep','','sesunit','w'); L = cell(4,5);
  for k = 1:5
    for lev = 0:2
      o.anonymize = lev; L{lev+1,k} = cat_io_dcm2bids_bids('sessionLabel',fullfile(Dl,'sub-A-1',S{k,1}),S{k,2},o,off);
    end
    o.anonymize = 1; L{4,k} = cat_io_dcm2bids_bids('sessionLabel',fullfile(Dl,'sub-B-1',T{k,1}),T{k,2},o,0);
  end
  e = {'ses-20250101','ses-20250101b','ses-20250103','ses-20250105','ses-20250105'; ...
       'ses-w0000','ses-w0000b','ses-w0000c','ses-w0001','ses-w0001'};
  % sessions without date (time only, a study in two directories) and with 
  % the date 1970-01-01 (placeholder of anonymized data, no shifted date)
  C = {'ses-1116','134431'; 'ses-1116-20180925','134431'; 'ses-1117','150000'; 'ses-1-19700101','19700101120000'};
  for k = 1:3, fakeSession(Dl,'sub-C-1',C(k,:)); end, fakeSession(Dl,'sub-D-1',C(4,:));
  o.anonymize = 1; Lc = cell(1,4);
  for k = 1:4, Lc{k} = cat_io_dcm2bids_bids('sessionLabel',fullfile(Dl,['sub-' char('C' + (k == 4)) '-1'],C{k,1}),C{k,2},o,0); end
  ec = {'ses-NA','ses-NA','ses-NAb','ses-19700101'};
  R(end+1) = struct('name','session-labels','ok',numel(unique(L(1,:))) == 4 && isequal(L(2,:),e(1,:)) && ...
    isequal(L(3,:),e(2,:)) && isequal(L(4,:),e(2,:)) && isequal(Lc,ec),'info',sprintf('level 1: %s; level 2: %s; import: %s; %s', ...
    strjoin(L(2,:),' '),strjoin(L(3,:),' '),strjoin(L(4,:),' '),strjoin(Lc,' ')));
end
% =========================================================================
function fakeSession(D,sub,s)
%fakeSession. Scan directory with an empty sidecar of a synthetic database
%  D with the subject sub and the session directory and time stamp s.
  P = fullfile(D,sub,s{1},'snr-0001-t1'); mkdir(P); fclose(fopen(fullfile(P,['f=t1=' s{2} '=1.json']),'w'));
end
% =========================================================================
function fakeScan(Ds,sc)
%fakeScan. Sidecar of a scan of a synthetic database Ds with the serial 
%  number, address, field strength, and manufacturer sc.
  P = fullfile(Ds,['sub-' sc{1} '-P1'],'ses-1-20250101','snr-0001-t1'); mkdir(P);
  cat_io_json(fullfile(P,'f=t1=20250101120000=1.json'),struct('DeviceSerialNumber',sc{1},'InstitutionAddress',sc{2}, ...
    'MagneticFieldStrength',sc{3},'Manufacturer',sc{4},'PatientID','P1','StudyID','1','SeriesNumber',1, ...
    'AcquisitionDateTime','2025-01-01T12:00:00'));
end
% =========================================================================
function R = testInputs(E)
%testInputs. Split and copied series, flat DICOM folder with dcm2niix name
%  conflicts (with the protocol tables), output 0, changed dcm2niix options
%  (derived ADC series), UTF-8 paths, and a BIDS dataset of DataLad.
  R = struct('name',{},'ok',{},'info',{});
  In = fullfile(E.Pout,'inputs'); fresh(In);
  mkdir(fullfile(In,'split','p1')); mkdir(fullfile(In,'split','p2'));
  mkdir(fullfile(In,'copy','A')); mkdir(fullfile(In,'copy','B'));
  Pu = fullfile(In,'umlaut',['Studie ' char(8211) ' M' char(252) 'ller']); mkdir(Pu);
  lnk = @(f,d) system(sprintf('ln -s "%s" "%s"',f,d));
  F = dir(fullfile(E.TU,'25_diff_standard_*_AP')); F = dir(fullfile(F(1).folder,F(1).name)); F = F(~[F.isdir]);
  for fi = 1:numel(F), lnk(fullfile(F(fi).folder,F(fi).name),fullfile(In,'split',sprintf('p%d',1 + (fi > 52)))); end
  F = dir(fullfile(E.TU,'05_gre_field_mapping*')); F = dir(fullfile(F(1).folder,F(1).name)); F = F(~[F.isdir]);
  for fi = 1:numel(F), lnk(fullfile(F(fi).folder,F(fi).name),fullfile(In,'copy','A')); lnk(fullfile(F(fi).folder,F(fi).name),fullfile(In,'copy','B')); end
  F = dir(fullfile(E.TU,'22_T1w_MPR')); F = dir(fullfile(F(1).folder,F(1).name)); F = F(~[F.isdir]);
  for fi = 1:numel(F), lnk(fullfile(F(fi).folder,F(fi).name),Pu); end

  % (1) split series and copy
  T7 = fullfile(E.Pout,'run7'); fresh(T7); job = E.job0; job.outdir = {T7};
  job.data = {fullfile(In,'split'); fullfile(In,'copy')}; job.opts.output = 1;
  txt = evalc('cat_io_dcm2bids(job);');
  side = sidecars(fullfile(T7,'study','catDCM2BIDSdb'));
  R(end+1) = struct('name','split-and-copy','ok',contains(txt,'split over several directories') && numel(side)==2 && ...
    ~any(contains(side,'diff')),'info',sprintf('%d sidecars',numel(side)));

  % (2) flat DICOM folder with dcm2niix name conflicts: output 1, then 2
  %     (sidecars of the full conversion by dcm2niix)
  Pf = fullfile(E.Pout,['fullconversion_' E.dcm2niixver]); % (reference of the dcm2niix version of the import)
  if ~exist(Pf,'dir') || isempty(dir(fullfile(Pf,'*.json'))) % (also an empty folder, e.g. of an interrupted conversion)
    if ~exist(Pf,'dir'), mkdir(Pf); end; IM = cat_vol_findfiles(E.MU,'IMAGES',struct('dirs',1));
    system(sprintf('"%s" -f "%%f=%%p=%%t=%%s" -p n -ba n -d 0 -i y -b y -z y -o "%s" "%s"',E.dcm2niix,Pf,IM{1}));
  end
  T8 = fullfile(E.Pout,'run8'); fresh(T8); job = E.job0; job.outdir = {T8}; job.data = {E.MU}; job.opts.output = 1;
  evalc('cat_io_dcm2bids(job);');
  DB8 = fullfile(T8,'study','catDCM2BIDSdb'); s1 = sidecars(DB8);
  Pt8 = fullfile(T8,'study','catDCM2BIDSdb-tables','protocols.csv');
  Tp1 = readcell(Pt8,'Delimiter',',','FileType','text'); Pv1 = cat_vol_findfiles(fullfile(T8,'study','catDCM2BIDSdb-tables','protocols'),'*.json');
  job.opts.output = 2; tic; txt = evalc('cat_io_dcm2bids(job);'); t = toc;
  s2  = sidecars(DB8);
  img = cat_vol_findfiles(DB8,'*=*=*=*.nii.gz');
  img = img(~contains(img,'anon_') & ~cellfun(@(f) any(regexp(spm_file(f,'basename'),'^(c0|m|wc|mwc|l0|wl0|y_)')),img));
  full = cat_vol_findfiles(Pf,'*.json');
  s2f  = s2(ismember(spm_file(s2,'filename'),spm_file(full,'filename')));
  same = cellfun(@(f) isequal(sidecarFields(f),sidecarFields(full{strcmp(spm_file(full,'filename'),spm_file(f,'filename'))})),s2f);
  R(end+1) = struct('name','flat-folder','ok',numel(s2)==numel(img)+1 && all(same) && numel(same)==numel(img) && ...
    contains(txt,'no image after DICOM conversion'),'info',sprintf('output 1: %d sidecars; output 2: %d sidecars, %d images, equal to dcm2niix: %d/%d, %0.0f s', ...
    numel(s1),numel(s2),numel(img),nnz(same),numel(same),t));

  % protocol variants: the conversion (output 2) of the JSON-only scans
  % (output 1) keeps the variants (with the image fields)
  Tp2 = readcell(Pt8,'Delimiter',',','FileType','text'); Pv2 = cat_vol_findfiles(fullfile(T8,'study','catDCM2BIDSdb-tables','protocols'),'*.json');
  hdr = Tp2(1,:); c = find(strcmp(hdr,'VoxelSize'));
  isna = @(x) isa(x,'missing') || (ischar(x) && strcmp(x,'n/a')) || (isnumeric(x) && isscalar(x) && isnan(x));
  na1 = cellfun(isna,Tp1(2:end,c)); na2 = cellfun(isna,Tp2(2:end,c));
  R(end+1) = struct('name','protocol-table-columns','ok',all(ismember({'Sites','SiteList','VoxelSize','NumberOfSlices','NumberOfVolumes'},hdr)) && ...
    all(na1) && nnz(na2) < numel(na2) && isequal(sort(Tp1(2:end,1)),sort(Tp2(2:end,1))) && numel(Pv1)==size(Tp1,1)-1 && numel(Pv2)==size(Tp2,1)-1, ...
    'info',sprintf('variants %d -> %d, VoxelSize n/a %d -> %d',size(Tp1,1)-1,size(Tp2,1)-1,nnz(na1),nnz(na2)));
  Pj = Pv2(contains(Pv2,[filesep 'anat' filesep])); J = cat_io_json(Pj{1});
  R(end+1) = struct('name','protocol-json-image','ok',all(isfield(J,{'VoxelSize','NumberOfSlices','NumberOfVolumes'})),'info','');
  B8 = listing(fullfile(T8,'study','BIDS-pseudonymized')); bj = B8(contains(B8,'.json') & contains(B8,'/anat/') & ~contains(B8,'derivatives'));
  J  = cat_io_json(fullfile(T8,'study','BIDS-pseudonymized',regexprep(bj{1},'^\S+\s+\./','')));
  R(end+1) = struct('name','bids-json-imagefields','ok',isfield(J,'Dimensions') && isfield(J,'VoxelSize'),'info','');

  % (3) output 0 (overview only)
  T9 = fullfile(E.Pout,'run9'); fresh(T9); job = E.job0; job.outdir = {T9}; job.data = {fullfile(In,'copy')}; job.opts.output = 0;
  try evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, id = e.message; end
  R(end+1) = struct('name','output0','ok',strcmp(id,'none') && exist(fullfile(T9,'study','catDCM2BIDSdb-tables','protocols.csv'),'file') && ...
    ~exist(fullfile(T9,'study','BIDS-pseudonymized','BIDS'),'dir'),'info',id);

  % (4) derived image: import without ignoreScouts (ADC kept in the 
  %     database), then with ignoreScouts and output 2 (ADC ignored as 
  %     derived image, see ignoredDerived; dcm2niix always with -i n, i.e. 
  %     the directory is not read again)
  T10 = fullfile(E.Pout,'run10'); fresh(T10); job = E.job0; job.outdir = {T10}; job.data = {E.ADC};
  job.opts.output = 1; job.opts.ignoreScouts = 0; evalc('cat_io_dcm2bids(job);');
  n0 = numel(cat_vol_findfiles(fullfile(T10,'study','catDCM2BIDSdb'),'catDCM2BIDSimport_*.json'));
  job.opts.ignoreScouts = 1; job.opts.output = 2; txt = evalc('cat_io_dcm2bids(job);');
  Lb = listing(fullfile(T10,'study','BIDS-pseudonymized'));
  R(end+1) = struct('name','options-changed','ok',n0==1 && ~contains(txt,'no image') && contains(txt,'1 directories imported before') && ...
    ~any(contains(Lb,'.nii')),'info',sprintf('ADC sidecars with -i n: %d, BIDS images: %d',n0,nnz(contains(Lb,'.nii'))));
  job.data = {fullfile(T10,'study','catDCM2BIDSdb')}; txt = evalc('cat_io_dcm2bids(job);');
  R(end+1) = struct('name','derived-db-input','ok',contains(txt,'ignore derived image') && ~contains(txt,'no image'),'info','');

  % (5) DICOM path with special characters (UTF-8): stored correctly and
  %     recognized in the next run (no dcm2niix call)
  T11 = fullfile(E.Pout,'run11'); fresh(T11); job = E.job0; job.outdir = {T11}; job.data = {fileparts(Pu)}; job.opts.output = 1;
  evalc('cat_io_dcm2bids(job);');
  Pi = cat_vol_findfiles(fullfile(T11,'study','catDCM2BIDSdb'),'catDCM2BIDSimport_*.json'); I = jsondecode(fileread(Pi{1},'Encoding','UTF-8'));
  cleanLog(E); job.opts.Pdcm2nii = E.fake; evalc('cat_io_dcm2bids(job);');
  R(end+1) = struct('name','utf8-path','ok',contains(I.DICOMdirectory,['Studie ' char(8211) ' M' char(252) 'ller']) && ~exist(E.flog,'file'),'info','');

  % (6) no temporary directories left
  tl = [dir(fullfile(T7,'study','catDCM2BIDSdb','catDCM2BIDStmp*')); dir(fullfile(T8,'study','catDCM2BIDSdb','catDCM2BIDStmp*'))];
  R(end+1) = struct('name','no-tmp-left','ok',isempty(tl),'info','');

  % (7) BIDS dataset of DataLad (simulated): an annexed image (link to the
  % read-only content in .git/annex), an image without content (link 
  % without file), and the output within a DataLad dataset: the available 
  % image is imported as writable copy, the missing image is reported (no
  % error of the session without images), and the output gets a warning
  Dl = fullfile(E.Pout,'datalad'); fresh(Dl); ds = fullfile(Dl,'ds'); 
  Pt = cat_vol_findfiles(fullfile(E.B,'BIDS'),'sub-*_T1w.nii.gz'); Pt = Pt(~contains(Pt,'derivatives')); 
  K  = fullfile(ds,'.git','annex','objects','aa','bb','KEY','KEY'); mkdir(fileparts(K)); copyfile(Pt{1},K); fileattrib(K,'-w','a');
  for si = 1:2
    Pa = fullfile(ds,sprintf('sub-%02d',si),'anat'); mkdir(Pa); 
    copyfile(spm_file(spm_file(Pt{1},'ext',''),'ext','.json'),fullfile(Pa,sprintf('sub-%02d_T1w.json',si)));
  end
  system(sprintf('ln -s ../../.git/annex/objects/aa/bb/KEY/KEY "%s"',fullfile(ds,'sub-01','anat','sub-01_T1w.nii.gz')));
  system(sprintf('ln -s ../../.git/annex/objects/cc/dd/MISSING/MISSING "%s"',fullfile(ds,'sub-02','anat','sub-02_T1w.nii.gz')));
  cat_io_dcm2bids_helper('writeText',fullfile(ds,'dataset_description.json'),'{"Name":"DataLad test","BIDSVersion":"1.9.0"}');
  mkdir(fullfile(ds,'.datalad')); mkdir(fullfile(Dl,'out','.datalad')); 
  job = E.job0; job.outdir = {fullfile(Dl,'out')}; job.data = {ds}; 
  try out = evalc('cat_io_dcm2bids(job);'); id = 'none'; catch e, out = ''; id = e.message; end
  Pdb = cat_vol_findfiles(fullfile(Dl,'out','study','catDCM2BIDSdb'),'ds=*.nii*'); wr = false;
  if isscalar(Pdb), [~,a] = fileattrib(Pdb{1}); wr = a.UserWrite == 1; end
  Pb = cat_vol_findfiles(fullfile(Dl,'out','study','BIDS-pseudonymized','BIDS'),'sub-*_T1w.nii.gz');
  fileattrib(K,'+w','a'); 
  R(end+1) = struct('name','datalad-bids','ok',strcmp(id,'none') && isscalar(Pdb) && wr && isscalar(Pb) && ...
    contains(out,'datalad get') && contains(out,'DataLad dataset'), ...
    'info',sprintf('%s, database copies %d (writable %d), BIDS images %d, hints: %d %d',id,numel(Pdb),wr,numel(Pb), ...
    contains(out,'datalad get'),contains(out,'DataLad dataset')));
end
% =========================================================================
function R = testPhenotypic(E)
%testPhenotypic. Phenotypic data tables (see cat_io_dcm2bids_phenotypic) of
%  the JE7 sessions (JSON-only output): session rows by name (German column
%  names, other order than the DICOM name), birth date (dd.mm.yyyy), and 
%  scan date with their data in the sessions files (also of a reserved 
%  column name), subject rows with an unknown name and a conflict, key 
%  columns by number with other names, no names and birth dates in the BIDS
%  directory, the data of the database without these tables, and a 
%  corrected table. The tables are created from the imported data.
  R  = struct('name',{},'ok',{},'info',{});
  T  = fullfile(E.Pout,'phenotypic'); fresh(T); mkdir(T);
  DB = fullfile(T,'study','catDCM2BIDSdb'); B = fullfile(T,'study','BIDS-pseudonymized','BIDS');
  job = E.job0; job.outdir = {T}; job.data = E.dcm; job.opts.output = 1;
  evalc('cat_io_dcm2bids(job);'); job.data = {DB};

  % sessions of the database (DICOM name, PatientID, birth date, scan date)
  S = struct('pn',{},'pid',{},'birth',{},'date',{});
  D = dir(fullfile(DB,'sub-*','ses-*')); D = D([D.isdir]);
  for di = 1:numel(D)
    P = cat_vol_findfiles(fullfile(D(di).folder,D(di).name),'*=*=*=*.json'); 
    P = P(~startsWith(spm_file(P,'filename'),'catDCM2BIDS'));
    J = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(P{1}));
    S(end+1) = struct('pn',char(J.PatientName),'pid',char(J.PatientID),'birth',char(J.PatientBirthDate), ...
      'date',char(datetime(J.ScanDate,'Format','yyyy-MM-dd'))); %#ok<AGROW>
  end
  n = arrayfun(@(s) nnz(strcmp({S.pid},s.pid)),S);
  [~,o] = sortrows([-n(:) posixtime(datetime({S.date},'InputFormat','yyyy-MM-dd'))']);
  S = S(o); % sessions of the subject with two sessions first (by date)
  given  = @(s) strjoin(subsref(strsplit(s.pn,'^'),struct('type','()','subs',{{2:numel(strsplit(s.pn,'^'))}})),' ');
  family = @(s) subsref(strsplit(s.pn,'^'),struct('type','{}','subs',{{1}}));
  birth  = @(s) germanDate(s.birth);
  Pp = fullfile(T,'tables'); mkdir(Pp);

  % (1) session rows (name, birth date, scan date) and subject rows (name):
  % the subject with one session, the subject with two sessions, an unknown
  % name, and a conflict (the first session of the subject with two sessions
  % by a subject row and a visit)
  C1 = {'Vorname','Nachname','Geburtsdatum','Scan date','MMSE','Group','age'};
  for k = 1:numel(S)
    C1(end+1,:) = {given(S(k)),family(S(k)),birth(S(k)),S(k).date,sprintf('%d',20+k),sprintf('g%d',k),sprintf('%d',60+k)}; %#ok<AGROW>
  end
  C2 = {'Name','Visit','Diagnosis'; [given(S(3)) ' ' family(S(3))],'','dx3'; [given(S(1)) ' ' family(S(1))],'','dx1'; ...
        'Erika Musterfrau','','dx0'; [family(S(1)) ' ' given(S(1))],'1','dx9'};
  P1 = fullfile(Pp,'visits.csv'); P2 = fullfile(Pp,'diagnoses.tsv');
  cat_io_csv(P1,C1); cat_io_csv(P2,C2,'','',struct('delimiter','\t'));
  job1 = job; job1.dataset.phenotypic = {P1; P2}; evalc('cat_io_dcm2bids(job1);');
  [V,ok1] = sessionValues(B,S,{'MMSE','Group','age_user','Diagnosis'});
  exp1 = arrayfun(@(k) {sprintf('%d',20+k),sprintf('g%d',k),sprintf('%d',60+k)},1:numel(S),'UniformOutput',false);
  st1 = reportStatus(fullfile(T,'study','catDCM2BIDSdb-tables','phenotypic_visits.csv'));
  st2 = reportStatus(fullfile(T,'study','catDCM2BIDSdb-tables','phenotypic_diagnoses.csv'));
  okv = ok1 && isequal(cellfun(@(v) v(1:3),V,'UniformOutput',false),exp1) && ...
    isequal(cellfun(@(v) v{4},V,'UniformOutput',false),{'n/a','dx1','dx3'});
  R(end+1) = struct('name','session-rows','ok',okv && all(startsWith(st1,'assigned (1 sessions)')), ...
    'info',sprintf('%d sessions, values in the sessions files: %d, visits: %s',numel(S),okv,strjoin(unique(st1),'; ')));
  R(end+1) = struct('name','subject-rows','ok',strcmp(st2{1},'assigned (1 sessions)') && ...
    startsWith(st2{2},'conflict (1 of 2') && strcmp(st2{3},'no session') && startsWith(st2{4},'conflict'), ...
    'info',strjoin(st2,'; '));

  % (2) key columns by number with other names (name, name, scan date as
  % yyyymmdd): the data of the other tables are kept (database)
  C3 = {'A','B','C','score'};
  for k = 1:numel(S), C3(end+1,:) = {family(S(k)),given(S(k)),strrep(S(k).date,'-',''),sprintf('%d',k)}; end %#ok<AGROW>
  P3 = fullfile(Pp,'scores.csv'); cat_io_csv(P3,C3);
  job3 = job; job3.dataset.phenotypic = {P3}; job3.dataset.phenotypickeys = 1:3; evalc('cat_io_dcm2bids(job3);');
  V3 = sessionValues(B,S,{'score','MMSE'});
  R(end+1) = struct('name','keys-by-number','ok',isequal(cellfun(@(v) v{1},V3,'UniformOutput',false),arrayfun(@(k) sprintf('%d',k),1:numel(S),'UniformOutput',false)) && ...
    isequal(cellfun(@(v) v{2},V3,'UniformOutput',false),cellfun(@(e) e{1},exp1,'UniformOutput',false)), ...
    'info',sprintf('score: %s, MMSE kept: %s',strjoin(cellfun(@(v) v{1},V3,'UniformOutput',false),' '),strjoin(cellfun(@(v) v{2},V3,'UniformOutput',false),' ')));

  % (3) no names and birth dates in the BIDS directory
  txt = lower(strjoin(cellfun(@fileread,[cat_vol_findfiles(B,'*.tsv'); cat_vol_findfiles(B,'*.json')],'UniformOutput',false),' '));
  w   = unique(regexp(lower(strjoin({S.pn},' ')),'[a-z]{3,}','match'));
  bd  = unique([{S.birth} cellfun(birth,num2cell(S),'UniformOutput',false)]);
  found = [w(cellfun(@(x) ~isempty(regexp(txt,['(?<![a-z])' x '(?![a-z])'],'once')),w)) bd(cellfun(@(x) contains(txt,x),bd))];
  R(end+1) = struct('name','no-names-in-BIDS','ok',isempty(found),'info',sprintf('%d name words and birth dates, %d found',numel(w) + numel(bd),numel(found)));

  % (4) corrected table: new value of the first session
  C1{2,5} = '99'; cat_io_csv(P1,C1);
  job4 = job; job4.dataset.phenotypic = {P1}; evalc('cat_io_dcm2bids(job4);');
  V4 = sessionValues(B,S,{'MMSE','Diagnosis','score'});
  R(end+1) = struct('name','corrected-table','ok',strcmp(V4{1}{1},'99') && isequal(cellfun(@(v) v{1},V4(2:end),'UniformOutput',false), ...
    cellfun(@(e) e{1},exp1(2:end),'UniformOutput',false)) && strcmp(V4{2}{2},'dx1') && strcmp(V4{3}{3},'3'), ...
    'info',sprintf('MMSE: %s, diagnoses and scores kept: %d',strjoin(cellfun(@(v) v{1},V4,'UniformOutput',false),' '), ...
    strcmp(V4{2}{2},'dx1') && strcmp(V4{3}{3},'3')));
end
% =========================================================================
function t = germanDate(d)
%germanDate. Date yyyy-mm-dd as dd.mm.yyyy ('' if not available).
  try
    t = char(datetime(d,'InputFormat','yyyy-MM-dd','Format','dd.MM.yyyy'));
  catch
    t = '';
  end
end
% =========================================================================
function [V,ok] = sessionValues(B,S,cols)
%sessionValues. Values of the columns cols of each session S (by the subject
%  of its PatientID and the session of its scan date) in the sessions files
%  of the BIDS directory B (n/a if not available).
  V = cell(1,numel(S)); ok = true;
  sub = dir(fullfile(B,'sub-*')); sub = sub([sub.isdir]);
  for k = 1:numel(S)
    V{k} = repmat({'n/a'},1,numel(cols));
    si = find(endsWith(lower(regexprep({sub.name},'[^A-Za-z0-9]','')),lower(regexprep(S(k).pid,'[^A-Za-z0-9]',''))),1);
    if isempty(si), ok = false; continue; end
    Ps = fullfile(B,sub(si).name,[sub(si).name '_sessions.tsv']);
    if ~exist(Ps,'file'), ok = false; continue; end
    C = cat_io_csv(Ps,'','',struct('delimiter','\t','convert2double',0));
    ri = find(strcmp(C(:,1),['ses-' strrep(S(k).date,'-','')]),1);
    if isempty(ri), ok = false; continue; end
    for ci = 1:numel(cols)
      c = find(strcmp(C(1,:),cols{ci}),1);
      if ~isempty(c), V{k}{ci} = C{ri,c}; end
    end
  end
end
% =========================================================================
function st = reportStatus(P)
%reportStatus. Status of the rows of a match report of the phenotypic data.
  st = {};
  if ~exist(P,'file'), return; end
  C  = cat_io_csv(P,'','',struct('convert2double',0));
  st = C(2:end,strcmp(C(1,:),'status'))';
end
% =========================================================================
function R = testExport(E)
%testExport. Anonymized DICOM export (see cat_io_dcm2bids_export) of a
%  classic (Berlin) and an enhanced (JE7) session: zip files and their
%  content tables (images of the database, series directories, and files
%  of the zip), export tables without personal data, reimport of the
%  export (only the cleaned fields of the sidecars differ, the same BIDS
%  names and participants, defaced images), repeated run (no new export),
%  the protocol scans of a JSON-only import (images converted only for the
%  defacing), the cleaning of a file with a name part in a description,
%  and the defacing of a Siemens mosaic (extensive anonymization).
  R  = struct('name',{},'ok',{},'info',{});
  T  = fullfile(E.Pout,'export'); fresh(T);
  DB = fullfile(T,'study','catDCM2BIDSdb'); X = [DB '-dicom'];
  job = E.job0; job.outdir = {T}; job.data = {E.BER; E.dcmses}; job.dicts = E.dicts; job.opts.dicomexport = 2;

  % (1) import with the export of all scans: one run directory with a zip
  % file and its content table per session
  tic; evalc('cat_io_dcm2bids(job);'); t = toc;
  Pz = cat_vol_findfiles(X,'*.zip'); L = load(fullfile(X,'exports.mat'));
  okc = ~isempty(Pz);
  for zi = 1:numel(Pz) % content table: images of the database with the series directory and files of the zip
    U = tempname; unzip(Pz{zi},U);
    F = strrep(strrep(cat_vol_findfiles(U,'*.dcm'),[U filesep],''),filesep,'/'); rmdir(U,'s');
    C = cat_io_csv(spm_file(Pz{zi},'ext','.csv'),'','',struct('convert2double',0));
    col = @(n) C(2:end,strcmp(C(1,:),n));
    [ser,si] = unique(col('Series')); nd = str2double(col('DICOMfiles')); 
    nf  = cellfun(@(s) nnz(startsWith(F,[s '/'])),ser);
    img = cellfun(@(p) ~isempty(dir(fullfile(DB,fileparts(p),['*=' spm_file(p,'filename')]))),col('Image'));
    okc = okc && sum(nf) == numel(F) && isequal(nf(:),nd(si)) && all(img) && size(C,1) > 1;
  end
  R(end+1) = struct('name','export','ok',isscalar(dir(fullfile(X,'2*'))) && numel(Pz) == 2 && okc && ...
    numel(L.Z) == 2 && sum([L.Z.Scans]) == numel(L.R) && any([L.R.Defaced]), ...
    'info',sprintf('%0.0f s, %d zip files, %d scans (%d defaced), %d files, content tables: %d', ...
    t,numel(Pz),numel(L.R),nnz([L.R.Defaced]),sum([L.R.Files]),okc));

  % (2) no personal data in the tables of the export (parts of the patient
  % names, birth dates, and the DICOM directories)
  txt  = lower(strjoin(cellfun(@fileread,[cat_vol_findfiles(X,'*.csv'); {fullfile(X,'readme.txt')}], ...
    'UniformOutput',false),' '));
  pats = {};
  for P = sidecars(DB)'
    J = cat_io_json(P{1}); I = cat_io_json(spm_file(P{1},'prefix','catDCM2BIDSimport_'));
    if isfield(J,'PatientName') && isfield(J,'PatientID')
      p = strsplit(lower(char(J.PatientName)),{'^',' ','-','_'});
      pats = [pats p(cellfun(@numel,p) >= 4 & ~cellfun(@(q) contains(lower(char(J.PatientID)),q),p))]; %#ok<AGROW>
    end
    if isfield(J,'PatientBirthDate')
      pats = [pats {char(J.PatientBirthDate), strrep(char(J.PatientBirthDate),'-','')}]; %#ok<AGROW>
    end
    pats{end+1} = lower(I.DICOMdirectory); %#ok<AGROW>
  end
  pats = unique(pats); found = pats(cellfun(@(q) contains(txt,q),pats));
  R(end+1) = struct('name','export-private','ok',isempty(found), ...
    'info',sprintf('%d personal strings, %d found',numel(pats),numel(found)));

  % (3) reimport of the export: only the cleaned fields of the sidecars
  % differ, the same BIDS names and participants, the defaced images as
  % the anonymized images of the first import, and other images identical
  T2 = fullfile(E.Pout,'export2'); fresh(T2); U = fullfile(T2,'dicom');
  for zi = 1:numel(Pz), unzip(Pz{zi},U); end
  job2 = E.job0; job2.outdir = {T2}; job2.data = {U}; job2.dicts = E.dicts;
  tic; evalc('cat_io_dcm2bids(job2);'); t = toc;
  DB2 = fullfile(T2,'study','catDCM2BIDSdb');
  S1  = sidecars(DB); S1 = S1(ismember(cellfun(@(P) sidecarKey(P,1),S1,'UniformOutput',false),{L.R.SeriesInstanceUID}));
  S2  = sidecars(DB2);
  [~,i1,i2] = intersect(cellfun(@(P) sidecarKey(P),S1,'UniformOutput',false),cellfun(@(P) sidecarKey(P),S2,'UniformOutput',false));
  dif = {}; dface = 0; nface = 0; same = true;
  for j = 1:numel(i1)
    A = cat_io_json(S1{i1(j)}); B = cat_io_json(S2{i2(j)});
    f = union(fieldnames(A),fieldnames(B));
    dif = [dif; f(cellfun(@(x) ~isfield(A,x) || ~isfield(B,x) || ~isequal(A.(x),B.(x)),f))]; %#ok<AGROW>
    Pa = spm_file(S1{i1(j)},'prefix','anon_','ext','.nii.gz'); P2 = spm_file(S2{i2(j)},'ext','.nii.gz');
    if exist(Pa,'file') % defaced: as the anonymized image of the first import
      Ya = readNiiGz(Pa); Yb = readNiiGz(P2); Ya(isnan(Ya)) = 0;
      dface = max(dface,nnz(Ya ~= Yb) / numel(Ya)); nface = nface + 1;
    elseif dir(P2).bytes < 50e6 % other images (without the large ones)
      same = same && isequal(readNiiGz(spm_file(S1{i1(j)},'ext','.nii.gz')),readNiiGz(P2));
    end
  end
  dif = unique(dif);
  B1  = fullfile(T,'study','BIDS-pseudonymized'); B2 = fullfile(T2,'study','BIDS-pseudonymized');
  pt1 = cat_vol_findfiles(B1,'participants.tsv'); pt2 = cat_vol_findfiles(B2,'participants.tsv');
  pt  = ~isempty(pt1) && ~isempty(pt2) && isequal(fileread(pt1{1}),fileread(pt2{1}));
  bn  = isequal(names(listing(B1)),names(listing(B2)));
  R(end+1) = struct('name','export-reimport','ok',numel(i1) == numel(S1) && ...
    all(ismember(dif,{'PatientName','PatientBirthDate','AccessionNumber','ReferringPhysicianName'})) && ...
    bn && pt && nface > 0 && dface < 1e-4 && same, ...
    'info',sprintf(['%0.0f s, %d of %d sidecars, differing fields: %s, BIDS names equal: %d, participants ' ...
    'equal: %d, defaced images: %d (max %0.4f%% voxels differ), other images identical: %d'], ...
    t,numel(i1),numel(S1),strjoin(dif',' '),bn,pt,nface,100*dface,same));

  % (4) the same import again: no new export
  tic; evalc('cat_io_dcm2bids(job);'); t = toc; L4 = load(fullfile(X,'exports.mat'));
  R(end+1) = struct('name','export-again','ok',isscalar(dir(fullfile(X,'2*'))) && numel(L4.R) == numel(L.R), ...
    'info',sprintf('%0.0f s, run directories: %d, scans: %d',t,numel(dir(fullfile(X,'2*'))),numel(L4.R)));

  % (5) protocol scans (conform or semiconform) of a JSON-only import with
  % only the T1w protocol (the other scans are nonconform): the images are
  % converted only for the defaced scans
  T3 = fullfile(E.Pout,'export3'); fresh(T3); DB3 = fullfile(T3,'study','catDCM2BIDSdb');
  Pp = fullfile(T3,'protocols','anat'); mkdir(Pp);
  copyfile(fullfile(E.dicts.Pprotocoldirs{1},'anat','*T1w*.json'),Pp);
  job3 = job; job3.outdir = {T3}; job3.data = {E.BER}; job3.dicts.Pprotocoldirs = {fileparts(Pp)};
  job3.opts.output = 1; job3.opts.dicomexport = 1;
  tic; evalc('cat_io_dcm2bids(job3);'); t = toc;
  L5 = load(fullfile([DB3 '-dicom'],'exports.mat'));
  Ts = cat_io_csv(fullfile(T3,'study','catDCM2BIDSdb-tables','scans.csv'),'','',struct('convert2double',0));
  sid = Ts(2:end,strcmp(Ts(1,:),'ScanID')); sst = Ts(2:end,strcmp(Ts(1,:),'Status'));
  nexp = numel(unique(regexprep(sid(ismember(sst,{'conform','semiconform'})),'\|.*$','')));
  nimg = numel(cat_vol_findfiles(DB3,'*=*=*=*.nii.gz'));
  R(end+1) = struct('name','export-json-protocols','ok',numel(L5.R) == nexp && nexp < numel(unique(regexprep(sid,'\|.*$',''))) && ...
    any([L5.R.Defaced]) && isequal([L5.R.Defaced] == 1,strcmp({L5.R.Datatype},'anat')) && nimg == nnz([L5.R.Defaced]), ...
    'info',sprintf('%0.0f s, %d of %d series (conform or semiconform), %d defaced, %d images converted', ...
    t,numel(L5.R),numel(unique(regexprep(sid,'\|.*$',''))),nnz([L5.R.Defaced]),nimg));

  % (6) cleaning of a file (JE7 T1w) with a part of the patient name in the
  % study description: the name part is found, the patient name is the
  % PatientID, the birth date is the 15th, the cleaned elements are empty,
  % and all other top-level elements are unchanged (besides sequences with
  % cleaned elements)
  F  = dir(E.dcmseries); F = F(~[F.isdir] & ~startsWith({F.name},'.'));
  fid = fopen(fullfile(E.dcmseries,F(1).name)); b0 = fread(fid,inf,'*uint8')'; fclose(fid);
  E0 = cat_io_dcm2bids_export('dicomElements',b0);
  pn = elementText(b0,E0,'00100010'); pid = elementText(b0,E0,'00100020');
  p  = strsplit(lower(pn),{'^',' ','-','_'}); p = p(cellfun(@numel,p) >= 3 & ~cellfun(@(q) contains(lower(pid),q),p));
  sd = find(E0.depth == 0 & E0.tag == hex2dec('00081030'),1);
  if ~isempty(p) && ~isempty(sd) && E0.len(sd) >= numel(p{1}) + 2
    b0(E0.voff(sd) + (0:numel(p{1}) + 1)) = uint8([' ' upper(p{1}) ' ']);
  end
  [b1,info] = cat_io_dcm2bids_export('anonymizeFile',b0,[],struct('anonymize',1));
  E1  = cat_io_dcm2bids_export('dicomElements',b1);
  cln = hex2dec([cat_io_dcm2bids_export('cleanedElements'); {'00100010'; '00100030'}]);
  chg = false; emp = true;
  for ei = find(E0.depth == 0)'
    ej = find(E1.depth == 0 & E1.tag == E0.tag(ei),1);
    if ismember(E0.tag(ei),cln(1:end-2)), emp = emp && E1.len(ej) == 0; continue; end
    if ~ismember(E0.tag(ei),cln) && ~strcmp(E0.vr(ei,:),'SQ')
      chg = chg || ~isequal(b0(E0.voff(ei) + (0:E0.len(ei) - 1)),b1(E1.voff(ej) + (0:E1.len(ej) - 1)));
    end
  end
  bd = elementText(b1,E1,'00100030');
  R(end+1) = struct('name','export-clean','ok',~isempty(p) && any(strcmp(info.names,'(0008,1030)')) && ...
    strcmp(elementText(b1,E1,'00100010'),pid) && endsWith(bd,'15') && emp && ~chg, ...
    'info',sprintf('name parts: %d, found in: %s, patient name = PatientID: %d, birth date 15th: %d, cleaned empty: %d, others unchanged: %d', ...
    numel(p),strjoin(info.names,' '),strcmp(elementText(b1,E1,'00100010'),pid),endsWith(bd,'15'),emp,~chg));

  % (7) extensive anonymization: the defaced Siemens mosaic (field map) of
  % the export as the anonymized image of the BIDS export
  T4 = fullfile(E.Pout,'export4'); fresh(T4);
  Ds = dir(fullfile(E.BER,'**','*SpinEchoFieldMap*')); Ds = Ds([Ds.isdir]);
  job4 = E.job0; job4.outdir = {T4}; job4.data = {fullfile(Ds(1).folder,Ds(1).name)};
  job4.opts.anonymize = 2; job4.opts.dicomexport = 2;
  evalc('cat_io_dcm2bids(job4);');
  Pz4 = cat_vol_findfiles(fullfile(T4,'study','catDCM2BIDSdb-dicom'),'*.zip');
  Pb  = cat_vol_findfiles(fullfile(T4,'study','BIDS-anonymized'),'*.nii.gz');
  Pr  = cat_vol_findfiles(fullfile(T4,'study','catDCM2BIDSdb'),'*=*=*=*.nii.gz');
  Pr  = Pr(~startsWith(spm_file(Pr,'filename'),'anon_'));
  if ~isempty(Pz4) && ~isempty(Pb) && ~isempty(Pr)
    U = fullfile(T4,'unzip'); unzip(Pz4{1},U); C = fullfile(T4,'conv'); mkdir(C);
    system(sprintf('"%s" -b n -z n -f out -o "%s" "%s" > /dev/null 2>&1',E.dcm2niix,C,U));
    Ye = spm_read_vols(spm_vol(fullfile(C,'out.nii'))); Ya = readNiiGz(Pb{1}); Ya(isnan(Ya)) = 0;
    dm = nnz(abs(Ye - Ya) > 0.5) / numel(Ye); z = nnz(Ye == 0 & readNiiGz(Pr{1}) ~= 0) / numel(Ye);
  else
    dm = inf; z = 0;
  end
  R(end+1) = struct('name','export-mosaic','ok',dm < 1e-3 && z > 0.01, ...
    'info',sprintf('voxels set to 0 by the defacing: %0.1f%%, voxels differing from the BIDS image: %0.3f%%',100*z,100*dm));

  % (8) anonymization level 2 (T1w, field map, and DWI reference scan of 
  % Berlin): subject label with the site (without siteID only manufacturer
  % and field strength) and number and session label of the time since the
  % first session (weeks), without personal data and location (city, 
  % institution) in the text files of the BIDS and report directories, 
  % without scanner model and coils and with the type of the software 
  % version in the sidecars and report protocols, and without acquisition 
  % times in the NIfTI headers
  T5 = fullfile(E.Pout,'export5'); fresh(T5);
  D5 = [dir(fullfile(E.BER,'**','*mprage*')); dir(fullfile(E.BER,'**','*SpinEchoFieldMap_PA*')); dir(fullfile(E.BER,'**','*refscan*'))];
  D5 = D5([D5.isdir]);
  job5 = E.job0; job5.outdir = {T5}; job5.data = fullfile({D5.folder},{D5.name})';
  job5.opts.anonymize = 2; job5.opts.dicomexport = 2;
  evalc('cat_io_dcm2bids(job5);');
  A5  = fullfile(T5,'study','BIDS-anonymized'); B5 = fullfile(A5,'BIDS');
  sub = dir(fullfile(B5,'sub-*')); sub = sub([sub.isdir]); ses = dir(fullfile(B5,'sub-*','ses-*'));
  S5  = sidecars(fullfile(T5,'study','catDCM2BIDSdb')); J5 = cat_io_json(S5{1});
  dat = regexprep(char(J5.AcquisitionDateTime),'T.*$','');
  pats = {char(J5.PatientID),char(J5.DeviceSerialNumber),char(J5.PatientBirthDate),strrep(char(J5.PatientBirthDate),'-',''), ...
    dat,strrep(dat,'-','')};
  p = strsplit(lower(char(J5.PatientName)),{'^',' ','-','_'}); pats = [pats p(cellfun(@numel,p) >= 4)];
  addr = strtrim(strsplit(char(fieldText5(J5,'InstitutionAddress')),',')); % location and scanner (level 2)
  loc  = {fieldText5(J5,'InstitutionName'),addr{min(2,end)},fieldText5(J5,'ManufacturersModelName'), ...
    fieldText5(J5,'ReceiveCoilName'),fieldText5(J5,'CoilString')};
  pats = [pats loc(cellfun(@numel,loc) >= 4)];
  pats = unique(pats(~cellfun('isempty',pats)));
  sw   = fieldText5(J5,'SoftwareVersions'); % (only in the DICOM export)
  Pt  = {}; % text files of the BIDS and report directories (shareable)
  for Pd = {A5,[A5 '-report']}, for ext = {'*.tsv','*.json','*.csv','*.txt'}, Pt = [Pt; cat_vol_findfiles(Pd{1},ext{1})]; end, end %#ok<AGROW>
  txt = lower(strjoin(cellfun(@fileread,Pt,'UniformOutput',false),' '));
  found = [pats {sw}]; found = found(cellfun(@(q) numel(q) >= 4 && contains(txt,lower(q)),found));
  Pj5 = [cat_vol_findfiles(B5,'sub-*.json'); cat_vol_findfiles(fullfile([A5 '-report'],'protocols'),'*.json')]; nscn = 0; swt = {};
  for k = 1:numel(Pj5)
    J = cat_io_json(Pj5{k}); 
    nscn = nscn + any(isfield(J,{'ManufacturersModelName','ReceiveCoilName','CoilString','DeviceSerialNumber','InstitutionName'}));
    if isfield(J,'SoftwareVersions'), swt{end+1} = char(J.SoftwareVersions); end %#ok<AGROW>
  end
  swt = unique(swt); 
  acq5 = 'n/a'; % acquisition time of the session (date of the label, 3-hour interval of the day)
  if isscalar(sub)
    St5 = cat_io_csv(fullfile(B5,sub(1).name,[sub(1).name '_sessions.tsv']),'','',struct('delimiter','\t'));
    if any(strcmp(St5(1,:),'acq_time')), acq5 = char(string(St5{2,strcmp(St5(1,:),'acq_time')})); end
  end
  Pn = cat_vol_findfiles(B5,'*.nii.gz'); tim = sum(cellfun(@(P) contains(niiDescrip(P),'Time='),Pn));
  okb = isscalar(sub) && ~isempty(regexp(sub(1).name,'^sub-[A-Z]{2}\d{2}T\d{5}$','once')) && ...
    ~isempty(ses) && all(strcmp({ses.name},'ses-w0000')) && isempty(found) && tim == 0 && ~isempty(Pn) && ...
    nscn == 0 && ~isempty(swt) && all(cellfun(@(x) ~isempty(regexp(x,'^[A-Z]{1,3}\d*$','once')),swt)) && ...
    ~isempty(regexp(acq5,'^1900-01-01T(01|04|07|10|13|16|19|22):30:00$','once'));
  if isempty(sub), sname = 'none'; else, sname = sub(1).name; end
  R(end+1) = struct('name','anonymized-bids','ok',okb,'info',sprintf(['%s/%s, personal strings found: %d of %d, ' ...
    'NIfTI with time: %d of %d, files with scanner fields: %d of %d, software %s, acq_time %s'],sname,strjoin({ses.name},' '), ...
    numel(found),numel(pats),tim,numel(Pn),nscn,numel(Pj5),strjoin(swt,'/'),acq5));

  % (8b) siteID of the site dictionary: site label with the siteID (e.g.
  % 07SI30T) in the subject labels and sessions files of the level 2
  T7 = fullfile(E.Pout,'export7'); fresh(T7); mkdir(T7); 
  Pd = fullfile(T7,'sites.json'); 
  cat_io_json(Pd,struct('DeviceSerialNumber',char(J5.DeviceSerialNumber),'InstitutionAbbreviation','Tst','siteID','07')); 
  job7 = job5; job7.outdir = {T7}; job7.opts.dicomexport = 0; 
  job7.dicts = struct('Pprotocoldirs',{{''}},'Pcenterdict',{{Pd}}); 
  evalc('cat_io_dcm2bids(job7);'); 
  sub7 = dir(fullfile(T7,'study','BIDS-anonymized','BIDS','sub-*')); sub7 = sub7([sub7.isdir]); site7 = '';
  if isscalar(sub7)
    St = cat_io_csv(fullfile(sub7.folder,sub7.name,[sub7.name '_sessions.tsv']),'','',struct('delimiter','\t'));
    site7 = char(string(St{2,strcmp(St(1,:),'site')}));
  end
  if isempty(sub7), sname7 = 'none'; else, sname7 = sub7(1).name; end
  R(end+1) = struct('name','anonymized-siteID','ok',isscalar(sub7) && ~isempty(regexp(sname7,'^sub-07[A-Z]{2}\d{2}T\d{5}$','once')) && ...
    ~isempty(regexp(site7,'^07[A-Z]{2}\d{2}T$','once')),'info',sprintf('%s, site %s',sname7,site7));

  % (9) DICOM export of the level 2: subject label, the date of the session
  % label (1900-01-01 for the first session), times since the start of the
  % session at the center of its 3-hour interval of the day (as acq_time of
  % the BIDS export), the shifted birth date on the 15th of the month, and 
  % site in the header, no original IDs, dates, or
  % UIDs in the files, and its import (level 1) gives the subject and 
  % session labels, age, and slice timing of (8)
  Pz5 = cat_vol_findfiles(fullfile(T5,'study','catDCM2BIDSdb-dicom'),'*.zip');
  U5  = fullfile(T5,'unzip'); for zi = 1:numel(Pz5), unzip(Pz5{zi},U5); end
  F5  = cat_vol_findfiles(U5,'*.dcm');
  orig = [pats {char(J5.StudyInstanceUID),char(J5.SeriesInstanceUID)}]; nf = 0; % (also location, model, and coils)
  for k = 1:numel(F5)
    fid = fopen(F5{k}); t = lower(char(fread(fid,inf,'*uint8')')); fclose(fid);
    nf = nf + any(cellfun(@(q) contains(t,lower(q)),orig));
  end
  hp = ''; hd = ''; hs = ''; ht = ''; hb = ''; at = nan(1,numel(F5));
  for k = 1:numel(F5) % times since the start of the session (acquisition times of all files)
    fid = fopen(F5{k}); b = fread(fid,inf,'*uint8')'; fclose(fid); El = cat_io_dcm2bids_export('dicomElements',b);
    at(k) = str2double(regexprep(elementText(b,El,'00080032'),'\..*$',''));
    if k == 1
      hp = elementText(b,El,'00100020'); hd = elementText(b,El,'00080020'); hs = elementText(b,El,'00181000');
      ht = elementText(b,El,'00080030'); hb = elementText(b,El,'00100030');
    end
  end
  T6 = fullfile(E.Pout,'export6'); fresh(T6);
  job6 = E.job0; job6.outdir = {T6}; job6.data = {U5}; evalc('cat_io_dcm2bids(job6);');
  B6 = fullfile(T6,'study','BIDS-pseudonymized','BIDS');
  sub6 = dir(fullfile(B6,'sub-*')); sub6 = sub6([sub6.isdir]); ses6 = dir(fullfile(B6,'sub-*','ses-*'));
  same = isequal({sub6.name},{sub.name}) && isequal(sort({ses6.name}),sort({ses.name}));
  a5 = cat_io_csv(fullfile(B5,'participants.tsv'),'','',struct('delimiter','\t'));
  a6 = cat_io_csv(fullfile(B6,'participants.tsv'),'','',struct('delimiter','\t'));
  sage = isequal(round(a5{2,strcmp(a5(1,:),'age')}),round(a6{2,strcmp(a6(1,:),'age')}));
  % slice timing of the import equal to the anonymized BIDS export (sidecars with the same names)
  Pj = cat_vol_findfiles(B5,'*.json'); nst = 0; dst = 0;
  for k = 1:numel(Pj)
    J1 = cat_io_json(Pj{k}); P2 = strrep(Pj{k},B5,B6);
    if ~isfield(J1,'SliceTiming') || ~exist(P2,'file'), continue; end
    J2 = cat_io_json(P2); nst = nst + 1;
    if isfield(J2,'SliceTiming') && isequal(size(J1.SliceTiming),size(J2.SliceTiming))
      dst = max(dst,max(abs(J1.SliceTiming(:) - J2.SliceTiming(:))));
    else
      dst = inf;
    end
  end
  R(end+1) = struct('name','anonymized-dicom','ok',~isempty(F5) && nf == 0 && strcmp(['sub-' hp],sname) && ...
    strcmp(hd,'19000101') && numel(ht) >= 6 && numel(acq5) >= 19 && strcmp(ht(1:6),acq5([12 13 15 16 18 19])) && ...
    min(at) >= floor(str2double(ht)) && numel(hb) == 8 && strcmp(hb(7:8),'15') && str2double(hb(1:4)) < 1900 && ...
    ~isempty(regexp(hs,'^[A-Z]{2}\d{2}T$','once')) && same && sage && nst > 0 && dst < 1e-4, ...
    'info',sprintf(['%d files (%d with original strings), StudyDate %s, StudyTime %s, last acquisition %06.0f, scanner %s, ' ...
    'birth date on the 15th before 1900: %d, labels of the import equal: %d, age equal: %d, slice timing of %d scans equal: %d'], ...
    numel(F5),nf,hd,ht,max(at),hs,numel(hb) == 8 && strcmp(hb(7:8),'15') && str2double(hb(1:4)) < 1900,same,sage,nst,dst < 1e-4));
end
% =========================================================================
function v = fieldText5(S,fn)
%fieldText5. Text of the field fn of the structure S ('' if not available).
  if isfield(S,fn) && ~isempty(S.(fn)) && ~isnumeric(S.(fn)), v = char(S.(fn)); else, v = ''; end
end
% =========================================================================
function d = niiDescrip(P)
%niiDescrip. Description of the header of a NIfTI-1 file (also gzipped).
  if endsWith(P,'.gz')
    gz = java.util.zip.GZIPInputStream(java.io.FileInputStream(P)); h = zeros(1,348,'uint8');
    for k = 1:348, x = gz.read(); if x < 0, break; end; h(k) = x; end
    gz.close();
  else
    fid = fopen(P); h = fread(fid,348,'*uint8')'; fclose(fid);
  end
  d = char(h(149:228)); d = d(d ~= char(0));
end
% =========================================================================
function k = sidecarKey(P,~)
%sidecarKey. SeriesInstanceUID and image suffix (last part of the name) of
%  a sidecar of the database (only the UID with a second argument).
  J = cat_io_json(P); k = J.SeriesInstanceUID;
  if nargin < 2, n = strsplit(spm_file(P,'basename'),'='); k = [k '|' n{end}]; end
end
% =========================================================================
function Y = readNiiGz(P)
%readNiiGz. Image of a NIfTI (also gzipped, unzipped into a temporary
%  directory).
  T = tempname; mkdir(T); c = onCleanup(@() rmdir(T,'s'));
  if endsWith(P,'.gz'), Pu = gunzip(P,T); P = Pu{1}; end
  Y = spm_read_vols(spm_vol(P));
end
% =========================================================================
function v = elementText(b,E,tag)
%elementText. Text value of a top-level DICOM element (see dicomElements
%  of cat_io_dcm2bids_export).
  ei = find(E.depth == 0 & E.tag == hex2dec(tag),1);
  if isempty(ei) || E.len(ei) == 0, v = ''; return; end
  v = char(b(E.voff(ei) + (0:E.len(ei) - 1))); v = strtrim(v(v ~= char(0)));
end
% =========================================================================
function R = testParallel(E)
%testParallel. Subject packages by two parallel processes (nproc 2, see 
%  cat_io_dcm2bids_parallel) and by one process (nproc 0) of the JE7 
%  sessions (two subjects) with QC, preprocessing, and rendering (real 
%  dcm2niix): the same BIDS files (raw data and derivatives), tables of the
%  database (scans, subjects, sites, protocol variants), participants and 
%  sessions files, quality ratings, and rendered pages, without remaining
%  queue directory, and the processing time of both runs. 
  R = struct('name',{},'ok',{},'info',{});
  O = fullfile(E.Pout,'parallel'); fresh(O); 
  run = {'serial',0; 'parallel',2}; t = zeros(1,2); txt = cell(1,2); 
  for ri = 1:2
    job = struct('outdir',{{fullfile(O,run{ri,1})}},'subdir','study','dicts',E.dicts,'nproc',run{ri,2}); 
    job.data = E.dcm; 
    job.opts = struct('runqc',1,'preprocessing',1,'Pdcm2nii',E.dcm2niix,'reportfig',struct('run',0)); %#ok<STRNU>
    tic; txt{ri} = evalc('cat_io_dcm2bids(job);'); t(ri) = toc; 
  end
  S = fullfile(O,'serial','study'); P = fullfile(O,'parallel','study'); 
  rel = @(L,root) strrep(L,root,'<out>'); 

  % BIDS files (names of raw data and derivatives, contents of the raw 
  % sidecars and TSV files)
  Ls = listing(fullfile(S,'BIDS-pseudonymized')); Lp = listing(fullfile(P,'BIDS-pseudonymized')); 
  same = isequal(names(Ls),names(Lp)); 
  tx  = @(L) sel(L(cellfun('isempty',regexp(L,'derivatives','once'))),'\.(json|tsv)$'); 
  ndiff = numel(setxor(tx(Ls),tx(Lp))); 
  R(end+1) = struct('name','BIDS files','ok',same && ndiff == 0 && ~isempty(Ls), ...
    'info',sprintf('%d files (%d parallel), %d different raw sidecars/TSV',numel(Ls),numel(Lp),ndiff/2)); 

  % tables of the database and participants (with the output directory)
  Ts = cat_vol_findfiles(fullfile(S,'catDCM2BIDSdb-tables'),'*.*'); nd = 0; 
  Ts = Ts(~endsWith(Ts,'.mat')); % (MAT files with their creation time)
  Ts = [Ts; cat_vol_findfiles(fullfile(S,'BIDS-pseudonymized'),'participants.tsv')]; 
  for ti = 1:numel(Ts)
    Tp = strrep(Ts{ti},S,P); 
    nd = nd + (~exist(Tp,'file') || ~strcmp(rel(fileread(Ts{ti}),S),rel(fileread(Tp),P))); 
  end
  np = nnz(~endsWith(cat_vol_findfiles(fullfile(P,'catDCM2BIDSdb-tables'),'*.*'),'.mat')) + ...
    numel(cat_vol_findfiles(fullfile(P,'BIDS-pseudonymized'),'participants.tsv')); 
  R(end+1) = struct('name','tables','ok',nd == 0 && np == numel(Ts) && numel(Ts) > 3, ...
    'info',sprintf('%d tables, %d different',numel(Ts),nd)); 

  % quality ratings of the scans (QC JSON)
  Qs = cat_vol_findfiles(fullfile(S,'BIDS-pseudonymized'),'catDCM2BIDSqc_*.json'); dq = nan(1,numel(Qs)); 
  for qi = 1:numel(Qs)
    Qp = strrep(Qs{qi},S,P); 
    if exist(Qp,'file')
      a = cat_io_json(Qs{qi}); b = cat_io_json(Qp); 
      if isfield(a,'qualityrating') && isfield(b,'qualityrating') && isnumeric(a.qualityrating.SQR) && isnumeric(b.qualityrating.SQR)
        dq(qi) = abs(double(a.qualityrating.SQR) - double(b.qualityrating.SQR)); 
      end
    end
  end
  R(end+1) = struct('name','quality ratings','ok',~isempty(Qs) && all(dq(~isnan(dq)) < 0.05) && nnz(isnan(dq)) == 0, ...
    'info',sprintf('%d scans, max SQR difference %0.3f',numel(Qs),max([0 dq]))); 

  % rendered pages (overview, rows, subjects, sessions)
  Gs = names(listing(fullfile(S,'BIDS-pseudonymized-report','render'))); 
  Gp = names(listing(fullfile(P,'BIDS-pseudonymized-report','render'))); 
  pg = @(L) L(cellfun('isempty',regexp(L,'/(cache|rows)/','once'))); 
  R(end+1) = struct('name','rendering','ok',isequal(pg(Gs),pg(Gp)) && ~isempty(sel(Gp,'/subjects/p01_')), ...
    'info',sprintf('%d pages (%d parallel), %d subject pages',numel(pg(Gs)),numel(pg(Gp)),numel(sel(Gp,'/subjects/p')))); 

  % parallel processing: packages of the workers, no queue directory, time
  npk = regexp(txt{2},'(\d+) of (\d+) packages processed','tokens'); 
  R(end+1) = struct('name','parallel run','ok',~isempty(npk) && strcmp(npk{end}{1},npk{end}{2}) && ...
    isempty(dir(fullfile(P,'catDCM2BIDSdb','catDCM2BIDSpar_*'))) && ~contains(txt{2},'failed'), ...
    'info',sprintf('%0.0f s (one process %0.0f s, ratio %0.2f)',t(2),t(1),t(2)/t(1))); 
end
% =========================================================================
function R = testLevels(E,levels)
%testLevels. Preprocessing levels (see levelParameters in 
%  cat_io_dcm2bids_pp_anat) of the Berlin pilot session: one import (with
%  QC) and the preprocessing of each level (rerun) with the processing time
%  of the segmentation, DTI, and rs-fMRI (duration in their JSON files) and
%  the deviation of their measures from level 1 (relative difference of the
%  volumes, Dice of the label map, absolute difference of the DTI and 
%  rs-fMRI measures, and the correlation of the FA maps and connectivity 
%  matrices). The results of each level are saved in 
%  Pout/pplevels/level<L>.mat (level 1 is only processed once for each
%  version of the results, see levelParameters) and the
%  table in Pout/pplevels/pplevels.csv. Moreover, the coregistration to the
%  T1w of a synthetic EPI with a known offset is tested (see coregAnat) and
%  that the DTI includes only voxels with signal in all volumes (see
%  estimateDTI in cat_io_dcm2bids_pp_dwi).
  R = struct('name',{},'ok',{},'info',{});
  O = fullfile(E.Pout,'pplevels'); if ~exist(O,'dir'), mkdir(O); end
  job0 = struct('outdir',{{O}},'subdir','study','dicts',E.dicts);
  job0.opts = struct('runqc',1,'preprocessing',0,'render',struct('run',0),'reportfig',struct('run',0));
  DB = fullfile(O,'study','catDCM2BIDSdb');
  if ~exist(DB,'dir') % import with QC (once)
    job = job0; job.data = {E.BER}; evalc('cat_io_dcm2bids(job);'); %#ok<STRNU>
  end
  levels = sort(unique([1 levels(:)']),'descend'); % fast levels first (early errors), the reference last
  for li = 1:numel(levels)
    Pl = fullfile(O,sprintf('level%s.mat',strrep(sprintf('%g',levels(li)),'.','p')));
    p  = cat_io_dcm2bids_pp_anat('levelParameters',levels(li)); S = [];
    if levels(li) == 1 && exist(Pl,'file'), load(Pl,'S'); end
    if ~isfield(S,'version') || S.version < p.version % (level 1 only once for each version of the results)
      job = job0; job.data = {DB}; job.opts.preprocessing = 2; job.opts.rerun = 1; job.opts.pplevel = levels(li); %#ok<STRNU>
      t0 = datetime('now'); t = tic; evalc('cat_io_dcm2bids(job);'); ttotal = toc(t);
      S = levelResults(DB,O,levels(li),t0); S.ttotal = ttotal; S.version = p.version;
      save(Pl,'S');
    end
    L(li) = S; %#ok<AGROW>
  end

  % table of the deviations from level 1 (sorted by level)
  [~,so] = sort([L.level]); L = L(so); ref = L(1);
  hdr = {'level','time_s','time_ratio','total_s','TIV_pct','GMV_pct','WMV_pct','DiceGM','DiceWM', ...
    'FA_WM','MD_WM','r_FA','tSNR','GCOR','within','between','segregation','r_conn'};
  C = cell(numel(L),numel(hdr));
  for li = 1:numel(L)
    S = L(li);
    C(li,:) = {S.level, S.time, S.time / ref.time, S.ttotal, ...
      100 * (S.anat.TIV / ref.anat.TIV - 1), 100 * (S.anat.GMV / ref.anat.GMV - 1), 100 * (S.anat.WMV / ref.anat.WMV - 1), ...
      dice(S.anat.c0,ref.anat.c0,2), dice(S.anat.c0,ref.anat.c0,3), ...
      S.dwi.FA_WM - ref.dwi.FA_WM, 1000 * (S.dwi.MD_WM - ref.dwi.MD_WM), mapcorr(S.dwi.FA,ref.dwi.FA), ...
      S.func.tSNR - ref.func.tSNR, S.func.GCOR - ref.func.GCOR, S.func.within - ref.func.within, ...
      S.func.between - ref.func.between, S.func.segregation - ref.func.segregation, mapcorr(S.func.conn,ref.func.conn)};
    okl = all(isfinite([S.time S.anat.TIV S.dwi.FA_WM S.func.tSNR]));
    R(end+1) = struct('name',sprintf('level %g',S.level),'ok',okl, ...
      'info',sprintf('%5.0f s (%5.3f), TIV %+5.1f%%, Dice GM %0.3f, FA %+0.3f, r(FA) %0.3f, tSNR %+5.1f, segregation %+0.3f, r(conn) %0.3f', ...
      C{li,[2 3 5 8 10 12 13 17 18]})); %#ok<AGROW>
  end
  cat_io_dcm2bids_helper('writeTable',fullfile(O,'pplevels.csv'),[hdr; C],'csv');

  % coregistration to the T1w (see coregAnat in cat_io_dcm2bids_pp_anat) of
  % a synthetic EPI (T2-like image of the segmentation with 2 mm voxels) 
  % with a known offset of 40 mm and 10 degree (e.g. a repositioning of the
  % subject between the scans): recovered within 2 mm and 2 degree with a
  % Dice of the brain masks above 0.9 (and a median displacement of the 
  % non-linear coregistration below 1 mm), whereas the scanner alignment
  % fails. Moreover, the non-linear coregistration (see nonlinCoreg) of the
  % synthetic EPI with a smooth distortion along y (up to 4 mm in opposite
  % directions of both hemispheres, i.e. not rigid, e.g. of an EPI without
  % distortion correction) recovers the true positions of the displaced 
  % brain voxels (median residual of the voxels with more than 1 mm below 
  % 1 mm and below 60% of the residual of the rigid coregistration).
  Pn = cat_vol_findfiles(DB,'*=*=*=*.nii.gz');
  A  = cat_io_dcm2bids_pp_anat('sessionAnat',Pn,repmat({'anat'},size(Pn)),repmat({'T1w'},size(Pn)));
  ok = false; info = 'no T1w segmentation'; okn = false; infon = info;
  if ~isempty(A)
    o  = struct('gzipi',1); Pc = cat_io_dcm2bids_helper('prepNii',{A.Pc0},o,0); Vc = spm_vol(Pc{1}); Yc0 = spm_read_vols(Vc);
    cat_io_dcm2bids_helper('prepNiigz',Pc,o);
    Vs = struct('dim',ceil(Vc.dim .* sqrt(sum(Vc.mat(1:3,1:3).^2)) / 2),'mat',[]); 
    Vs.mat = Vc.mat * diag([2 ./ sqrt(sum(Vc.mat(1:3,1:3).^2)) 1]); Vs.mat(1:3,4) = Vc.mat(1:3,4);
    [x,y,z] = ndgrid(1:Vs.dim(1),1:Vs.dim(2),1:Vs.dim(3)); W0 = Vs.mat * [x(:)'; y(:)'; z(:)'; ones(1,numel(x))];
    lab = @(W) reshape(spm_sample_vol(Yc0,subsref(Vc.mat \ W,substruct('()',{1,':'})),subsref(Vc.mat \ W,substruct('()',{2,':'})), ...
      subsref(Vc.mat \ W,substruct('()',{3,':'})),1),Vs.dim);
    syn = @(L) interp1(0:3,[0 1 0.55 0.35],min(3,max(0,L))); % CSF bright, GM medium, WM dark (b0-like)
    L  = lab(W0); Ys = syn(L);
    Mt = spm_matrix([5 -10 40 10*pi/180 3*pi/180 -4*pi/180]); % true offset of the synthetic scan
    Vs.mat = Mt \ Vs.mat;                                    % scan world coordinates (shifted)
    evalc('[M,Q] = cat_io_dcm2bids_pp_anat(''coregAnat'',Ys,Vs,A,o);');
    q  = spm_imatrix(M / Mt); e = [norm(q(1:3)) max(abs(q(4:6))) * 180/pi]; % residual error
    dsp = nan; if isfield(Q,'displacement'), dsp = Q.displacement(1); end
    ok = e(1) < 2 && e(2) < 2 && Q.Dice > 0.9 && Q.Dice_scanner < 0.8 && dsp < 1;
    info = sprintf('residual %0.1f mm, %0.1f deg, Dice %0.3f (scanner %0.3f), non-linear displacement %0.1f mm',e,Q.Dice,Q.Dice_scanner,dsp);
    % smooth distortion along y (T1w world) of the synthetic scan
    u  = 4 * sin(2*pi * x(:)' / Vs.dim(1)) .* sin(pi * z(:)' / Vs.dim(3));
    Wt = W0 + [zeros(1,numel(u)); u; zeros(2,numel(u))];  % true T1w world positions of the scan voxels
    Ld = lab(Wt); b = Ld(:) > 0.5 & abs(u(:)) > 1;           % displaced brain voxels
    evalc('[M2,Q2,D2] = cat_io_dcm2bids_pp_anat(''coregAnat'',syn(Ld),Vs,A,o);');
    Wr = M2 * Vs.mat * [x(b)'; y(b)'; z(b)'; ones(1,nnz(b))]; rr = median(sqrt(sum((Wr(1:3,:) - Wt(1:3,b)).^2,1)));
    rn = nan; dn = nan; if isfield(Q2,'Dice_nonlinear'), dn = Q2.Dice_nonlinear; end
    if ~isempty(D2), Yi = reshape(D2.Yi,[],3)'; rn = median(sqrt(sum((double(Yi(:,b)) - Wt(1:3,b)).^2,1))); end
    okn  = rn < 1 && rn < 0.6 * rr;
    infon = sprintf('distortion up to 4 mm: residual median %0.2f mm (rigid %0.2f mm), Dice %0.3f (rigid %0.3f)',rn,rr,dn,Q2.Dice);
  end
  R(end+1) = struct('name','coregistration','ok',ok,'info',info);
  R(end+1) = struct('name','coregistration-nonlinear','ok',okn,'info',infon);

  % DTI only of the voxels with signal in all volumes (see estimateDTI in
  % cat_io_dcm2bids_pp_dwi), i.e. not of the edge slices of the field of
  % view that are missing in some volumes after the realignment (with high
  % FA in the direction of the slices; the realigned volumes are resliced 
  % in memory by the realignment parameters of the QC)
  ok = false; info = 'no FA map and realignment of the DWI';
  Pf = cat_vol_findfiles(DB,'catDCM2BIDSdwi_FA_0*.nii.gz');
  if ~isempty(Pf)
    nam = regexprep(spm_file(Pf{1},'filename'),'^catDCM2BIDSdwi_FA_|\.nii\.gz$','');
    Pd  = fullfile(fileparts(Pf{1}),[nam '.nii.gz']); Prp = fullfile(fileparts(Pf{1}),['rp_' nam '.txt']);
    if exist(Pd,'file') && exist(Prp,'file')
      FA = readNii(Pf(1),O,0); Pt = fullfile(O,'tmp_r.nii.gz'); copyfile(Pd,Pt); Pu = gunzip(Pt); delete(Pt);
      Vr = cat_io_dcm2bids_helper('realignedVolumes',Pu{1},Prp,struct('hrrealign',0)); z = false(Vr(1).dim);
      for vi = 1:numel(Vr), z = z | ~(cat_io_dcm2bids_helper('realignedVolume',Vr,vi) > 0); end
      delete(Pu{1});
      J  = cat_io_json(fullfile(fileparts(Pf{1}),['catDCM2BIDSdwi_' nam '.json'])); nx = nan;
      if isfield(J,'qc') && isfield(J.qc,'voxels_without_signal'), nx = J.qc.voxels_without_signal; end
      nz = nnz(isfinite(FA) & z); ok = nz <= 0.001 * nnz(isfinite(FA)) && isfinite(nx);
      info = sprintf('FA of %d voxels without signal in some volumes (%g brain voxels excluded)',nz,nx);
    end
  end
  R(end+1) = struct('name','dti-signal','ok',ok,'info',info);
end
% =========================================================================
function S = levelResults(DB,O,level,t0)
%levelResults. Processing time (sum of the durations of all results that 
%  were written after t0, i.e. by the run of the level) and measures of the
%  first scan of each type of the session of the level test (i.e. of the 
%  scans of all levels).
  S  = struct('level',level,'time',0,'ttotal',nan,'anat',struct(),'dwi',struct(),'func',struct());
  Pj = [cat_vol_findfiles(DB,'catDCM2BIDSsegus_*.json'); cat_vol_findfiles(DB,'catDCM2BIDSdwi_0*.json'); ...
    cat_vol_findfiles(DB,'catDCM2BIDSrsfmri_0*.json')];
  for ji = 1:numel(Pj)
    D = dir(Pj{ji}); 
    if datetime(D.datenum,'ConvertFrom','datenum') >= t0
      J = cat_io_json(Pj{ji}); S.time = S.time + J.method.duration; 
    end
  end
  % segmentation of the (first) anatomical scan
  Pj = cat_vol_findfiles(DB,'catDCM2BIDSsegus_*.json'); J = cat_io_json(Pj{1});
  S.anat = struct('TIV',J.TIV,'GMV',J.aGMV,'WMV',J.aWMV,'c0',readNii(cat_vol_findfiles(fileparts(Pj{1}),'c0anon_*.nii.gz'),O,level));
  % DTI and rs-fMRI of the first scans
  Pd = cat_vol_findfiles(DB,'catDCM2BIDSdwi_0*.json'); J = cat_io_json(Pd{1});
  S.dwi  = struct('FA_WM',J.global.FA_WM,'MD_WM',J.global.MD_WM, ...
    'FA',readNii({strrep(strrep(Pd{1},'catDCM2BIDSdwi_','catDCM2BIDSdwi_FA_'),'.json','.nii.gz')},O,level));
  Pf = cat_vol_findfiles(DB,'catDCM2BIDSrsfmri_0*.json'); J = cat_io_json(Pf{1});
  Pc = cat_vol_findfiles(fileparts(Pf{1}),'catDCM2BIDSrsfmri_conn_*.tsv'); C = readcell(Pc{1},'FileType','text','Delimiter','\t');
  C  = C(2:end,2:end); C(~cellfun(@(x) isnumeric(x) && isscalar(x),C)) = {nan}; % (n/a or empty)
  S.func = struct('tSNR',J.qc.tSNR,'GCOR',J.connectivity.GCOR,'within',J.connectivity.within, ...
    'between',J.connectivity.between,'segregation',J.connectivity.segregation,'conn',cell2mat(C));
end
% =========================================================================
function Y = readNii(P,O,level)
%readNii. Image of the (zipped) NIfTI P{1} (temporary copy in O).
  Pt = fullfile(O,sprintf('tmp_%g.nii.gz',level)); copyfile(P{1},Pt);
  Pu = gunzip(Pt); Y = single(spm_read_vols(spm_vol(Pu{1}))); delete(Pu{1}); delete(Pt);
end
% =========================================================================
function d = dice(A,B,c)
%dice. Dice coefficient of the label c of the label maps A and B (rounded).
  if ~isequal(size(A),size(B)), d = nan; return; end
  A = round(A) == c; B = round(B) == c; d = 2 * nnz(A & B) / (nnz(A) + nnz(B));
end
% =========================================================================
function r = mapcorr(A,B)
%mapcorr. Correlation of two maps/matrices of the same size (finite values).
  if ~isequal(size(A),size(B)), r = nan; return; end
  m = isfinite(A) & isfinite(B); if nnz(m) < 3, r = nan; return; end
  c = corrcoef(double(A(m)),double(B(m))); r = c(1,2);
end
% =========================================================================
function [sA,sB,sesA] = dbSubjects(E)
%dbSubjects. Database subjects of the import (sA with two sessions, sB
%  with one) and the sessions of sA sorted by their date.
  D  = dir(fullfile(E.DB,'sub-*')); D = D([D.isdir]);
  ns = arrayfun(@(d) numel(dir(fullfile(E.DB,d.name,'ses-*'))),D);
  sA = D(find(ns >= 2,1)).name; sB = D(find(ns == 1,1)).name;
  S  = dir(fullfile(E.DB,sA,'ses-*')); sesA = {S([S.isdir]).name};
  [~,so] = sort(regexprep(sesA,'^.*-','')); sesA = sesA(so);
end
% =========================================================================
function b = bidsSub(s,DB)
%bidsSub. BIDS subject of a database subject sub-<serial>-<PatientID> 
%  (subIDform 2), i.e. with the site abbreviation of the site registry of 
%  the database (without site dictionary, see siteRegistry in 
%  cat_io_dcm2bids_protocols) or the serial number.
  tok  = regexp(s,'^sub-(.+?)-(.*)$','tokens','once'); site = tok{1};
  Preg = fullfile(DB,'catDCM2BIDSsites.json');
  if exist(Preg,'file')
    S = cat_io_json(Preg); i = find(strcmp(cellfun(@(x) char(string(x)),{S.DeviceSerialNumber},'UniformOutput',false),tok{1}),1);
    if ~isempty(i), site = S(i).InstitutionAbbreviation; end
  end
  b = ['sub-' cat_io_dcm2bids_bids('bidsLabel',site,1) regexprep(tok{2},'[^A-Za-z0-9]','')];
end
function b = bidsSes(s), b = ['ses-' regexprep(s,'^.*-','')]; end    % BIDS session of a database session (date)
% =========================================================================
function L = listing(B)
%listing. Relative files of the BIDS directory B with md5 (without mat files).
  if ismac, md5 = 'md5 -r'; else, md5 = 'md5sum'; end
  [~,txt] = system(sprintf('cd "%s" 2>/dev/null && find . -type f ! -name "*.mat" -exec %s {} \\; | sort -k2',B,md5));
  L = strsplit(strtrim(txt),newline)'; L = L(~cellfun('isempty',L));
end
function L = scans(L), L = L(~cellfun('isempty',regexp(L,'/(anat|func|dwi|fmap)/|catDCM2BIDSaffine','once'))); end
function L = names(L), L = regexprep(L,'^\S+\s+',''); end
function L = sel(L,pat), L = L(~cellfun('isempty',regexp(L,pat,'once'))); end
function S = sidecars(DB), S = cat_vol_findfiles(DB,'*=*=*=*.json'); S = S(~startsWith(spm_file(S,'filename'),'catDCM2BIDS')); end
function V = sidecarFields(P) % (fields of a sidecar of dcm2niix, without the warnings added by the import, see convertImages)
  V = jsondecode(fileread(P)); if isfield(V,'ConversionWarnings'), V = rmfield(V,'ConversionWarnings'); end
end
% =========================================================================
function resetBIDS(T)
%resetBIDS. Remove the BIDS outputs of a run (keep the database and its tables).
  D = dir(fullfile(T,'study','BIDS*')); D = D([D.isdir]); % (BIDS directories of all anonymization levels)
  for di = 1:numel(D), rmdir(fullfile(D(di).folder,D(di).name),'s'); end
end
function fresh(T), if exist(T,'dir'), rmdir(T,'s'); end, end
function cleanLog(E), if exist(E.flog,'file'), delete(E.flog); end, end
