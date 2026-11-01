function varargout = cat_io_dcm2bids_phenotypic(action,varargin)
%cat_io_dcm2bids_phenotypic. Phenotypic data of the sessions of cat_io_dcm2bids.
%  User tables (CSV or TSV, e.g. clinical or behavioral data of the scan
%  sessions, see cat_io_csv) are matched to the sessions of the internal
%  database by their key columns. The other columns (data) are stored in
%  the database (catDCM2BIDSphenotypic.json of the session directory) and
%  exported to the BIDS sessions files (sub-<label>_sessions.tsv, see
%  writeSubjectTSV in cat_io_dcm2bids_bids). The key columns are not
%  exported, as they usually contain personal data (e.g. names).
%
%  Key columns are given by their numbers (job.dataset.phenotypickeys) or
%  detected by their names (ignoring case and special characters):
%    participant participantid, participant, subject, subjectid, sub, proband
%                (ID of the participant as stored in the DICOM PatientID by
%                the scanner, or the subject of the database, sub-<serial>-<ID>,
%                but not the BIDS participant label of the export)
%    PatientID   patientid, pid, dicomid (only the DICOM PatientID)
%    name        patientname, name, fullname, lastname, surname, familyname,
%                firstname, givenname, forename, nachname, vorname
%                (the words of all name columns, e.g. "Max" and "Mustermann",
%                and of the DICOM name, e.g. "Mustermann^Max", in any order
%                and with at least two words, where further words of one
%                of them are ignored, e.g. titles or second names)
%    birth date  birthdate, dateofbirth, dob, patientbirthdate, birthday,
%                geburtsdatum (year and month)
%    scan date   scandate, date, datetime, acquisitiondate(time), studydate,
%                sessiondate, mrdate, datum, messdatum (day, and the closest
%                time for sessions of the same day)
%    visit       visit, timepoint, tp, wave, session, sessionid, ses, besuch
%                (number of the session of the participant by date, also of
%                several scanners, or the StudyID)
%  Key columns given by number with another name are scan or birth dates
%  (if their values are dates), visits (numbers up to 99), or names. Dates
%  are read as yyyy-mm-dd, yyyymmdd, dd.mm.yyyy, or dd/mm/yyyy (with an
%  optional time). Empty key values are ignored. See
%  DCM2BIDS/study/cat_study_phenotype.csv as template (with the column
%  description cat_study_phenotype.json).
%
%  Key values can be patterns: "*" (all values), wildcards "*" and "?"
%  (e.g. 002_S_* for the participants of the ADNI site 002, or Muster*), a
%  year or month of dates (2023, 2023-05), and ranges a..b of dates or
%  visits (e.g. 2023-09-01..2026-12-31, 2023.., or 2..3).
%
%  Rows with an exact scan date or visit (session rows) have to match one
%  session of each subject of the database of one participant (the
%  PatientID, i.e. also sessions of several scanners), rows without
%  (participant rows) the sessions of one participant. Rows with patterns
%  can match any sessions. The data of a session are combined from all its
%  rows by their specificity (level): session rows replace the values of
%  participant rows, these the values of pattern rows, and pattern rows the
%  values of former pattern rows (e.g. study defaults in a first row
%  "*,*,*,*,DZPG,JE1,NC" and site or participant values below). Empty
%  cells replace nothing. Different values of participant or session rows
%  of a session are a conflict, i.e. the session gets no value of this
%  column. The data of a table (by its file name) are replaced by each
%  import, i.e. corrected tables can be imported again. The data of several
%  tables are combined (by the order of the table names for equal columns).
%
%  The column names are the table headers without a description in
%  brackets at their end (e.g. "SEX (male = 1, female = 2)" gives SEX).
%  The phenotype dictionary (JSON file, e.g. DCM2BIDS/study/
%  cat_study_phenotype_dictionary.json, see readDictionary) renames columns
%  with a variant of a field name (e.g. gender -> sex) and recodes their
%  values with a variant of a level (e.g. weiblich -> F). Numbers are only
%  recoded by the coding of each table, i.e. the description of the header
%  (e.g. 1 -> male -> M) or the Levels of the JSON sidecar of the table
%  (table.json of table.csv). Other values are kept (with a warning).
%  Columns of the sessions file (session_id, site, age, weight, SOR) get
%  the suffix "_user". The match of each row (with its key values, i.e.
%  personal data) is given in the private tables of the database
%  (catDCM2BIDSdb-tables/phenotypic_<table>.csv).
%
%  varargout = cat_io_dcm2bids_phenotypic(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_phenotypic('importTables',Ptables,keys,Pmdbdir,Ptbldir,Pdict)
%    ok = cat_io_dcm2bids_phenotypic('checkTables',Ptables,keys,Pdict)
%    data = cat_io_dcm2bids_phenotypic('sessionData',Pdbses)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_bids, cat_io_csv.

  switch action
    case {'importTables','checkTables','sessionData'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_phenotypic:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function importTables(Ptables,keys,Pmdbdir,Ptbldir,Pdict)
%importTables. Import of the user tables Ptables (cell of CSV/TSV files)
%  with the key columns keys (column numbers, empty for the detection by
%  their names) into the sessions of the database Pmdbdir with the report
%  of each table in Ptbldir and the phenotype dictionary Pdict (JSON file,
%  see the help of this file).
  if ~exist('Pdict','var'), Pdict = ''; end
  Ptables = cellstr(Ptables); Ptables = Ptables(~cellfun('isempty',strtrim(Ptables)));
  if isempty(Ptables), return; end
  fprintf('\nDCM2BIDS - phenotypic data\n');
  Dict = readDictionary(Pdict,1);
  S = dbSessions(Pmdbdir);
  for ti = 1:numel(Ptables)
    try
      T = readTable(Ptables{ti},keys,Dict);
    catch e
      cat_io_cprintf('err','  %s: %s\n',spm_file(Ptables{ti},'filename'),e.message);
      continue
    end
    [cand,status,level] = matchRows(T,S);
    [D,src,status] = mergeRows(T,S,cand,status,level);
    storeTable(T,S,D,src);
    lname = {'all','pattern','participant','session'};
    R = [{'row','status','level','sessions'} T.hdr(T.kc)];
    for ri = 1:size(T.rows,1)
      R(end+1,:) = [{ri + 1, status{ri}, lname{level(ri) + 1}, strjoin(cellfun(@(d) dbName(d,Pmdbdir), ...
        {S(cand{ri}).dir},'UniformOutput',false),' ')} T.rows(ri,T.kc)]; %#ok<AGROW>
    end
    cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['phenotypic_' T.name '.csv']),R,'csv');
    n = @(s) nnz(strncmp(status,s,numel(s)));
    fprintf('  %s: %d rows, %d assigned (%d sessions), %d without session, %d ambiguous, %d with conflicts (keys: %s)\n', ...
      T.name,size(T.rows,1),n('assigned'),nnz(~cellfun('isempty',src)),n('no session'),n('ambiguous'), ...
      nnz(contains(status,'conflict')), ...
      strjoin(cellfun(@(h,t) sprintf('%s = %s',h,t),T.hdr(T.kc),T.types,'UniformOutput',false),', '));
    if ~isempty(T.renamed)
      fprintf('    dictionary columns: %s (%d values recoded)\n',strjoin(T.renamed,', '),T.recoded);
    end
    if ~isempty(T.unknown)
      cat_io_cprintf('warn','    values not in the dictionary (kept): %s\n',strjoin(T.unknown,', '));
    end
  end
  tableConflicts(S);
end
% =========================================================================
function tableConflicts(S)
%tableConflicts. Warning for the columns of several tables with different
%  values in the same sessions S (e.g. the studyName of two studies with
%  the same participant IDs, whose participant rows match the sessions of
%  both studies), where the value of the last table (by name) is used (see
%  sessionData). Session rows (e.g. with the scan date) avoid this.
  C = containers.Map;
  for si = 1:numel(S)
    P = fullfile(S(si).dir,'catDCM2BIDSphenotypic.json'); if ~exist(P,'file'), continue; end
    J = cat_io_json(P); if ~isfield(J,'tables') || ~isstruct(J.tables), continue; end
    t = sort(fieldnames(J.tables));
    for a = 1:numel(t)
      for b = a + 1:numel(t)
        Da = J.tables.(t{a}); Db = J.tables.(t{b});
        if ~isfield(Da,'data') || ~isfield(Db,'data') || ~isstruct(Da.data) || ~isstruct(Db.data), continue; end
        for f = reshape(intersect(fieldnames(Da.data),fieldnames(Db.data)),1,[])
          if ~strcmp(valueText(Da.data.(f{1})),valueText(Db.data.(f{1})))
            k = sprintf('%s and %s: %s',t{a},t{b},f{1});
            if isKey(C,k), C(k) = C(k) + 1; else, C(k) = 1; end
          end
        end
      end
    end
  end
  for k = keys(C)
    cat_io_cprintf('warn','  different values of the tables %s in %d session(s) (the value of the last table is used)\n',k{1},C(k{1}));
  end
end
% =========================================================================
function s = valueText(v)
%valueText. Text of a value of the phenotypic data for the comparison of
%  tables (e.g. 35.1 and '35.1' are equal).
  if ischar(v) || (isstring(v) && isscalar(v)), s = char(v);
  elseif (isnumeric(v) || islogical(v)) && isscalar(v), s = sprintf('%g',v);
  else, s = jsonencode(v);
  end
end
% =========================================================================
function ok = checkTables(Ptables,keys,Pdict)
%checkTables. Check of the user tables Ptables (readable, with key
%  columns) and of the phenotype dictionary Pdict at the start of the
%  batch, i.e. before the (time-consuming) DICOM import. ok is false if a
%  table cannot be imported.
  if ~exist('Pdict','var'), Pdict = ''; end
  ok = true;
  Ptables = cellstr(Ptables); Ptables = Ptables(~cellfun('isempty',strtrim(Ptables)));
  if isempty(Ptables), return; end
  Dict = readDictionary(Pdict,0);
  for ti = 1:numel(Ptables)
    if ~exist(Ptables{ti},'file')
      cat_io_cprintf('err','  Phenotypic table "%s" not found.\n',Ptables{ti}); ok = false; continue
    end
    try
      readTable(Ptables{ti},keys,Dict);
    catch e
      cat_io_cprintf('err','  Phenotypic table "%s": %s\n',spm_file(Ptables{ti},'filename'),e.message); ok = false;
    end
  end
end
% =========================================================================
function data = sessionData(Pdbses)
%sessionData. Phenotypic data of the session directory Pdbses of the
%  database (structure with the data columns of all tables in the order of
%  the table names, empty structure without data).
  data = struct(); P = fullfile(Pdbses,'catDCM2BIDSphenotypic.json');
  if ~exist(P,'file'), return; end
  J = cat_io_json(P);
  if ~isfield(J,'tables') || ~isstruct(J.tables), return; end
  for t = sort(fieldnames(J.tables))'
    if ~isfield(J.tables.(t{1}),'data') || ~isstruct(J.tables.(t{1}).data), continue; end
    D = J.tables.(t{1}).data;
    for f = fieldnames(D)', data.(f{1}) = D.(f{1}); end
  end
end
% =========================================================================
function S = dbSessions(Pmdbdir)
%dbSessions. Sessions of the database with the fields for the matching:
%  directory, subject (sub-<serial>-<PatientID>), PatientID, participant
%  (normalized PatientID with the birth month, if available, as different
%  sites may use the same IDs, e.g. 001; the subject without PatientID),
%  words of the
%  PatientName, birth date, scan date and time, StudyID, and the visit
%  (number of the session of the participant by date and time, i.e. also
%  of several subjects of the database of one PatientID, e.g. of different
%  scanners), by the first scan sidecar of each session.
  S = struct('dir',{},'sub',{},'pid',{},'part',{},'name',{},'birth',{},'date',{},'time',{},'study',{},'visit',{});
  D = dir(fullfile(Pmdbdir,'sub-*','ses-*')); D = D([D.isdir]);
  for di = 1:numel(D)
    V = []; Psnr = dir(fullfile(D(di).folder,D(di).name,'snr-*'));
    for si = 1:numel(Psnr)
      P = cat_io_dcm2bids_db('getDBscanJSONs',fullfile(Psnr(si).folder,Psnr(si).name));
      if ~isempty(P), V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(P{1})); break; end
    end
    if isempty(V), continue; end
    [~,sub] = fileparts(D(di).folder);
    pid = normID(V.PatientID); birth = parseDate(V.PatientBirthDate);
    if isempty(pid), part = sub; elseif isnat(birth), part = pid; else, part = [pid '-' char(birth,'yyyyMM')]; end
    S(end+1) = struct('dir',fullfile(D(di).folder,D(di).name),'sub',sub,'pid',pid,'part',part, ...
      'name',{nameWords(V.PatientName)},'birth',birth,'date',dateshift(parseDate(char(V.ScanDate)),'start','day'), ...
      'time',parseDate(fieldText(V,'AcquisitionDateTime')),'study',normID(V.StudyID),'visit',0); %#ok<AGROW>
  end
  for p = unique({S.part})
    i = find(strcmp({S.part},p{1}));
    t = datenum([S(i).date]); t(isnan(t)) = inf; %#ok<DATNM>
    h = datenum([S(i).time]); h(isnan(h)) = inf; %#ok<DATNM>
    [~,o] = sortrows([t(:) h(:) (1:numel(i))']);
    for k = 1:numel(o), S(i(o(k))).visit = k; end
  end
end
% =========================================================================
function T = readTable(P,keys,Dict)
%readTable. User table P (CSV or TSV) with its header hdr, rows (text),
%  key columns kc with their types, data columns dc with their names cols
%  (valid names without the columns of the sessions file) and their data
%  (values recoded by the phenotype dictionary Dict, see readDictionary),
%  the renamed columns of the dictionary, the number of recoded values,
%  and the values that are not in the dictionary (unknown).
  if ~exist('Dict','var'), Dict = readDictionary(''); end
  C = cat_io_csv(P,'','',struct('convert2double',0));
  if isempty(C), error('cat_io_dcm2bids_phenotypic:empty','empty table'); end
  C(cellfun(@(c) ~ischar(c),C)) = {''};
  C = strtrim(C);
  C{1,1} = regexprep(C{1,1},['^(' char(65279) '|' char([239 187 191]) ')'],''); % byte order mark
  C = C([true; any(~cellfun('isempty',C(2:end,:)),2)],:);   % without empty rows
  if size(C,1) < 2, error('cat_io_dcm2bids_phenotypic:empty','no rows'); end
  hdr = C(1,:); rows = C(2:end,:);
  types = cellfun(@keyType,hdr,'UniformOutput',false);
  if isempty(keys)
    kc = find(~cellfun('isempty',types));
  else
    kc = unique(keys(keys >= 1 & keys <= numel(hdr)));
    for c = kc(cellfun('isempty',types(kc))), types{c} = contentType(rows(:,c)); end
  end
  if isempty(kc), error('cat_io_dcm2bids_phenotypic:noKeys','no key columns (see cat_io_dcm2bids_phenotypic)'); end
  dc     = setdiff(1:numel(hdr),kc);
  titles = cellfun(@columnTitle,hdr(dc),'UniformOutput',false);
  cols   = cellfun(@(h) matlab.lang.makeValidName(regexprep(transliterate(h),'\s+','_')),titles,'UniformOutput',false);

  % dictionary: names of the columns (variants) and their values (levels)
  data = rows(:,dc); renamed = {}; recoded = 0; unknown = {};
  side = tableSidecar(P);
  for k = 1:numel(dc)
    di = find(arrayfun(@(d) any(strcmp(d.variants,normID(transliterate(titles{k})))),Dict),1);
    if isempty(di), continue; end
    if ~strcmp(cols{k},Dict(di).name), renamed{end+1} = sprintf('%s -> %s',titles{k},Dict(di).name); end %#ok<AGROW>
    cols{k} = Dict(di).name;
    if isempty(Dict(di).codes), continue; end
    coding = tableCoding(hdr{dc(k)},titles{k},side);
    unk = {};
    for ri = 1:size(data,1)
      [v,known] = dictValue(Dict(di),data{ri,k},coding);
      if ~known, unk{end+1} = data{ri,k}; end %#ok<AGROW>
      recoded = recoded + ~strcmp(v,data{ri,k}); data{ri,k} = v;
    end
    [u,~,ui] = unique(unk);
    for ui2 = 1:numel(u), unknown{end+1} = sprintf('%s = %s (%d)',cols{k},u{ui2},nnz(ui == ui2)); end %#ok<AGROW>
  end

  res  = ismember(lower(cols),{'session_id','site','age','weight','sor'});
  cols(res) = strcat(cols(res),'_user');
  cols = matlab.lang.makeUniqueStrings(cols);
  T = struct('file',P,'name',matlab.lang.makeValidName(spm_file(P,'basename')),'hdr',{hdr},'rows',{rows}, ...
    'kc',kc,'types',{types(kc)},'dc',dc,'cols',{cols},'data',{data},'renamed',{renamed},'recoded',recoded, ...
    'unknown',{unknown});
end
% =========================================================================
function Dict = readDictionary(P,verb)
%readDictionary. Phenotype dictionary of the JSON file P (e.g.
%  DCM2BIDS/study/cat_study_phenotype_dictionary.json) as structure array
%  with the column name, its normalized variants (incl. the name), and the
%  codes of its levels with their normalized variants (incl. the code and
%  its description). Variants with only digits are ignored, as numeric
%  codes are defined by each table (see tableCoding). Warnings for missing
%  or invalid files if verb (default).
  if ~exist('verb','var'), verb = 1; end
  Dict = struct('name',{},'variants',{},'codes',{},'levels',{});
  P = strtrim(char(P)); if isempty(P), return; end
  if ~exist(P,'file')
    cat_io_cprintf('warn','  Phenotype dictionary "%s" not found - no recoding.\n',P); return
  end
  try
    J = jsondecode(cat_io_dcm2bids_helper('readText',P));
  catch e
    cat_io_cprintf('warn','  Phenotype dictionary "%s" is no valid JSON file (%s) - no recoding.\n',P,e.message); return
  end
  nv = @(c) unique(cellfun(@(x) normID(transliterate(x)),cellstr(c(:)'),'UniformOutput',false));
  nonum = @(c) c(cellfun('isempty',regexp(c,'^\d*$','once')));
  for f = fieldnames(J)'
    if strncmp(f{1},'x_',2) || ~isstruct(J.(f{1})), continue; end  % comments ("_..." fields)
    E = J.(f{1}); vars = f(1);
    if isfield(E,'Variants') && ~isempty(E.Variants), vars = [vars cellstr(E.Variants(:)')]; end %#ok<AGROW>
    codes = {}; levels = {};
    if isfield(E,'Levels') && ~isempty(E.Levels)
      L = E.Levels; if isstruct(L), L = num2cell(L); end
      for li = 1:numel(L)
        if ~isfield(L{li},'Code') || isempty(L{li}.Code), continue; end
        lv = {char(L{li}.Code)};
        if isfield(L{li},'Description') && ~isempty(L{li}.Description), lv{end+1} = char(L{li}.Description); end %#ok<AGROW>
        if isfield(L{li},'Variants') && ~isempty(L{li}.Variants), lv = [lv cellstr(L{li}.Variants(:)')]; end %#ok<AGROW>
        codes{end+1} = char(L{li}.Code); levels{end+1} = nonum(nv(lv)); %#ok<AGROW>
      end
    end
    Dict(end+1) = struct('name',f{1},'variants',{nonum(nv(vars))},'codes',{codes},'levels',{levels}); %#ok<AGROW>
  end
  if verb, fprintf('  dictionary: %s (%s)\n',spm_file(P,'filename'),strjoin({Dict.name},', ')); end
end
% =========================================================================
function [v,known] = dictValue(E,v,coding)
%dictValue. Value v of a column with the dictionary entry E: the code of
%  the level with the variant v (e.g. "F" for "female" or "weiblich"),
%  otherwise the code of the label of v in the coding of the table (map of
%  normalized codes, e.g. "2" -> "female", see tableCoding), otherwise v
%  (known = false, empty values are known).
  known = true; if isempty(v), return; end
  li = find(cellfun(@(l) any(strcmp(l,normID(transliterate(v)))),E.levels),1);
  if isempty(li) && isKey(coding,normID(v))
    v  = coding(normID(v));
    li = find(cellfun(@(l) any(strcmp(l,normID(transliterate(v)))),E.levels),1);
  end
  if isempty(li), known = false; else, v = E.codes{li}; end
end
% =========================================================================
function side = tableSidecar(P)
%tableSidecar. JSON sidecar of the table P (same name, e.g. table.json of
%  table.csv) with the descriptions of the columns (BIDS-like, e.g. with
%  "Levels": {"1": "male", "2": "female"}), empty structure without.
  side = struct(); Pj = spm_file(P,'ext','.json');
  if ~exist(Pj,'file'), return; end
  try
    side = jsondecode(cat_io_dcm2bids_helper('readText',Pj));
  catch e
    cat_io_cprintf('warn','  Sidecar "%s" is no valid JSON file (%s) - ignored.\n',Pj,e.message);
  end
  if ~isstruct(side) || ~isscalar(side), side = struct(); end
end
% =========================================================================
function coding = tableCoding(h,title,side)
%tableCoding. Coding of the values of a column of a table (map of the
%  normalized codes to their labels) by the description of the header h,
%  e.g. "SEX (male = 1, female = 2)" or "hand (1: right, 2: left)", or by
%  the levels of the column title in the sidecar of the table (see
%  tableSidecar). Codes of numbers are only defined in this way, i.e. for
%  each table, never by the dictionary.
  coding = containers.Map('KeyType','char','ValueType','char');
  d = regexp(char(h),'[\(\[]([^\(\)\[\]]*)[\)\]]\s*$','tokens','once');
  if ~isempty(d)
    for p = strsplit(d{1},{',',';'})
      ab = regexp(p{1},'^\s*(.+?)\s*[=:]\s*(.+?)\s*$','tokens','once');
      if isempty(ab), continue; end
      if ~isempty(regexp(ab{1},'^\d+$','once')) || (isempty(regexp(ab{2},'^\d+$','once')) && numel(ab{1}) <= numel(ab{2}))
        coding(normID(ab{1})) = ab{2};
      else
        coding(normID(ab{2})) = ab{1};
      end
    end
  end
  for f = fieldnames(side)'
    if ~strcmp(normID(f{1}),normID(transliterate(title))) || ~isstruct(side.(f{1})) || ...
       ~isfield(side.(f{1}),'Levels') || ~isstruct(side.(f{1}).Levels), continue; end
    L = side.(f{1}).Levels;
    for c = fieldnames(L)'
      if ischar(L.(c{1})), coding(normID(regexprep(c{1},'^x(?=\d)',''))) = L.(c{1}); end  % jsondecode: "1" -> x1
    end
  end
end
% =========================================================================
function h = columnTitle(h)
%columnTitle. Column header without a description in brackets at its end
%  (e.g. "SEX (male = 1, female = 2)" gives "SEX"), the header if nothing
%  else remains.
  t = strtrim(regexprep(char(h),'\s*[\(\[][^\(\)\[\]]*[\)\]]\s*$',''));
  if ~isempty(t), h = t; end
end
% =========================================================================
function t = keyType(h)
%keyType. Type of a key column by its name h ('' for other columns).
  K = {
    'pid',   {'patientid','pid','dicomid','dicompatientid'}
    'sub',   {'participantid','participant','subject','subjectid','subid','sub','proband','probandid'}
    'name',  {'patientname','name','fullname','lastname','surname','familyname','nachname', ...
              'firstname','givenname','forename','vorname'}
    'birth', {'birthdate','dateofbirth','dob','patientbirthdate','birthday','geburtsdatum','gebdatum'}
    'date',  {'scandate','date','datetime','scandatetime','acquisitiondate','acquisitiondatetime','studydate', ...
              'sessiondate','mrdate','mridate','scantime','datum','messdatum','untersuchungsdatum'}
    'visit', {'visit','visitid','timepoint','tp','wave','session','sessionid','ses','besuch'}};
  n = regexprep(h,'[^A-Za-z0-9]',''); t = '';
  for ki = 1:size(K,1)
    if any(strcmpi(n,K{ki,2})), t = K{ki,1}; return; end
  end
end
% =========================================================================
function t = contentType(v)
%contentType. Type of a key column without a key name by its values v:
%  dates (anydate, scan or birth date), visits (numbers up to 99), or names.
  v = v(~cellfun('isempty',v));
  if ~isempty(v) && mean(cellfun(@(x) ~isnat(parseDate(x)),v)) > 0.8
    t = 'anydate';
  elseif ~isempty(v) && all(~cellfun('isempty',regexp(v,'^\d{1,2}$','once')))
    t = 'visit';
  else
    t = 'name';
  end
end
% =========================================================================
function [cand,status,level] = matchRows(T,S)
%matchRows. Sessions of each row of the table T (cand, indices of the
%  sessions S), the status of each row, and its level (specificity: 0 -
%  only wildcards "*", 1 - patterns, 2 - participant, 3 - session; see the
%  help of this file). Rows with exact key values (without patterns) have
%  to match one participant (the PatientID, also of several subjects of
%  the database, e.g. of different scanners) and session rows (scan date
%  or visit) one session of each subject of the database.
  nr = size(T.rows,1); cand = cell(nr,1); status = repmat({''},nr,1); level = zeros(nr,1); part = {S.part};
  for ri = 1:nr
    ok = true(1,numel(S)); w = {}; d0 = NaT; nk = 0; nstar = 0; pat = false; exsub = false; exses = false;
    for k = 1:numel(T.kc)
      v = T.rows{ri,T.kc(k)}; t = T.types{k};
      if isempty(v), continue; end
      nk = nk + 1;
      if strcmp(v,'*'), nstar = nstar + 1; continue; end
      [m,isp,d] = keyMatch(t,v,S);
      ok  = ok & m; pat = pat || isp;
      if isp, continue; end
      if any(strcmp(t,{'date','anydate','visit'})), exses = true; else, exsub = true; end
      if strcmp(t,'name'), w = [w nameWords(v)]; end %#ok<AGROW>
      if strcmp(t,'date'), d0 = d; end
    end
    if ~isempty(w)
      w  = unique(w);
      ok = ok & arrayfun(@(s) nameMatch(w,s.name),S);
    end
    c = find(ok);
    if nk == 0, status{ri} = 'no key values'; continue; end
    if     nstar == nk,                level(ri) = 0;
    elseif exses && (exsub || ~pat),   level(ri) = 3;
    elseif exsub,                      level(ri) = 2;
    else,                              level(ri) = 1;
    end
    exact = ~pat && nstar < nk;
    if exact && exses && ~isnat(d0) && any([hour(d0) minute(d0)]) % closest time of sessions of the same day (of each subject)
      cs = [];
      for sub = unique({S(c).sub})
        ci = c(strcmp({S(c).sub},sub{1}));
        [~,mi] = min(abs([S(ci).time] - d0)); cs(end+1) = ci(mi); %#ok<AGROW>
      end
      c = sort(cs);
    end
    cand{ri} = c;
    if isempty(c)
      status{ri} = 'no session';
    elseif exact && numel(unique(part(c))) > 1
      status{ri} = sprintf('ambiguous (%d participants)',numel(unique(part(c))));
    elseif exact && exses && numel(c) > numel(unique({S(c).sub}))
      status{ri} = sprintf('ambiguous (%d sessions)',numel(c));
    else
      status{ri} = sprintf('assigned (%d sessions)',numel(c));
    end
  end
end
% =========================================================================
function [m,isp,d] = keyMatch(t,v,S)
%keyMatch. Sessions S that match the value v of a key column of the type
%  t (m), whether v is a pattern (isp), and the date d of a scan date.
%  Patterns are the wildcards "*" and "?" (e.g. 002_S_* or Muster*),
%  ranges a..b of dates or visits (open without a or b), and a year or
%  month of dates (e.g. 2023 or 2023-05). Name values without wildcards
%  are compared by matchRows (all name words of a row).
  d = NaT; m = true(1,numel(S));
  isdate = any(strcmp(t,{'date','birth','anydate'}));
  isp = ~isempty(regexp(v,'[*?]','once')) || (any(strcmp(t,{'date','birth','anydate','visit'})) && contains(v,'..')) || ...
        (isdate && ~isempty(regexp(v,'^\d{4}([-/]\d{1,2})?$','once')));
  switch t
    case {'pid','sub'}
      x = regexprep(v,'^sub-','');
      if isp, eq = @(a) ~cellfun('isempty',regexp(a,globRegex(x),'once')); else, eq = @(a) strcmp(a,normID(x)); end
      m = eq({S.pid});
      if strcmp(t,'sub'), m = m | eq(cellfun(@(s) normID(regexprep(s,'^sub-','')),{S.sub},'UniformOutput',false)); end
    case 'name'
      if isp
        for p = regexp(lower(transliterate(v)),'[a-z0-9*?]+','match')
          m = m & arrayfun(@(s) any(~cellfun('isempty',regexp(s.name,globRegex(p{1}),'once'))),S);
        end
      end
    case {'date','birth','anydate'}
      if isp
        [lo,hi] = dateRange(v);
        in = @(x) ~isnat(x) && x >= lo && x < hi;
        switch t
          case 'date',    m = arrayfun(@(s) in(s.date),S);
          case 'birth',   m = arrayfun(@(s) in(s.birth),S);
          case 'anydate', m = arrayfun(@(s) in(s.date) || in(s.birth),S);
        end
      else
        d = parseDate(v); d1 = dateshift(d,'start','day');
        switch t
          case 'date',    m = arrayfun(@(s) isequal(s.date,d1),S);
          case 'birth',   m = arrayfun(@(s) sameMonth(s.birth,d),S);
          case 'anydate', m = arrayfun(@(s) isequal(s.date,d1) || sameMonth(s.birth,d),S);
        end
      end
    case 'visit'
      x = regexprep(v,'^ses-','');
      r = regexp(x,'^\s*(\d*)\s*\.\.\s*(\d*)\s*$','tokens','once');
      if ~isempty(r)
        lo = str2double(r{1}); if isnan(lo), lo = -inf; end
        hi = str2double(r{2}); if isnan(hi), hi = inf; end
        m  = [S.visit] >= lo & [S.visit] <= hi;
      elseif isp
        m  = ~cellfun('isempty',regexp({S.study},globRegex(x),'once'));
      else
        m  = [S.visit] == str2double(x) | strcmp({S.study},normID(x));
      end
  end
end
% =========================================================================
function re = globRegex(v)
%globRegex. Regular expression of a pattern v with the wildcards "*" (any
%  characters) and "?" (one character) for normalized identifiers (see
%  normID), e.g. '^002s.*$' for 002_S_*.
  p = regexp(v,'[*?]','split'); w = regexp(v,'[*?]','match');
  re = ['^' normID(p{1})];
  for i = 1:numel(w), re = [re strrep(strrep(w{i},'*','.*'),'?','.') normID(p{i+1})]; end %#ok<AGROW>
  re = [re '$'];
end
% =========================================================================
function [lo,hi] = dateRange(v)
%dateRange. Interval [lo,hi) of a date pattern: a date, year (2023), or
%  month (2023-05), or a range a..b of them (open without a or b).
  p = strtrim(strsplit(v,'..')); if isscalar(p), p = [p p]; end
  lo = periodBound(p{1},0); hi = periodBound(p{end},1);
end
% =========================================================================
function b = periodBound(x,isend)
%periodBound. Start (or end if isend) of the period of a date, year, or
%  month x (no bound for empty x, NaT for invalid dates).
  if isempty(x)
    if isend, b = datetime(9999,12,31); else, b = datetime(1,1,1); end
  elseif ~isempty(regexp(x,'^\d{4}$','once'))
    b = datetime(str2double(x) + isend,1,1);
  elseif ~isempty(regexp(x,'^\d{4}[-/]\d{1,2}$','once'))
    ym = str2double(regexp(x,'\d+','match')); b = datetime(ym(1),ym(2),1) + calmonths(isend);
  else
    b = dateshift(parseDate(x),'start','day') + days(isend);
  end
end
% =========================================================================
function [D,src,status] = mergeRows(T,S,cand,status,level)
%mergeRows. Data of each session S (D, structures of the data columns)
%  and its rows (src) by the assigned rows of the table T (cand, status)
%  in the order of their level (specificity) and position: rows of the
%  levels 0 and 1 (wildcards and patterns) replace the values of former
%  rows (e.g. study defaults in the first row and site values below),
%  rows of the levels 2 and 3 (participant and session) the values of
%  lower levels, whereas different values of rows of the same level are
%  a conflict (the column gets no value, given in the status of the rows).
%  Empty values do not replace values.
  ok = strncmp(status,'assigned',8); nr = numel(status);
  D = repmat({struct()},1,numel(S)); src = cell(1,numel(S));
  conf = repmat({{}},nr,1); nconf = zeros(nr,1);
  for si = 1:numel(S)
    r = find(ok & cellfun(@(c) any(c == si),cand));
    if isempty(r), continue; end
    [~,o] = sortrows([level(r) r]); r = r(o); data = struct(); rc = [];
    for l = 0:3
      rl = r(level(r) == l); if isempty(rl), continue; end
      for k = 1:numel(T.cols)
        v = T.data(rl,k); e = cellfun('isempty',v);
        if all(e), continue; end
        if l < 2 || isscalar(unique(v(~e)))
          data.(T.cols{k}) = v{find(~e,1,'last')};
        else
          if isfield(data,T.cols{k}), data = rmfield(data,T.cols{k}); end
          for ri = rl(~e)', conf{ri} = union(conf{ri},T.cols(k)); end
          rc = union(rc,rl(~e));
        end
      end
    end
    nconf(rc) = nconf(rc) + 1;
    D{si} = data; src{si} = r;
  end
  for ri = find(nconf > 0)'
    status{ri} = sprintf('%s, conflict of %s in %d sessions',status{ri},strjoin(conf{ri},' '),nconf(ri));
  end
end
% =========================================================================
function storeTable(T,S,D,src)
%storeTable. Data of the table T in the database: the data of this table
%  (by its name) are removed from all sessions S and the data D of the
%  sessions with rows (src) are added (catDCM2BIDSphenotypic.json of the
%  session directories).
  for si = 1:numel(S)
    P = fullfile(S(si).dir,'catDCM2BIDSphenotypic.json');
    if exist(P,'file'), J = cat_io_json(P); else, J = struct('tables',struct()); end
    if ~isfield(J,'tables') || ~isstruct(J.tables), J.tables = struct(); end
    had = isfield(J.tables,T.name);
    if had, J.tables = rmfield(J.tables,T.name); end
    if ~isempty(src{si})
      J.tables.(T.name) = struct('table',spm_file(T.file,'filename'),'rows',src{si}(:)' + 1, ...
        'imported',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')),'data',D{si});
    elseif ~had
      continue
    end
    if isempty(fieldnames(J.tables)), delete(P); else, cat_io_json(P,J); end
  end
end
% =========================================================================
function ok = nameMatch(a,b)
%nameMatch. True if the name words a are in b or b in a, with at least
%  two words (one word only if both have only one word).
  na = numel(a); nb = numel(b);
  ok = na > 0 && nb > 0 && min(na,nb) >= min(2,max(na,nb)) && (all(ismember(a,b)) || all(ismember(b,a)));
end
% =========================================================================
function ok = sameMonth(a,b)
%sameMonth. True if the dates a and b have the same year and month.
  ok = ~isnat(a) && ~isnat(b) && year(a) == year(b) && month(a) == month(b);
end
% =========================================================================
function w = nameWords(v)
%nameWords. Words of a name in lower case (special characters
%  transliterated) with at least two letters or digits, e.g. {'max',
%  'mustermann'} of "Mustermann^Max" (or {'pilot01'} of "Pilot01").
  w = regexp(lower(transliterate(char(v))),'[a-z0-9]+','match');
  w = unique(w(cellfun(@numel,w) > 1));
end
% =========================================================================
function s = transliterate(s)
%transliterate. Umlauts and other special characters of names (as
%  assurePatientDCMfields of cat_io_dcm2bids_db).
  s = cat_io_strrep(s,{'ä','ö','ü','Ä','Ö','Ü','ß','é','è','ê','á','à','â','ó','ò','ô','í','ì','ú','ù','ç','ñ'}, ...
                      {'ae','oe','ue','Ae','Oe','Ue','ss','e','e','e','a','a','a','o','o','o','i','i','u','u','c','n'});
end
% =========================================================================
function x = normID(v)
%normID. Identifier in lower case without special characters.
  if isnumeric(v), v = sprintf('%d',v); end
  x = lower(regexprep(char(v),'[^A-Za-z0-9]',''));
end
% =========================================================================
function v = fieldText(S,fn)
%fieldText. Text of the field fn of the structure S ('' if not available).
  if isfield(S,fn) && ~isempty(S.(fn)), v = char(S.(fn)); else, v = ''; end
end
% =========================================================================
function d = parseDate(v)
%parseDate. Date (and time) of a value in the formats yyyy-mm-dd,
%  yyyymmdd, dd.mm.yyyy, or dd/mm/yyyy with an optional time (HH:mm or
%  HH:mm:ss, also of yyyy-mm-ddTHH:mm:ss.ssssss), NaT otherwise.
  d = NaT;
  if isdatetime(v), d = v; return; end
  if isnumeric(v), v = sprintf('%d',v); end
  v = regexprep(strtrim(char(v)),'(\d)T(\d)','$1 $2'); v = regexprep(v,'(:\d\d)\.\d+$','$1');
  if isempty(v) || isempty(regexp(v,'\d{4}','once')), return; end
  F = {'yyyy-M-d H:m:s','yyyy-M-d H:m','yyyy-M-d','yyyyMMdd HHmmss','yyyyMMddHHmmss','yyyyMMdd', ...
       'd.M.yyyy H:m:s','d.M.yyyy H:m','d.M.yyyy','d/M/yyyy H:m:s','d/M/yyyy H:m','d/M/yyyy','yyyy/M/d'};
  for fi = 1:numel(F)
    try
      d = datetime(v,'InputFormat',F{fi});
      if ~isnat(d), return; end
    catch
      d = NaT;
    end
  end
end
% =========================================================================
function n = dbName(P,Pmdbdir)
%dbName. Session directory P relative to the database (sub-*/ses-*).
  n = strrep(strrep(P,[Pmdbdir filesep],''),filesep,'/');
end
