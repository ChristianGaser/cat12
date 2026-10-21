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
%    PatientID   patientid, pid, dicomid
%    subject     participantid, participant, subject, subjectid, sub, proband
%                (PatientID or the subject of the database, sub-<serial>-<ID>)
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
%                (number of the session of the subject in the database by
%                date, or the StudyID)
%  Key columns given by number with another name are scan or birth dates
%  (if their values are dates), visits (numbers up to 99), or names. Dates
%  are read as yyyy-mm-dd, yyyymmdd, dd.mm.yyyy, or dd/mm/yyyy (with an
%  optional time). Empty key values are ignored.
%
%  Rows with a scan date or visit (session rows) have to match exactly one
%  session, rows without (subject rows) the sessions of one subject (i.e.
%  their data are given to all sessions of the subject). Sessions with
%  several rows of a table get no data of this table (conflict). The data
%  of a table (by its file name) are replaced by each import, i.e.
%  corrected tables can be imported again. The data of several tables are
%  combined (by the order of the table names for equal columns). The 
%  column names are the table headers without a description in brackets at
%  their end (e.g. "SEX (male = 1, female = 2)" gives SEX). Columns of the
%  sessions file (session_id, site, age, weight, SOR) get the suffix
%  "_user". The match of each row (with its key values, i.e. personal data)
%  is given in the private tables of the database
%  (catDCM2BIDSdb-tables/phenotypic_<table>.csv).
%
%  varargout = cat_io_dcm2bids_phenotypic(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_phenotypic('importTables',Ptables,keys,Pmdbdir,Ptbldir)
%    data = cat_io_dcm2bids_phenotypic('sessionData',Pdbses)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_bids, cat_io_csv.

  switch action
    case {'importTables','sessionData'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_phenotypic:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function importTables(Ptables,keys,Pmdbdir,Ptbldir)
%importTables. Import of the user tables Ptables (cell of CSV/TSV files)
%  with the key columns keys (column numbers, empty for the detection by
%  their names) into the sessions of the database Pmdbdir with the report
%  of each table in Ptbldir (see the help of this file).
  Ptables = cellstr(Ptables); Ptables = Ptables(~cellfun('isempty',strtrim(Ptables)));
  if isempty(Ptables), return; end
  fprintf('\nDCM2BIDS - phenotypic data\n');
  S = dbSessions(Pmdbdir);
  for ti = 1:numel(Ptables)
    try
      T = readTable(Ptables{ti},keys);
    catch e
      cat_io_cprintf('err','  %s: %s\n',spm_file(Ptables{ti},'filename'),e.message);
      continue
    end
    [assign,status,cand] = matchRows(T,S);
    storeTable(T,S,assign);
    R = [{'row','status','sessions'} T.hdr(T.kc)];
    for ri = 1:size(T.rows,1)
      R(end+1,:) = [{ri + 1, status{ri}, strjoin(cellfun(@(d) dbName(d,Pmdbdir),{S(cand{ri}).dir}, ...
        'UniformOutput',false),' ')} T.rows(ri,T.kc)]; %#ok<AGROW>
    end
    cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['phenotypic_' T.name '.csv']),R,'csv');
    n = @(s) nnz(strncmp(status,s,numel(s)));
    fprintf('  %s: %d rows, %d assigned (%d sessions), %d without session, %d ambiguous, %d conflicts (keys: %s)\n', ...
      T.name,size(T.rows,1),n('assigned'),nnz(assign > 0),n('no session'),n('ambiguous'),n('conflict'), ...
      strjoin(cellfun(@(h,t) sprintf('%s = %s',h,t),T.hdr(T.kc),T.types,'UniformOutput',false),', '));
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
%  directory, subject (sub-<serial>-<PatientID>), PatientID, words of the
%  PatientName, birth date, scan date and time, StudyID, and the visit
%  (number of the session of the subject by date), by the first scan
%  sidecar of each session.
  S = struct('dir',{},'sub',{},'pid',{},'name',{},'birth',{},'date',{},'time',{},'study',{},'visit',{});
  D = dir(fullfile(Pmdbdir,'sub-*','ses-*')); D = D([D.isdir]);
  for di = 1:numel(D)
    V = []; Psnr = dir(fullfile(D(di).folder,D(di).name,'snr-*'));
    for si = 1:numel(Psnr)
      P = cat_io_dcm2bids_db('getDBscanJSONs',fullfile(Psnr(si).folder,Psnr(si).name));
      if ~isempty(P), V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(P{1})); break; end
    end
    if isempty(V), continue; end
    [~,sub] = fileparts(D(di).folder);
    S(end+1) = struct('dir',fullfile(D(di).folder,D(di).name),'sub',sub,'pid',normID(V.PatientID), ...
      'name',{nameWords(V.PatientName)},'birth',parseDate(V.PatientBirthDate),'date',dateshift(parseDate(char(V.ScanDate)),'start','day'), ...
      'time',parseDate(fieldText(V,'AcquisitionDateTime')),'study',normID(V.StudyID),'visit',0); %#ok<AGROW>
  end
  for sub = unique({S.sub})
    i = find(strcmp({S.sub},sub{1}));
    t = datenum([S(i).date]); t(isnan(t)) = inf; %#ok<DATNM>
    [~,o] = sortrows([t(:) (1:numel(i))']);
    for k = 1:numel(o), S(i(o(k))).visit = k; end
  end
end
% =========================================================================
function T = readTable(P,keys)
%readTable. User table P (CSV or TSV) with its header hdr, rows (text),
%  key columns kc with their types, and data columns dc with their names
%  cols (valid names without the columns of the sessions file).
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
  dc   = setdiff(1:numel(hdr),kc);
  cols = cellfun(@(h) matlab.lang.makeValidName(regexprep(transliterate(columnTitle(h)),'\s+','_')),hdr(dc),'UniformOutput',false);
  res  = ismember(lower(cols),{'session_id','site','age','weight','sor'});
  cols(res) = strcat(cols(res),'_user');
  cols = matlab.lang.makeUniqueStrings(cols);
  T = struct('file',P,'name',matlab.lang.makeValidName(spm_file(P,'basename')),'hdr',{hdr},'rows',{rows}, ...
    'kc',kc,'types',{types(kc)},'dc',dc,'cols',{cols});
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
function [assign,status,cand] = matchRows(T,S)
%matchRows. Sessions of each row of the table T (cand, indices of the
%  sessions S), the row assigned to each session (assign, 0 for none, see
%  the help of this file), and the status of each row.
  nr = size(T.rows,1); cand = cell(nr,1); status = cell(nr,1); sesrow = false(nr,1);
  for ri = 1:nr
    ok = true(1,numel(S)); w = {}; d0 = NaT; nk = 0;
    for k = 1:numel(T.kc)
      v = T.rows{ri,T.kc(k)};
      if isempty(v), continue; end
      nk = nk + 1;
      switch T.types{k}
        case 'pid'
          ok = ok & strcmp({S.pid},normID(v));
        case 'sub'
          x  = normID(regexprep(v,'^sub-',''));
          ok = ok & (strcmp({S.pid},x) | strcmp(cellfun(@(s) normID(regexprep(s,'^sub-','')),{S.sub},'UniformOutput',false),x));
        case 'name'
          w  = [w nameWords(v)]; %#ok<AGROW>
        case 'birth'
          d  = parseDate(v);
          ok = ok & arrayfun(@(s) sameMonth(s.birth,d),S);
        case 'date'
          d0 = parseDate(v); sesrow(ri) = true;
          ok = ok & arrayfun(@(s) isequal(s.date,dateshift(d0,'start','day')),S);
        case 'anydate'
          d  = parseDate(v); sesrow(ri) = true;
          ok = ok & arrayfun(@(s) isequal(s.date,dateshift(d,'start','day')) || sameMonth(s.birth,d),S);
        case 'visit'
          x  = regexprep(v,'^ses-',''); sesrow(ri) = true;
          ok = ok & ([S.visit] == str2double(x) | strcmp({S.study},normID(x)));
      end
    end
    if ~isempty(w)
      w  = unique(w);
      ok = ok & arrayfun(@(s) nameMatch(w,s.name),S);
    end
    c = find(ok);
    if numel(c) > 1 && ~isnat(d0) && any([hour(d0) minute(d0)]) % closest time of sessions of the same day
      [~,ci] = min(abs([S(c).time] - d0)); c = c(ci);
    end
    if nk == 0, c = []; status{ri} = 'no key values'; end
    cand{ri} = c;
  end

  % assignment of the rows to the sessions
  assign = zeros(1,numel(S));
  for ri = 1:nr
    c = cand{ri};
    if ~isempty(status{ri}), continue; end
    if isempty(c)
      status{ri} = 'no session';
    elseif sesrow(ri) && numel(c) > 1
      status{ri} = sprintf('ambiguous (%d sessions)',numel(c));
    elseif ~sesrow(ri) && numel(unique({S(c).sub})) > 1
      status{ri} = sprintf('ambiguous (%d subjects)',numel(unique({S(c).sub})));
    else
      assign(c(assign(c) == 0)) = ri; assign(c(assign(c) > 0 & assign(c) ~= ri)) = -1;
      status{ri} = sprintf('assigned (%d sessions)',numel(c));
    end
  end
  for ri = find(strncmp(status,'assigned',8))'
    if any(assign(cand{ri}) == -1)
      status{ri} = sprintf('conflict (%d of %d sessions with other rows)',nnz(assign(cand{ri}) == -1),numel(cand{ri}));
    end
  end
  assign(assign < 0) = 0;
end
% =========================================================================
function storeTable(T,S,assign)
%storeTable. Data of the table T in the database: the data of this table
%  (by its name) are removed from all sessions S and the data of the rows
%  assigned to the sessions (assign) are added (catDCM2BIDSphenotypic.json
%  of the session directories).
  for si = 1:numel(S)
    P = fullfile(S(si).dir,'catDCM2BIDSphenotypic.json');
    if exist(P,'file'), J = cat_io_json(P); else, J = struct('tables',struct()); end
    if ~isfield(J,'tables') || ~isstruct(J.tables), J.tables = struct(); end
    had = isfield(J.tables,T.name);
    if had, J.tables = rmfield(J.tables,T.name); end
    if assign(si) > 0
      r = assign(si); data = struct();
      for k = 1:numel(T.dc)
        if ~isempty(T.rows{r,T.dc(k)}), data.(T.cols{k}) = T.rows{r,T.dc(k)}; end
      end
      J.tables.(T.name) = struct('table',spm_file(T.file,'filename'),'row',r + 1, ...
        'imported',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')),'data',data);
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
