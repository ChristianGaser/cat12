function varargout = cat_io_dcm2bids_anonymize(action,varargin)
%cat_io_dcm2bids_anonymize. Subject numbers, session labels, and date 
%  shifts of the anonymization level 2 (opts.anonymize) of cat_io_dcm2bids.
%  Each subject of the database gets a number (five digits, in the order of
%  the registration and the first session) and a date offset in days that
%  shifts its first session to 1900-01-01. Both are stored in the registry
%  of the database (catDCM2BIDSanonymization.json, i.e. the key of the 
%  anonymized data that stays at the site that creates it) and are never
%  changed. The BIDS export uses the site and number as subject label (e.g.
%  sub-07SI30T00001) and the time since the first session in the unit 
%  opts.sesunit (d-days, w-weeks, m-months, y-years) as session label 
%  (e.g. ses-w0000 and ses-w0052 one year later, see session). The DICOM
%  export uses the subject label as PatientID and PatientName, shifts the 
%  dates of a session to the date of its label (1900-01-01 plus the time
%  since the first session, i.e. with the resolution of the unit) and its
%  times so that the session starts at the center of its 3-hour interval 
%  of the day (01:30, 04:30, ..., 22:30, i.e. the time of the day is only 
%  kept roughly, but the timing of the scans and slices is kept, see 
%  shiftSession and dayInterval), the birth date by the date offset of the
%  subject (on the 15th of the month, i.e. the age with an accuracy of 
%  about a month), and replaces the UIDs by UIDs of the
%  same length (hash with the secret of the registry, besides the standard
%  UIDs 1.2.840.10008.*), also in private elements (e.g. the Siemens CSA 
%  header, with their dates). 
%  The site label of the level 2 contains no location, only the manually 
%  defined siteID of the scanner (empty by default) with the manufacturer 
%  and field strength (see setupSites in cat_io_dcm2bids_protocols), and 
%  the BIDS sidecars contain no scanner model and coil and only the type of
%  the software version (see sidecar). 
%
%  varargout = cat_io_dcm2bids_anonymize(action,varargin)
%
%  Actions (see the help of the local functions):
%    A = cat_io_dcm2bids_anonymize('registry',Pmdbdir)
%    [number,offset] = cat_io_dcm2bids_anonymize('subject',A,dbsub)
%    [label,offset] = cat_io_dcm2bids_anonymize('session',date,offset,unit)
%    t = cat_io_dcm2bids_anonymize('acqTime',stamp,offset,unit,level)
%    s = cat_io_dcm2bids_anonymize('dayInterval',s)
%    s = cat_io_dcm2bids_anonymize('shiftDate',s,offset)
%    s = cat_io_dcm2bids_anonymize('shiftSession',s,vr,X)
%    u = cat_io_dcm2bids_anonymize('hashUID',u,salt)
%    V = cat_io_dcm2bids_anonymize('sidecar',V,level)
%    t = cat_io_dcm2bids_anonymize('softwareType',s)
%    b = cat_io_dcm2bids_anonymize('hashUIDsInBytes',b,salt)
%    b = cat_io_dcm2bids_anonymize('shiftDatesInBytes',b,X,known)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_export.

  switch action
    case {'registry','subject','session','acqTime','dayInterval','shiftDate','shiftSession','hashUID','hashUIDsInBytes', ...
          'shiftDatesInBytes','sidecar','softwareType'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_anonymize:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function A = registry(Pmdbdir)
%registry. Numbers and date offsets of all subjects of the database 
%  Pmdbdir (subject directories sub-<serial>-<PatientID>), with new 
%  subjects added to the registry (numbered by their first session): 
%  A.salt (secret of the UID hashes) and A.map (containers.Map of the 
%  subject directories with the number and offset).
  P = fullfile(Pmdbdir,'catDCM2BIDSanonymization.json');
  S = struct('subject',{},'number',{},'offset',{},'registered',{});
  rs = RandStream('mt19937ar','Seed','shuffle'); % independent of the global random numbers
  if exist(P,'file')
    J = cat_io_json(P); salt = J.salt;
    if isfield(J,'subjects') && ~isempty(J.subjects)
      T = J.subjects(:)';
      if ~isfield(T,'number'), n = num2cell(1:numel(T)); [T.number] = n{:}; end % registry with codes
      S = struct('subject',{T.subject},'number',{T.number},'offset',{T.offset},'registered',{T.registered});
    end
  else
    salt = char(randi(rs,[0 15],1,32) + '0'); salt(salt > '9') = salt(salt > '9') + 39; % hexadecimal
  end
  D = dir(fullfile(Pmdbdir,'sub-*')); D = D([D.isdir]);
  new = setdiff({D.name},{S.subject});
  first = repmat({'99999999'},size(new)); % first session (first valid study date of the file names)
  for ni = 1:numel(new)
    F = dir(fullfile(Pmdbdir,new{ni},'ses-*','snr-*','*=*=*=*.json')); dat = {};
    for fi = 1:numel(F)
      p = cat_io_dcm2bids_db('nameParts',spm_file(F(fi).name,'basename'));
      dat = [dat regexp(p{min(3,end)},'^(\d{8})','tokens','once')]; %#ok<AGROW>
    end
    for d = unique(dat)
      try
        datetime(d{1},'InputFormat','yyyyMMdd'); first{ni} = d{1}; break
      catch
        % no valid date (e.g. 10101000 of anonymized data)
      end
    end
  end
  [~,o] = sortrows([first(:) new(:)]); new = new(o); first = first(o);
  for ni = 1:numel(new)
    if strcmp(first{ni},'99999999')
      offset = 0;
    else
      offset = days(datetime(1900,1,1) - datetime(first{ni},'InputFormat','yyyyMMdd'));
    end
    S(end+1) = struct('subject',new{ni},'number',max([0 S.number]) + 1,'offset',offset, ...
      'registered',char(datetime('now','Format','yyyy-MM-dd'))); %#ok<AGROW>
  end
  if ~isempty(new) || ~exist(P,'file')
    cat_io_json(P,struct('salt',salt,'subjects',S));
  end
  A.salt = salt;
  A.map  = containers.Map('KeyType','char','ValueType','any');
  for si = 1:numel(S), A.map(S(si).subject) = struct('number',S(si).number,'offset',S(si).offset); end
end
% =========================================================================
function [number,offset] = subject(A,dbsub)
%subject. Number (five digits, e.g. 00001) and date offset (days) of the 
%  database subject dbsub (directory name sub-<serial>-<PatientID>) of the
%  registry A.
  if ~isKey(A.map,dbsub)
    error('cat_io_dcm2bids_anonymize:unknownSubject','Subject "%s" is not registered.',dbsub);
  end
  e = A.map(dbsub); number = sprintf('%05d',e.number); offset = e.offset;
end
% =========================================================================
function [label,offset] = session(date,offset,unit)
%session. Session label of the scan date (yyyymmdd) of a subject with the 
%  date offset of the subject (days, its first session at 1900-01-01, 0 for
%  shifted dates, e.g. of an anonymized DICOM export): the time since the 
%  first session in the unit (d-days, w-weeks, m-months, y-years, rounded)
%  with four digits, e.g. w0000 and w0052 one year later (n for sessions
%  before the first session, e.g. wn0002), and the date offset of the 
%  session that shifts its date to the date of the label (1900-01-01 plus
%  the time of the label, e.g. 1901-01-01 for y0001). 
  len = struct('d',1,'w',7,'m',365.25/12,'y',365.25);
  if ~isfield(len,unit)
    error('cat_io_dcm2bids_anonymize:unit','Unknown unit "%s" (d, w, m, or y).',unit);
  end
  d0 = datetime(1900,1,1); d = datetime(date(1:8),'InputFormat','yyyyMMdd');
  n  = round(days(d + days(offset) - d0) / len.(unit));
  if n < 0, label = sprintf('%sn%04d',unit,-n); else, label = sprintf('%s%04d',unit,n); end
  switch unit
    case 'd', dq = d0 + days(n);
    case 'w', dq = d0 + days(7 * n);
    case 'm', dq = d0 + calmonths(n);
    case 'y', dq = d0 + calyears(n);
  end
  offset = days(dq - d);
end
% =========================================================================
function t = acqTime(stamp,offset,unit,level)
%acqTime. Acquisition time of a session (acq_time of the BIDS sessions 
%  file) by its study date and time stamp (yyyymmddHHMMSS) for the 
%  anonymization level (default 2): 
%    0 - the date and time of the study (e.g. 2025-03-10T11:50:05)
%    1 - the date and the full hour (rounded down, i.e. 10:00 for 10:00 to 
%        10:59, e.g. 2025-03-10T11:00:00)
%    2 - the date of the session label (1900-01-01 plus the time since the 
%        first session in the unit with the date offset of the subject, see
%        session) and the center of the 3-hour interval of the time of the 
%        day (see dayInterval), e.g. 1900-12-31T10:30:00
%  ('n/a' without date or time). Shifted dates (e.g. of an anonymized DICOM
%  export) give the same time again for the level 2. 
  if ~exist('level','var'), level = 2; end
  t = 'n/a';
  if isempty(regexp(char(stamp),'^\d{12}','once')), return; end
  try
    d   = datetime(stamp(1:8),'InputFormat','yyyyMMdd');
    hms = [str2double(stamp(9:10)) str2double(stamp(11:12)) 0];
    if numel(stamp) >= 14, hms(3) = str2double(stamp(13:14)); end
    tod = [3600 60 1] * hms(:);
    if ~(isfinite(tod) && tod >= 0 && tod < 86400), return; end
    switch level
      case 0, s = tod;
      case 1, s = 3600 * floor(tod / 3600);
      otherwise
        [~,so] = session(stamp(1:8),offset,unit); d = d + days(so); s = dayInterval(tod);
    end
    t = sprintf('%sT%02d:%02d:%02d',char(d,'yyyy-MM-dd'),floor(s / 3600),floor(mod(s,3600) / 60),mod(s,60));
  catch
    % no valid date
  end
end
% =========================================================================
function c = dayInterval(s)
%dayInterval. Center of the 3-hour interval of the time of the day s 
%  (seconds), i.e. 01:30, 04:30, ..., 22:30 (as seconds, NaN for times 
%  that are not within a day), the time of the day of the anonymization 
%  level 2 (e.g. for daytime effects without the exact time).
  if ~(isfinite(s) && s >= 0 && s < 86400), c = nan; return; end
  c = (3 * floor(s / 10800) + 1.5) * 3600;
end
% =========================================================================
function s = shiftDate(s,offset)
%shiftDate. Dates yyyymmdd (8 digits, e.g. DICOM dates DA and the date of
%  datetimes DT) of the text s shifted by offset days (other text kept).
  if isempty(s) || offset == 0, return; end
  [i0,i1] = regexp(s,'(?<!\d)\d{8}(?!\d)','start','end');
  if isempty(i0) && numel(s) >= 8 && all(isstrprop(s(1:8),'digit')), i0 = 1; i1 = 8; end % datetime
  for k = 1:numel(i0)
    try
      d = datetime(s(i0(k):i1(k)),'InputFormat','yyyyMMdd') + days(offset);
      s(i0(k):i1(k)) = char(d,'yyyyMMdd');
    catch
      % no date
    end
  end
end
% =========================================================================
function s = shiftSession(s,vr,X)
%shiftSession. Dates and times of the DICOM value s (VR DA, TM, or DT, 
%  several values separated by backslashes) of a session with the study 
%  date X.date (yyyymmdd) and start time X.time (seconds of the day), 
%  shifted to the date of the session label (X.date plus X.offset days) 
%  with the start of the session at X.start (seconds of the day, e.g. the 
%  center of its 3-hour interval, see dayInterval, default 00:00:00), i.e.
%  the times keep the time since the start of the session (e.g. the slice
%  timing is kept, also for sessions over midnight, whereas dates without 
%  time keep the date of the label). Dates without time that are not on the day
%  of the session or the next day (e.g. of a calibration) are shifted by 
%  X.offset days. Without (valid) study date, the dates are shifted by
%  X.offset days and the times as above. Values that are no dates or times
%  are kept.
  if isfield(X,'start') && isfinite(X.start), st = X.start; else, st = 0; end
  if isempty(X.date)
    if ~strcmp(vr,'TM'), s = shiftDate(s,X.offset); return; end
  else
    d0 = datetime(X.date,'InputFormat','yyyyMMdd'); dq = d0 + days(X.offset);
  end
  pat = struct('DA','^(\s*)(\d{8})()()()','TM','^(\s*)()(\d{2})(:?\d{2})?(:?\d{2})?', ...
    'DT','^(\s*)(\d{8})(\d{2})?(\d{2})?(\d{2})?');
  v = strsplit(s,'\');
  for k = 1:numel(v)
    [t,m] = regexp(v{k},pat.(vr),'tokens','match','once');
    if isempty(t), continue; end
    [ws,da,hh,mm,ss] = t{:};
    try
      if isempty(hh) % date: the session (also over midnight) or shifted
        d = datetime(da,'InputFormat','yyyyMMdd');
        if any(days(d - d0) == [0 1]), d = dq; else, d = d + days(X.offset); end
        n = [ws char(d,'yyyyMMdd')];
      else % time (and date): time since the start of the session
        tod = 3600 * str2double(hh) + 60 * str2double(['0' strrep(mm,':','')]) + str2double(['0' strrep(ss,':','')]);
        if isempty(da)
          tod = mod(tod - X.time + st,86400); n = ws;
        else
          d = dq + seconds(days(datetime(da,'InputFormat','yyyyMMdd') - d0) * 86400 + tod - X.time + st);
          tod = round(seconds(d - dateshift(d,'start','day'))); n = [ws char(dateshift(d,'start','day'),'yyyyMMdd')];
        end
        n = [n sprintf('%02d',floor(tod / 3600))];
        if ~isempty(mm), n = [n mm(1:end-2) sprintf('%02d',floor(mod(tod,3600) / 60))]; end
        if ~isempty(ss), n = [n ss(1:end-2) sprintf('%02d',mod(tod,60))]; end
      end
      v{k} = [n v{k}(numel(m) + 1:end)];
    catch
      % no date or time
    end
  end
  s = strjoin(v,'\');
end
% =========================================================================
function u = hashUID(u,salt)
%hashUID. UID of the same length with the digits of a hash (SHA-256) of
%  the secret salt and the UID u (the same UID gives the same new UID) and
%  without leading zeros of its components, besides the standard UIDs of
%  DICOM (1.2.840.10008.*, e.g. transfer syntax and SOP class).
  u = strtrim(u(u ~= char(0)));
  if numel(u) < 8 || startsWith(u,'1.2.840.10008.') || ~all(ismember(u,'0123456789.')), return; end
  d  = find(u ~= '.'); dig = '';
  md = java.security.MessageDigest.getInstance('SHA-256'); k = 0;
  while numel(dig) < numel(d)
    h   = typecast(int8(md.digest(uint8(sprintf('%s|%s|%d',salt,u,k)))),'uint8');
    dig = [dig char('0' + mod(double(h(:)'),10))]; k = k + 1; %#ok<AGROW>
  end
  u(d) = dig(1:numel(d));
  for p = [1 find(u == '.') + 1] % no leading zeros of components
    if p < numel(u) && u(p) == '0' && u(p + 1) ~= '.', u(p) = '1'; end
  end
end
% =========================================================================
function b = hashUIDsInBytes(b,salt)
%hashUIDsInBytes. UIDs in the bytes b (e.g. the text of the Siemens CSA
%  header with the UIDs of reference images) replaced by UIDs of the same
%  length (see hashUID), i.e. numbers with at least three dots and 16
%  characters (also after a dot, e.g. in a list).
  t = char(b(:)');
  [i0,i1] = regexp(t,'(?<!\d)\d+(\.\d+){3,}','start','end');
  for k = 1:numel(i0)
    if i1(k) - i0(k) + 1 >= 16
      b(i0(k):i1(k)) = uint8(hashUID(t(i0(k):i1(k)),salt));
    end
  end
end
% =========================================================================
function b = shiftDatesInBytes(b,X,known)
%shiftDatesInBytes. Dates in the bytes b (e.g. the text of private 
%  elements) shifted as the dates and times of the session X (see 
%  shiftSession): yyyymmdd (with the time HHMMSS of a datetime) and 
%  yyyy-mm-dd of the years 1900-2099 at the beginning of a number, and the
%  known dates (cell of yyyymmdd, e.g. the dates of the header) also after
%  a dot (e.g. in parts of UIDs, but not in decimals by chance).
  if ~exist('known','var'), known = {}; end
  t = char(b(:)');
  [i0,i1] = regexp(t,'(?<![\d.])(19|20)\d\d(-?)(0[1-9]|1[0-2])\2(0[1-9]|[12]\d|3[01])','start','end');
  for di = 1:numel(known)
    [k0,k1] = regexp(t,['(?<!\d)' known{di}],'start','end');
    i0 = [i0 k0]; i1 = [i1 k1]; %#ok<AGROW>
  end
  [i0,ui] = unique(i0); i1 = i1(ui);
  for k = 1:numel(i0)
    s = t(i0(k):i1(k));
    if numel(s) == 10 % yyyy-mm-dd
      d = shiftSession(strrep(s,'-',''),'DA',X); b(i0(k):i1(k)) = uint8([d(1:4) '-' d(5:6) '-' d(7:8)]);
    elseif i1(k) + 6 <= numel(t) && ~isempty(regexp(t(i1(k)+1:i1(k)+6),'^([01]\d|2[0-3])[0-5]\d[0-5]\d$','once'))
      b(i0(k):i1(k)+6) = uint8(shiftSession(t(i0(k):i1(k)+6),'DT',X)); % datetime
    else
      b(i0(k):i1(k)) = uint8(shiftSession(s,'DA',X));
    end
  end
end
% =========================================================================
function V = sidecar(V,level)
%sidecar. Scanner fields of the BIDS sidecar V (or a protocol variant) of 
%  the anonymization level 2 (level > 1), as a scanner (and therefore its
%  site) can be identified by its model, coils, and software version in a
%  study with few scanners: the model (ManufacturersModelName), the coil 
%  fields (ReceiveCoilName, CoilString, ReceiveCoilActiveElements), and 
%  scanner identifiers (DeviceSerialNumber, StationName, institution) are
%  removed and the SoftwareVersions are replaced by their type (see 
%  softwareType, e.g. XA for syngo MR XA60), whereas the manufacturer and 
%  field strength are kept (e.g. for the harmonization). 
  if level < 2 || ~isstruct(V), return; end
  rm = {'ManufacturersModelName','ReceiveCoilName','CoilString','ReceiveCoilActiveElements', ...
    'DeviceSerialNumber','StationName','InstitutionName','InstitutionalDepartmentName','InstitutionAddress'};
  V  = rmfield(V,intersect(fieldnames(V),rm));
  if isfield(V,'SoftwareVersions')
    t = softwareType(V.SoftwareVersions);
    if isempty(t), V = rmfield(V,'SoftwareVersions'); else, V.SoftwareVersions = t; end
  end
end
% =========================================================================
function t = softwareType(s)
%softwareType. Type (major version) of the scanner software version s:
%  the software line of Siemens (e.g. XA for syngo MR XA60, VE for syngo MR
%  E11, VB for syngo_MR_B17 or N4_VB17A), the platform of GE (e.g. DV for 
%  27\LX\MR Software release:DV26.0), the major version of Philips and other
%  manufacturers with R (e.g. R5 for 5.6.1 or R3 for R3.2.2), '' if not 
%  available. 
  t = ''; s = char(string(s)); 
  if isempty(s), return; end
  tok = regexp(s,'release:\s*([A-Za-z]+)\d','tokens','once');                       % GE
  if ~isempty(tok), t = upper(tok{1}); return; end
  tok = regexp(s,'(?<![A-Za-z])(X[A-Z]|V[A-E]|[A-E])\d{2}','tokens','once');          % Siemens
  if ~isempty(tok)
    t = upper(tok{1}); if isscalar(t), t = ['V' t]; end
    return
  end
  tok = regexp(s,'(?<![A-Za-z0-9.])R?(\d+)(?=[.\\_]|$)','tokens','once');            % Philips and others
  if ~isempty(tok), t = ['R' tok{1}]; end
end
