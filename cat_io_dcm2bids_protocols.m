function varargout = cat_io_dcm2bids_protocols(action,varargin)
%cat_io_dcm2bids_protocols. Sites, protocols, and overview tables of cat_io_dcm2bids.
%  The site dictionary (getSites, setupSites, or without dictionary the 
%  site registry of the database, see siteRegistry and siteProposal) and
%  the protocol definitions (getProtocols, JSON files with the relevant MR
%  parameters) are read once.
%  Each scan is tested against the protocols (testProtocols, with fields 
%  that are normalized for the scanner software and dcm2niix version, see 
%  normalizeFields, and the protocol fit rating, see protocolRating). Its 
%  protocol variant (protocolVariant, normalized JSON files that can be 
%  used to define new protocols) and its patient (updateDBtables) are added
%  to the tables of the database (catDCM2BIDSdb-tables, see also 
%  cat_io_dcm2bids_tables). cleanupVjson removes the patient-specific 
%  fields of a sidecar or keeps only the protocol (or BIDS) relevant fields.
%
%  varargout = cat_io_dcm2bids_protocols(action,varargin)
%
%  Actions (see the help of the local functions):
%    sites = cat_io_dcm2bids_protocols('getSites',Pcenterdict)
%    protocols = cat_io_dcm2bids_protocols('getProtocols',Pprodictdirs)
%    site = cat_io_dcm2bids_protocols('setupSites',Vjson,sites[,opts])
%    siteID = cat_io_dcm2bids_protocols('siteID',serial,sites)
%    sites = cat_io_dcm2bids_protocols('siteRegistry',Pmdbdir,opts,dict)
%    [abbr,city,rest] = cat_io_dcm2bids_protocols('siteProposal',address,field,vendor,opts)
%    C = cat_io_dcm2bids_protocols('siteColors',Pcenterdict,Pmdbdir,labels)
%    d = cat_io_dcm2bids_protocols('acquisitionDuration',V)
%    [match,pname,mismatchstr,Pnfailedid,Pfname,PT] = cat_io_dcm2bids_protocols('testProtocols',Pjson,Vjson,protocols,tol,job,sites)
%    variant = cat_io_dcm2bids_protocols('protocolVariant',Ptbldir,PT,ProtocolName,datatype,type,tol,opts)
%    Vjson = cat_io_dcm2bids_protocols('cleanupVjson',Vjson,type)
%    cat_io_dcm2bids_protocols('updateDBtables',Vjson,Ptbldir,sites,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'getSites','getProtocols','setupSites','siteID','siteRegistry','siteProposal','siteColors','acquisitionDuration','testProtocols', ...
          'protocolVariant','cleanupVjson','updateDBtables'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_protocols:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function sites = getSites(Pcenterdict)
%getSites. Site dictionary of the JSON file Pcenterdict{1}, i.e. a list of
%  scanners with the fields DeviceSerialNumber and InstitutionAbbreviation
%  (further fields are optional and can differ between the scanners, e.g.
%  the Color of the site in the report figure, see siteColors, and the 
%  siteID of the anonymization level 2, see setupSites) as cell array with
%  the serial number, abbreviation, and siteID ('' if not defined) of each
%  scanner. Color entries (Cities and Color, see siteColors) are skipped 
%  and other entries without these two fields are ignored (warning).
  sites = cell(0,3);
  if isempty(Pcenterdict) || isempty(Pcenterdict{1}), return; end
  if ~exist(Pcenterdict{1},'file')
    error('cat_io_dcm2bids:Pcenterdict','Center dictionary file "%s" is not existing.', Pcenterdict{1});
  end
  S = jsonList(cat_io_json(Pcenterdict{1}));
  for si = 1:numel(S)
    if isstruct(S{si}) && all(isfield(S{si},{'DeviceSerialNumber','InstitutionAbbreviation'}))
      sites(end+1,:) = {char(string(S{si}.DeviceSerialNumber)), char(string(S{si}.InstitutionAbbreviation)), ...
        char(string(fieldOr(S{si},'siteID')))}; %#ok<AGROW>
    elseif ~(isstruct(S{si}) && all(isfield(S{si},{'Cities','Color'})))
      cat_io_cprintf('warn','  Site dictionary: entry %d without DeviceSerialNumber or InstitutionAbbreviation ignored.\n',si);
    end
  end
end
% =========================================================================
function C = siteColors(Pcenterdict,Pmdbdir,labels)
%siteColors. Colors of the site labels (RGB hex, e.g. #A82117, '' if not
%  defined) by the site dictionary Pcenterdict{1} (see getSites): the
%  Color of the scanner entry, otherwise the Color of the color entry of
%  the city of the scanner (address of the dictionary or of the site
%  registry of the database Pmdbdir, see siteRegistry), i.e. of the 
%  entries with the fields Cities (list of cities, compared transliterated
%  and without case as beginning of the city of the address, "*" for all
%  other cities) and Color, e.g. 
%    {"Name": "DZPG Halle-Jena-Magdeburg", "Cities": ["Halle","Jena","Magdeburg"], "Color": "#A82117"}
%  for all scanners of a site, also for scanners that are added later.
%  Sites with the same color get a slightly different brightness (in the 
%  order of their labels: the color, lighter, darker, ...).
  C = repmat({''},1,numel(labels));
  if isempty(Pcenterdict) || isempty(Pcenterdict{1}) || ~exist(Pcenterdict{1},'file'), return; end
  J    = jsonList(cat_io_json(Pcenterdict{1}));
  isen = @(e,f) isstruct(e) && all(isfield(e,f));
  scn  = J(cellfun(@(e) isen(e,{'DeviceSerialNumber','InstitutionAbbreviation'}),J)); % scanners
  col  = J(cellfun(@(e) isen(e,{'Cities','Color'}),J));                                % colors of cities
  Preg = fullfile(Pmdbdir,'catDCM2BIDSsites.json'); reg = {};
  if exist(Preg,'file'), reg = jsonList(cat_io_json(Preg)); reg = reg(cellfun(@(e) isen(e,{'InstitutionAbbreviation'}),reg)); end
  low  = @(s) lower(cat_io_dcm2bids_bids('bidsLabel',char(string(s))));
  abbr = @(E) cellfun(@(e) cat_io_dcm2bids_bids('bidsLabel',char(string(e.InstitutionAbbreviation)),1),E, ...
    'UniformOutput',false); % (labels as setupSites)
  for li = 1:numel(labels)
    e = [scn(strcmpi(abbr(scn),labels{li})); reg(strcmpi(abbr(reg),labels{li}))];
    if ~isempty(e) && isfield(e{1},'Color'), C{li} = hexColor(e{1}.Color); end
    if ~isempty(C{li}), continue; end
    city = ''; 
    if ~isempty(e) && isfield(e{1},'InstitutionAddress'), city = low(addressCity(e{1}.InstitutionAddress)); end
    for ci = 1:numel(col) % city of a color entry, otherwise all other cities ("*")
      cities = cellstr(col{ci}.Cities);
      if ~isempty(city) && any(cellfun(@(c) ~strcmp(c,'*') && startsWith(city,low(c)),cities))
        C{li} = hexColor(col{ci}.Color); break
      end
    end
    if isempty(C{li})
      ci = find(cellfun(@(c) any(strcmp(cellstr(c.Cities),'*')),col),1);
      if ~isempty(ci), C{li} = hexColor(col{ci}.Color); end
    end
  end
  % brightness of the sites with the same color (in the order of their labels)
  [~,so] = sort(lower(cellfun(@char,labels(:)','UniformOutput',false)));
  [uc,~,ui] = unique(C(so)); ui = ui(:)';
  for ci = find(~cellfun('isempty',uc(:)'))
    k = so(ui == ci);
    for ki = 2:numel(k)
      rgb = sscanf(C{k(ki)}(2:end),'%2x')' / 255; s = 0.18 * ceil((ki-1)/2);
      if mod(ki,2) == 0, rgb = rgb + s * (1 - rgb); else, rgb = rgb * (1 - s); end % lighter, darker
      C{k(ki)} = hexColor(rgb);
    end
  end
end
% =========================================================================
function h = hexColor(c)
%hexColor. Color c (RGB hex text with or without #, or RGB values from 0
%  to 1) as text #RRGGBB ('' if not valid).
  h = '';
  if isnumeric(c) && numel(c) == 3
    h = sprintf('#%02X%02X%02X',round(255 * min(1,max(0,double(c(:)')))));
  elseif ischar(c) || isstring(c)
    c = regexprep(char(c),'^#','');
    if ~isempty(regexp(c,'^[0-9a-fA-F]{6}$','once')), h = ['#' upper(c)]; end
  end
end
% =========================================================================
function city = addressCity(address)
%addressCity. City of the address of a scanner (second entry of the
%  Siemens address, e.g. "Street,City,District,...", '' if not available).
  parts = strtrim(strsplit(char(string(address)),','));
  if numel(parts) < 2, city = ''; else, city = parts{2}; end
end
% =========================================================================
function C = jsonList(S)
%jsonList. Entries of a JSON list of objects as cell array of structures, 
%  i.e. of a structure array (objects with the same fields) or a cell array
%  (objects with different fields) of jsondecode.
  if isstruct(S), C = num2cell(S(:)); elseif iscell(S), C = S(:); else, C = {}; end
end
% =========================================================================
function protocols = getProtocols(Pprodictdirs)
  protocols = cell(numel(Pprodictdirs),3); pdi = 0; 
  if isempty(Pprodictdirs) || isempty(Pprodictdirs{1}), return; end
  for di = 1:numel(Pprodictdirs)
    if isempty(Pprodictdirs{di}), continue; end  
    if ~exist(Pprodictdirs{di},'dir')
      error('cat_io_dcm2bids:protocoldir','Protocol directory %d  "%s" does not exist.\n',di,Pprodictdirs{di}); 
    end
    Pprodictsubdirs = cat_vol_findfiles(Pprodictdirs{di},'*.json');
    Pprodictsubdirs(cat_io_contains(Pprodictsubdirs,[filesep 'qc'])) = []; 
    if isempty(Pprodictsubdirs) % 'cat_io_dcm2bids:noProtocols',
      cat_io_cprintf('warn','No protocols found in:\n  %s\n',Pprodictdirs{di});
    end
    for fi = 1:numel(Pprodictsubdirs)
      pdi = pdi + 1; 
      protocols{pdi,1} = spm_file(Pprodictdirs{di},'basename'); 
      protocols{pdi,2} = Pprodictsubdirs{fi}; 
      protocols{pdi,3} = normalizeFields(cat_io_json(Pprodictsubdirs{fi})); % see testProtocols
    end
  end
  protocols = protocols(1:max(1,pdi),:); % remove empty rows of directories without protocols
  if isempty(protocols{1}) % 'cat_io_dcm2bids:noProtocols',
    cat_io_cprintf('warn','No protocols found in any protocol directory:\n  %s\n', ...
      strjoin(Pprodictdirs(:)',sprintf('\n  ')));
  end
end
% =========================================================================
function site = setupSites(Vjson,sites,opts)
%setupSites. Site label of a scan by its DeviceSerialNumber and the site
%  dictionary sites (serial number, abbreviation, and siteID of each 
%  scanner, first the same serial number, then a serial number that 
%  contains it), otherwise the serial number. 
%  For the anonymization level 2 (opts.anonymize > 1), the label contains 
%  no location (city or abbreviation of the dictionary): the siteID of the
%  scanner (a manually defined code, e.g. given centrally for all sites of
%  a federated study, empty by default) with the code of the manufacturer 
%  and the field strength (see siteRest), e.g. 07SI30T, or without siteID 
%  only the manufacturer and field strength, e.g. SI30T, i.e. the site is 
%  not available (and subject labels of other databases are not unique).
%  Scans of an anonymized DICOM export keep their site label (serial number
%  with the form of a site label, see isSiteLabel).
  if ~exist('opts','var'), opts = struct(); end
  serial = char(string(Vjson.DeviceSerialNumber)); 
  [sid,k] = siteID(serial,sites); 
  if isfield(opts,'anonymize') && opts.anonymize > 1
    if isSiteLabel(serial)
      site = serial; % (anonymized DICOM export)
    else
      site = [cat_io_dcm2bids_bids('bidsLabel',sid,1) siteRest(numField(Vjson,'MagneticFieldStrength'), ...
        textField(Vjson,'Manufacturer'),opts)]; 
    end
  elseif ~isempty(k)
    site = sites{k,2};
  else
    site = serial; 
  end
  site = cat_io_dcm2bids_bids('bidsLabel',site,1); 
end
% =========================================================================
function [sid,k] = siteID(serial,sites)
%siteID. siteID (anonymization level 2, '' if not defined) and row k of 
%  the scanner with the serial number in the site dictionary sites (see 
%  getSites and siteRegistry: first the same serial number, then a serial
%  number that contains it, empty if not available).
  sid = ''; k = []; 
  if isempty(sites), return; end
  k = find(strcmp(cellfun(@(x) char(string(x)),sites(:,1),'UniformOutput',false),char(serial)),1); 
  if isempty(k), k = find(cat_io_contains(sites(:,1),char(serial)),1); end
  if ~isempty(k) && size(sites,2) > 2, sid = char(string(sites{k,3})); end
end
% =========================================================================
function ok = isSiteLabel(s)
%isSiteLabel. Serial number with the form of a site label of an anonymized
%  DICOM export, i.e. a code (abbreviation or siteID, optional) and the 
%  manufacturer (letters) and field strength (two digits and T), e.g. 
%  BerSI30T, 07SI30T, or SI30T.
  ok = ~isempty(regexp(char(s),'^[A-Za-z0-9]{0,16}[A-Z]\d{2}T$','once'));
end
% =========================================================================
function sites = siteRegistry(Pmdbdir,opts,dict)
%siteRegistry. Site abbreviations of the scanners of the database Pmdbdir
%  (as cell array of getSites): the abbreviations of the site dictionary
%  dict (see getSites) and for all other scanners the abbreviations of the
%  registry of the database (catDCM2BIDSsites.json), i.e. the BIDS labels
%  contain no serial numbers. New scanners get the abbreviation of their 
%  city, manufacturer, and field strength (see siteProposal, e.g. 
%  HeiSI30T), where scanners with the same abbreviation (e.g. several 
%  scanners of a site, Heidelberg and Heidenheim, or an abbreviation of the
%  dictionary) get a capital letter A to Z instead of the last letter of 
%  the city (e.g. HeASI30T and HeBSI30T, ordered by their serial number, 
%  with a warning and the serial number after Z). Serial numbers with the
%  form of a site label (e.g. of the anonymized DICOM export) keep it. The 
%  abbreviations of the registry are never changed, as they are part of
%  the BIDS labels, i.e. a new scanner with the abbreviation (proposal) of
%  a registered scanner gets the next free letter. Abbreviations are 
%  compared without case (e.g. directories of macOS). The entries have an
%  empty siteID (the code of the anonymization level 2 that has to be 
%  defined manually, see setupSites), which is added to the entries of 
%  former versions. sites contains the serial number, abbreviation, and 
%  siteID of each scanner (dictionary first).
  Preg = fullfile(Pmdbdir,'catDCM2BIDSsites.json');
  R = struct('DeviceSerialNumber',{},'InstitutionAbbreviation',{},'siteID',{},'Proposal',{},'InstitutionName',{}, ...
    'InstitutionAddress',{},'MagneticFieldStrength',{},'Manufacturer',{},'Registered',{});
  nosid = false; 
  if exist(Preg,'file') % (entries with other fields, e.g. without manufacturer of older versions)
    J = jsonList(cat_io_json(Preg)); fn = fieldnames(R)'; 
    for ji = 1:numel(J)
      e = cell2struct(cellfun(@(f) fieldOr(J{ji},f),fn,'UniformOutput',false),fn,2); R(end+1) = e; %#ok<AGROW>
      nosid = nosid || ~isfield(J{ji},'siteID'); 
    end
  end
  if ~exist('dict','var') || isempty(dict), dict = cell(0,3); end
  if size(dict,2) < 3, dict(:,3) = {''}; end
  known = [cellfun(@(x) char(string(x)),{R.DeviceSerialNumber},'UniformOutput',false) ...
           cellfun(@(x) char(string(x)),dict(:,1)','UniformOutput',false)];

  % new scanners of the database (first sidecar of each subject)
  N = struct('serial',{},'name',{},'address',{},'field',{},'vendor',{});
  S = dir(fullfile(Pmdbdir,'sub-*')); S = S([S.isdir]);
  for si = 1:numel(S)
    tok = regexp(S(si).name,'^sub-(.+?)-','tokens','once'); % serial number of the subject directory
    if ~isempty(tok) && any(strcmp(known,tok{1})), continue; end
    P = cat_vol_findfiles(fullfile(S(si).folder,S(si).name),'*=*=*=*.json',struct('depth',3)); % (sub-*/ses-*/snr-*)
    P = P(~startsWith(spm_file(P,'filename'),'catDCM2BIDS'));
    if isempty(P), continue; end
    V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(P{1}));
    serial = char(V.DeviceSerialNumber);
    if any(strcmp([known {N.serial}],serial)), continue; end
    N(end+1) = struct('serial',serial,'name',textField(V,'InstitutionName'), ...
      'address',textField(V,'InstitutionAddress'),'field',numField(V,'MagneticFieldStrength'), ...
      'vendor',textField(V,'Manufacturer')); %#ok<AGROW>
  end
  if ~isempty(N)
    [~,o] = sort({N.serial}); N = N(o);
    prop = cell(1,numel(N)); city = prop; rest = prop;
    for ni = 1:numel(N)
      [prop{ni},city{ni},rest{ni}] = siteProposal(N(ni).address,N(ni).field,N(ni).vendor,opts);
      if isSiteLabel(N(ni).serial) % site label of an anonymized DICOM export
        prop{ni} = N(ni).serial; 
      end
    end
    used = [cellfun(@char,{R.InstitutionAbbreviation},'UniformOutput',false) ...
            cellfun(@(x) char(string(x)),dict(:,2)','UniformOutput',false)]; % (case-insensitive)
    for ni = 1:numel(N)
      abbr = prop{ni};
      if any(strcmpi(used,abbr)) || any(strcmpi({R.Proposal},abbr)) || nnz(strcmpi(prop,abbr)) > 1
        abbr = ''; % another scanner with this abbreviation (also before): capital letter
        for L = 'A':'Z'
          if ~any(strcmpi(used,[city{ni}(1:end-1) L rest{ni}])), abbr = [city{ni}(1:end-1) L rest{ni}]; break; end
        end
        if isempty(abbr)
          abbr = cat_io_dcm2bids_bids('bidsLabel',N(ni).serial);
          cat_io_cprintf('warn','  No free site abbreviation %s for scanner %s (A to Z used): serial number used.\n', ...
            prop{ni},N(ni).serial);
        end
      end
      used{end+1} = abbr; %#ok<AGROW>
      R(end+1) = struct('DeviceSerialNumber',N(ni).serial,'InstitutionAbbreviation',abbr,'siteID','','Proposal',prop{ni}, ...
        'InstitutionName',N(ni).name,'InstitutionAddress',N(ni).address,'MagneticFieldStrength',N(ni).field, ...
        'Manufacturer',N(ni).vendor,'Registered',char(datetime('now','Format','yyyy-MM-dd'))); %#ok<AGROW>
      fprintf('  New site abbreviation %s of scanner %s (%s)\n',abbr,N(ni).serial,N(ni).address);
    end
  end
  if ~isempty(N) || nosid % (new scanners or the siteID field of former versions)
    cat_io_json(Preg,R);
  end
  sid   = cellfun(@(x) char(string(x)),{R.siteID},'UniformOutput',false); 
  sites = [dict(:,1:3); {R.DeviceSerialNumber}' {R.InstitutionAbbreviation}' sid']; % dictionary first
end
% =========================================================================
function [abbr,city,rest] = siteProposal(address,field,vendor,opts)
%siteProposal. Abbreviation of a scanner with a constant length: the first
%  opts.sitelength letters of the city of its address (second entry of the
%  Siemens address, e.g. "Street,City,District,...", transliterated, first
%  letter as capital, Unk if not available), the code of the manufacturer
%  vendor with opts.sitevendor letters (see vendorCode, e.g. SI or S for 
%  Siemens), and the field strength in tenths of Tesla with two digits 
%  (rounded, e.g. 15T, 30T, or 70T for 1.5, 3, and 6.98 T, 00T if not 
%  available), e.g. HeiSI30T, with the city and the rest of the 
%  abbreviation.
  n = opts.sitelength; 
  parts = strtrim(strsplit(char(address),','));
  city  = regexprep(cat_io_dcm2bids_bids('bidsLabel',char(parts(min(2,numel(parts))))),'[^a-zA-Z]','');
  if numel(parts) < 2 || numel(city) < max(2,n), city = 'Unknown'; end
  city = [upper(city(1)) lower(city(2:n))]; 
  rest = siteRest(field,vendor,opts);
  abbr = [city rest];
end
% =========================================================================
function rest = siteRest(field,vendor,opts)
%siteRest. Manufacturer and field strength part of the site label: the 
%  code of the manufacturer vendor with opts.sitevendor letters (see 
%  vendorCode, default 2) and the field strength in tenths of Tesla with 
%  two digits (rounded, e.g. 15T, 30T, or 70T for 1.5, 3, and 6.98 T, 00T 
%  if not available), e.g. SI30T.
  if ~isfield(opts,'sitevendor'), opts.sitevendor = 2; end
  if ~isempty(field) && isfinite(field), fstr = sprintf('%02dT',round(10*field - 0.01)); else, fstr = '00T'; end
  rest = [vendorCode(vendor,opts.sitevendor) fstr];
end
% =========================================================================
function v = vendorCode(vendor,n)
%vendorCode. Code of the manufacturer vendor with n letters (0 to 2): SI
%  (Siemens), GE, PH (Philips), CA (Canon), TO (Toshiba), HI (Hitachi), FU
%  (Fujifilm), UI (United Imaging), BR (Bruker), ES (Esaote), HY 
%  (Hyperfine), NE (Neusoft), otherwise the first letters of the 
%  manufacturer (XX if not available), or for one letter the first letter
%  of the code (e.g. S for Siemens).
  codes = {'siemens','SI'; '^ge( |$)|general electric','GE'; 'philips','PH'; 'canon','CA'; ...
    'toshiba','TO'; 'hitachi','HI'; 'fujifilm','FU'; 'united imaging|^uih','UI'; 'bruker','BR'; ...
    'esaote','ES'; 'hyperfine','HY'; 'neusoft','NE'};
  m = lower(strtrim(char(vendor))); 
  k = find(~cellfun('isempty',regexp(m,codes(:,1),'once')),1); 
  if isempty(k), v = [upper(regexprep(m,'[^a-z]','')) 'XX']; else, v = codes{k,2}; end
  v = v(1:n);
end
% =========================================================================
function d = acquisitionDuration(V)
%acquisitionDuration. Acquisition duration of a scan in seconds by its 
%  sidecar V: the AcquisitionDuration of dcm2niix (e.g. the total scan time
%  of Siemens or of enhanced DICOM), otherwise the number of volumes times
%  the TR of time series (NumberOfVolumes or the 4th image dimension of 
%  Dimensions), NaN if not available (e.g. single volumes without scan 
%  time). (Not in the BIDS sidecars, where AcquisitionDuration is the 
%  duration of a volume and excludes the RepetitionTime.)
  num = @(fn) isfield(V,fn) && isnumeric(V.(fn)) && ~isempty(V.(fn));
  d   = nan;
  if num('AcquisitionDuration')
    d = double(V.AcquisitionDuration(1));
  elseif num('RepetitionTime')
    nv = nan;
    if num('NumberOfVolumes'), nv = double(V.NumberOfVolumes(1)); 
    elseif num('Dimensions') && numel(V.Dimensions) > 3, nv = double(V.Dimensions(4)); end
    if nv > 1, d = nv * double(V.RepetitionTime(1)); end
  end
  d = round(d,1);
end
% =========================================================================
function v = fieldOr(S,fn)
%fieldOr. Field fn of the structure S ('' if not available).
  if isstruct(S) && isfield(S,fn), v = S.(fn); else, v = ''; end
end
% =========================================================================
function v = textField(S,fn)
%textField. Text of the field fn of the structure S ('' if not available).
  if isfield(S,fn) && ~isempty(S.(fn)) && ~isnumeric(S.(fn)), v = char(S.(fn)); else, v = ''; end
end
% =========================================================================
function x = numField(S,fn)
%numField. Number of the field fn of the structure S (NaN if not available).
  if isfield(S,fn) && isnumeric(S.(fn)) && ~isempty(S.(fn)), x = double(S.(fn)(1)); else, x = nan; end
end

% =========================================================================
function [match,pname,mismatchstr,Pnfailedid,Pfname,PT] = testProtocols(Pjson,Vjson,protocols,tol,job,sites)
%testProtocols. Test a scan (sidecar Vjson) against the protocols.
%  The fields of the scan and of the protocols are normalized, i.e., the
%  same protocol results in the same values for different scanner software
%  versions (VE/XA) or dcm2niix versions (see normalizeFields). Numeric
%  fields have to match within the tolerance tol (in percent, only for the
%  rounding of the values), text fields and the entries of lists (order
%  independent) have to be equal (see compareProtocol). Protocol fields
%  that are not available for the scan are listed as 'n/a' without counting.
%  Deviations of key parameters (contrast and resolution) count at least 1
%  and of critical parameters (e.g. phase-encoding direction, sequence, 
%  image component, b-value shells) 10, i.e., such scans are nonconform 
%  (see compareProtocol and job.opts.protocol), as well as a direction in
%  the series name (e.g. PA of "SpinEchoFieldMap_PA") that does not fit to
%  the direction in the protocol name (e.g. a PA field map that was 
%  acquired with AP phase encoding is not conform to the AP protocol).
%  The protocol is selected in two steps to avoid the assignment to another
%  protocol with more equal fields: the protocols with the same sequence
%  (PulseSequenceDetails) and key and critical parameters are preferred, 
%  and of these the protocol with the fewest mismatches is used.
%  The protocol fit measure PFM is the (weighted) number of deviations from
%  the closest protocol and PFR its rating (marks 0.5-10.5), which defines 
%  the protocol status by job.opts.protocol (see protocolRating): conform 
%  scans are named by the protocol, semiconform scans are rated by the QC 
%  definition of the closest protocol, and nonconform scans are unknown. 
%  The scan line shows the (closest) protocol and the PFR in the colors of
%  the marks. The shortened (normalized) protocol of the scan is saved next
%  to its sidecar in the database as template for protocol definitions.
%  match is the index of the conform protocol (0 if not conform), pname its
%  name (otherwise the DICOM protocol name), mismatchstr the mismatches of
%  each protocol {field, scan value, protocol value}, Pnfailedid the index
%  of the closest protocol(s), Pfname the protocol file of conform and 
%  semiconform scans, and PT the result of the test for the tables of the 
%  database (see cat_io_dcm2bids_tables): the protocol fit measure and 
%  rating (PFM, PFR, NaN without protocols), the status, the closest 
%  protocol (name, file, and tested fields), its deviations ({field, scan
%  value, protocol value}, without fields that are not available), and the
%  normalized sidecar (Vn).

  Vjsonn = normalizeFields(Vjson);
  if isempty(protocols) || isempty(protocols{1,1})
  % no given protocols
    mismatchstr{1}  = cell(0,3);
    mismatchcnt(1)  = 0;
    rankcnt(1)      = 0;
  else
    mismatchstr = cell(1,size(protocols,1));
    mismatchcnt = zeros(1,size(protocols,1));
    rankcnt     = zeros(1,size(protocols,1));
    sdesc = ''; if isfield(Vjson,'SeriesDescription'), sdesc = char(Vjson.SeriesDescription); end
    sdir  = nameDirection([char(Vjson.ProtocolName) ' ' sdesc]); 
    for pri = 1:size(protocols,1)
      [mismatchcnt(pri),mismatchstr{pri},seqdiff,keydiff] = ...
        compareProtocol(Vjsonn,protocols{pri,3},tol,job.opts);
      % direction of the series name vs. direction of the protocol name (critical)
      pdir = nameDirection(spm_file(protocols{pri,2},'basename')); 
      if ~isempty(sdir) && ~isempty(pdir) && ~strcmp(sdir,pdir)
        mismatchcnt(pri) = mismatchcnt(pri) + 10; 
        mismatchstr{pri}(end+1,:) = {'NameDirection', sdir, pdir}; 
        keydiff = true; 
      end
      % protocol selection: same sequence and key parameters first
      rankcnt(pri) = mismatchcnt(pri) + 100*seqdiff + 10*keydiff;
    end
  end
  [~,bestid] = min(rankcnt);
  mincnt     = mismatchcnt(bestid);

  % protocol fit measure and rating with the protocol status
  PT = struct('PFM',nan,'PFR',nan,'status','','protocol','','file','','fields',{{}}, ...
    'deviations',{cell(0,3)},'Vn',Vjsonn); 
  if ~isempty(protocols) && ~isempty(protocols{1,1})
    [PT.PFR,PT.status] = protocolRating(mincnt,job.opts); 
    PT.PFM      = mincnt; 
    PT.protocol = spm_file(protocols{bestid,2},'basename'); 
    PT.file     = protocols{bestid,2}; 
    FNb         = fieldnames(protocols{bestid,3}); 
    PT.fields   = FNb(~strncmp(FNb,'x_',2))'; % tested fields (without disabled fields)
    PT.deviations = mismatchstr{bestid}(~strcmp(mismatchstr{bestid}(:,2),'n/a'),:); 
    proSubdir   = PT.status; 
  end
  PFR = PT.PFR; 

  % display progress: (closest) protocol and rating in the colors of the 
  % marks (as the QC ratings, separated by "|")
  % =======================================================================
  fname1     = sprintf('%3d) %s', Vjson.SeriesNumber, spm_str_manip(strrep(sprintf('%s_%s_%s', ...
                strrep(Vjson.PatientID, strrep(char(Vjson.ScanDate),'-',''),''), ...
                strrep(char(Vjson.ScanDate),'-',''), Vjson.ProtocolName),'__','_'),'l55'));
  fprintf('%60s : ',fname1);
  if isnan(PFR)
  % no given protocols: no match but use the DICOM protocol name
    match      = 0;
    Pnfailedid = [];
    pname      = Vjson.ProtocolName;
    Pfname     = '';
    datatype = cat_io_dcm2bids_bids('setupDatatype',pname,Vjson);
    cat_io_cprintf([0 0 0.5],sprintf('%-50s%8s |', [datatype filesep pname], '-'));
    return
  end
  Pnfailedid = find(rankcnt == rankcnt(bestid)); % closest protocols (with the same rank)
  switch proSubdir
    case 'conform'
      match  = bestid; 
      pname  = spm_file(protocols{bestid,2},'basename');
      Pfname = protocols{bestid,2};
      pname0 = spm_str_manip(pname,'l50'); 
    case 'semiconform'
      match  = 0; 
      pname  = Vjson.ProtocolName;
      Pfname = protocols{bestid,2};
      pname0 = spm_str_manip(spm_file(protocols{bestid,2},'basename'),'l50'); 
    otherwise % nonconform
      match  = 0; 
      pname  = Vjson.ProtocolName;
      Pfname = '';
      pname0 = sprintf('Unknown %s protocol',cat_io_dcm2bids_bids('setupDatatype',pname,Vjson)); 
  end
  MarkColor = cat_io_colormaps('marks+',40); % as the QC ratings (see cat_io_dcm2bids_qc)
  col2mark  = @(val) MarkColor(min(size(MarkColor,1)-3,max(1,floor( val/9.5 * size(MarkColor,1)))),:); 
  cat_io_cprintf(col2mark(PFR),sprintf('%-50s%8.1f',pname0,PFR)); 
  fprintf(' |'); 

  % shortened version of the protocol of the scan (normalized values of all
  % fields of the protocols) next to its sidecar as template for protocol 
  % definitions
  proname = strrep(strrep(Vjson.ProtocolName,'_','-'),' ','-');
  site    = setupSites(Vjson,sites);
  FN      = cellfun(@(x) fieldnames(x), protocols(:,3) , 'UniformOutput', false );
  FNM     = intersect(unique(vertcat(FN{:})),fieldnames(Vjsonn));
  if isempty(FNM), return; end
  Vjson4  = struct();
  for fni = 1:numel(FNM)
    Vjson4.(FNM{fni}) = Vjsonn.(FNM{fni});
  end
  cat_io_json(fullfile(fileparts(Pjson),sprintf('%s-%s.json',proname,site)), cleanupVjson(Vjson4));
end
% =========================================================================
function [PFR,status] = protocolRating(PFM,opts)
%protocolRating. Rating PFR (marks 0.5-10.5) of the protocol fit measure
%  PFM, i.e., the weighted number of deviations from the closest protocol 
%  (see compareProtocol), with the best mark 0.5 for opts.protocol.PFM(1) 
%  and the mark 4.5 for opts.protocol.PFM(2) (linear), and the protocol 
%  status by the limits opts.protocol.status: conform (PFR < status(1)), 
%  semiconform, or nonconform (PFR > status(2)). 
  s   = opts.protocol.PFM; 
  PFR = max(0.5,min(10.5, 0.5 + (PFM - s(1)) / (s(2) - s(1)) * 4)); 
  if PFR < opts.protocol.status(1)
    status = 'conform'; 
  elseif PFR > opts.protocol.status(2)
    status = 'nonconform'; 
  else
    status = 'semiconform'; 
  end
end
% =========================================================================
function variant = protocolVariant(Ptbldir,PT,ProtocolName,datatype,type,tol,opts)
%protocolVariant. Protocol variant of a scan for the tables of the database
%  (see cat_io_dcm2bids_tables): the normalized fields of the scan (PT.Vn of
%  testProtocols) that are tested by the closest protocol of conform and 
%  semiconform scans (i.e. without the disabled fields, e.g. with 
%  representation differences of the scanner software), otherwise the 
%  protocol fields (see cleanupVjson, without coil fields), separated by the
%  DICOM protocol name and the image type (BIDS suffix and entities, e.g. 
%  sbref or epi_dir-AP_part-phase). Scans with the same fields (within the
%  tolerance tol) share a variant, which is saved as protocol definition in
%    catDCM2BIDSdb-tables/protocols/<datatype>/<protocol>_<type>[_<n>].json
%  variant is the name of this file. Further fields of a scan (e.g. the
%  image fields after the conversion of a JSON-only scan) are added to its
%  variant, as well as the acquisition duration (see acquisitionDuration),
%  i.e. protocols of these variants also test the duration (within 5%).
  if any(strcmp(PT.status,{'conform','semiconform'}))
    FN = PT.fields; 
  else
    FN = [fieldnames(cleanupVjson(PT.Vn,'protocol'))' ...
      {'NumberOfSlices','NumberOfVolumes','VoxelSize','SliceOrder','BValueShells'}]; 
    FN = setdiff(FN,{'CoilString','ReceiveCoilName','ReceiveCoilActiveElements', ...
      'PhaseEncodingSteps'},'stable'); 
  end
  V = struct(); 
  for fni = 1:numel(FN)
    if isfield(PT.Vn,FN{fni}), V.(FN{fni}) = PT.Vn.(FN{fni}); end
  end
  d = acquisitionDuration(PT.Vn); 
  if ~isnan(d), V.AcquisitionDuration = d; end % (informative, compared within 5%)

  % existing variants of the protocol and image type 
  name = [getFileString(strrep(strrep(ProtocolName,'_','-'),' ','-')) '_' type]; 
  Pdir = fullfile(Ptbldir,'protocols',datatype); 
  if ~exist(Pdir,'dir'), mkdir(Pdir); end
  D    = dir(fullfile(Pdir,[name '*.json'])); 
  D    = D( ~cellfun('isempty',regexp({D.name},['^' regexptranslate('escape',name) '(_\d+)?\.json$'],'once')) ); 
  for di = 1:numel(D)
    Pv  = fullfile(Pdir,D(di).name);
    Vv  = normalizeFields(cat_io_json(Pv));
    FNv = fieldnames(Vv); FNs = fieldnames(V);
    % the same fields or more fields of the scan or of the variant (e.g. the
    % image fields of a scan that is converted after a JSON-only run) with
    % equal values of the common fields (compareProtocol ignores fields of
    % the variant that the scan does not have)
    if (all(ismember(FNv,FNs)) || all(ismember(FNs,FNv))) && compareProtocol(V,Vv,tol,opts) == 0
      variant = spm_file(D(di).name,'basename');
      FNn = setdiff(FNs,FNv,'stable'); % add the further fields of the scan
      if ~isempty(FNn)
        for fni = 1:numel(FNn), Vv.(FNn{fni}) = V.(FNn{fni}); end
        cat_io_json(Pv,Vv);
      end
      return
    end
  end

  % new variant
  variant = name; ci = 0; 
  while exist(fullfile(Pdir,[variant '.json']),'file')
    ci = ci + 1; variant = sprintf('%s_%d',name,ci); 
  end
  cat_io_json(fullfile(Pdir,[variant '.json']),V); 
end
% =========================================================================
function V = normalizeFields(V)
%normalizeFields. Software and converter independent representation of the
%  protocol relevant fields of a sidecar or protocol, i.e., the same
%  protocol results in the same values for Siemens VE and XA software and
%  different dcm2niix versions:
%    ImageType     merged with ImageTypeText (XA) as sorted list of upper
%                  case entries without format entries (MOSAIC, NONE,
%                  FMRI), i.e., e.g. ND/DIS2D/DIS3D (distortion correction),
%                  M/P, or NORM (prescan normalize) remain
%    ScanOptions, ScanningSequence, SequenceVariant
%                  sorted entries without the external trigger EXT (that
%                  does not change the images)
%    SequenceName  also from PulseSequenceName (XA) without '*' and numbers,
%                  e.g. *tfl3d1_16ns (VE) and *tfl3d1_256ns (XA) as tfl3d1ns
%    PulseSequenceDetails    sequence without directory, e.g. tfl
%    MultibandAccelerationFactor   from the ImageComments (e.g. 'Unaliased
%                  MB8/PE3/LB') if not given (XA60 data of dcm2niix 2022)
%    NumberOfVolumes   4th image dimension (Dimensions, if the image exists)
%    NumberOfSlices    number of slices (SliceTiming), otherwise the 3rd
%                  image dimension (e.g. of 3D images)
%    SliceOrder    acquisition order of the slices (fMRI), i.e., the index
%                  of the acquisition time of each slice (SliceTiming in ms)
%                  that is independent of the TR and of the timing precision
%                  of the software, with equal indices for the simultaneously
%                  acquired slices of multiband
%    BValueShells  b-value shells of DWI (from the b-values of the bval file
%                  of the image, see cat_io_dcm2bids): clusters of similar
%                  b-values (gaps below 5% + 10 s/mm2, e.g. 1190-1205 of VE 
%                  or the b-values of XA that are rounded to 50) with at 
%                  least 2% of the volumes (i.e. without single directions of
%                  other b-values) as mean value rounded to 50 (i.e. b=5 of 
%                  the b0 images of VE is 0)
%  The protocols (getProtocols) are normalized too, i.e., they can use
%  both representations.

  % sorted lists
  if isfield(V,'ImageTypeText')
    if isfield(V,'ImageType')
      V.ImageType = [reshape(cellstr(V.ImageType),[],1); reshape(cellstr(V.ImageTypeText),[],1)];
    else
      V.ImageType = V.ImageTypeText;
    end
  end
  if isfield(V,'ImageType') && (iscellstr(V.ImageType) || ischar(V.ImageType))
    V.ImageType = setdiff(upper(reshape(cellstr(V.ImageType),[],1)),{'MOSAIC','NONE','FMRI'});
  end
  for fn = {'ScanOptions','ScanningSequence','SequenceVariant'}
    if isfield(V,fn{1}) && ischar(V.(fn{1}))
      V.(fn{1}) = strjoin(setdiff(upper(strsplit(V.(fn{1}),'\\')),{'','EXT'}),'\\'); % escaped backslash
    end
  end

  % sequence names
  if ~isfield(V,'SequenceName') && isfield(V,'PulseSequenceName')
    V.SequenceName = V.PulseSequenceName;
  end
  if isfield(V,'SequenceName') && ischar(V.SequenceName)
    V.SequenceName = regexprep(V.SequenceName,'\*|_\d+','');
  end
  if isfield(V,'PulseSequenceDetails') && ischar(V.PulseSequenceDetails)
    V.PulseSequenceDetails = lower(regexprep(V.PulseSequenceDetails,'^.*[\\/]',''));
  end

  % multiband factor of older dcm2niix versions for XA60
  if ~isfield(V,'MultibandAccelerationFactor') && isfield(V,'ImageComments') && ischar(V.ImageComments)
    mb = regexp(V.ImageComments,'MB(\d+)','tokens','once');
    if ~isempty(mb), V.MultibandAccelerationFactor = str2double(mb{1}); end
  end

  % number of volumes and slices (image fields without image, e.g. NaN, are missing)
  for fn = {'Dimensions','VoxelSize'}
    if isfield(V,fn{1}) && (~isnumeric(V.(fn{1})) || any(isnan(V.(fn{1})(:)))), V = rmfield(V,fn{1}); end
  end
  if isfield(V,'Dimensions') && numel(V.Dimensions) >= 3
    dm = [reshape(V.Dimensions,1,[]) 1];
    if ~isfield(V,'NumberOfVolumes'), V.NumberOfVolumes = dm(4); end
    if ~isfield(V,'NumberOfSlices'),  V.NumberOfSlices  = dm(3); end
  end
  if isfield(V,'SliceTiming') && isnumeric(V.SliceTiming) && numel(V.SliceTiming) > 1
    V.NumberOfSlices = numel(V.SliceTiming);
    [~,~,so] = unique(round(1000*double(V.SliceTiming(:)))); 
    V.SliceOrder = so'; 
  end

  % b-value shells of DWI 
  if isfield(V,'BValues') && isnumeric(V.BValues) && ~isempty(V.BValues)
    b  = sort(double(V.BValues(:))); 
    id = cumsum([1; diff(b) > 0.05*b(2:end) + 10]); 
    V.BValueShells = []; 
    for ci = 1:id(end)
      if sum(id==ci) >= max(1,0.02*numel(b))
        V.BValueShells(end+1) = 50*round(mean(b(id==ci))/50); 
      end
    end
  end
end
% =========================================================================
function [cnt,str,seqdiff,keydiff,nfields] = compareProtocol(V,P,tol,opts)
%compareProtocol. Mismatches of the (normalized) scan fields V with the
%  (normalized) protocol P: numeric fields count 0.5 if they differ by more
%  than tol percent (b-value shells by more than 10%), text fields 1, and 
%  lists the number of different entries (independent of their order). 
%  The key parameters opts.protocol.key (contrast and resolution, e.g. TE, 
%  TR, voxel size) count at least 1, the critical parameters 
%  opts.protocol.critical (e.g. phase-encoding direction, sequence) and the
%  image component of the ImageType (magnitude, phase, real, imaginary) 
%  count 10, i.e., such scans are nonconform (see protocolRating). The 
%  AcquisitionDuration counts 0.5 if it differs by more than 5%. str 
%  lists the mismatches as {field, scan value, protocol value}, also for 
%  the protocol fields that are not available for the scan (as 'n/a', 
%  without counting). seqdiff and keydiff indicate a different sequence 
%  (PulseSequenceDetails) or key or critical parameter for the selection of
%  the protocol, and nfields is the number of tested protocol fields.

  crit = opts.protocol.critical; key = opts.protocol.key; 
  cnt  = 0; str = cell(0,3); seqdiff = false; keydiff = false;
  FN   = fieldnames(P);
  FN   = FN(~strncmp(FN,'x_',2)); % comments and disabled fields (e.g. "_ImageType")
  nfields = numel(FN);
  for fni = 1:numel(FN)
    v = []; p = P.(FN{fni});
    if isfield(V,FN{fni}), v = V.(FN{fni}); end
    if isempty(v) && ~isempty(p)
      str(end+1,:) = {FN{fni}, 'n/a', val2str(p)}; %#ok<AGROW>
      continue
    end
    if strcmp(FN{fni},'BValueShells') && isnumeric(v) && isnumeric(p)
      % shells within 10% (e.g. the rounded b-values of XA)
      v = double(v(:)); p = double(p(:));
      d = 0.5 * ~(numel(v) == numel(p) && all(abs(v - p) <= 0.1*p + 25));
    elseif strcmp(FN{fni},'AcquisitionDuration') && isnumeric(v) && isnumeric(p)
      % duration within 5% (e.g. the scan time of dcm2niix or volumes times TR)
      d = 0.5 * ~(abs(double(v(1)) - double(p(1))) <= 0.05 * double(p(1)));
    elseif (isnumeric(v) || islogical(v)) && (isnumeric(p) || islogical(p))
      v = double(v(:)); p = double(p(:));
      d = 0.5 * ~(numel(v) == numel(p) && all(abs(v - p) <= tol/100 * abs(p) + eps));
    elseif (iscellstr(v) || ischar(v)) && (iscellstr(p) || ischar(p)) && (iscell(v) || iscell(p)) %#ok<ISCLSTR> JSON lists
      d = numel(setxor(upper(cellstr(v)),upper(cellstr(p))));
      if strcmp(FN{fni},'ImageType') && ~strcmp(imageComponent(v),imageComponent(p))
        d = 10; % critical: e.g. phase vs. magnitude image
      end
    elseif ischar(v) && ischar(p)
      d = double(~strcmp(v,p));
    else
      d = 1;
    end
    if d > 0
      if any(strcmp(FN{fni},crit))
        d = 10;           % critical parameters (nonconform)
      elseif any(strcmp(FN{fni},key))
        d = max(1,d);     % contrast and resolution
      end
      cnt = cnt + d;
      str(end+1,:) = {FN{fni}, val2str(V.(FN{fni})), val2str(P.(FN{fni}))}; %#ok<AGROW>
      seqdiff = seqdiff || strcmp(FN{fni},'PulseSequenceDetails');
      keydiff = keydiff || any(strcmp(FN{fni},[crit key])) || d >= 10;
    end
  end
end
% =========================================================================
function c = imageComponent(it)
%imageComponent. Image component of the ImageType entries: phase (P or
%  PHASE), real, imaginary, otherwise magnitude.
  it = upper(cellstr(it));
  if any(ismember(it,{'P','PHASE'}))
    c = 'phase';
  elseif any(strcmp(it,'REAL'))
    c = 'real';
  elseif any(strcmp(it,'IMAGINARY'))
    c = 'imaginary';
  else
    c = 'magnitude';
  end
end
% =========================================================================
function d = nameDirection(str)
%nameDirection. Phase-encoding direction token (AP, PA, LR, RL, SI, IS) of 
%  a series or protocol name (e.g. PA of "SpinEchoFieldMap_PA"), otherwise
%  empty. 
  t = regexpi(str,'(?:^|[_ -])(AP|PA|LR|RL|SI|IS)(?=$|[_ -])','tokens','once');
  if isempty(t), d = ''; else, d = upper(t{1}); end
end
% =========================================================================
function s = val2str(v)
%val2str. Value as text for the mismatch tables.
  if ischar(v)
    s = v;
  elseif iscell(v)
    s = strjoin(cellfun(@val2str,reshape(v,1,[]),'UniformOutput',false),' ');
  elseif islogical(v)
    s = strjoin(cellstr(char(string(reshape(v,1,[])))),' ');
  elseif isnumeric(v)
    s = strjoin(arrayfun(@(x) sprintf('%g',x),reshape(v,1,[]),'UniformOutput',false),' ');
  else
    s = sprintf('%dx%d %s',size(v,1),size(v,2),class(v));
  end
end
% =========================================================================
function Vjson = cleanupVjson(Vjson,type)
% remove Patient specific entries. type={'basic'*|'protocol'}.

  if ~exist('type','var'), type = 'basic'; end

  FN    = fieldnames(Vjson); 
  Vjson = rmfield(Vjson, FN(cat_io_contains(FN,'Patient')));  
  
  % maybe 2-3 levels of data reduction ?
  % (0) all DCM2NII, (1) light cleanup, (2) severe cleanup

  switch type
    case 'basic'
      % negative list, i.e. remove these entries
      RFn = {
        'SeriesInstanceUID'; 'StudyInstanceUID'; 'StudyID'; 
        'InstitutionName'; 'InstitutionalDepartmentName'; 'InstitutionAddress'; 
        'StationName'; 
        ... 'ProcedureStepDescription'; 
        'BodyPartExamined';
        'AcquisitionTime'; 
        ... 'AcquisitionDateTime'
        'ShimSetting'; 'TxRefAmp'; 
        }; 
    case {'protocol','bids'}
      % positive list, i.e. keep only these entries
      % this should be for an overview of various files
      RFn = {
        'Modality'; 
        'MagneticFieldStrength'; 
        'Manufacturer';
        'PatientPosition';
        'MRAcquisitionType';
        'ScanningSequence';
        'SequenceVariant';
        'ScanOptions';
        'ImageType';
        'NonlinearGradientCorrection';
        'SliceThickness';
        'EchoTime';
        'RepetitionTime';
        'SpoilingState';
        'InversionTime';
        'FlipAngle';
        'PartialFourier';
        'BaseResolution';
        'DiffusionScheme';
        'PhaseResolution';
        'CoilString';
        'RefLinesPE';
        'CoilCombinationMethod';
        'MatrixCoilMode';
        'PercentPhaseFOV';
        'PercentSampling';
        'PhaseEncodingSteps';
        'AcquisitionMatrixPE';
        'DwellTime';
        'ReconMatrixPE';
        'InPlanePhaseEncodingDirectionDICOM';
        'DerivedVendorReportedEchoSpacing';
        ... 
        ... 'PulseSequenceName'; % not always and operator specific
        ... 'ProcedureStepDescription'; % free text of the site (e.g. names)
        'PhaseEncodingDirection';
        ...'SliceTiming';
        'InPlanePhaseEncodingDirectionDICOM';
        'PixelBandwidth';
        'TotalReadoutTime';
        'EffectiveEchoSpacing';
        'DerivedVendorReportedEchoSpacing';
        'ParallelReductionFactorInPlane';
        'BandwidthPerPixelPhaseEncode';
        'EchoTrainLength';
        'MultibandAccelerationFactor';
        }; 
      if strcmp(type,'bids')
        % further technical (not identifying) fields of the BIDS sidecars, 
        % e.g. required by BIDS (TaskName for func, EchoTime1/2 for phasediff,
        % ASL and qMRI parameters)
        RFn = [RFn; {
          'TaskName'; 'SliceTiming'; 'SliceEncodingDirection'; 'SpacingBetweenSlices'; 
          'EchoTime1'; 'EchoTime2'; 'RepetitionTimeExcitation'; 'RepetitionTimePreparation'; 
          'NumberOfVolumesDiscardedByScanner'; 'NumberOfVolumesDiscardedByUser'; 'DelayTime'; 
          'ManufacturersModelName'; 'SoftwareVersions'; 'ReceiveCoilName'; 'PulseSequenceType'; 
          'SequenceName'; 'B0FieldIdentifier'; 'B0FieldSource'; 'MTState'; 'NumberShots'; 
          'ArterialSpinLabelingType'; 'PostLabelingDelay'; 'BackgroundSuppression'; 'M0Type'; 
          'TotalAcquiredPairs'; 'LabelingDuration'; 'PCASLType'; 'VascularCrushing'; 
          'AcquisitionVoxelSize'; 'BolusCutOffFlag'; 'LookLocker'; 'Units'; 
          'ConversionWarnings'; ... warnings of dcm2niix without paths and file names (see conversionWarnings in cat_io_dcm2bids_importDCM)
          'IntendedFor'}]; % remapped to the new file names (see remapIntendedFor in cat_io_dcm2bids_bids)
      end
      RFn = setdiff( fieldnames(Vjson), RFn );   
  end
  RFn   = intersect( RFn , fieldnames(Vjson) ); 
  Vjson = rmfield(Vjson,RFn); 
end
% =========================================================================
function updateDBtables(Vjson, Ptbldir, sites, opts) 
%updateDBtables. Update the subject table of the database (subjects.csv in
%  catDCM2BIDSdb-tables, one row per patient) for a scan with its JSON 
%  sidecar Vjson. The tables of the scans, sites, and protocols are created 
%  at the end of each run (see cat_io_dcm2bids_tables). 
  site    = setupSites(Vjson,sites);
  Thdr    = {'PatientID','PatientName', 'PatientSex','PatientBirthDate','SiteID'}; 
  Tnewrow = {Vjson.PatientID, Vjson.PatientName, ...
             Vjson.PatientSex, Vjson.PatientBirthDate, ...
             site}; 
  cat_io_dcm2bids_helper('updateTable', ...
    spm_file(fullfile(Ptbldir,'subjects'),'ext',opts.tableformat),Thdr,Tnewrow,1,opts.rerun);
end
% =========================================================================
function str = getFileString(str)
  str = cat_io_strrep(str, ...
    {'ä','ü','ö',' '}, {'ae','ue','oe','_'});
  strd = double(str);
  str(strd<48 & strd~=45) = '_';
  str(strd>57 & strd<65)  = '_';
  str(strd>90 & strd<97 & strd~=95)  = '_';
  str(strd>122) = 'X';
end
