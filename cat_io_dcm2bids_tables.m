function varargout = cat_io_dcm2bids_tables(action,varargin)
%cat_io_dcm2bids_tables. Tables of the database of cat_io_dcm2bids.
%  The internal tables in catDCM2BIDSdb-tables give an overview of all
%  scans of the database that were processed (in this or previous runs)
%  and contain personal data (e.g. PatientID, dates), i.e., they are not
%  for sharing (see cat_io_dcm2bids_report for the report tables):
%    scans.<fmt>       one row per scan with the site, session, BIDS type,
%                      protocol test (closest protocol, status, PFM, PFR,
%                      deviations), protocol variant, quality rating (SQR),
%                      acquisition duration (s, see acquisitionDuration in
%                      cat_io_dcm2bids_protocols), and scanner fields (the
%                      cumulative list is kept in scans.mat)
%    sites.<fmt>       one row per scanner with its identity, software and
%                      coils, scan period, number of subjects, sessions, and
%                      scans (also of the typical scan types opts.scantypes),
%                      protocol status, and median ratings
%    sites.json        template for the site dictionary (protocols/sites.json)
%                      with the abbreviation of the dictionary or a proposal
%    protocols.<fmt>   one row per protocol variant (see protocolVariant in
%                      cat_io_dcm2bids_protocols) with the sites using it
%                      and its main parameters (with the acquisition 
%                      duration)
%    protocols_sites.<fmt>  one row per site and protocol variant
%    protocols-details/<datatype>/<protocol>/<site>-deviations.<fmt>
%                      deviations of conform and semiconform scans of a
%                      site from the (closest) protocol
%    protocols/<datatype>/<variant>.json
%                      protocol definition of each variant (see 
%                      protocolVariant, variants without scans are removed)
%    subjects.<fmt>    one row per patient (see updateDBtables in
%                      cat_io_dcm2bids_protocols)
%  where <fmt> is opts.tableformat (csv or tsv). The tables (besides the
%  subjects and protocol variants) are created again at the end of each run
%  from the cumulative scan list.
%
%  varargout = cat_io_dcm2bids_tables(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_tables('initTables',Ptbldir)
%    R = cat_io_dcm2bids_tables('scanRecord',Vjson,scanid,site,session,datatype,suffix,ent,variant,PT,BIDSfile)
%    cat_io_dcm2bids_tables('writeTables',Ptbldir,R,sites,opts)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_protocols, cat_io_dcm2bids_report.

  switch action
    case {'initTables','scanRecord','writeTables'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_tables:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function initTables(Ptbldir)
%initTables. Remove the tables of older versions (studies, protocol lists
%  by DICOM protocol name, and protocol details by status) once, i.e., if
%  the cumulative scan list (scans.mat) does not exist yet.
  if exist(fullfile(Ptbldir,'scans.mat'),'file') || ~exist(Ptbldir,'dir'), return; end
  for f = {'studies.csv','studies.tsv','protocols.csv','protocols.tsv'}
    if exist(fullfile(Ptbldir,f{1}),'file'), delete(fullfile(Ptbldir,f{1})); end
  end
  for d = {'study_subject_lists','protocols','protocols-details'}
    if exist(fullfile(Ptbldir,d{1}),'dir'), rmdir(fullfile(Ptbldir,d{1}),'s'); end
  end
end
% =========================================================================
function R = scanRecord(Vjson,scanid,site,session,datatype,suffix,ent,variant,PT,BIDSfile)
%scanRecord. Record of a scan for the tables of the database: its identity
%  (scan ID, patient, site, session of the database, date, series), BIDS
%  type and file, protocol test (see testProtocols in
%  cat_io_dcm2bids_protocols) and variant, and scanner fields of its
%  sidecar Vjson. The quality rating SQR is added after the QC.
  str = @(fn) fieldString(Vjson,fn);
  num = @(fn) fieldNumber(Vjson,fn);
  R = struct('ScanID',scanid, 'PatientID',str('PatientID'), 'Site',site, ...
    'DeviceSerialNumber',str('DeviceSerialNumber'), 'Session',session, ...
    'AcquisitionDateTime',str('AcquisitionDateTime'), 'SeriesNumber',num('SeriesNumber'), ...
    'ProtocolName',str('ProtocolName'), 'SeriesDescription',str('SeriesDescription'), ...
    'Datatype',datatype, 'Suffix',suffix, 'Entities',ent, 'BIDSfile',BIDSfile, ...
    'Variant',variant, 'Protocol',PT.protocol, 'Status',PT.status, 'PFM',PT.PFM, 'PFR',PT.PFR, ...
    'Deviations',{PT.deviations}, 'SQR',nan, ...
    'Manufacturer',str('Manufacturer'), 'ManufacturersModelName',str('ManufacturersModelName'), ...
    'MagneticFieldStrength',num('MagneticFieldStrength'), 'SoftwareVersions',str('SoftwareVersions'), ...
    'ReceiveCoilName',str('ReceiveCoilName'), 'ReceiveCoilActiveElements',str('ReceiveCoilActiveElements'), ...
    'AcquisitionDuration',cat_io_dcm2bids_protocols('acquisitionDuration',Vjson), ...
    'TxRefAmp',num('TxRefAmp'), 'ProcedureStepDescription',str('ProcedureStepDescription'), ...
    'InstitutionName',str('InstitutionName'), 'InstitutionalDepartmentName',str('InstitutionalDepartmentName'), ...
    'InstitutionAddress',str('InstitutionAddress'), 'StationName',str('StationName'));
end
% =========================================================================
function writeTables(Ptbldir,R,sites,opts)
%writeTables. Add the scan records R of this run to the cumulative scan
%  list (scans.mat, records of the same scan are replaced) and create the
%  scan, site, and protocol tables of the database (see the help of this
%  file).
  if isempty(R), return; end
  Pmat = fullfile(Ptbldir,'scans.mat');
  if exist(Pmat,'file') % (records of older versions with other fields)
    S = load(Pmat,'R'); S = S.R;
    if ~isempty(S), R = cat_io_dcm2bids_helper('catRecords',S(~ismember({S.ScanID},{R.ScanID})),R); end
  end
  [~,so] = sort(strcat({R.Session},'|',cellfun(@(x) sprintf('%06d',x),{R.SeriesNumber},'UniformOutput',false),'|',{R.ScanID}));
  R = R(so);
  save(Pmat,'R');
  ext = ['.' opts.tableformat];

  writeScans(fullfile(Ptbldir,['scans' ext]),R,opts);
  writeSites(Ptbldir,R,sites,opts);
  writeProtocols(Ptbldir,R,opts);
end
% =========================================================================
function writeScans(Ptable,R,opts)
%writeScans. Table with one row per scan of the cumulative scan list.
  FN = setdiff(fieldnames(R),{'Deviations'},'stable');
  T  = [FN' {'Deviations'}];
  for ri = 1:numel(R)
    row = cellfun(@(fn) R(ri).(fn),FN','UniformOutput',false);
    T(end+1,:) = [row {deviationString(R(ri))}]; %#ok<AGROW>
  end
  cat_io_dcm2bids_helper('writeTable',Ptable,T,opts.tableformat);
end
% =========================================================================
function writeSites(Ptbldir,R,sites,opts)
%writeSites. Site table with one row per scanner (DeviceSerialNumber) and
%  the template of the site dictionary (sites.json, see getSites in
%  cat_io_dcm2bids_protocols) with the abbreviation of the dictionary (or
%  of the site registry of the database without dictionary, see 
%  siteRegistry in cat_io_dcm2bids_protocols) or a proposal for new 
%  scanners (city, manufacturer, and field strength, e.g. MagSI30T), and an
%  empty Color of the site in the report figure (e.g. #A82117, see 
%  siteColors in cat_io_dcm2bids_protocols).
  dev = {R.DeviceSerialNumber}; dev(cellfun('isempty',dev)) = {R(cellfun('isempty',dev)).Site};
  [udev,~,di] = unique(dev,'stable');
  [types,tnam] = scanTypes(opts);
  rtype = arrayfun(@recordType,R,'UniformOutput',false);
  T = [{'Site','DeviceSerialNumber','InstitutionName','InstitutionalDepartmentName','InstitutionAddress', ...
    'StationName','Manufacturer','ManufacturersModelName','MagneticFieldStrength','SoftwareVersions', ...
    'ReceiveCoilName','FirstScan','LastScan','Subjects','Sessions','Scans'}, ...
    reshape([strcat(tnam,'_scans'); strcat(tnam,'_subjects')],1,[]), ...
    {'conform','semiconform','nonconform','PFR','SQR','TxRefAmp','ProcedureStepDescription'}];
  J = struct('DeviceSerialNumber',{},'InstitutionAbbreviation',{},'Color',{},'InstitutionName',{}, ...
    'InstitutionalDepartmentName',{},'InstitutionAddress',{},'MagneticFieldStrength',{}, ...
    'Manufacturer',{},'ManufacturersModelName',{},'StationName',{});
  for si = 1:numel(udev)
    Rs   = R(di==si);
    pid  = {Rs.PatientID};
    date = sort({Rs.AcquisitionDateTime}); date = date(~cellfun('isempty',date));
    if isempty(date), date = {'',''}; end
    row  = {Rs(1).Site, udev{si}, mostFrequent({Rs.InstitutionName}), mostFrequent({Rs.InstitutionalDepartmentName}), ...
      mostFrequent({Rs.InstitutionAddress}), mostFrequent({Rs.StationName}), mostFrequent({Rs.Manufacturer}), ...
      mostFrequent({Rs.ManufacturersModelName}), median([Rs.MagneticFieldStrength],'omitnan'), ...
      valueCounts({Rs.SoftwareVersions}), valueCounts({Rs.ReceiveCoilName}), ...
      date{1}(1:min(10,end)), date{end}(1:min(10,end)), numel(unique(pid)), numel(unique({Rs.Session})), numel(Rs)};
    for ti = 1:numel(types)
      it  = strcmp(rtype(di==si),types{ti});
      row = [row {sum(it), numel(unique(pid(it)))}]; %#ok<AGROW>
    end
    row = [row {sum(strcmp({Rs.Status},'conform')), sum(strcmp({Rs.Status},'semiconform')), ...
      sum(strcmp({Rs.Status},'nonconform')), median([Rs.PFR],'omitnan'), median([Rs.SQR],'omitnan'), ...
      median([Rs.TxRefAmp],'omitnan'), valueCounts({Rs.ProcedureStepDescription})}]; %#ok<AGROW>
    T(end+1,:) = row; %#ok<AGROW>

    % site dictionary entry (abbreviation of the dictionary or proposal)
    abbr = '';
    if ~isempty(sites)
      id = find(cat_io_contains(sites(:,1),udev{si}),1); % as setupSites in cat_io_dcm2bids_protocols
      if ~isempty(id), abbr = char(sites{id,2}); end
    end
    if isempty(abbr), abbr = cat_io_dcm2bids_protocols('siteProposal',row{5},row{9},row{7},opts); end
    J(end+1) = struct('DeviceSerialNumber',udev{si},'InstitutionAbbreviation',abbr,'Color','', ...
      'InstitutionName',row{3},'InstitutionalDepartmentName',row{4},'InstitutionAddress',row{5}, ...
      'MagneticFieldStrength',row{9},'Manufacturer',row{7},'ManufacturersModelName',row{8}, ...
      'StationName',row{6}); %#ok<AGROW>
  end
  cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['sites.' opts.tableformat]),T,opts.tableformat);
  cat_io_dcm2bids_helper('writeJSON',fullfile(Ptbldir,'sites.json'),J);
end
% =========================================================================
function writeProtocols(Ptbldir,R,opts)
%writeProtocols. Tables of the protocol variants (protocols.<fmt>, one row
%  per variant with the sites using it, and protocols_sites.<fmt>, one row
%  per site and variant) and the deviations of conform and semiconform
%  scans of each site from the (closest) protocol (protocols-details).
  R = R(~cellfun('isempty',{R.Variant}));

  % remove the variants that are no longer used by any scan (e.g. after the
  % protocol test with other protocols)
  Pv = cat_vol_findfiles(fullfile(Ptbldir,'protocols'),'*.json');
  if ~isempty(Pv)
    used = strcat({R.Datatype},'/',{R.Variant});
    vnam = strcat(spm_file(spm_file(Pv,'path'),'basename'),'/',spm_file(Pv,'basename'));
    for vi = reshape(find(~ismember(vnam,used)),1,[]), delete(Pv{vi}); end
  end
  if isempty(R), return; end
  ext  = ['.' opts.tableformat];
  par  = {'MRAcquisitionType','PulseSequenceDetails','EchoTime','RepetitionTime','InversionTime', ...
    'FlipAngle','VoxelSize','NumberOfSlices','NumberOfVolumes','MultibandAccelerationFactor', ...
    'PhaseEncodingDirection','AcquisitionDuration'};
  hdr  = {'Variant','Datatype','Type','ProtocolName','Protocol','Status','PFR','Deviations', ...
    'Scans','Sessions','Subjects','FirstScan','LastScan'};

  % variants (also of all sites) with their parameters
  [uvar,~,vi] = unique(strcat({R.Datatype},'/',{R.Variant}),'stable');
  T  = [hdr(1) {'Sites','SiteList'} hdr(2:end) par];
  Ts = [{'Site'} hdr par];
  for ui = 1:numel(uvar)
    Rv  = R(vi==ui);
    Pv  = fullfile(Ptbldir,'protocols',Rv(1).Datatype,[Rv(1).Variant '.json']);
    if exist(Pv,'file'), V = cat_io_json(Pv); else, V = struct(); end
    pv  = cellfun(@(fn) fieldValue(V,fn),par,'UniformOutput',false);
    T(end+1,:) = [{Rv(1).Variant, numel(unique({Rv.Site})), valueCounts({Rv.Site})}, variantRow(Rv), pv]; %#ok<AGROW>
    [usite,~,si] = unique({Rv.Site},'stable');
    for sj = 1:numel(usite)
      Ts(end+1,:) = [usite(sj) {Rv(1).Variant} variantRow(Rv(si==sj)) pv]; %#ok<AGROW>
    end
  end
  cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['protocols' ext]),T,opts.tableformat);
  cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['protocols_sites' ext]),Ts,opts.tableformat);

  % deviations of each site from the (closest) protocol
  Pdet = fullfile(Ptbldir,'protocols-details');
  if exist(Pdet,'dir'), rmdir(Pdet,'s'); end
  R = R(ismember({R.Status},{'conform','semiconform'}) & ~cellfun('isempty',{R.Protocol}));
  [ukey,~,ki] = unique(strcat({R.Datatype},'|',{R.Protocol},'|',{R.Site}),'stable');
  for ui = 1:numel(ukey)
    Rk  = R(ki==ui);
    D   = vertcat(Rk.Deviations);
    if isempty(D), continue; end
    T   = {'Field','ProtocolValue','SiteValues','DeviatingScans','Scans'};
    [ufn,~,fi] = unique(D(:,1),'stable');
    for fj = 1:numel(ufn)
      T(end+1,:) = {ufn{fj}, mostFrequent(D(fi==fj,3)), valueCounts(D(fi==fj,2)), sum(fi==fj), numel(Rk)}; %#ok<AGROW>
    end
    cat_io_dcm2bids_helper('writeTable',fullfile(Pdet,Rk(1).Datatype,Rk(1).Protocol, ...
      [Rk(1).Site '-deviations' ext]),T,opts.tableformat);
  end
end
% =========================================================================
function row = variantRow(Rv)
%variantRow. Columns of the scans Rv of a protocol variant (Datatype to
%  LastScan of the header in writeProtocols).
  date = sort({Rv.AcquisitionDateTime}); date = date(~cellfun('isempty',date));
  if isempty(date), date = {'',''}; end
  row = {Rv(1).Datatype, [Rv(1).Suffix Rv(1).Entities], mostFrequent({Rv.ProtocolName}), ...
    mostFrequent({Rv.Protocol}), mostFrequent({Rv.Status}), median([Rv.PFR],'omitnan'), ...
    mostFrequent(arrayfun(@deviationString,Rv,'UniformOutput',false)), ...
    numel(Rv), numel(unique({Rv.Session})), numel(unique(strcat({Rv.Site},'|',{Rv.PatientID}))), ...
    date{1}(1:min(10,end)), date{end}(1:min(10,end))};
end
% =========================================================================
function [types,names] = scanTypes(opts)
%scanTypes. Typical scan types (datatype/type, e.g. anat/T1w, see
%  opts.scantypes) and their column names (the type, or the datatype and
%  type if the type is not unique, e.g. func_sbref and dwi_sbref).
  types = cellstr(opts.scantypes);
  names = regexprep(types,'^[^/]*/','');
  dup   = cellfun(@(n) sum(strcmp(names,n)) > 1, names);
  names(dup) = strrep(types(dup),'/','_');
  names = regexprep(names,'[^a-zA-Z0-9_]','_');
end
% =========================================================================
function t = recordType(r)
%recordType. Scan type of a record (datatype/suffix with the image part and
%  inversion, as the image types of cat_io_dcm2bids_render).
  tok = regexp(r.Entities,'_(part-(?:phase|real|imag)|inv-\d+)(?=_|$)','tokens');
  t   = [r.Datatype '/' r.Suffix sprintf('_%s',strjoin(cellfun(@(x) x{1},tok,'UniformOutput',false),'_'))];
  if isempty(tok), t = [r.Datatype '/' r.Suffix]; end
end
% =========================================================================
function s = deviationString(r)
%deviationString. Deviations {field, scan value, protocol value} of the
%  record r from its (closest) protocol as text, e.g. "EchoTime: 0.03
%  (0.037); ImageType: DIS2D M (M ND)", 'none' without deviations, and
%  empty without protocol test.
  D = r.Deviations;
  if isempty(r.Protocol)
    s = '';
  elseif isempty(D)
    s = 'none';
  else
    s = strjoin(cellfun(@(f,v,p) sprintf('%s: %s (%s)',f,v,p),D(:,1)',D(:,2)',D(:,3)', ...
      'UniformOutput',false),'; ');
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
%  (most frequent first), e.g. "syngo MR XA60 (120); syngo MR E11 (3)".
  c = cellfun(@char,c,'UniformOutput',false);
  c = c(~cellfun('isempty',c));
  if isempty(c), v = ''; return; end
  [u,~,ui] = unique(c);
  n = accumarray(ui(:),1);
  [n,so] = sort(n,'descend'); u = u(so);
  v = strjoin(arrayfun(@(i) sprintf('%s (%d)',u{i},n(i)),1:numel(u),'UniformOutput',false),'; ');
end
% =========================================================================
function v = fieldValue(V,fn)
%fieldValue. Value of the field fn of V (empty if not available).
  if isfield(V,fn), v = V.(fn); else, v = []; end
end
% =========================================================================
function s = fieldString(V,fn)
%fieldString. Field fn of a sidecar as text (empty if not available).
  s = '';
  if isfield(V,fn) && ~isempty(V.(fn))
    v = V.(fn);
    if iscell(v), v = strjoin(cellfun(@char,v(:)','UniformOutput',false),' '); end
    if isnumeric(v), v = strjoin(arrayfun(@(x) sprintf('%g',x),v(:)','UniformOutput',false),' '); end
    s = char(v);
  end
end
% =========================================================================
function x = fieldNumber(V,fn)
%fieldNumber. Field fn of a sidecar as number (NaN if not available).
  x = nan;
  if isfield(V,fn) && isnumeric(V.(fn)) && ~isempty(V.(fn)), x = double(V.(fn)(1)); end
  if isfield(V,fn) && ischar(V.(fn)), x = str2double(V.(fn)); end
end
