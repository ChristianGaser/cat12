function varargout = cat_io_dcm2bids_bids(action,varargin)
%cat_io_dcm2bids_bids. BIDS naming and output files of cat_io_dcm2bids.
%  Naming: BIDS conform labels (bidsLabel) and the datatype, suffix, and 
%  task of a scan by its protocol name (setupDatatype, setupSuffix, 
%  setupTask). BIDSoutputs lists the datatypes and suffixes of the output 
%  structure (batch dependencies), getBIDSoutputs collects the exported 
%  images. 
%  Output files: sessions.tsv of each subject (writeSubjectTSV, with SOR), 
%  participants.tsv and the private participants table (writeParticipantTSV,
%  writePrivateTSV), the IntendedFor references of BIDS inputs 
%  (remapIntendedFor), dataset_description.json, README.md, and 
%  participants.json (writeDatasetFiles with the information of the imported
%  BIDS datasets, see getBIDSsources), and the gzip status of the output 
%  (gzipunzipOutputdata).
%
%  varargout = cat_io_dcm2bids_bids(action,varargin)
%
%  Actions (see the help of the local functions):
%    label = cat_io_dcm2bids_bids('bidsLabel',str,camelCase)
%    datatype = cat_io_dcm2bids_bids('setupDatatype',pro)
%    [suffix,pro] = cat_io_dcm2bids_bids('setupSuffix',datatype,pro,ser)
%    task = cat_io_dcm2bids_bids('setupTask',datatype,pro)
%    out = cat_io_dcm2bids_bids('BIDSoutputs')
%    out = cat_io_dcm2bids_bids('getBIDSoutputs',PBIDS,out)
%    cat_io_dcm2bids_bids('remapIntendedFor',intended,srcmap)
%    cat_io_dcm2bids_bids('writeSubjectTSV',Vjson,sub,ses,Poutdir,BIDSsubdir,anon,SOR)
%    cat_io_dcm2bids_bids('writeParticipantTSV',P,Pparticipants,anon)
%    cat_io_dcm2bids_bids('writePrivateTSV',P,Pprivate)
%    sources = cat_io_dcm2bids_bids('getBIDSsources',roots)
%    cat_io_dcm2bids_bids('writeDatasetFiles',Pbids,job,sources,convver)
%    cat_io_dcm2bids_bids('gzipunzipOutputdata',Poutdir,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'bidsLabel','setupDatatype','setupSuffix','setupTask','BIDSoutputs', ...
          'getBIDSoutputs','remapIntendedFor','writeSubjectTSV', ...
          'writeParticipantTSV','writePrivateTSV','getBIDSsources', ...
          'writeDatasetFiles','gzipunzipOutputdata'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_bids:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function label = bidsLabel(str,camelCase)
%bidsLabel. Assure BIDS conform labels.
  if ~exist('camelCase','var'), camelCase = 0; end
  str = char(str);
  
  % convert german characters
  map = { char(228),'ae'; char(246),'oe'; char(252),'ue'; ...   % ä ö ü
          char(196),'Ae'; char(214),'Oe'; char(220),'Ue'; ...   % Ä Ö Ü
          char(223),'ss' };                                      % ß
  for k = 1:size(map,1)
      str = strrep(str, map{k,1}, map{k,2});
  end
  
  % other chars
  if usejava('jvm')
    form = javaMethod('valueOf', 'java.text.Normalizer$Form', 'NFD');
    str  = char(java.text.Normalizer.normalize(java.lang.String(str), form));
    % remove combining diacritical marks (U+0300 to U+036F)
    c = double(str);
    str = char(c(~(c >= 768 & c <= 879)));
  end
  
  % 
  if camelCase
    str = regexprep(str, '[\s_-]+(\w)', '${upper($1)}');
  end

  % remove unsupported characters
  label = regexprep(str, '[^a-zA-Z0-9]', '');
  
  if isempty(label) & ~isempty(str)
    error('bidsLabel:empty', 'No BIDS-conform characters in label "%s".', str);
  end
end
% =========================================================================
function datatype = setupDatatype(pro)
% first level
  if cat_io_contains( lower(pro) , {'mpr','mp2r','t1w','t2w','pdw','flair','inv','uni','t1','t2','tse'} ) && ...
    ~cat_io_contains( lower(pro) , {'fmri','rest','bold','func','rmri','dti','dwi','diff','dmri','field',} )
    datatype  = 'anat'; 
  elseif cat_io_contains( lower(pro) , {'fmri','rest','bold','func','rmri','rs','tb'} )
    datatype  = 'func'; 
  elseif cat_io_contains( lower(pro) , {'dti','dwi','diff','dmri'} )
    datatype  = 'dwi';
  elseif cat_io_contains( lower(pro) , {'mrs'} )
    datatype  = 'mrs';
  elseif cat_io_contains( lower(pro) , {'field','fieldmap','mag','b0','b1','phase'} )
    datatype  = 'fmap';
  else
    % unknown protocols
    cat_io_cprintf('red', sprintf('\n  Unknown BIDS datatype (anat/func/...) for protocol "%s"\n', lower(pro)) ) 
    datatype = 'other';
  end
end
% =========================================================================
function [suffix,pro] = setupSuffix(datatype,pro,ser)
  suffix = '';
  if strcmp( datatype , 'anat')
    if cat_io_contains( lower(pro) , {'t1map'})
      suffix = 'T1map'; 
    elseif cat_io_contains( lower(pro) , {'t2map'})
      suffix = 'T2map'; 
    elseif cat_io_contains( lower(pro) , {'t2starmap'})
      suffix = 'T2starmap'; 
    elseif cat_io_contains( lower(pro) , {'pdt2'})
      suffix = 'PDT2'; 
    elseif cat_io_contains( lower(pro) , {'mtr'})
      suffix = 'MTR'; 
    elseif cat_io_contains( lower(pro) , {'mt'})
      suffix = 'MT'; 
    elseif cat_io_contains( lower(pro) , {'mpr','mp2r','t1'})
      suffix = 'T1w'; 
    elseif cat_io_contains( lower(pro) , {'t2'}) && cat_io_contains( lower(pro) , {'star'})
      suffix = 'T2star';
    elseif cat_io_contains( lower(pro) , {'t2'})
      suffix = 'T2w';
    elseif cat_io_contains( lower(pro) , {'pd'})
      suffix = 'PDw';
    elseif cat_io_contains( lower(pro) , {'flair'})
      suffix = 'FLAIR';
    end
  
  elseif strcmp( datatype , 'dwi')
    if cat_io_contains( lower(pro) , {'sbref'}) || cat_io_contains( lower(ser) , {'sbref'})
      suffix = 'sbref';
    else
      suffix = 'dwi';
    end
    %pro = cat_io_strrep(pro,{'sbref'},{''});
  
  elseif strcmp( datatype , 'fmap')
    if cat_io_contains( lower(pro) , {'fieldmap'})
      suffix = 'fieldmap'; 
    elseif cat_io_contains( lower(pro) , {'phase'})
      suffix = 'phase'; 
    else
      suffix = 'epi'; 
    end
  
  elseif strcmp( datatype , 'func')  
    if cat_io_contains( lower(pro) , {'sbref'}) || cat_io_contains( lower(ser) , {'sbref'})
      suffix = 'sbref'; 
    else
      suffix = 'bold'; 
    end
  end
end
% =========================================================================
function task = setupTask(datatype,pro)
  if strcmp( datatype , 'func')
    if cat_io_contains( lower(pro) , {'rs','rest'})
      task = sprintf('_task-rest');
    elseif cat_io_contains( lower(pro) , {'motor'})
      task = sprintf('_task-motor');
    elseif cat_io_contains( lower(pro) , {'lang'})
      task = sprintf('_task-language');
    elseif cat_io_contains( lower(pro) , {'stroop'})
      task = sprintf('_task-stroop');
    elseif cat_io_contains( lower(pro) , {'nback'})
      task = sprintf('_task-nback');
    elseif cat_io_contains( lower(pro) , {'memory'})
      task = sprintf('_task-memory');
    elseif cat_io_contains( lower(pro) , {'faces'})
      task = sprintf('_task-faces');
    elseif cat_io_contains( lower(pro) , {'reward'})
      task = sprintf('_task-reward');
    elseif cat_io_contains( lower(pro) , {'gambling'})
      task = sprintf('_task-gambling');
    elseif cat_io_contains( lower(pro) , {'oddball'})
      task = sprintf('_task-oddball');
    elseif cat_io_contains( lower(pro) , {'attention'})
      task = sprintf('_task-attention');
    elseif cat_io_contains( lower(pro) , {'inhibition'})
      task = sprintf('_task-inhibition');
    elseif cat_io_contains( lower(pro) , {'go'})
      task = sprintf('_task-gonogo');
    elseif cat_io_contains( lower(pro) , {'flanker'})
      task = sprintf('_task-flanker');
    elseif cat_io_contains( lower(pro) , {'social'})
      task = sprintf('_task-social');
    elseif cat_io_contains( lower(pro) , {'pain'})
      task = sprintf('_task-pain');
    elseif cat_io_contains( lower(pro) , {'auditory'})
      task = sprintf('_task-auditory');
    elseif cat_io_contains( lower(pro) , {'visual','picture'})
      task = sprintf('_task-visual');
    elseif cat_io_contains( lower(pro) , {'movie'})
      task = sprintf('_task-movie');
    elseif cat_io_contains( lower(pro) , {'audio'})
      task = sprintf('_task-audio');
    elseif cat_io_contains( lower(pro) , {'working'})
      task = sprintf('_task-workingmemory');
    else
      task = sprintf('_task-other');
    end
  else
    task = '';
  end
end
% =========================================================================
function out = BIDSoutputs
%BIDSoutputs. Fixed list of BIDS datatypes and (lower case) suffixes that is
%  used for the output structure and the batch dependencies (defined in
%  cat_io_dcm2bids_defaults). It covers the suffixes defined by
%  setupDatatype and setupSuffix.

  def  = cat_io_dcm2bids_defaults;
  list = def.BIDSoutputs;
  out  = struct();
  for di = 1:size(list,1)
    for si = 1:numel(list{di,2})
      out.(list{di,1}).(list{di,2}{si}) = {};
    end
  end
end
% =========================================================================
function out = getBIDSoutputs(PBIDS,out)
%getBIDSoutputs. Add the raw BIDS images (all sessions, also from protocol 
%  and mismatch subdirectories but no derivatives) to out.datatype.suffix. 
%  Suffixes that are not in the fixed list are added as further fields. 
  if ~exist(PBIDS,'dir'), return; end
  P = cat_vol_findfiles(PBIDS,'sub-*.nii*'); 
  P = P( ~cellfun('isempty',regexp(P,'\.nii(\.gz)?$','once')) ); 
  P = P( ~cat_io_contains(P,[filesep 'derivatives' filesep]) ); 
  for fi = 1:numel(P)
    [pp,ff] = fileparts(P{fi}); 
    ff      = regexprep(ff,'\.nii$',''); 
    parts   = strsplit(ff,'_'); 
    dt      = spm_file(pp,'basename'); 
    sx      = lower(parts{end}); 
    if isvarname(dt) && isvarname(sx) 
      if ~isfield(out,dt) || ~isfield(out.(dt),sx), out.(dt).(sx) = {}; end
      out.(dt).(sx){end+1,1} = P{fi}; 
    end
  end
end
% =========================================================================
function remapIntendedFor(intended, srcmap)
%remapIntendedFor. Replace the IntendedFor references of exported BIDS 
%  inputs (e.g. fieldmaps) by the new file names (relative to the subject 
%  directory). The references of the source dataset are relative to the 
%  subject directory or BIDS URIs ("bids::sub-..."). References to files that
%  are not exported (in the same BIDS directory and subject) are removed. 
  nmiss = 0; 
  for ii = 1:numel(intended)
    I   = intended{ii}; 
    new = {}; 
    for ei = 1:numel(I.entries)
      e = strrep(I.entries{ei},'\','/'); 
      if startsWith(e,'bids::')
        src = cat_io_dcm2bids_helper('normPath',fullfile(I.srcroot, e(7:end))); 
      else
        src = cat_io_dcm2bids_helper('normPath',fullfile(I.srcroot, I.srcsub, e)); 
      end
      if isKey(srcmap,src) && strcmp(srcmap(src).root,I.root) && strcmp(srcmap(src).sub,I.sub)
        new{end+1} = strrep(srcmap(src).rel,'\','/'); %#ok<AGROW>
      else
        nmiss = nmiss + 1; 
      end
    end
    J = jsondecode(cat_io_dcm2bids_helper('readText',I.json)); 
    if isempty(new)
      J = rmfield(J,'IntendedFor'); 
    else
      J.IntendedFor = new(:); 
    end
    cat_io_dcm2bids_helper('writeJSON',I.json, J); 
  end
  if nmiss > 0
    cat_io_cprintf('warn','  %d IntendedFor references to files that were not exported were removed.\n', nmiss); 
  end
end
% =========================================================================
function writeSubjectTSV(Vjson,sub,ses,Poutdir,BIDSsubdir,anon,SOR)
%writeSubjectTSV. BIDS sessions file of a subject (sub-<label>_sessions.tsv) 
%  with one row per session: age, weight, and the session orientation 
%  rating SOR (see writeSessionAffines in cat_io_dcm2bids_pp). The rows of
%  other sessions are kept (by column name, i.e., also from files with other
%  columns). Missing values are n/a. 
  if ~exist('SOR','var') || isempty(SOR), SOR = nan; end
  Psubject = fullfile(Poutdir,BIDSsubdir,sub,sprintf('%s_sessions.tsv',sub)); % BIDS sessions file
  hdr      = {'session_id','age','weight','SOR'}; 
  Tsubject = hdr; 
  if exist(Psubject,'file')
    T0 = cat_io_csv(Psubject, '','', struct('delimiter','\t','convert2double',-1)); 
    [isc,ci] = ismember(hdr,T0(1,:)); 
    for ri = 2:size(T0,1)
      if strcmp(T0{ri,1},ses), continue; end % replaced by the current session
      Tsubject(end+1,:) = {'n/a'}; Tsubject(end,isc) = T0(ri,ci(isc)); %#ok<AGROW>
    end
  end
  Tsubject(end+1,:) = {ses, round(Vjson.PatientAge,2-anon), round(Vjson.PatientWeight), round(SOR,2)}; 
  Tsubject(cellfun(@(v) isnumeric(v) && isscalar(v) && isnan(v), Tsubject)) = {'n/a'}; 
  [~,so] = sort(Tsubject(2:end,1)); Tsubject(2:end,:) = Tsubject(so+1,:); 
  
  % write
  cat_io_csv(Psubject,Tsubject,'','',struct('delimiter','\t')); 
end
% =========================================================================
function writeParticipantTSV(P,Pparticipants,anon)
% Create participant file (eg. OpenNeuro) with sex and age (sex first as it
% is more robust, the session-specific age is in the subject files) and 
% further columns of imported BIDS datasets (extra.columns/extra.values). 
% P is a structure array with the subject (sub), its data (V), and further 
% columns (extra) that is added to (or replaced in) the table Pparticipants.
  if exist(Pparticipants,'file')
    Tparticipants = cat_io_csv(Pparticipants, '','', struct('delimiter','\t','convert2double',-2)); 
  else
    Tparticipants = {'participant_id','sex','age'}; 
  end
  for pi = 1:numel(P)
    pidpa = find( strcmp( Tparticipants(2:end,1) , P(pi).sub ) ) + 1;
    if isempty(pidpa), pidpa = size(Tparticipants,1) + 1; end
    Tparticipants(pidpa,1:3) = {P(pi).sub, P(pi).V.PatientSex, round(P(pi).V.PatientAge,2-anon) };
    extra = P(pi).extra; 
    if isstruct(extra) && isfield(extra,'columns') && ~isempty(extra.columns)
      cols = cellstr(extra.columns); vals = cellstr(string(extra.values)); 
      for ci = 1:numel(cols)
        cid = find(strcmp(Tparticipants(1,:),cols{ci}),1); 
        if isempty(cid), cid = size(Tparticipants,2) + 1; Tparticipants(:,cid) = {'n/a'}; Tparticipants{1,cid} = cols{ci}; end
        Tparticipants{pidpa,cid} = vals{ci}; 
      end
    end
  end
  Tparticipants(cellfun('isempty',Tparticipants)) = {'n/a'}; 
  [~,so] = sort(cellfun(@(x) char(string(x)), Tparticipants(2:end,1), 'UniformOutput', false)); % by ID (columns can have mixed types)
  Tparticipants(2:end,:) = Tparticipants(so+1,:); 

  cat_io_csv(Pparticipants,Tparticipants,'','',struct('delimiter','\t')); 
end
% =========================================================================
function writePrivateTSV(P,Pprivate)
% private and participant data
% The private.tsv should contain fields that are removed in the BIDS
% processing such as the real Patient name and his birth data etc. 
% It might be saved in another directory to avoid unwanted uploading?
% P is a structure array with the subject (sub) and its data (V) that is 
% added to (or replaced in) the table Pprivate.
  if exist(Pprivate,'file')
    Tprivate = cat_io_csv(Pprivate, '','', struct('delimiter','\t','convert2double',-2)); 
    pidnum   = find(cellfun(@isnumeric,Tprivate(:,2))); 
    Tprivate(pidnum,2) = cellfun(@num2str,Tprivate(pidnum,2),'UniformOutput',false); 
  else
    if ~exist(fileparts(Pprivate),'dir'), mkdir(fileparts(Pprivate)); end
    Tprivate = {'participant_id','PatientID','PatientName', ...
                'PatientSex','PatientAge','PatientWeight', ...
                'PatientBirthDate','AcquisitionDateTime'}; 
  end
  for pi = 1:numel(P)
    pidpa = find( strcmp( Tprivate(2:end,1) , P(pi).sub ) ) + 1;
    if isempty(pidpa), pidpa = size(Tprivate,1) + 1; end
    V = P(pi).V; 
    % private table with possible subject name 
    Tprivate(pidpa,:) = {P(pi).sub, V.PatientID, V.PatientName, ...
                         V.PatientSex, V.PatientAge, V.PatientWeight', ...
                         V.PatientBirthDate, V.AcquisitionDateTime};
  end
  [~,so] = sort(cellfun(@(x) char(string(x)), Tprivate(2:end,1), 'UniformOutput', false)); 
  Tprivate(2:end,:) = Tprivate(so+1,:); 
  cat_io_csv(Pprivate,Tprivate,'','',struct('delimiter','\t')); 

  % subreports
end
% =========================================================================
function sources = getBIDSsources(roots)
%getBIDSsources. Name, license, README text, and participants.json of the 
%  imported BIDS datasets (for the dataset files of the output). 
  sources = struct('Name',{},'License',{},'README',{},'participants',{}); 
  for ri = 1:numel(roots)
    if ~exist(roots{ri},'dir'), continue; end
    src = struct('Name',spm_file(roots{ri},'basename'),'License','','README','','participants',[]); 
    try
      D = jsondecode(cat_io_dcm2bids_helper('readText',fullfile(roots{ri},'dataset_description.json'))); 
      if isfield(D,'Name') && ~isempty(D.Name), src.Name = sprintf('%s (%s)',char(D.Name),spm_file(roots{ri},'basename')); end
      if isfield(D,'License'), src.License = char(D.License); end
      if isfield(D,'Licence') && isempty(src.License), src.License = char(D.Licence); end
    end
    Pr = dir(fullfile(roots{ri},'README*')); 
    if ~isempty(Pr), src.README = cat_io_dcm2bids_helper('readText',fullfile(roots{ri},Pr(1).name)); end
    if exist(fullfile(roots{ri},'participants.json'),'file')
      try %#ok<TRYNC>
        src.participants = jsondecode(cat_io_dcm2bids_helper('readText',fullfile(roots{ri},'participants.json'))); 
      end
    end
    sources(end+1) = src; %#ok<AGROW>
  end
end
% =========================================================================
function writeDatasetFiles(Pbids, job, sources, convver)
%writeDatasetFiles. BIDS dataset files of the BIDS directory Pbids: 
%  dataset_description.json (from job.dataset, with automatic BIDSVersion, 
%  DatasetType, and GeneratedBy), README.md (head from the job.dataset.README
%  text file, extended by the README files of imported BIDS datasets), and 
%  participants.json (description of the columns of participants.tsv). 
%  sources is an optional structure array of imported BIDS datasets (fields
%  Name, License, README, participants) for the license, README, and 
%  column descriptions. convver are the dcm2niix versions of the imports 
%  (from the import information of the scans). 
  if ~exist('sources','var'), sources = struct('Name',{},'License',{},'README',{},'participants',{}); end
  if ~exist('convver','var'), convver = {}; end
  ds    = job.dataset; 
  clean = @(c) reshape(c(~cellfun('isempty',strtrim(cellstr(c)))),1,[]); 

  % dataset_description.json (only non-empty fields)
  D = struct(); 
  if isempty(strtrim(ds.Name)), D.Name = job.subdir; else, D.Name = strtrim(ds.Name); end
  D.BIDSVersion = '1.10.0'; 
  D.DatasetType = 'raw'; 
  lic = mostRestrictiveLicense([{ds.License}, {sources.License}], D.Name, {sources.Name}); 
  if ~isempty(lic),                             D.License            = lic; end
  if ~isempty(clean(ds.Authors)),               D.Authors            = clean(ds.Authors); end
  if ~isempty(strtrim(ds.HowToAcknowledge)),    D.HowToAcknowledge   = strtrim(ds.HowToAcknowledge); end
  if ~isempty(strtrim(ds.Acknowledgements)),    D.Acknowledgements   = strtrim(ds.Acknowledgements); end
  if ~isempty(clean(ds.Funding)),               D.Funding            = clean(ds.Funding); end
  if ~isempty(clean(ds.EthicsApprovals)),       D.EthicsApprovals    = clean(ds.EthicsApprovals); end
  if ~isempty(clean(ds.ReferencesAndLinks)),    D.ReferencesAndLinks = clean(ds.ReferencesAndLinks); end
  if ~isempty(strtrim(ds.DatasetDOI)),          D.DatasetDOI         = strtrim(ds.DatasetDOI); end
  try
    [~,catver] = cat_version; 
  catch
    catver = ''; 
  end
  D.GeneratedBy = {struct('Name','catDCM2BIDS','Version',char(catver), ...
    'Description','DICOM/BIDS import, protocol evaluation, anonymization, and QC of CAT12', ...
    'CodeURL','https://github.com/ChristianGaser/cat12')}; 
  if ~isempty(convver) % DICOM imports
    D.GeneratedBy{end+1} = struct('Name','dcm2niix','Version',strjoin(cellstr(convver),', '), ...
      'CodeURL','https://github.com/rordenlab/dcm2niix'); 
  end
  cat_io_dcm2bids_helper('writeJSON',fullfile(Pbids,'dataset_description.json'), D); 

  % README.md: head and the README files of imported BIDS datasets
  Preadme = char(ds.README); 
  if ~isempty(Preadme) && exist(Preadme,'file')
    txt = cat_io_dcm2bids_helper('readText',Preadme); 
  else
    txt = sprintf(['# %s\n\nBIDS dataset created by catDCM2BIDS (CAT12), including the evaluation ' ...
      'of the MRI protocols, anonymization, and quality control.\n'], D.Name); 
  end
  for si = 1:numel(sources)
    if ~isempty(sources(si).README)
      txt = sprintf('%s\n\n## %s\n\n%s\n', strtrim(txt), sources(si).Name, strtrim(sources(si).README)); 
    end
  end
  cat_io_dcm2bids_helper('writeText',fullfile(Pbids,'README.md'), txt); 

  % participants.json: description of the columns of participants.tsv
  Pparticipants = fullfile(Pbids,'participants.tsv'); 
  if exist(Pparticipants,'file')
    fid  = fopen(Pparticipants,'r','n','UTF-8'); hdr = fgetl(fid); fclose(fid); 
    cols = strsplit(strtrim(hdr),'\t'); 
    P = struct(); 
    for ci = 2:numel(cols)
      col = cols{ci}; 
      switch col
        case 'sex'
          P.sex = struct('LongName','sex', 'Description','Phenotypical sex of the participant.', ...
            'Levels',struct('M','male','F','female','O','other')); 
        case 'age'
          P.age = struct('LongName','age', 'Description', ...
            'Age of the participant at the (first) scan (session-specific ages are in the sub-*_sessions.tsv files).', ...
            'Units','year'); 
        otherwise
          % description of the source datasets (if available)
          desc = struct('LongName',col,'Description',''); 
          for si = 1:numel(sources)
            if isstruct(sources(si).participants) && isfield(sources(si).participants,col)
              desc = sources(si).participants.(col); break
            end
          end
          if isvarname(col), P.(col) = desc; end
      end
    end
    cat_io_dcm2bids_helper('writeJSON',fullfile(Pbids,'participants.json'), P); 
  end
end
% =========================================================================
function lic = mostRestrictiveLicense(lics, name, srcnames)
%mostRestrictiveLicense. The most restrictive license of the dataset 
%  setting lics{1} and the licenses of the source datasets lics{2:end} (as 
%  SPDX identifiers, free text is mapped), with a warning if it differs from
%  the dataset setting. Unknown licenses of source datasets are more 
%  restrictive than all known ones. 
  order = {'CC0-1.0','PDDL-1.0','ODC-BY-1.0','CC-BY-4.0','CC-BY-SA-4.0','CC-BY-NC-4.0'}; 
  spdx  = cellfun(@licenseSPDX, lics, 'UniformOutput', false); 
  rank  = zeros(size(spdx)); 
  for li = 1:numel(spdx)
    if isempty(lics{li})
      rank(li) = 0;                          % no license entry
    elseif any(strcmp(order,spdx{li}))
      rank(li) = find(strcmp(order,spdx{li})); 
    else
      rank(li) = numel(order) + 1;           % unknown/other license
    end
  end
  lic = spdx{1}; 
  if numel(rank) > 1 && max(rank(2:end)) > rank(1)
    [~,mi] = max(rank); 
    lic = spdx{mi}; 
    cat_io_cprintf('warn',['  The license "%s" of "%s" is less restrictive than the license "%s" of the ' ...
      'imported dataset "%s" - the more restrictive license is used.\n'], lics{1}, name, lics{mi}, srcnames{mi-1}); 
  end
end
% =========================================================================
function spdx = licenseSPDX(lic)
%licenseSPDX. SPDX identifier of a (free text) license. 
  l = lower(char(lic)); spdx = strtrim(char(lic)); 
  if     ~isempty(regexp(l,'by[- ]?nc|non[- ]?commercial','once')),  spdx = 'CC-BY-NC-4.0'; 
  elseif ~isempty(regexp(l,'by[- ]?sa|share[- ]?alike','once')),     spdx = 'CC-BY-SA-4.0'; 
  elseif ~isempty(regexp(l,'odc[- ]?by','once')),                     spdx = 'ODC-BY-1.0'; 
  elseif ~isempty(regexp(l,'cc[ -]?by|creative commons attribution','once')), spdx = 'CC-BY-4.0'; 
  elseif ~isempty(regexp(l,'cc[- ]?0|cc[- ]?o\>|public domain dedication|creative commons zero','once')), spdx = 'CC0-1.0'; 
  elseif contains(l,'pddl'),                                      spdx = 'PDDL-1.0'; 
  end
end
% =========================================================================
function gzipunzipOutputdata(Poutdir,opts)
%gzipunzipOutputdata. Assure the internal .nii.gz status defined by gzipi
  BIDSsubdirs = cat_vol_findfiles( Poutdir , 'BIDS*', struct('dirs',1)); 
  if opts.gzipi == 1  &&  opts.gzipe == 0  &&  ~isempty(BIDSsubdirs)
    % gunzip all files in the result directory
    fprintf('DCM2BIDS - gunzip BIDS output data')
    for bi = 1:numel(BIDSsubdirs)
      Pniigz = cat_vol_findfiles( BIDSsubdirs{bi} , '*.nii.gz' ); 
      for fi = 1:numel(Pniigz), gunzip(Pniigz{fi}); delete(Pniigz{fi}); end
    end
    fprintf('done.\n')
  elseif opts.gzipi == 0  &&  opts.gzipe == 1  &&  ~isempty(BIDSsubdirs)
    % gzip all files in the result directory
    fprintf('DCM2BIDS - gzip BIDS output data')
    for bi = 1:numel(BIDSsubdirs)
      Pniigz = cat_vol_findfiles( BIDSsubdirs{bi} , '*.nii' ); 
      for fi = 1:numel(Pniigz), gzip(Pniigz{fi}); delete(Pniigz{fi}); end
    end
    fprintf('done.\n')
  end
end
