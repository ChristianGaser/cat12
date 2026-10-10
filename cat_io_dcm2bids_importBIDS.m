function Psnr = cat_io_dcm2bids_importBIDS(Pbids, Pmdbdir, opts)
%cat_io_dcm2bids_importBIDS. Import the MRI scans (NIfTI and JSON) of BIDS datasets. 
%  For each NIfTI of the MRI datatypes (anat, func, dwi, fmap, perf) of the 
%  selection Pbids(i).sel of the dataset Pbids(i).root (without derivatives,
%  sourcedata, and code), the metadata is merged from all JSON files that 
%  apply by the BIDS inheritance principle (sidecar and general JSON files 
%  of the dataset, subject, and session level with the same suffix and a 
%  subset of the entities). The scan is stored in the database like DICOM 
%  scans: 
%    sub-<dataset>-<subject>/ses-<session>[-<date>]/snr-<series>-<protocol>/
%      <dataset>=<protocol>=<time or session>=<series>.json
%      catDCM2BIDSimport_<...>.json (source file, datatype, suffix, entities, 
%                                    participant data)
%  where the dataset (folder) name avoids collisions of subject IDs of 
%  different datasets, the session is 1 without session directory, the 
%  scan date is taken from the scans.tsv/sessions.tsv (acq_time) or JSON 
%  (AcquisitionDateTime) if available, and the series number is taken from
%  the JSON or numbered by the file order of the subject/session. Sex and 
%  age (and further columns) are taken from the participants.tsv (and the 
%  sessions.tsv). The images are copied only if required (see 
%  convertImages in cat_io_dcm2bids_importDCM). Previously imported files 
%  are not imported again (only with opts.rerun). 
%  Psnr is the sorted list of the scan directories of the selection. 
%
%  Psnr = cat_io_dcm2bids_importBIDS(Pbids,Pmdbdir,opts)
%
%  Pbids   .. BIDS selections (structure array with the fields root and sel, 
%             see getInputSources in cat_io_dcm2bids_import)
%  Pmdbdir .. database directory (catDCM2BIDSdb)
%  opts    .. options of cat_io_dcm2bids (e.g. rerun)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_db.

  Psnr = {}; 
  if isempty(Pbids), return; end
  fprintf('DCM2BIDS - import %d BIDS dataset selection(s)\n', numel(Pbids)); 
  mri = {'anat','func','dwi','fmap','perf'}; 
  [~,~,~,Pimpsrc,Pimpsnr] = cat_io_dcm2bids_db('getImportIndex',Pmdbdir); 
  fs = regexptranslate('escape',filesep); 

  for bi = 1:numel(Pbids)
    R      = Pbids(bi).root; 
    dsname = spm_file(R,'basename'); 
    dslab  = cat_io_dcm2bids_bids('bidsLabel',dsname); 
    nnew   = 0; nold = 0; nside = 0; ninh = 0; nnone = 0; 

    % NIfTI files of the selection (without derivatives, sourcedata, code)
    P   = cat_vol_findfiles(Pbids(bi).sel,'*.nii*'); 
    P   = P( ~cellfun('isempty',regexp(P,'\.nii(\.gz)?$','once')) ); 
    rel = cellfun(@(p) p(numel(R)+2:end), P, 'UniformOutput', false); 
    ok  = ~cellfun('isempty',regexp(rel,['^sub-[^' fs ']+' fs '(ses-[^' fs ']+' fs ')?[^' fs ']+' fs '[^' fs ']+$'],'once')); 
    P = sort(P(ok)); 
    dt  = cellfun(@(p) spm_file(fileparts(p),'basename'), P, 'UniformOutput', false); 
    nonmri = ~ismember(dt,mri); 
    if any(nonmri)
      cat_io_cprintf([.5 .5 .5],'  %s: ignore %d images of other modalities (%s)\n', dsname, sum(nonmri), ...
        strjoin(unique(dt(nonmri)),', ')); 
    end
    P = P(~nonmri); dt = dt(~nonmri); 

    % participant data
    Tpart = readTSV(fullfile(R,'participants.tsv')); 
    sn    = containers.Map('KeyType','char','ValueType','double'); % series counter per session

    for fi = 1:numel(P)
      % previously imported file
      imported = strcmp(Pimpsrc, P{fi}); 
      if any(imported) && ~opts.rerun
        Psnr = [Psnr; Pimpsnr(imported)]; %#ok<AGROW>
        nold = nold + 1; 
        continue
      end

      % entities and metadata (BIDS inheritance)
      [~,name] = fileparts(regexprep(P{fi},'\.gz$','')); 
      ent      = parseBIDSname(name); 
      [J,nj,own] = bidsMetadata(R, P{fi}, ent, '.json'); 
      if own, nside = nside + 1; elseif nj, ninh = ninh + 1; else, nnone = nnone + 1; end
      sub = ent.entities.sub; 
      if isfield(ent.entities,'ses'), ses = ent.entities.ses; else, ses = '1'; end

      % acquisition time of the scans.tsv or sessions.tsv
      acqtime = ''; 
      subdir  = fullfile(R,['sub-' sub]); 
      if isfield(ent.entities,'ses')
        Tscans = readTSV(fullfile(subdir,['ses-' ses],sprintf('sub-%s_ses-%s_scans.tsv',sub,ses))); 
        frel   = P{fi}(numel(fullfile(subdir,['ses-' ses]))+2:end); 
      else
        Tscans = readTSV(fullfile(subdir,sprintf('sub-%s_scans.tsv',sub))); 
        frel   = P{fi}(numel(subdir)+2:end); 
      end
      v = tsvValue(Tscans,'filename',strrep(frel,filesep,'/'),'acq_time'); 
      if ~isempty(v), acqtime = v; end
      Tses = readTSV(fullfile(subdir,sprintf('sub-%s_sessions.tsv',sub))); 
      if isempty(acqtime)
        v = tsvValue(Tses,'session_id',['ses-' ses],'acq_time'); 
        if ~isempty(v), acqtime = v; end
      end
      if isempty(acqtime) && isfield(J,'AcquisitionDateTime'), acqtime = char(J.AcquisitionDateTime); end
      try
        acqdt = datetime(regexprep(acqtime,'Z$',''),'InputFormat','yyyy-MM-dd''T''HH:mm:ss'); 
      catch
        try
          acqdt = datetime(acqtime(1:min(10,end)),'InputFormat','yyyy-MM-dd'); 
        catch
          acqdt = NaT; 
        end
      end

      % series number of the JSON or by the file order of the session
      skey = [sub '|' ses]; 
      if isKey(sn,skey), sn(skey) = sn(skey) + 1; else, sn(skey) = 1; end
      if isfield(J,'SeriesNumber') && isnumeric(J.SeriesNumber) && isscalar(J.SeriesNumber)
        series = J.SeriesNumber; 
      else
        series = sn(skey); 
      end

      % protocol name of the JSON or by the entities (without sub, ses, run)
      if isfield(J,'ProtocolName') && ~isempty(J.ProtocolName)
        prot = char(J.ProtocolName); 
      else
        keys = setdiff(ent.keys,{'sub','ses','run'},'stable'); 
        prot = strjoin([cellfun(@(k) [k '-' ent.entities.(k)], keys, 'UniformOutput', false) {ent.suffix}],'_'); 
      end
      prot = regexprep(prot,'[=/\\:*?"<>| ]','_'); 

      % participant data (sex, age, further columns) and session age
      [psex, page, pcols, pvals] = participantData(Tpart, ['sub-' sub]); 
      v = str2double(tsvValue(Tses,'session_id',['ses-' ses],'age')); 
      if ~isnan(v), page = v; end

      % sidecar of the database (metadata with the identity fields)
      V = J; 
      V.PatientID          = sub; 
      V.DeviceSerialNumber = dslab; 
      V.StudyID            = ses; 
      V.SeriesNumber       = series; 
      V.ProtocolName       = prot; 
      if ~isfield(V,'SeriesDescription') || isempty(V.SeriesDescription), V.SeriesDescription = prot; end
      V.PatientSex         = psex; 
      V.PatientAge         = page; 
      if ~isnat(acqdt), V.AcquisitionDateTime = char(acqdt,'yyyy-MM-dd''T''HH:mm:ss'); end
      if isnat(acqdt), tpart = cat_io_dcm2bids_bids('bidsLabel',ses); else, tpart = char(acqdt,'yyyyMMddHHmmss'); end
      sidecar    = sprintf('%s=%s=%s=%d', dslab, prot, tpart, series); 
      V0         = cat_io_dcm2bids_db('assurePatientDCMfields',V); 
      Pdbdirpath = fullfile(Pmdbdir,cat_io_dcm2bids_db('getDBdir',V0)); 
      Pjsondb    = fullfile(Pdbdirpath,[sidecar '.json']); 
      if ~exist(Pdbdirpath,'dir'), mkdir(Pdbdirpath); end
      if ~exist(Pjsondb,'file') || opts.rerun
        cat_io_dcm2bids_helper('writeJSON',Pjsondb, V); 
        nnew = nnew + 1; 
      end

      % import information
      I = struct('SourceType','BIDS', 'SourceFile',P{fi}, 'SourceDataset',R, 'SourceDatasetName',dsname, ...
        'Datatype',dt{fi}, 'Suffix',ent.suffix, 'Entities',ent.entities, ...
        'Participant',struct('columns',{pcols},'values',{pvals}), ...
        'ImportDate',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
        'ScanID',cat_io_dcm2bids_db('getScanID',V0,strsplit(sidecar,'='))); 
      if strcmp(dt{fi},'dwi')
        [~,~,~,I.SourceBval] = bidsMetadata(R, P{fi}, ent, '.bval'); 
        [~,~,~,I.SourceBvec] = bidsMetadata(R, P{fi}, ent, '.bvec'); 
      end
      cat_io_dcm2bids_db('writeImportJSON',spm_file(Pjsondb,'prefix','catDCM2BIDSimport_'), I); 

      % files of the scan (e.g. events, physio, aslcontext) 
      A = bidsAssociated(R, P{fi}, ent); 
      for ai = 1:size(A,1)
        copyfile(A{ai,1}, fullfile(Pdbdirpath,sprintf('catDCM2BIDSassoc_%s_%s',sidecar,A{ai,2}))); 
      end
      Psnr{end+1,1} = Pdbdirpath; %#ok<AGROW>
    end
    fprintf('  %s: %d new scans, %d imported before (JSON: %d sidecar, %d general, %d none)\n', ...
      dsname, nnew, nold, nside, ninh, nnone); 
    if nnone > 0
      cat_io_cprintf('warn','  %s: %d images without JSON metadata (no protocol information).\n', dsname, nnone); 
    end
  end
  Psnr = unique(Psnr); 
  fprintf('\n'); 
end
% =========================================================================
function [M,n,own,Pfile] = bidsMetadata(R, Pdata, ent, ext)
%bidsMetadata. Metadata of a data file by the BIDS inheritance principle. 
%  All files with the extension ext (e.g. .json, .bval) in the directory of
%  the data file and its parent directories up to the dataset root R apply,
%  if they have the same suffix and all their entities are entities of the
%  data file with the same values. JSON files are merged from the top level
%  to the data file (more specific files overwrite fields). n is the number 
%  of applying files, own is true for a sidecar with the name of the data 
%  file, and Pfile is the most specific file. 
  M = struct(); n = 0; own = false; Pfile = ''; 
  d = fileparts(Pdata); dirs = {d}; 
  while numel(d) > numel(R)
    d = fileparts(d); dirs = [{d} dirs]; %#ok<AGROW>
  end
  for di = 1:numel(dirs)
    D = dir(fullfile(dirs{di},['*' ext])); 
    cand = {}; nent = []; 
    for fi = 1:numel(D)
      e = parseBIDSname(regexprep(D(fi).name,[regexptranslate('escape',ext) '$'],'')); 
      if ~strcmp(e.suffix,ent.suffix), continue; end
      fit = all(cellfun(@(k) isfield(ent.entities,k) && strcmp(ent.entities.(k),e.entities.(k)), e.keys)); 
      if fit, cand{end+1} = fullfile(dirs{di},D(fi).name); nent(end+1) = numel(e.keys); end %#ok<AGROW>
    end
    [~,so] = sort(nent); 
    for ci = so
      n = n + 1; Pfile = cand{ci}; 
      if strcmp(ext,'.json')
        try
          J = jsondecode(cat_io_dcm2bids_helper('readText',cand{ci})); 
          fn = fieldnames(J); 
          for fni = 1:numel(fn), M.(fn{fni}) = J.(fn{fni}); end
        catch
          cat_io_cprintf('warn','  Cannot read "%s".\n',cand{ci}); 
        end
      end
      [~,cn] = fileparts(cand{ci}); [~,dn] = fileparts(regexprep(Pdata,'\.gz$','')); 
      own = own || strcmp(cn,dn); 
    end
  end
end
% =========================================================================
function A = bidsAssociated(R, Pdata, ent)
%bidsAssociated. Files that belong to a BIDS data file, i.e., events, 
%  physio, stim, aslcontext, and asllabeling files of the data file (same 
%  entities, also with further entities such as recording-cardiac) or 
%  inherited from higher levels (BIDS inheritance principle, subset of the 
%  entities). A is a cell with the files and their name part after the 
%  entities of the data file (e.g. "events.tsv", "recording-cardiac_physio.tsv.gz"). 
  A    = cell(0,2); 
  asuf = {'events','physio','stim','aslcontext','asllabeling'}; 
  aext = '\.(tsv|tsv\.gz|json|jpg|png)$'; 
  [pp,name] = fileparts(regexprep(Pdata,'\.gz$','')); 
  base = regexprep(name,['_' ent.suffix '$'],''); 

  % files of the data file (same directory and entities)
  D = dir(fullfile(pp,[base '_*'])); 
  for di = 1:numel(D)
    rest = D(di).name(numel(base)+2:end); 
    suf  = regexp(rest,['(^|_)([a-zA-Z]+)' aext],'tokens','once'); 
    if ~isempty(suf) && any(strcmp(suf{2},asuf)), A(end+1,:) = {fullfile(pp,D(di).name), rest}; end %#ok<AGROW>
  end

  % inherited files (only if there is no own file with the same name part)
  d = pp; dirs = {d}; 
  while numel(d) > numel(R)
    d = fileparts(d); dirs = [{d} dirs]; %#ok<AGROW>
  end
  for di = 1:numel(dirs)
    D = dir(dirs{di}); 
    for fi = 1:numel(D)
      if D(fi).isdir, continue; end
      tok = regexp(D(fi).name,['^(.*)' aext],'tokens','once'); 
      if isempty(tok), continue; end
      e = parseBIDSname(tok{1}); 
      if ~any(strcmp(e.suffix,asuf)) || (isfield(e.entities,'sub') && ~strcmp(e.entities.sub,ent.entities.sub)), continue; end
      % entities of the inherited file that are not entities of the data file 
      % (e.g. recording) are kept, all others have to fit
      keys = e.keys(cellfun(@(k) isfield(ent.entities,k), e.keys)); 
      if ~all(cellfun(@(k) strcmp(ent.entities.(k),e.entities.(k)), keys)), continue; end
      extra = setdiff(e.keys,keys,'stable'); 
      rest  = strjoin([cellfun(@(k) [k '-' e.entities.(k)], extra, 'UniformOutput', false) ...
        {[e.suffix regexp(D(fi).name,aext,'match','once')]}],'_'); 
      if strcmp(fullfile(dirs{di},D(fi).name), Pdata) || any(strcmp(A(:,2),rest)), continue; end
      A(end+1,:) = {fullfile(dirs{di},D(fi).name), rest}; %#ok<AGROW>
    end
  end
end
% =========================================================================
function ent = parseBIDSname(name)
%parseBIDSname. Entities (key-value pairs) and suffix of a BIDS file name 
%  (without extension). 
  parts = strsplit(name,'_'); 
  ent   = struct('keys',{{}},'entities',struct(),'suffix',''); 
  for k = 1:numel(parts)
    kv = regexp(parts{k},'^([a-zA-Z0-9]+)-(.+)$','tokens','once'); 
    if ~isempty(kv) && isvarname(kv{1})
      ent.keys{end+1} = kv{1}; ent.entities.(kv{1}) = kv{2}; 
    elseif k == numel(parts)
      ent.suffix = parts{k}; 
    end
  end
end
% =========================================================================
function T = readTSV(P)
%readTSV. Table of a TSV file as cell (first row with the header), empty 
%  if the file does not exist. 
  T = {}; 
  if ~exist(P,'file'), return; end
  txt   = regexprep(cat_io_dcm2bids_helper('readText',P),['^' char(65279)],''); % without byte order mark
  lines = regexp(txt,'\r?\n','split'); 
  lines = lines(~cellfun('isempty',strtrim(lines))); 
  if isempty(lines), return; end
  hdr = strtrim(strsplit(lines{1},'\t')); 
  T   = repmat({''},numel(lines),numel(hdr)); T(1,:) = hdr; 
  for li = 2:numel(lines)
    v = strtrim(strsplit(lines{li},'\t')); 
    T(li,1:min(end,numel(v))) = v(1:min(end,numel(hdr))); 
  end
end
% =========================================================================
function v = tsvValue(T, keycol, key, col)
%tsvValue. Value of the column col in the row with key in the column keycol
%  (empty if not available or n/a). 
  v = ''; 
  if isempty(T), return; end
  kc = find(strcmp(T(1,:),keycol),1); vc = find(strcmp(T(1,:),col),1); 
  if isempty(kc) || isempty(vc), return; end
  ri = find(strcmp(T(2:end,kc),key),1) + 1; 
  if ~isempty(ri) && ~any(strcmpi(T{ri,vc},{'n/a','na',''})), v = T{ri,vc}; end
end
% =========================================================================
function [sex, age, cols, vals] = participantData(T, sub)
%participantData. Sex (M/F/O or NA), age (years or NaN), and the further 
%  columns of a participant of a participants.tsv table T. 
  sex = 'NA'; age = NaN; cols = {}; vals = {}; 
  if isempty(T), return; end
  idc = find(strcmp(T(1,:),'participant_id'),1); 
  if isempty(idc), return; end
  ri  = find(strcmp(T(2:end,idc),sub),1) + 1; 
  if isempty(ri), return; end
  hdr  = lower(T(1,:)); 
  sexc = find(ismember(hdr,{'sex','gender'}),1); 
  agec = find(strcmp(hdr,'age'),1); 
  if isempty(agec), agec = find(startsWith(hdr,'age'),1); end
  if ~isempty(sexc)
    s = lower(strtrim(T{ri,sexc})); 
    if     startsWith(s,'m'), sex = 'M'; 
    elseif startsWith(s,'f') || startsWith(s,'w'), sex = 'F'; 
    elseif startsWith(s,'o') || startsWith(s,'d'), sex = 'O'; 
    end
  end
  if ~isempty(agec)
    a = regexp(T{ri,agec},'[0-9]+(\.[0-9]+)?','match','once'); 
    if ~isempty(a), age = str2double(a); end
  end
  other = setdiff(1:size(T,2),[idc sexc agec]); 
  cols  = T(1,other); vals = T(ri,other); 
end
