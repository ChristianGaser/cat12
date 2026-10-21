function varargout = cat_io_dcm2bids_bids(action,varargin)
%cat_io_dcm2bids_bids. BIDS naming and output files of cat_io_dcm2bids.
%  Naming: BIDS conform labels (bidsLabel), the session label of a 
%  database session (sessionLabel, by the anonymization level), and the 
%  datatype, suffix, entities, and task of a scan by its protocol name and
%  sidecar (setupDatatype, setupSuffix, entityString, setupTask). BIDSoutputs lists
%  the datatypes and suffixes of the output structure (batch dependencies),
%  getBIDSoutputs collects the exported images. 
%  Output files: sessions.tsv of each subject (writeSubjectTSV, with SOR and
%  the phenotypic data of cat_io_dcm2bids_phenotypic), 
%  participants.tsv and the private participants table (writeParticipantTSV,
%  writePrivateTSV), the IntendedFor references of BIDS inputs 
%  (remapIntendedFor), dataset_description.json, README.md, 
%  participants.json, and .bidsignore (writeDatasetFiles with the information
%  of the imported BIDS datasets, see getBIDSsources), the generalized MRI
%  protocols of the scans (writeProtocols), and the gzip status of the 
%  output (gzipunzipOutputdata).
%
%  varargout = cat_io_dcm2bids_bids(action,varargin)
%
%  Actions (see the help of the local functions):
%    label = cat_io_dcm2bids_bids('bidsLabel',str,camelCase)
%    ses = cat_io_dcm2bids_bids('sessionLabel',Pdbses,stamp,opts,offset)
%    datatype = cat_io_dcm2bids_bids('setupDatatype',pro,V)
%    [suffix,ent,label] = cat_io_dcm2bids_bids('setupSuffix',datatype,pro,V)
%    str = cat_io_dcm2bids_bids('entityString',E)
%    task = cat_io_dcm2bids_bids('setupTask',datatype,pro)
%    out = cat_io_dcm2bids_bids('BIDSoutputs')
%    out = cat_io_dcm2bids_bids('getBIDSoutputs',PBIDS,out)
%    cat_io_dcm2bids_bids('remapIntendedFor',intended,srcmap)
%    cat_io_dcm2bids_bids('writeSubjectTSV',Vjson,sub,ses,Poutdir,BIDSsubdir,anon,SOR,site,extra)
%    cat_io_dcm2bids_bids('writeParticipantTSV',P,Pparticipants,anon)
%    cat_io_dcm2bids_bids('writePrivateTSV',P,Pprivate)
%    sources = cat_io_dcm2bids_bids('getBIDSsources',roots)
%    cat_io_dcm2bids_bids('writeDatasetFiles',Pbids,job,sources,convver)
%    cat_io_dcm2bids_bids('writeProtocols',Pbids,Pp,Ptbldir)
%    cat_io_dcm2bids_bids('gzipunzipOutputdata',Poutdir,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'bidsLabel','sessionLabel','setupDatatype','setupSuffix','entityString','setupTask','BIDSoutputs', ...
          'getBIDSoutputs','remapIntendedFor','writeSubjectTSV', ...
          'writeParticipantTSV','writePrivateTSV','getBIDSsources', ...
          'writeDatasetFiles','writeProtocols','gzipunzipOutputdata'}
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
function ses = sessionLabel(Pdbses,stamp,opts,offset)
%sessionLabel. BIDS session label of the database session directory Pdbses
%  (ses-<StudyID>-<yyyymmdd>) with the study date and time stamp of its
%  files (yyyymmddHHMMSS): the label of a session without date or time 
%  (e.g. BIDS input), the date and time (opts.anonymize = 0), the date (1),
%  or the time since the first session of the subject (2, with the date 
%  offset of the subject, see session in cat_io_dcm2bids_anonymize, also 
%  for imports of anonymized data with dates from 1850 to 1949, invalid 
%  dates as date), and NA for a time without date (opts.anonymize > 0). 
%  Studies of a subject with the same label (e.g. two sessions of a day 
%  for opts.anonymize = 1 or of a week for 2) get a letter (b, c, ...) in
%  the order of their date, StudyID, and time.
  isdate = ~isempty(regexp(stamp,'^[0-9]{8}','once'));
  if ~isdate && isempty(regexp(stamp,'^[0-9]{1,7}$','once')) % no date or time (e.g. BIDS input)
    ses = ['ses-' bidsLabel(stamp)]; return
  end
  shifted = isdate && str2double(stamp(1:4)) >= 1850 && str2double(stamp(1:4)) < 1950; % dates of an anonymized export
  if opts.anonymize == 0 && ~shifted % date and time
    if isdate, ses = ['ses-' stamp(1:8)]; else, ses = ['ses-' stamp]; end
    if isdate && numel(stamp) > 8, ses = [ses opts.BIDSsep stamp(9:end)]; end
    return
  end
  if ~isdate
    label = @(t) 'NA'; % no date (and no time of the day)
  elseif opts.anonymize > 1 || shifted
    label = @(t) timeLabel(t(1:8),offset * (opts.anonymize > 1),opts.sesunit);
  else
    label = @(t) t(1:8); % date (opts.anonymize = 1)
  end
  ses = ['ses-' label(stamp)];

  % studies of the subject with this label (time stamp of their files and
  % StudyID of their directories, e.g. two directories of a study over 
  % midnight have the same key)
  [~,name] = fileparts(Pdbses); key = {}; own = '';
  D = dir(fullfile(fileparts(Pdbses),'ses-*')); D = D([D.isdir]);
  for di = 1:numel(D)
    F = dir(fullfile(D(di).folder,D(di).name,'snr-*','*=*=*=*.json')); 
    F = F(~startsWith({F.name},'catDCM2BIDS')); 
    if isempty(F), continue; end
    t = cat_io_dcm2bids_db('nameParts',spm_file(F(1).name,'basename')); t = t{min(3,end)}; 
    if isdate ~= ~isempty(regexp(t,'^[0-9]{8}','once')) || isempty(regexp(t,'^[0-9]','once')) || ...
       ~strcmp(['ses-' label(t)],ses), continue; end
    id = regexp(D(di).name,'^ses-(.*?)(-\d{8})?$','tokens','once'); if isempty(id), id = {''}; end
    if ~isempty(regexp(id{1},'^\d+$','once')), id{1} = sprintf('%012d',str2double(id{1})); end % numeric order
    if isdate, d = t(1:8); t = t(9:end); else, d = ''; end
    key{end+1} = [d char(1) id{1} char(1) t]; %#ok<AGROW> (separator before all characters)
    if strcmp(D(di).name,name), own = key{end}; end
  end
  r = find(strcmp(unique(key),own),1); 
  if r > 1, ses = [ses char('a' + r - 1)]; end
end
% =========================================================================
function l = timeLabel(d,offset,unit)
%timeLabel. Label of the time since the first session of the date d (see
%  session in cat_io_dcm2bids_anonymize) or d for an invalid date (e.g. 
%  10101000 of anonymized data).
  try
    l = cat_io_dcm2bids_anonymize('session',d,offset,unit);
  catch
    l = d;
  end
end
% =========================================================================
function datatype = setupDatatype(pro,V)
%setupDatatype. BIDS datatype of a scan by its protocol name pro, or for
%  protocol names without this information (e.g. of other vendors or sites)
%  by its sidecar V (see sidecarDatatype), otherwise other.
  if ~exist('V','var') || ~isstruct(V), V = struct(); end
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
  elseif cat_io_contains( lower(pro) , {'asl'} ) % perfusion (not supported)
    datatype  = 'other'; 
  else
    datatype = sidecarDatatype(V); 
    if isempty(datatype) % unknown protocols
      cat_io_cprintf('red', sprintf('\n  Unknown BIDS datatype (anat/func/...) for protocol "%s"\n', lower(pro)) ) 
      datatype = 'other';
    end
  end
end
% =========================================================================
function datatype = sidecarDatatype(V)
%sidecarDatatype. BIDS datatype of a scan by its sidecar V (with the image
%  fields Dimensions and BValues of cat_io_dcm2bids), i.e. for protocol 
%  names without datatype information: dwi (b-values above 50 s/mm2), fmap
%  (GRE field maps with EchoTime1/2 or Hz units, EPI with up to 3 volumes,
%  e.g. spin-echo field maps or b0 images), func (EPI time series or SBRef
%  in the series description), anat (other images with one volume), or ''
%  (e.g. ASL with PERFUSION in the ImageType, or other time series).
  has = @(fn) isfield(V,fn) && ~isempty(V.(fn));
  txt = @(fn) char(strjoin(string(V.(fn)),'\'));
  seq = ''; for fn = {'ScanningSequence','SequenceName','PulseSequenceName'}, if has(fn{1}), seq = [seq ' ' txt(fn{1})]; end, end %#ok<AGROW>
  epi = ~isempty(regexpi(seq,'(^|[^a-z])ep([^a-z]|fid|se|i)|epi','once')); 
  nv  = 1; if has('Dimensions') && isnumeric(V.Dimensions) && numel(V.Dimensions) > 3, nv = V.Dimensions(4); end
  ser = ''; if has('SeriesDescription'), ser = lower(txt('SeriesDescription')); end
  itype = ''; if has('ImageType'), itype = upper(txt('ImageType')); end
  if has('ArterialSpinLabelingType') || ~isempty(regexpi(seq,'asl','once')) || contains(itype,{'PERFUSION','ASL'})
    datatype = ''; % perfusion (not supported)
  elseif has('BValues') && isnumeric(V.BValues) && any(V.BValues(:) > 50)
    datatype = 'dwi';
  elseif (has('EchoTime1') && has('EchoTime2')) || (has('Units') && strcmpi(txt('Units'),'Hz'))
    datatype = 'fmap';
  elseif epi && contains(ser,'sbref')
    datatype = 'func';
  elseif epi
    if nv > 3, datatype = 'func'; else, datatype = 'fmap'; end
  elseif nv == 1 && has('RepetitionTime') && has('EchoTime')
    datatype = 'anat';
  else
    datatype = '';
  end
end
% =========================================================================
function suffix = anatSuffix(V)
%anatSuffix. BIDS suffix of an anatomical scan by the timing of its sidecar
%  V (seconds) for protocol names without this information: FLAIR (TI 
%  above 1.5 and TR above 4), T1w (TI from 0.3 to 1.5, TR below 1 with TE 
%  below 0.03, or single-echo gradient echo with TR below 0.1 and TE below
%  0.01), T2star (gradient echo with TE from 0.015 or several echoes), T2w
%  (TR above 2 and TE above 0.06), or PDw (TR above 2 and TE below 0.03), 
%  otherwise ''.
  val = @(fn) isfield(V,fn) && isnumeric(V.(fn)) && ~isempty(V.(fn));
  TR = nan; TE = nan; TI = nan;
  if val('RepetitionTime'), TR = V.RepetitionTime(1); end
  if val('EchoTime'),       TE = V.EchoTime(1); end
  if val('InversionTime'),  TI = V.InversionTime(1); end
  gre = isfield(V,'ScanningSequence') && contains(upper(char(strjoin(string(V.ScanningSequence),'\'))),'GR');
  me  = isfield(V,'EchoNumber') && ~isempty(V.EchoNumber); % multi-echo (dcm2niix)
  if TI > 1.5 && TR > 4
    suffix = 'FLAIR';
  elseif (TI >= 0.3 && TI <= 1.5) || (TR < 1 && TE < 0.03 && ~gre) || (gre && ~me && TR < 0.1 && TE < 0.01)
    suffix = 'T1w';
  elseif gre && (TE >= 0.015 || me)
    suffix = 'T2star';
  elseif TR > 2 && TE > 0.06
    suffix = 'T2w';
  elseif TR > 2 && TE < 0.03
    suffix = 'PDw';
  else
    suffix = '';
  end
end
% =========================================================================
function [suffix,ent,label] = setupSuffix(datatype,pro,V)
%setupSuffix. BIDS suffix, further entities, and acquisition label of a scan
%  by its datatype, protocol name pro, and JSON sidecar V:
%    anat  T1w, T2w, ... (by the protocol name, otherwise by the timing of
%          the sidecar, see anatSuffix), MP2RAGE with inv-<n> for 
%          the inversion images (e.g. INV1 or inv-1 in the protocol or 
%          series name), UNIT1 for the uniform image (UNI, or UNIFORM/UNI
%          in the ImageType), and T1map for the T1 map (T1 MAP in the 
%          ImageType, e.g. the T1_Images series of a Siemens MP2RAGE)
%    func  bold or sbref (SBRef in the protocol or series name)
%    dwi   dwi or sbref
%    fmap  epi for EPI field maps (e.g. spin-echo AP/PA) with dir-<label>
%          (see directionLabel), otherwise (e.g. GRE) magnitude1/magnitude2
%          (EchoNumber), phasediff, phase1/phase2 (phase images of the
%          echoes without EchoTime1/EchoTime2), or fieldmap (Units Hz)
%  Phase images (ImageType P or PHASE) get part-phase (real and imaginary
%  images part-real and part-imag) and multi-echo images (EchoNumber) of
%  anat and func echo-<n>. ent is the string of the entities in BIDS order
%  (e.g. "_dir-AP_part-phase"). The label is the protocol name without the
%  tokens of the information of the suffix and entities (e.g. SBRef, phase,
%  e1, AP, inv-1) and without trailing magnitude or phase tokens (e.g. 
%  "_phase_mag"), i.e. the images of one acquisition (SBRef and bold, 
%  field-map magnitude and phase, AP and PA, MP2RAGE inversions) get the 
%  same label (see cat_io_dcm2bids).
  if ~exist('V','var') || ~isstruct(V), V = struct(); end
  ser   = ''; if isfield(V,'SeriesDescription'), ser = char(V.SeriesDescription); end
  itype = {};
  for fn = {'ImageType','ImageTypeText'}
    if isfield(V,fn{1}), itype = [itype; upper(reshape(cellstr(V.(fn{1})),[],1))]; end %#ok<AGROW>
  end
  echo  = [];
  if isfield(V,'EchoNumber') && isnumeric(V.EchoNumber) && isscalar(V.EchoNumber), echo = V.EchoNumber; end
  if any(ismember(itype,{'P','PHASE'})),  part = 'phase';
  elseif any(strcmp(itype,'REAL')),       part = 'real';
  elseif any(strcmp(itype,'IMAGINARY')),  part = 'imag';
  else,                                   part = '';
  end
  partok = struct('phase',{{'phase','pha','ph'}},'real',{{'real'}},'imag',{{'imag','imaginary'}});
  issbref = cat_io_contains(lower(pro),{'sbref'}) || cat_io_contains(lower(ser),{'sbref'});
  suffix  = ''; ent = struct(); tok = {};

  if strcmp( datatype , 'anat')
    inv = regexpi([pro ' ' ser],'(?:^|[-_ ])inv[-_ ]?(\d)(?=$|[-_ ])','tokens','once'); 
    its = [' ' strjoin(reshape(itype,1,[]),' ') ' ']; % (e.g. T1, MAP and T1 MAP of dcm2niix)
    if ~isempty(inv) % inversion images of MP2RAGE
      suffix  = 'MP2RAGE'; 
      ent.inv = str2double(inv{1}); 
      tok     = {sprintf('inv%d',ent.inv)}; 
    elseif contains(its,{' T1 MAP ',' T1_MAP '}) % T1 map (e.g. of MP2RAGE)
      suffix = 'T1map'; 
    elseif contains(its,{' UNIFORM ',' UNI '}) || ...
        (~isempty(regexpi([pro ' ' ser],'(?:^|[-_ ])uni(?=$|[-_ ])','once')) && ...
        ~cat_io_contains( lower(pro) , {'t1w'})) % uniform image of MP2RAGE
      suffix = 'UNIT1'; 
    elseif cat_io_contains( lower(pro) , {'t1map'})
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
    else
      suffix = anatSuffix(V); % by the timing of the sidecar
    end
    if any(strcmp(suffix,{'T1w','T2w','PDw','T2star','FLAIR','PDT2','MP2RAGE'})) % no maps
      if ~isempty(echo), ent.echo = echo; end
      if ~isempty(part), ent.part = part; end
    end

  elseif any(strcmp( datatype , {'func','dwi'}))
    if issbref
      suffix = 'sbref'; tok = {'sbref'};
    elseif strcmp( datatype , 'func')
      suffix = 'bold';
    else
      suffix = 'dwi';
    end
    if ~isempty(echo) && strcmp( datatype , 'func'), ent.echo = echo; end
    if ~isempty(part), ent.part = part; end

  elseif strcmp( datatype , 'fmap')
    seq = '';
    if isfield(V,'ScanningSequence'), seq = upper(strjoin(cellstr(V.ScanningSequence),'\')); end
    if any(strcmp(strsplit(seq,{'\','_'}),'EP')) || (isempty(seq) && isfield(V,'TotalReadoutTime'))
      % EPI field maps with opposite phase-encoding directions
      suffix  = 'epi';
      ent.dir = directionLabel(pro,ser,V);
      tok     = {ent.dir};
      if ~isempty(part), ent.part = part; end
    elseif isfield(V,'Units') && strcmpi(char(V.Units),'Hz')
      suffix = 'fieldmap';
    elseif strcmp(part,'phase')
      % GRE phase difference (with EchoTime1 and EchoTime2) or phase images
      if ~isempty(echo) && ~(isfield(V,'EchoTime1') && isfield(V,'EchoTime2'))
        suffix = sprintf('phase%d',echo);
      else
        suffix = 'phasediff';
      end
    else
      % GRE magnitude images of the echoes
      if isempty(echo), echo = 1; end
      suffix = sprintf('magnitude%d',echo);
      tok    = {'mag','magnitude'};
    end
    if ~strcmp(suffix,'epi') && ~isempty(echo), tok = [tok {sprintf('e%d',echo)}]; end
  end
  if isfield(ent,'echo'), tok = [tok {sprintf('e%d',ent.echo)}]; end
  if isfield(ent,'part') || any(strcmp(suffix,{'phasediff','phase1','phase2'}))
    tok = [tok partok.(part)];
  end

  % entities in BIDS order and the label without the tokens of their 
  % information and without trailing tokens of the image part (e.g. 
  % "_phase_mag" of protocols with magnitude and phase images)
  ent   = entityString(ent);
  parts = regexp(regexprep(pro,'(?i)(inv)[-_ ](\d)','$1$2'),'[-_ ]+','split');
  parts = parts(~ismember(lower(parts),lower(tok)));
  while numel(parts) > 1 && any(strcmpi(parts{end},{'mag','magnitude','phase','pha','ph'}))
    parts(end) = [];
  end
  label = strjoin(parts,'-');
end
% =========================================================================
function d = directionLabel(pro,ser,V)
%directionLabel. Label of the dir entity of an EPI field map, i.e., the
%  direction token of the protocol or series name (AP, PA, LR, RL, SI, IS),
%  otherwise the direction of the PhaseEncodingDirection for the axes of
%  the dcm2niix output of axial or coronal slices (i: R>L, j: P>A, k: I>S).
%  BIDS requires the PhaseEncodingDirection for the dir entity, i.e., the
%  label is empty without it.
  d = '';
  if ~isfield(V,'PhaseEncodingDirection') || isempty(V.PhaseEncodingDirection), return; end
  t = regexpi([pro ' ' ser],'(?:^|[_ -])(AP|PA|LR|RL|SI|IS)(?=$|[_ -])','tokens','once');
  if ~isempty(t)
    d = upper(t{1});
  else
    map = {'i','RL'; 'i-','LR'; 'j','PA'; 'j-','AP'; 'k','IS'; 'k-','SI'};
    d   = map(strcmp(map(:,1),char(V.PhaseEncodingDirection)),2);
    if isempty(d), d = ''; else, d = d{1}; end
  end
end
% =========================================================================
function str = entityString(E)
%entityString. File name part of the entities in E (structure with the
%  keys of the entities as fields, e.g. imported BIDS names) that follow the
%  acq entity, in BIDS order, e.g. "_dir-AP_echo-1_part-phase". Other keys
%  (sub, ses, task, acq, run) are set by cat_io_dcm2bids.
  keys = {'ce','rec','dir','mod','echo','flip','inv','mt','part','chunk'};
  str  = '';
  for ki = 1:numel(keys)
    if isfield(E,keys{ki}) && ~isempty(E.(keys{ki}))
      v = E.(keys{ki});
      if isnumeric(v), v = sprintf('%d',v); end
      str = sprintf('%s_%s-%s', str, keys{ki}, char(v));
    end
  end
end
% =========================================================================
function task = setupTask(datatype,pro)
  if strcmp( datatype , 'func')
    if cat_io_contains( lower(pro) , {'rs','rest'})
      task = sprintf('_task-rest');
    elseif cat_io_contains( lower(pro) , {'checker'})                        % visual checkerboard stimulation (e.g. NKI-RS)
      task = sprintf('_task-checkerboard');
    elseif cat_io_contains( lower(pro) , {'breath'})                         % breath holding (cerebrovascular reactivity)
      task = sprintf('_task-breathhold');
    elseif cat_io_contains( lower(pro) , {'eye_mov','eyemov','eye-mov','eye mov','saccad','eyecal'}) % eye movement calibration
      task = sprintf('_task-eyemovement');
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
function writeSubjectTSV(Vjson,sub,ses,Poutdir,BIDSsubdir,anon,SOR,site,extra)
%writeSubjectTSV. BIDS sessions file of a subject (sub-<label>_sessions.tsv) 
%  with one row per session: site (abbreviation, see setupSites in 
%  cat_io_dcm2bids_protocols), age, weight, the session orientation rating
%  SOR (see writeSessionAffines in cat_io_dcm2bids_pp_anat), and the 
%  phenotypic data of the session (structure extra, see 
%  cat_io_dcm2bids_phenotypic). The rows of other sessions are kept (by 
%  column name, i.e., also from files with other columns, including their
%  phenotypic data). Missing values are n/a. 
  if ~exist('SOR','var') || isempty(SOR), SOR = nan; end
  if ~exist('site','var') || isempty(site), site = 'n/a'; end
  if ~exist('extra','var') || ~isstruct(extra), extra = struct(); end
  Psubject = fullfile(Poutdir,BIDSsubdir,sub,sprintf('%s_sessions.tsv',sub)); % BIDS sessions file
  hdr = {'session_id','site','age','weight','SOR'}; T0 = {};
  if exist(Psubject,'file')
    T0 = cat_io_csv(Psubject, '','', struct('delimiter','\t','convert2double',-1)); 
  end
  ext = fieldnames(extra)'; old = {};
  if ~isempty(T0), old = T0(1,~ismember(T0(1,:),hdr)); end
  hdr = [hdr old setdiff(ext,old,'stable')]; % phenotypic columns of all sessions
  Tsubject = hdr; 
  for ri = 2:size(T0,1)
    if strcmp(T0{ri,1},ses), continue; end % replaced by the current session
    [isc,ci] = ismember(hdr,T0(1,:)); 
    Tsubject(end+1,:) = {'n/a'}; Tsubject(end,isc) = T0(ri,ci(isc)); %#ok<AGROW>
  end
  row = repmat({'n/a'},1,numel(hdr)); 
  row(1:5) = {ses, site, round(Vjson.PatientAge,2-anon), round(Vjson.PatientWeight), round(SOR,2)}; 
  for e = ext, row{strcmp(hdr,e{1})} = extra.(e{1}); end
  Tsubject(end+1,:) = row; 
  Tsubject(cellfun(@(v) (isnumeric(v) && isscalar(v) && isnan(v)) || (ischar(v) && isempty(v)), Tsubject)) = {'n/a'}; 
  [~,so] = sort(Tsubject(2:end,1)); Tsubject(2:end,:) = Tsubject(so+1,:); 
  
  % write
  cat_io_csv(Psubject,Tsubject,'','',struct('delimiter','\t')); 
end
% =========================================================================
function writeParticipantTSV(P,Pparticipants,anon)
% Create participant file (eg. OpenNeuro) with sex and age (sex first as it
% is more robust, the session-specific age is in the subject files) and 
% further columns of imported BIDS datasets (extra.columns/extra.values). 
% The age is the minimum age of the subject, i.e. also of the sessions of 
% previous runs (the age of the existing row is kept if it is lower). 
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
    age = round(P(pi).V.PatientAge,2-anon); 
    if pidpa <= size(Tparticipants,1) && isnumeric(Tparticipants{pidpa,3}) && isscalar(Tparticipants{pidpa,3}) && ...
        ~isnan(Tparticipants{pidpa,3}) && ~(Tparticipants{pidpa,3} >= age) % (also for a new age NaN)
      age = Tparticipants{pidpa,3}; % minimum age of all runs
    end
    Tparticipants(pidpa,1:3) = {P(pi).sub, P(pi).V.PatientSex, age };
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
%  text file, extended by the README files of imported BIDS datasets), 
%  participants.json (description of the columns of participants.tsv), 
%  .bidsignore (files of catDCM2BIDS that are not part of BIDS), and the
%  dataset_description.json of the derivatives of catDCM2BIDS (BIDS 
%  derivatives dataset, see derivativesDir in cat_io_dcm2bids_helper). 
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
    'Description','DICOM/BIDS import, protocol evaluation, anonymization, and QC of CAT', ...
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
    txt = sprintf(['# %s\n\nBIDS dataset created by catDCM2BIDS (CAT), including the evaluation ' ...
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

  % .bidsignore: settings and readme of the output and the protocol 
  % differences of non-conform scans (matchfiles, see cat_io_dcm2bids)
  cat_io_dcm2bids_helper('writeText',fullfile(Pbids,'.bidsignore'), ...
    sprintf('catDCM2BIDSdb.mat\nreadme.txt\n**/matchfiles/\n'));   

  % dataset_description.json of the derivatives (QC, segmentation, session
  % registration, DTI, and rs-fMRI measures and maps of catDCM2BIDS)
  Pd = cat_io_dcm2bids_helper('derivativesDir',Pbids); 
  if exist(Pd,'dir')
    Dd = struct('Name',[D.Name ' (catDCM2BIDS derivatives)'],'BIDSVersion',D.BIDSVersion,'DatasetType','derivative', ...
      'GeneratedBy',{{struct('Name','catDCM2BIDS','Version',char(catver), ...
      'Description','QC, segmentation, session registration, DTI, and rs-fMRI measures and maps of CAT', ...
      'CodeURL','https://github.com/ChristianGaser/cat12')}}); 
    if ~isempty(lic), Dd.License = lic; end
    cat_io_dcm2bids_helper('writeJSON',fullfile(Pd,'dataset_description.json'), Dd); 
  end
end
% =========================================================================
function writeProtocols(Pbids,Pp,Ptbldir)
%writeProtocols. Generalized MRI protocols of the scans of the BIDS 
%  directory Pbids (with its protocol subdirectories), i.e. the protocol
%  variants of the database tables Ptbldir (normalized protocol fields 
%  without personal fields, see protocolVariant in 
%  cat_io_dcm2bids_protocols) of the scans of this BIDS directory (by the 
%  datatype and BIDS name without subject and session, i.e. also for other
%  subject labels of the anonymization) in the report directory Pp, e.g.
%    <BIDSdir>-report/protocols/<BIDSsubdir>/<datatype>_<variant>.json
%  with the number of scans of each variant (protocols.tsv), created again
%  in each run (the old protocols are kept without the database tables).
  Pold = fullfile(Pbids,'protocols'); % (location of a former version)
  if exist(fullfile(Pold,'protocols.tsv'),'file'), rmdir(Pold,'s'); end
  Pmat = fullfile(Ptbldir,'scans.mat'); if ~exist(Pmat,'file'), return; end
  S = load(Pmat,'R'); R = S.R; if isempty(R) || ~all(isfield(R,{'BIDSfile','Variant','Datatype'})), return; end
  if exist(Pp,'dir'), rmdir(Pp,'s'); end
  nos = @(f) regexprep(f,'^sub-[^_]+_(ses-[^_]+_)?|\.json$','');  % BIDS name without subject, session, and extension
  P   = cat_vol_findfiles(Pbids,'sub-*.json'); P = P(~contains(P,[filesep 'derivatives' filesep]));
  if isempty(P), return; end
  [~,dt] = cellfun(@(p) fileparts(fileparts(p)),P,'UniformOutput',false);    % datatype directories
  kb  = strcat(dt(:),'/',cellfun(@(p) nos(spm_file(p,'filename')),P(:),'UniformOutput',false));
  kr  = strcat({R.Datatype}','/',cellfun(nos,{R.BIDSfile}','UniformOutput',false));
  use = ismember(kr,kb) & ~cellfun('isempty',{R.Variant}');
  if ~any(use), return; end
  [uv,~,ui] = unique(strcat({R(use).Datatype}','/',{R(use).Variant}'));
  T = {'file','datatype','variant','scans'};
  for vi = 1:numel(uv)
    d = fileparts(uv{vi}); v = uv{vi}(numel(d)+2:end);
    Ps = fullfile(Ptbldir,'protocols',d,[v '.json']);
    if ~exist(Ps,'file'), continue; end
    if ~exist(Pp,'dir'), mkdir(Pp); end
    copyfile(Ps,fullfile(Pp,[d '_' v '.json']));
    T(end+1,:) = {[d '_' v '.json'], d, v, nnz(ui == vi)}; %#ok<AGROW>
  end
  cat_io_dcm2bids_helper('writeTable',fullfile(Pp,'protocols.tsv'),T,'tsv');
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
%  in the BIDS directories (BIDS*) of the output directory.
  D = dir(fullfile(Poutdir,'BIDS*')); D = D([D.isdir]); 
  BIDSsubdirs = fullfile(Poutdir,{D.name})'; 
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
