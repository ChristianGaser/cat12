function out = cat_io_BIDSquery(bidsdir, opt)
%cat_io_BIDSquery. Query a raw BIDS dataset for anatomical input files.
%
% ADD-ON (not part of upstream CAT12).
% ----------------------------------------------------------------------
% Upstream CAT12 is BIDS-aware only on the *output* side: cat_io_BIDS.m
% reverse-engineers sub-/ses-/anat entities from the paths of files that
% the user selected by hand, and uses them to place the derivatives.
% There is no way to point CAT12 at a BIDS dataset root and let it collect
% the input itself (cat_batch_bids.sh does this for exactly one subject
% with one hard-coded glob).
%
% This function adds that missing input layer. It is completely
% self-contained: it does not modify or depend on any upstream CAT12
% preprocessing code, and it emits ordinary absolute file paths, so the
% existing output-side BIDS logic in cat_io_BIDS.m keeps working unchanged.
% ----------------------------------------------------------------------
%
%   out = cat_io_BIDSquery(bidsdir, opt)
%
%   bidsdir           .. root of the raw BIDS dataset (the directory that
%                        contains the sub-* directories)
%
%   opt               .. query options (all optional)
%    .suffix          .. BIDS suffix / modality, char or cellstr
%                        (default 'T1w'; e.g. {'T1w','T2w'})
%    .sub             .. participant filter, cellstr of labels with or
%                        without the 'sub-' prefix ({} = all)
%    .ses             .. session filter, cellstr of labels with or without
%                        the 'ses-' prefix ({} = all)
%    .acq             .. regexp on the acq- entity ('' = any)
%    .run             .. regexp on the run- entity ('' = any)
%    .tsvfilter       .. participants.tsv filter, see below ('' = off)
%    .verb            .. 0-quiet, 1-summary (default), 2-detailed
%
%   .tsvfilter is a char expression of the form  "column op value", e.g.
%       'group == patient'    'age >= 18'    'sex ~= M'
%   Supported operators: == ~= != > >= < <= . Numeric comparison is used
%   when both sides parse as numbers, otherwise a case-insensitive string
%   comparison. Subjects missing from participants.tsv are reported and
%   excluded.
%
%   out               .. result structure
%    .bidsdir         .. resolved dataset root
%    .files           .. all matching files (cellstr, absolute paths)
%    .sub / .ses      .. per-file subject / session label
%    .subjects        .. unique subject labels (in file order)
%    .nses            .. number of distinct sessions per subject
%    .cross           .. files of subjects with <= 1 session
%                        -> feed into the cross-sectional pipeline
%    .long            .. files of subjects with >= 2 sessions
%    .longsubjects    .. 1xS cell, one cellstr of files per longitudinal
%                        subject -> feed into the longitudinal pipeline
%    .longsublabels   .. labels belonging to .longsubjects
%    .missing         .. report of things worth telling the user, see
%                        cat_io_BIDSquery_report
%
% Example:
%   out = cat_io_BIDSquery('/data/ds001', struct('suffix','T1w'));
%   matlabbatch{1}.spm.tools.cat.estwrite.data = out.cross;
%
% See also cat_io_BIDS, cat_conf_BIDSinput, cat_vol_findfiles
% ______________________________________________________________________
%
% BIDS input add-on for CAT12 (fork).
% Based on the CAT12 toolbox by Christian Gaser and Robert Dahnke,
% Structural Brain Mapping Group (https://neuro-jena.github.io).
% ______________________________________________________________________

  %#ok<*AGROW>

  if nargin < 1 || isempty(bidsdir)
    bidsdir = spm_select(1,'dir','Select the BIDS dataset root (contains sub-*)');
  end
  if nargin < 2, opt = struct(); end

  % -- defaults -------------------------------------------------------
  def.suffix    = {'T1w'};
  def.sub       = {};
  def.ses       = {};
  def.acq       = '';
  def.run       = '';
  def.tsvfilter = '';
  def.verb      = 1;
  opt = cat_io_checkinopt(opt,def);

  opt.suffix = cellstr(opt.suffix);
  opt.sub    = normalizeLabels(opt.sub,'sub-');
  opt.ses    = normalizeLabels(opt.ses,'ses-');

  bidsdir = char(bidsdir);
  if ~isempty(bidsdir) && bidsdir(end) == filesep, bidsdir(end) = []; end

  out          = struct();
  out.bidsdir  = bidsdir;
  out.opt      = opt;
  out.missing  = emptyReport();

  if ~exist(bidsdir,'dir')
    error('cat_io_BIDSquery:noDir','BIDS directory "%s" does not exist.',bidsdir);
  end

  % A raw BIDS root must contain sub-* directories. We only warn rather
  % than fail, because valid-but-unusual layouts should still be usable.
  if ~exist(fullfile(bidsdir,'dataset_description.json'),'file')
    out.missing.no_dataset_description = true;
  end

  subdirs = cat_vol_findfiles(bidsdir,'sub-*',struct('dirs',1,'depth',1));
  if isempty(subdirs)
    out.missing.no_subject_dirs = true;
    out = finalize(out,{},opt);
    cat_io_BIDSquery_report(out);
    return
  end

  % -- participants.tsv filter ----------------------------------------
  keepsub = {}; usetsv = false;
  if ~isempty(opt.tsvfilter)
    ptsv = fullfile(bidsdir,'participants.tsv');
    if ~exist(ptsv,'file')
      out.missing.no_participants_tsv = true;
    else
      [keepsub, tsvmiss] = filterParticipants(ptsv, opt.tsvfilter);
      out.missing.tsv_unknown_column = tsvmiss.unknown_column;
      out.missing.tsv_excluded       = tsvmiss.excluded;
      % An unreadable table or a missing column must not silently exclude
      % every subject - in that case the filter is genuinely ignored, which
      % is what the report tells the user.
      out.missing.tsv_unreadable     = tsvmiss.unreadable;
      usetsv = isempty(tsvmiss.unknown_column) && ~tsvmiss.unreadable;
    end
  end

  % -- collect files ---------------------------------------------------
  files = {};
  for si = 1:numel(subdirs)
    [~,sublabel] = fileparts(subdirs{si});

    if ~isempty(opt.sub) && ~any(strcmpi(sublabel,opt.sub)), continue; end
    if usetsv && ~any(strcmpi(sublabel,keepsub))
      out.missing.tsv_filtered_out{end+1} = sublabel; continue;
    end

    % anat dirs live either in sub-*/anat or sub-*/ses-*/anat
    anatdirs = cat_vol_findfiles(subdirs{si},'anat',struct('dirs',1,'maxdepth',2));
    if isempty(anatdirs)
      out.missing.no_anat_dir{end+1} = sublabel; continue;
    end

    subfiles = {};
    for ai = 1:numel(anatdirs)
      for xi = 1:numel(opt.suffix)
        % .nii and .nii.gz are both valid BIDS; CAT12 handles both.
        subfiles = [subfiles; cat_vol_findfiles(anatdirs{ai}, ...
          sprintf('*_%s.nii',opt.suffix{xi}), struct('depth',1))];
        subfiles = [subfiles; cat_vol_findfiles(anatdirs{ai}, ...
          sprintf('*_%s.nii.gz',opt.suffix{xi}), struct('depth',1))];
      end
    end

    if isempty(subfiles)
      out.missing.no_matching_file{end+1} = sublabel; continue;
    end

    files = [files; subfiles];
  end

  % -- entity filters (ses/acq/run) ------------------------------------
  keep = true(numel(files),1);
  for fi = 1:numel(files)
    ent = parseEntities(files{fi});
    if ~isempty(opt.ses)
      if isempty(ent.ses) || ~any(strcmpi(['ses-' ent.ses],opt.ses)), keep(fi) = false; end
    end
    if ~isempty(opt.acq) && isempty(regexpi(ent.acq,opt.acq,'once')), keep(fi) = false; end
    if ~isempty(opt.run) && isempty(regexpi(ent.run,opt.run,'once')), keep(fi) = false; end
  end
  files = files(keep);

  out = finalize(out,files,opt);

  if opt.verb, cat_io_BIDSquery_report(out,opt.verb); end
end
%==========================================================================
function out = finalize(out,files,opt)
%finalize. Derive per-file entities and the cross/long split.

  % Natural sort: the longitudinal pipeline relies on the timepoint order,
  % and a plain alphabetical sort would put ses-10 before ses-2.
  files = natsortunique(files(:));
  out.files = files;
  out.sub   = cell(numel(files),1);
  out.ses   = cell(numel(files),1);

  for fi = 1:numel(files)
    ent = parseEntities(files{fi});
    % keep the full BIDS labels ('sub-01', not '01') so that reports and
    % dependency names match what the user sees on disk
    out.sub{fi} = ['sub-' ent.sub];
    out.ses{fi} = ent.ses;
  end

  % Session-less data counts as a single session so that a dataset without
  % any ses-* level is treated as cross-sectional rather than as an error.
  ses = out.ses;
  ses(cellfun(@isempty,ses)) = {'_nosession_'};

  [out.subjects,~,subidx] = unique(out.sub,'stable');
  out.nses          = zeros(numel(out.subjects),1);
  out.cross         = {};
  out.long          = {};
  out.longsubjects  = {};
  out.longsublabels = {};

  for si = 1:numel(out.subjects)
    sel  = subidx == si;
    nses = numel(unique(ses(sel)));
    out.nses(si) = nses;

    if nses <= 1
      % <=1 session -> cross-sectional. Note that several runs within one
      % session are still cross-sectional; the longitudinal pipeline needs
      % genuinely separate sessions.
      out.cross = [out.cross; out.files(sel)];
      if sum(sel) > 1
        out.missing.multiple_runs_single_session{end+1} = out.subjects{si};
      end
    else
      out.long                    = [out.long; out.files(sel)];
      out.longsubjects{end+1}     = out.files(sel);
      out.longsublabels{end+1}    = out.subjects{si};
    end
  end

  out.ncross = numel(out.cross);
  out.nlong  = numel(out.long);
end
%==========================================================================
function files = natsortunique(files)
%natsortunique. unique() with a natural (human) ordering of digit groups.
% Every digit group is zero-padded in a temporary sort key only, so that
% ses-2 sorts before ses-10 as a user would expect.

  files = unique(files);
  if numel(files) < 2, return; end

  % Built without the ${...} dynamic regexprep expression, which MATLAB
  % supports but Octave does not.
  keys = cell(size(files));
  for i = 1:numel(files)
    s = files{i};
    [tok,startix,endix] = regexp(s,'\d+','match','start','end');
    k = ''; last = 0;
    for t = 1:numel(tok)
      k = [k s(last+1:startix(t)-1) sprintf('%012d',str2double(tok{t}))]; %#ok<AGROW>
      last = endix(t);
    end
    keys{i} = [k s(last+1:end)];
  end
  [~,idx] = sort(keys);
  files = files(idx);
end
%==========================================================================
function ent = parseEntities(file)
%parseEntities. Extract BIDS entities from a filename and its path.
% The filename is authoritative; the path is only a fallback, because a
% BIDS filename must carry sub- and ses- itself.

  ent = struct('sub','','ses','','acq','','run','','suffix','');

  [pth,nam] = fileparts(file);
  [~,nam2]  = fileparts(nam);           % strip the second ext of .nii.gz
  if ~isempty(nam2), nam = nam2; end

  parts = strsplit(nam,'_');
  for pi = 1:numel(parts)
    kv = strsplit(parts{pi},'-');
    if numel(kv) >= 2
      key = lower(kv{1});
      val = strjoin(kv(2:end),'-');
      if isfield(ent,key), ent.(key) = val; end
    elseif pi == numel(parts)
      ent.suffix = parts{pi};            % trailing token without '-' = suffix
    end
  end

  % fall back to the directory names if the filename is incomplete
  if isempty(ent.sub) || isempty(ent.ses)
    sdirs = strsplit(pth,filesep);
    for di = numel(sdirs):-1:1
      d = sdirs{di};
      if isempty(ent.ses) && strncmpi(d,'ses-',4), ent.ses = d(5:end); end
      if isempty(ent.sub) && strncmpi(d,'sub-',4), ent.sub = d(5:end); end
    end
  end
end
%==========================================================================
function [keep, miss] = filterParticipants(ptsv, expr)
%filterParticipants. Apply a "column op value" filter to participants.tsv.

  keep = {};
  miss = struct('unknown_column','','excluded',{{}},'unreadable',false);

  tok = regexp(strtrim(expr),'^(\S+)\s*(==|~=|!=|>=|<=|>|<)\s*(.+)$','tokens','once');
  if isempty(tok)
    error('cat_io_BIDSquery:badFilter', ...
      ['Cannot parse participants filter "%s". ' ...
       'Expected "column op value", e.g. "group == patient".'], expr);
  end
  [col, op, val] = deal(tok{1}, tok{2}, strtrim(tok{3}));
  if strcmp(op,'!='), op = '~='; end

  T = readTSV(ptsv);
  if isempty(T), miss.unreadable = true; return; end

  hdr    = T(1,:);
  idcol  = find(strcmpi(hdr,'participant_id'),1);
  valcol = find(strcmpi(hdr,col),1);
  if isempty(idcol),  miss.unknown_column = 'participant_id'; return; end
  if isempty(valcol), miss.unknown_column = col; return; end

  % == and ~= are always literal string comparisons, even when both sides
  % happen to parse as numbers. BIDS labels (session_id, run, site, etc.)
  % are strings by definition and a leading zero is part of the label, e.g.
  % session_id "01" must NOT equal "1". Only the ordering operators
  % (>,>=,<,<=), which have no string meaning, use numeric parsing.
  isOrder = any(strcmp(op,{'>','>=','<','<='}));
  refnum  = str2double(val);
  for ri = 2:size(T,1)
    id   = T{ri,idcol};
    cell_= T{ri,valcol};
    if isempty(id), continue; end

    if isOrder
      cellnum = str2double(cell_);
      if isnan(refnum) || isnan(cellnum)
        ok = false;   % ordering is undefined for non-numeric values
      else
        ok = compareNum(cellnum, op, refnum);
      end
    else
      ok = compareStr(cell_, op, val);
    end

    if ok
      keep{end+1} = id;
    else
      miss.excluded{end+1} = id;
    end
  end
end
%==========================================================================
function ok = compareNum(a,op,b)
  switch op
    case '==', ok = a == b;   case '~=', ok = a ~= b;
    case '>',  ok = a >  b;   case '>=', ok = a >= b;
    case '<',  ok = a <  b;   case '<=', ok = a <= b;
    otherwise, ok = false;
  end
end
%==========================================================================
function ok = compareStr(a,op,b)
% Ordering operators are meaningless for strings -> treated as no match.
  switch op
    case '==', ok =  strcmpi(strtrim(a),b);
    case '~=', ok = ~strcmpi(strtrim(a),b);
    otherwise, ok = false;
  end
end
%==========================================================================
function T = readTSV(fname)
%readTSV. Minimal tab-separated reader returning a cell matrix.
% Deliberately dependency-free (cat_io_csv is comma-oriented) and tolerant
% of ragged rows, which are common in hand-edited participants.tsv files.

  T = {};
  fid = fopen(fname,'r');
  if fid < 0, return; end
  c = textscan(fid,'%s','Delimiter','\n','Whitespace','');
  fclose(fid);
  lines = c{1};
  lines = lines(~cellfun(@isempty,lines));
  if isempty(lines), return; end

  % Strip a UTF-8 BOM, which spreadsheet tools (e.g. Excel) commonly add
  % and which otherwise corrupts the "participant_id" header match.
  bom = char([239 187 191]);
  if strncmp(lines{1}, bom, numel(bom))
    lines{1} = lines{1}(numel(bom)+1:end);
  end

  rows = cellfun(@(l) strsplit(l,'\t','CollapseDelimiters',false), lines, ...
    'UniformOutput', false);
  ncol = max(cellfun(@numel,rows));
  T    = repmat({''},numel(rows),ncol);
  for ri = 1:numel(rows)
    T(ri,1:numel(rows{ri})) = rows{ri};
  end
end
%==========================================================================
function L = normalizeLabels(L,prefix)
%normalizeLabels. Accept '01', 'sub-01' and comma separated lists alike.

  if isempty(L), L = {}; return; end
  if ischar(L), L = strtrim(strsplit(L,',')); end
  L = cellstr(L);
  L = L(~cellfun(@isempty,L));
  for i = 1:numel(L)
    if ~strncmpi(L{i},prefix,numel(prefix))
      L{i} = [prefix L{i}];
    end
  end
end
%==========================================================================
function r = emptyReport()
  r = struct( ...
    'no_dataset_description',      false, ...
    'no_subject_dirs',             false, ...
    'no_participants_tsv',         false, ...
    'tsv_unknown_column',          '', ...
    'tsv_unreadable',              false, ...
    'tsv_excluded',                {{}}, ...
    'tsv_filtered_out',            {{}}, ...
    'no_anat_dir',                 {{}}, ...
    'no_matching_file',            {{}}, ...
    'multiple_runs_single_session',{{}} );
end
