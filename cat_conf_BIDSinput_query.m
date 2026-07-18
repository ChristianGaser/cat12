function out = cat_conf_BIDSinput_query(job, quiet)
%cat_conf_BIDSinput_query. Translate a batch job into a cat_io_BIDSquery call.
%
% ADD-ON (not part of upstream CAT12). See cat_conf_BIDSinput.
%
%   out = cat_conf_BIDSinput_query(job [, quiet])
%
% Shared by the batch prog and the batch vout. Because cfg_util calls vout
% on every user interaction, the result is cached per job configuration -
% otherwise editing a batch would rescan the whole dataset on every
% keystroke. Results are never cached across a changed configuration, so
% the cache cannot hide an edit.
%
% On any error an empty but complete result structure is returned, since
% vout must not throw while the batch is being edited (the dataset root can
% legitimately be undefined at that moment).
%
% See also cat_conf_BIDSinput, cat_io_BIDSquery
% ______________________________________________________________________
%
% BIDS input add-on for CAT12 (fork).
% ______________________________________________________________________

  persistent cachekey cacheval

  if nargin < 2, quiet = 0; end

  opt = struct( ...
    'suffix',    {getfieldd(job,'suffix',{'T1w'})}, ...
    'sub',       getfieldd(job,'sub',''), ...
    'ses',       getfieldd(job,'ses',''), ...
    'acq',       getfieldd(job,'acq',''), ...
    'run',       getfieldd(job,'run',''), ...
    'tsvfilter', getfieldd(job,'tsvfilter',''), ...
    'verb',      getfieldd(job,'verb',1) );

  bidsdir = getfieldd(job,'bidsdir','');
  if iscell(bidsdir)
    if isempty(bidsdir), bidsdir = ''; else, bidsdir = bidsdir{1}; end
  end

  if quiet, opt.verb = 0; end

  key = jobkey(bidsdir, opt);
  if quiet && ~isempty(cachekey) && strcmp(key, cachekey)
    out = cacheval;
    return
  end

  try
    out = cat_io_BIDSquery(bidsdir, opt);
  catch err
    if ~quiet, rethrow(err); end
    % vout must stay silent and non-fatal during batch editing
    out = struct('bidsdir',bidsdir,'opt',opt,'files',{{}},'sub',{{}},'ses',{{}}, ...
      'subjects',{{}},'nses',[],'cross',{{}},'long',{{}},'longsubjects',{{}}, ...
      'longsublabels',{{}},'ncross',0,'nlong',0,'missing',struct());
  end

  cachekey = key;
  cacheval = out;
end
%==========================================================================
function v = getfieldd(s, f, d)
  if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
%==========================================================================
function key = jobkey(bidsdir, opt)
%jobkey. Cache key covering every input that can change the result.
  key = sprintf('%s|%s|%s|%s|%s|%s|%s', bidsdir, ...
    strjoin(cellstr(opt.suffix),','), tostr(opt.sub), tostr(opt.ses), ...
    tostr(opt.acq), tostr(opt.run), tostr(opt.tsvfilter));
end
%==========================================================================
function s = tostr(v)
  if iscell(v), s = strjoin(cellstr(v),','); elseif ischar(v), s = v;
  else, s = mat2str(v); end
end
