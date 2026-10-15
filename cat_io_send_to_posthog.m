function cat_io_send_to_posthog(event, props)
% ______________________________________________________________________
%
% Send an anonymous event to PostHog (EU cloud) for CAT user statistics.
% This is sent in addition to the Piwik information of
% cat_io_send_to_server and is also disabled by the flag
% cat.extopts.send_info in cat_defaults.m.
%
%   cat_io_send_to_posthog(event[, props]);
%
%   event  .. event name, e.g. 'cat12_start', 'cat12_run', 'cat12_error'
%   props  .. struct with event properties (char, numeric or logical);
%             a field 'stack' with an MException stack is converted to
%             'function:line | ...' (innermost first, without paths) and
%             its first entry is added as 'error_location'
%
%   Events with a stack are additionally sent as '$exception' event for
%   PostHog error tracking, with error_id and error_message as exception
%   type and message and the stack frames (function, file name without
%   path, line).
%
%   All events contain the properties cat_release, cat_revision, system,
%   standalone and matlab_release. A random ID per MATLAB session is used
%   as distinct_id and no person profiles are created. Any error (e.g. no
%   network, old MATLAB or Octave without webwrite/jsonencode) is ignored
%   to never interrupt CAT.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  persistent session_id

  % public write-only token of PostHog project 306100 (EU cloud)
  token = 'phc_qnZB2BmwhjeFk3c6yyERFtQHYBfRkkjrp42jAy36VpxX';
  url   = 'https://eu.i.posthog.com/batch/';

  try
    if ~cat_get_defaults('extopts.send_info'), return; end

    % anonymous ID per MATLAB session; tempname is used instead of rand/randi
    % to avoid changing the user's random number generator state
    if isempty(session_id)
      [~, session_id] = fileparts(tempname);
    end

    if nargin < 2 || isempty(props), props = struct; end

    % only send function names, file names and line numbers of the stack
    % (no paths) and add the filename for subfunctions (e.g. 'cat_run>run_job:277')
    exception = [];
    if isfield(props,'stack') && isstruct(props.stack)
      catdir = [fileparts(mfilename('fullpath')) filesep];
      stack  = props.stack;
      str    = cell(1,numel(stack));
      frames = cell(1,numel(stack));
      for si = 1:numel(stack)
        [~, fname, ext] = fileparts(stack(si).file);
        if isempty(fname) || strcmp(fname,stack(si).name)
          str{si} = sprintf('%s:%d',stack(si).name,stack(si).line);
        else
          str{si} = sprintf('%s>%s:%d',fname,stack(si).name,stack(si).line);
        end
        frames{si} = struct('platform', 'custom', 'lang', 'matlab', ...
          'function_', stack(si).name, 'filename', [fname ext], 'lineno', stack(si).line, ...
          'resolved', true, 'in_app', strncmp(stack(si).file,catdir,numel(catdir)));
      end
      if isempty(str), str = {''}; end
      if ~isfield(props,'error_location'), props.error_location = str{1}; end
      props.stack = strjoin(str,' | ');

      % for error tracking only use the frames up to the CAT error handler,
      % because the calling frames depend on how CAT was started (GUI, batch,
      % standalone) and would split the same error into different issues;
      % PostHog expects the frames ordered from the entry point to the error
      ind = find(strcmp({stack.name},'cat_run_newcatch'),1);
      if isempty(ind), ind = numel(stack); end
      exception = struct('type', 'MException', 'value', '', ...
        'mechanism', struct('handled', true, 'synthetic', false), ...
        'stacktrace', struct('type', 'raw', 'frames', {frames(ind:-1:1)}));
      if isfield(props,'error_id') && ~isempty(props.error_id), exception.type = props.error_id; end
      if isfield(props,'error_message'), exception.value = props.error_message; end
    end

    % CAT version and system to differentiate all events
    [CATrel, CATver]     = cat_version;
    props.cat_release    = CATrel;
    props.cat_revision   = ['r' CATver];
    props.system         = computer;
    props.standalone     = isdeployed;
    props.matlab_release = version('-release');

    % anonymous event without person profile
    props.x_process_person_profile = false;

    batch = {struct('event', event, 'distinct_id', session_id, 'properties', props)};
    if ~isempty(exception)
      props.x_exception_list = {exception};
      batch{2} = struct('event', '$exception', 'distinct_id', session_id, 'properties', props);
    end
    data = struct('api_key', token, 'batch', {batch});

    % MATLAB field names cannot start with '$' (x_name -> $name) and we
    % avoid the keyword function as field name (function_ -> function)
    json = regexprep(jsonencode(data), {'"x_(\w+)":', '"function_":'}, {'"\$$1":', '"function":'});

    webwrite(url, json, weboptions('MediaType', 'application/json', 'Timeout', 2));
  catch %#ok<CTCH>
    % ignore all errors
  end
end
