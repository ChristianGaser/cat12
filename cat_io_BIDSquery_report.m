function cat_io_BIDSquery_report(out, verb)
%cat_io_BIDSquery_report. Print the result of a cat_io_BIDSquery call.
%
% ADD-ON (not part of upstream CAT12). See cat_io_BIDSquery.
%
%   cat_io_BIDSquery_report(out [, verb])
%
%   out   .. output structure of cat_io_BIDSquery
%   verb  .. 1-summary (default), 2-list every subject
%
% Reports which subjects were found, how the cross-sectional /
% longitudinal split was made, and - importantly - what was *not* found,
% so that a silently empty or half-empty selection cannot go unnoticed.
%
% See also cat_io_BIDSquery
% ______________________________________________________________________
%
% BIDS input add-on for CAT12 (fork).
% ______________________________________________________________________

  if nargin < 2, verb = 1; end
  m = out.missing;

  fprintf('\nCAT12 BIDS input query\n');
  fprintf('  dataset:      %s\n', out.bidsdir);
  fprintf('  suffix:       %s\n', strjoin(cellstr(out.opt.suffix),', '));
  fprintf('  files found:  %d  (%d subjects)\n', numel(out.files), numel(out.subjects));
  fprintf('  cross-sect.:  %d files  (%d subjects with <=1 session)\n', ...
    out.ncross, sum(out.nses <= 1));
  fprintf('  longitudinal: %d files  (%d subjects with >=2 sessions)\n', ...
    out.nlong, numel(out.longsubjects));

  if verb > 1 && ~isempty(out.subjects)
    fprintf('\n  subjects:\n');
    for si = 1:numel(out.subjects)
      if out.nses(si) <= 1, kind = 'cross'; else, kind = 'long '; end
      fprintf('    %-20s %s  %d session(s)\n', out.subjects{si}, kind, out.nses(si));
    end
  end

  % -- warnings ------------------------------------------------------
  warned = false;

  if m.no_subject_dirs
    warned = warn(['No sub-* directories found. Is "%s" really the root of a raw ' ...
      'BIDS dataset (the directory containing the sub-* folders)?\n'], out.bidsdir);
  end
  if m.no_dataset_description
    warned = warn(['No dataset_description.json in the dataset root. The directory ' ...
      'may still work, but it is not a valid BIDS dataset.\n']);
  end
  if m.no_participants_tsv
    warned = warn(['A participants filter was requested but no participants.tsv ' ...
      'exists in the dataset root - the filter was ignored.\n']);
  end
  if isfield(m,'tsv_unreadable') && m.tsv_unreadable
    warned = warn(['participants.tsv could not be read or is empty - the ' ...
      'filter was ignored.\n']);
  end
  if ~isempty(m.tsv_unknown_column)
    warned = warn(['Column "%s" not found in participants.tsv - the filter was ' ...
      'ignored.\n'], m.tsv_unknown_column);
  end
  warned = listwarn(warned, m.tsv_filtered_out, ...
    'excluded by the participants.tsv filter');
  warned = listwarn(warned, m.no_anat_dir, ...
    'skipped: no anat directory found');
  warned = listwarn(warned, m.no_matching_file, ...
    'skipped: anat directory contains no file matching the requested suffix');
  warned = listwarn(warned, m.multiple_runs_single_session, ...
    ['has several files within a single session (runs/acquisitions) - ' ...
     'treated as cross-sectional, NOT longitudinal']);

  if ~warned && ~isempty(out.files)
    fprintf('  no missing data detected.\n');
  end
  fprintf('\n');
end
%==========================================================================
function warned = warn(varargin)
  cwrite(['  Warning: ' varargin{1}], varargin{2:end});
  warned = true;
end
%==========================================================================
function cwrite(varargin)
%cwrite. Coloured output, falling back to fprintf.
% cat_io_cprintf needs a Java command window and fails e.g. under Octave or
% with -nodisplay. A warning must never be the thing that aborts the query.

  try
    cat_io_cprintf('warn', varargin{:});
  catch
    fprintf(varargin{:});
  end
end
%==========================================================================
function warned = listwarn(warned, list, msg)
%listwarn. Report a list of subjects, truncated to keep the log readable.

  if isempty(list), return; end
  warned = warn('%d subject(s) %s:\n', numel(list), msg);
  nshow  = min(numel(list),10);
  for i = 1:nshow
    cwrite('    %s\n', list{i});
  end
  if numel(list) > nshow
    cwrite('    ... and %d more\n', numel(list) - nshow);
  end
end
