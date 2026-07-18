function bidsinput = cat_conf_BIDSinput(expert)
%cat_conf_BIDSinput. Batch configuration for the BIDS input add-on.
%
% ADD-ON (not part of upstream CAT12).
% ----------------------------------------------------------------------
% Adds a batch module that takes a raw BIDS dataset root plus entity
% filters and produces the input file lists for the CAT12 pipelines, so
% that input selection no longer has to be done by hand.
%
% The module deliberately does NOT modify the data selectors of
% "CAT: Segmentation" or "CAT: Longitudinal". It is wired up through
% ordinary SPM batch dependencies instead, which keeps the diff against
% upstream CAT12 to a single registration line in tbx_cfg_cat.m.
%
% Dependencies produced:
%   "BIDS: cross-sectional files"     -> CAT: Segmentation / Volumes
%   "BIDS: longitudinal (Subj k)"     -> CAT: Longitudinal / Subject k
%   "BIDS: all files"                 -> anything else
%
% The cross-sectional / longitudinal split is automatic: subjects with
% zero or one session go to the cross-sectional list, subjects with two or
% more sessions go to the longitudinal lists.
%
% See also cat_io_BIDSquery, cat_io_BIDS, tbx_cfg_cat
% ______________________________________________________________________
%
% BIDS input add-on for CAT12 (fork).
% ______________________________________________________________________

  if nargin < 1, expert = 0; end

  % -- dataset root ---------------------------------------------------
  bidsdir            = cfg_files;
  bidsdir.tag        = 'bidsdir';
  bidsdir.name       = 'BIDS dataset root';
  bidsdir.filter     = 'dir';
  bidsdir.ufilter    = '.*';
  bidsdir.num        = [1 1];
  bidsdir.help       = {
    'Select the root directory of the raw BIDS dataset, i.e. the directory that directly contains the sub-* folders (and normally dataset_description.json and participants.tsv). '
    ''
    'Both sub-*/anat/ and sub-*/ses-*/anat/ layouts are supported, as are .nii and .nii.gz files. '
    ''};

  % -- entity filters -------------------------------------------------
  suffix             = cfg_menu;
  suffix.tag         = 'suffix';
  suffix.name        = 'Modality (BIDS suffix)';
  suffix.labels      = {'T1w','T2w','FLAIR','PDw','T1w + T2w'};
  suffix.values      = {{'T1w'},{'T2w'},{'FLAIR'},{'PDw'},{'T1w','T2w'}};
  suffix.val         = {{'T1w'}};
  suffix.help        = {
    'BIDS suffix of the images to process. CAT12 preprocessing is designed for T1w data; the other options are mainly useful for the tools modules. '
    ''};

  sub                = cfg_entry;
  sub.tag            = 'sub';
  sub.name           = 'Participants';
  sub.strtype        = 's';
  sub.num            = [0 Inf];
  sub.val            = {''};
  sub.help           = {
    'Comma separated list of participants to include, e.g. "01,02,05" or "sub-01,sub-02". Leave empty to use all participants. '
    ''};

  ses                = cfg_entry;
  ses.tag            = 'ses';
  ses.name           = 'Sessions';
  ses.strtype        = 's';
  ses.num            = [0 Inf];
  ses.val            = {''};
  ses.help           = {
    'Comma separated list of sessions to include, e.g. "01,02" or "ses-baseline". Leave empty to use all sessions. '
    ''
    'Note that restricting the sessions also changes the automatic cross-sectional/longitudinal decision, because that decision is based on the number of *selected* sessions per subject. '
    ''};

  acq                = cfg_entry;
  acq.tag            = 'acq';
  acq.name           = 'Acquisition filter (acq-)';
  acq.strtype        = 's';
  acq.num            = [0 Inf];
  acq.val            = {''};
  acq.hidden         = expert < 1;
  acq.help           = {
    'Regular expression applied to the acq- entity, e.g. "mprage". Leave empty to accept any (or no) acq- entity. '
    ''};

  run                = cfg_entry;
  run.tag            = 'run';
  run.name           = 'Run filter (run-)';
  run.strtype        = 's';
  run.num            = [0 Inf];
  run.val            = {''};
  run.hidden         = expert < 1;
  run.help           = {
    'Regular expression applied to the run- entity, e.g. "^1$" to keep only run-1. Leave empty to accept any (or no) run- entity. '
    ''};

  tsvfilter          = cfg_entry;
  tsvfilter.tag      = 'tsvfilter';
  tsvfilter.name     = 'participants.tsv filter';
  tsvfilter.strtype  = 's';
  tsvfilter.num      = [0 Inf];
  tsvfilter.val      = {''};
  tsvfilter.help     = {
    'Select participants by a column of participants.tsv, given as "column op value". Leave empty to disable. '
    ''
    'Examples:   group == patient   |   age >= 18   |   sex ~= M '
    ''
    'Supported operators are ==, ~= (or !=), >, >=, < and <=. Numeric comparison is used if both the column entry and the value are numbers, otherwise a case-insensitive string comparison is done. Ordering operators on non-numeric values never match. '
    ''
    'Participants that are missing from participants.tsv are excluded and reported. '
    ''};

  verb               = cfg_menu;
  verb.tag           = 'verb';
  verb.name          = 'Verbose output';
  verb.labels        = {'summary','summary + subject list'};
  verb.values        = {1,2};
  verb.val           = {1};
  verb.help          = {'Amount of information printed about found and missing data. ' ''};

  % -- main branch ----------------------------------------------------
  bidsinput          = cfg_exbranch;
  bidsinput.tag      = 'bidsinput';
  bidsinput.name     = 'BIDS: Data selection';
  bidsinput.val      = {bidsdir suffix sub ses acq run tsvfilter verb};
  bidsinput.prog     = @cat_conf_BIDSinput_run;
  bidsinput.vout     = @vout_BIDSinput;
  bidsinput.help     = {
    'Collect anatomical input files from a raw BIDS dataset and pass them on to the CAT12 pipelines via batch dependencies. '
    ''
    'This module replaces the manual selection of individual NIfTI files. It scans the dataset, applies the entity and participants.tsv filters, and then splits the result automatically: '
    ''
    '  * subjects with no or exactly one session become the "cross-sectional files" dependency, which is meant for "CAT: Segmentation"; '
    '  * subjects with two or more sessions become one "longitudinal (Subj k)" dependency each, which is meant for the "Subject" entries of "CAT: Longitudinal". '
    ''
    'Subjects that are skipped (no anat directory, no matching file, filtered out) are reported in the MATLAB command window so that an incomplete selection does not pass unnoticed. '
    ''
    'Note that the dependency list is built while the batch is being edited, which means the dataset is scanned at edit time. If the dataset changes on disk, re-select the dataset root to refresh the dependencies. '
    ''
    'The files handed over are ordinary absolute paths, so the existing CAT12 BIDS output logic (cat_io_BIDS.m, the "bids_folder" default) applies unchanged and derivatives are written to the usual derivatives/CAT* location. '
    ''};
end
%==========================================================================
function dep = vout_BIDSinput(job)
%vout_BIDSinput. Declare the batch dependencies.
% The number of longitudinal dependencies depends on the data, so the
% dataset has to be scanned already while the batch is edited. The result
% is cached because cfg_util calls vout on every UI update.

  out = cat_conf_BIDSinput_query(job, 1);

  dep             = cfg_dep;
  dep.sname       = sprintf('BIDS: cross-sectional files (%d)', numel(out.cross));
  dep.src_output  = substruct('.','cross');
  dep.tgt_spec    = cfg_findspec({{'filter','image','strtype','e'}});

  for k = 1:numel(out.longsubjects)
    cdep            = cfg_dep;
    cdep.sname      = sprintf('BIDS: longitudinal %s (%d files)', ...
                        out.longsublabels{k}, numel(out.longsubjects{k}));
    cdep.src_output = substruct('.','longsubjects','{}',{k});
    cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(end+1)      = cdep; %#ok<AGROW>
  end

  cdep            = cfg_dep;
  cdep.sname      = sprintf('BIDS: all files (%d)', numel(out.files));
  cdep.src_output = substruct('.','files');
  cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  dep(end+1)      = cdep;
end
