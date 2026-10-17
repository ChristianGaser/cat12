function varargout = cat_io_dcm2bids_readme(action,varargin)
%cat_io_dcm2bids_readme. Readme files of the output of cat_io_dcm2bids.
%  The readme.txt files describe the output directory (<outdir>/<subdir>)
%  and its subdirectories with their content, personal data, sharing, and
%  whether they can be removed (i.e. are created again by the batch):
%    readme.txt                 overview of all directories (main)
%    catDCM2BIDSdb/             internal database (DO NOT CHANGE)
%    catDCM2BIDSdb-tables/      internal tables
%    catDCM2BIDSdb-dicom/       anonymized DICOM export
%    <BIDSdir><level>/          BIDS datasets of the anonymization level,
%                               i.e. BIDS (no anonymization),
%                               BIDS-pseudonymized, or BIDS-anonymized
%    <BIDSdir><level>-private/  private participant tables
%    <BIDSdir><level>-report/   report tables, figures, and renderings
%  The files are written again in each run (for the existing directories),
%  i.e. own notes should not be added to them. The disclaimer of the
%  anonymization (best effort, the user is responsible) is given in the
%  readme files of the exports and as warning of each run.
%
%  varargout = cat_io_dcm2bids_readme(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_readme('writeReadmes',Poutdir,job)
%    cat_io_dcm2bids_readme('writeReadme',type,Pdir,job,level)
%    txt = cat_io_dcm2bids_readme('disclaimer',level)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_export, cat_io_dcm2bids_anonymize.

  switch action
    case {'writeReadmes','writeReadme','disclaimer'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_readme:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function writeReadmes(Poutdir,job)
%writeReadmes. Readme files of the output directory Poutdir and of its
%  existing subdirectories (see the help of this file).
  D = dirs(Poutdir,job);
  writeReadme('main',Poutdir,job,job.opts.anonymize);
  for di = 1:size(D,1)
    if ~isempty(D{di,2}) && exist(fullfile(Poutdir,D{di,1}),'dir')
      writeReadme(D{di,2},fullfile(Poutdir,D{di,1}),job,D{di,3});
    end
  end
  % readme of older versions in the BIDS datasets (see the readme of the BIDS directory)
  Pold = {};
  if exist(fullfile(Poutdir,job.BIDSdir),'dir'), Pold = cat_vol_findfiles(fullfile(Poutdir,job.BIDSdir),'readme.txt',struct('depth',1)); end
  for k = 1:numel(Pold)
    if startsWith(cat_io_dcm2bids_helper('readText',Pold{k}),'catDCM2BIDS BIDS export sub-directory'), delete(Pold{k}); end
  end
end
% =========================================================================
function txt = disclaimer(level)
%disclaimer. Disclaimer of the anonymization level (text with line breaks).
  switch level
    case 0
      txt = {
        '!!! catDCM2BIDS without anonymization (level 0): the output contains personal data (IDs, dates, times, !!!'
        '!!! and images with faces), i.e. it is not for sharing.                                                !!!'};
    otherwise
      name = {'pseudonymization','anonymization'};
      txt = {
        sprintf('!!! catDCM2BIDS %s (level %d) is a best effort of the software without any warranty. !!!',name{level},level)
        'The software removes personal data of the DICOM headers and the faces of the images (see the help of the'
        'anonymization), but it cannot guarantee that a person cannot be identified, e.g. by defaced (or also'
        'skull-stripped) images, rare phenotypes, small sites, or other information. The user is responsible for'
        'checking the data (e.g. the phenotypic data) and for the data protection rules of the study before sharing.'};
  end
  txt = strjoin(txt,newline);
end
% =========================================================================
function writeReadme(type,Pdir,job,level)
%writeReadme. Readme file of the directory Pdir of the type main, db,
%  tables, dicom, bids, private, or report of the anonymization level (the
%  job gives the names of the BIDS directories and the table format).
  if ~isfield(job,'BIDSdir'), job.BIDSdir = 'BIDS'; end
  if ~exist('level','var') || isempty(level), level = job.opts.anonymize; end
  fmt = job.opts.tableformat; b = baseName(job.BIDSdir);
  switch type
    case 'main'
      txt = mainText(Pdir,job);
    case 'db'
      txt = {
        '!!! DO NOT CHANGE THIS DIRECTORY !!!'
        ''
        'Internal database of catDCM2BIDS (CAT, cat_io_dcm2bids) with the imported scans, their processing, and the'
        'phenotypic data of the sessions. All other directories of the output are created from this database and'
        sprintf('can be removed (%s*, catDCM2BIDSdb-tables, catDCM2BIDSdb-dicom), i.e. they are created again by the',b)
        'batch (e.g. with this database as input), whereas changes of the database (e.g. moved, renamed, or deleted'
        'files) can corrupt the later runs. The database refers to the original DICOM directories (for the image'
        'conversion and the DICOM export): keep the original data!'
        ''
        'CONTAINS PERSONAL DATA (names, IDs, and birth dates of the patients, original images with faces, and the'
        'key of the anonymization): NOT FOR SHARING!'
        ''
        'Layout: sub-<DeviceSerialNumber>-<PatientID>/ses-<StudyID>-<ScanDate>/snr-<SeriesNumber>-<ProtocolName>/'
        '  <folder>=<protocol>=<time>=<series>.json   DICOM header (sidecar of dcm2niix)'
        '  <folder>=<protocol>=<time>=<series>.nii.gz image (and .bval/.bvec) if required (not for JSON-only outputs)'
        '  catDCM2BIDSimport_*.json                   import information (original DICOM directory, dcm2niix options)'
        '  <protocol>-<site>.json                     protocol of the scan (template of protocol definitions)'
        '  anon_*.nii.gz                              defaced image (anonymization)'
        '  catDCM2BIDSqc_*.json/.mat, *_denoisedavg.nii.gz'
        '                                             quality control (denoised average of 4D data)'
        '  catDCM2BIDSsegus_*.json, c0/m/wc0/wc1/mwc1/l0/wl0/y_anon_*.nii.gz, *_seg8.mat'
        '                                             segmentation of the anatomical scan of each session'
        '                                             (T1w, otherwise UNIT1, T2w, or FLAIR; volumes, regions)'
        '  catDCM2BIDSdwi_*, catDCM2BIDSrsfmri_*, r*.nii.gz, rp_*.txt'
        '                                             DTI and rs-fMRI measures and maps, realignment'
        '  ses-*/catDCM2BIDSmaff8.mat                 affine registration of the session to MNI space'
        '  ses-*/catDCM2BIDSphenotypic.json           phenotypic data of the session (user tables)'
        '  catDCM2BIDSsites.json                      site abbreviations of the scanners (BIDS labels, never changed)'
        '  catDCM2BIDSanonymization.json              numbers and date offsets of the anonymization level 2, i.e. the'
        '                                             key of the anonymized data (never changed)'
        '  catDCM2BIDStmp_*                           temporary directories (removed automatically)'
        ''
        'Directories of this database (the database itself, sub-*, ses-*, or snr-*) can be selected as input to'
        'export the imported data again (e.g. with other protocols or options).'
        };
    case 'tables'
      txt = {
        'Internal tables of the database catDCM2BIDSdb (see cat_io_dcm2bids_tables) with an overview of the scans'
        'that were processed in this and previous runs.'
        ''
        'CONTAINS PERSONAL DATA (PatientIDs, dates, names in the phenotypic reports): NOT FOR SHARING!'
        ''
        sprintf('  scans.%-20s one row per scan (site, session, BIDS type, protocol test, quality rating)',fmt)
        '  scans.mat                  cumulative list of the scans of all runs (basis of the tables)'
        sprintf('  sites.%-20s one row per scanner (sites.json: template of the site dictionary)',fmt)
        sprintf('  protocols.%-16s one row per protocol variant (protocols/<datatype>/<variant>.json)',fmt)
        sprintf('  protocols_sites.%-10s one row per site and protocol variant',fmt)
        '  protocols-details/         deviations of the scans of each site from the closest protocol'
        sprintf('  subjects.%-17s one row per patient',fmt)
        '  phenotypic_<table>.csv     match of the rows of the phenotypic data tables to the sessions'
        ''
        'Can be removed: the tables are created again by the next run for its scans (i.e. run the batch with the'
        'database as input to get the tables of all scans).'
        };
    case 'dicom'
      txt = [{
        'DICOM export of catDCM2BIDS (see cat_io_dcm2bids_export), e.g. to transfer the imported data to a central'
        'site that runs the batch again.'
        ''
        'THIS IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED. KEEP THE ORIGINAL DATA!'
        ''
        '  <yyyyMMdd-HHmmss>/      one directory per run with the new scans'
        '    sub-*_ses-*.zip       DICOM files of a session in the layout of the database (sub-*/ses-*/snr-*/<series>_<instance>.dcm)'
        sprintf('    %-22scontent of the zip file: the images of the series (database path, BIDS file for',['sub-*_ses-*.' fmt])
        '                          level 2), the series directory, number of DICOM files, MB, and acquisition duration'
        sprintf('  %-24sexported scans of all runs',['scans.' fmt])
        sprintf('  %-24sexported zip files of all runs',['sessions.' fmt])
        '  exports.mat             lists of the tables: each scan is exported only once (for an anonymization level),'
        '                          i.e. the run directories can be removed after the upload'
        ''
        'Pseudonymized (level 1): names of the patient (replaced by the PatientID), physicians, and operators,'
        'addresses, phone numbers, other IDs, and comments are cleaned and the birth date is set to the 15th of the'
        'month, whereas the fields used by catDCM2BIDS are kept (e.g. PatientID, sex, age, dates, UIDs, scanner,'
        'protocol, private fields). The PatientID has to be a pseudonym!'
        'Anonymized (level 2): furthermore, the PatientID and PatientName are replaced by the BIDS subject label'
        '(site and number, e.g. BerSI30T00001) and the StudyID by the session label (time since the first session,'
        'e.g. w0052), the dates of a session are shifted to the date of its label (first session of a subject at'
        '1900-01-01, e.g. 1901-01-01 for w0052) and its times to the time since its start (00:00:00, i.e. the timing'
        'of the scans and slices is kept), the UIDs are replaced by UIDs of the same length (also in private'
        'fields), the scanner and institution are replaced by the site abbreviation and city, the station,'
        'department, and the study and procedure descriptions are removed, parts of the patient name, the'
        'PatientID, and the serial number in other fields are replaced by X, and the export uses the BIDS labels'
        'as names. Its import gives the labels of the anonymized BIDS export (also with the anonymization level 1).'
        'Images: defaced (anatomical scans for level 1, all scans for level 2).'
        ''
        'Can be removed: all scans are exported again by the next run.'
        ''}; {disclaimer(max(1,level))}];
    case 'bids'
      txt = [{
        sprintf('BIDS export of catDCM2BIDS (https://bids.neuroimaging.io) with one BIDS dataset per subdirectory (e.g. %s/BIDS,',b)
        'or one dataset per set of protocols with protocol-specific subdirectories and "mismatch" for the scans of other'
        'protocols):'
        '  sub-<label>/ses-<label>/<datatype>/   images (NIfTI) and cleaned sidecars (JSON)'
        '  sub-<label>/sub-<label>_sessions.tsv  sessions (site, age, weight, orientation rating, phenotypic data)'
        '  participants.tsv/.json                participants (age, sex)'
        '  derivatives/                          quality control, segmentation, DTI, and rs-fMRI measures and maps'
        '  dataset_description.json, README.md   description of the dataset (dataset descriptor of the batch)'
        '  catDCM2BIDSdb.mat                     settings of the batch that define the BIDS structure'
        ''}; levelText(level); {
        ''
        'Can be removed: created again by the batch (with the database catDCM2BIDSdb as input).'
        ''}; {disclaimer(level)}];
    case 'private'
      txt = {
        '!!! CONTAINS PERSONAL DATA - NOT FOR SHARING !!!'
        ''
        'Private tables of the BIDS export: <BIDSsubdir>/private_participants.tsv with the BIDS subject label and'
        'the PatientID, PatientName, sex, age, weight, birth date, and acquisition date of each subject, i.e. the key'
        'to identify the subjects of the BIDS export (e.g. in case of findings).'
        ''
        'Can be removed: created again by the batch.'
        };
    case 'report'
      txt = [{
        'Reports of the BIDS export (see cat_io_dcm2bids_report, cat_io_dcm2bids_reportfig, and cat_io_dcm2bids_render)'
        'with the BIDS labels of the subjects:'
        '  tables/<BIDSsubdir>/      report tables (scans, scans_<type>, sessions, project) with JSON column descriptions'
        '  regions/<BIDSsubdir>/     regional results of each session (anat, dwi, func)'
        '  render/<BIDSsubdir>/      rendered slices of the BIDS images (PDF of each subject)'
        '  report_<BIDSsubdir>_<design>.pdf   report figure (one PDF per design, e.g. CAT_bright)'
        ''}; levelText(level); {
        ''
        'Can be removed: created again by the batch.'
        ''}; {disclaimer(level)}];
    otherwise
      error('cat_io_dcm2bids_readme:unknownType','Unknown readme type "%s".',type);
  end
  if ischar(txt), txt = {txt}; end
  cat_io_dcm2bids_helper('writeText',fullfile(Pdir,'readme.txt'),[strjoin(txt,newline) newline]);
end
% =========================================================================
function txt = levelText(level)
%levelText. Description of the BIDS labels and data of the anonymization
%  level (cell of lines).
  switch level
    case 0
      txt = {
        'NOT ANONYMIZED: subject labels with the PatientID, session labels with the scan date and time, images with'
        'faces (besides the defacing settings), i.e. personal data: NOT FOR SHARING!'};
    case 1
      txt = {
        'PSEUDONYMIZED: subject labels with the site and the PatientID (has to be a pseudonym), session labels with the'
        'scan date, age with one decimal, defaced anatomical images, no acquisition times in the NIfTI headers, i.e.'
        'intended for sharing according to the data protection rules of the study.'};
    otherwise
      txt = {
        'ANONYMIZED (for recipients without the key): subject labels with the site and a number (e.g.'
        'sub-BerSI30T00001), session labels with the time since the first session of the subject (e.g. ses-w0000 and'
        'ses-w0052 one year later, unit d-days, w-weeks, m-months, or y-years), age in years, defaced images of all'
        'scans, no acquisition times in the NIfTI headers. The key (numbers and date offsets) stays in the database'
        '(catDCM2BIDSanonymization.json) and the private tables, i.e. the creator can identify the subjects (e.g. in'
        'case of findings).'};
  end
end
% =========================================================================
function b = baseName(B)
%baseName. BIDS directory without the name of the anonymization level.
  b = regexprep(B,'-(pseudonymized|anonymized)$','');
end
% =========================================================================
function D = dirs(Poutdir,job)
%dirs. Subdirectories of the output directory with their readme type,
%  anonymization level, content, personal data, sharing, and whether they
%  are created again (rows without type continue the content of the
%  previous row): the database, its tables and DICOM export, and the BIDS
%  directories (BIDS, private, report) of the current level and of the
%  existing directories of other levels.
  D = {
    'catDCM2BIDSdb',        'db',     nan, 'internal database (DICOM headers, original images, processing,', 'YES',           'NO',       'NO'
    '',                     '',       nan, 'phenotypic data, key of the anonymization)',                     '',              '',         ''
    'catDCM2BIDSdb-tables', 'tables', nan, 'internal tables (scans, sites, protocols, subjects,',            'YES',           'NO',       'yes'
    '',                     '',       nan, 'phenotypic match reports)',                                      '',              '',         ''
    'catDCM2BIDSdb-dicom',  'dicom',  job.opts.anonymize, 'DICOM export (zip files of the sessions)',        'pseudo/anonym.', 'transfer', 'yes'
    };
  b = baseName(job.BIDSdir); lev = {'','-pseudonymized','-anonymized'};
  pd  = {'YES','pseudonymized','anonymized'}; sh = {'NO','yes','yes'};
  for L = 0:2
    B = [b lev{L+1}];
    if L ~= job.opts.anonymize && ~exist(fullfile(Poutdir,B),'dir'), continue; end
    D = [D; {
      B,              'bids',    L, 'BIDS datasets (images, cleaned sidecars, sessions and',       pd{L+1}, sh{L+1}, 'yes'
      '',             '',        L, 'participants files, derivatives)',                            '',      '',      ''
      [B '-private'], 'private', L, 'private participant tables (names, IDs, birth dates)',        'YES',   'NO',    'yes'
      [B '-report'],  'report',  L, 'report tables, figures, and renderings',                      pd{L+1}, sh{L+1}, 'yes'
      }]; %#ok<AGROW>
  end
end
% =========================================================================
function txt = mainText(Poutdir,job)
%mainText. Readme of the output directory: the batch, the settings of the
%  last run, the table of the subdirectories, and the disclaimer.
  try
    [catname,catrev] = cat_version; catver = sprintf('%s (r%s)',catname,catrev);
  catch
    catver = 'CAT';
  end
  o = job.opts;
  set = sprintf('output level %d, anonymization %d, preprocessing %d (level %g), DICOM export %d', ...
    o.output,o.anonymize,o.preprocessing,o.pplevel,o.dicomexport);
  if isfield(job,'dicts') && isfield(job.dicts,'Pprotocoldirs') && ~isempty(char(job.dicts.Pprotocoldirs))
    set = sprintf('%s, protocols %s',set,strjoin(spm_file(cellstr(job.dicts.Pprotocoldirs),'basename'),', '));
  end
  txt = {
    'catDCM2BIDS output directory'
    '============================'
    ''
    'Output of the catDCM2BIDS batch of CAT (Computational Anatomy Toolbox, cat_io_dcm2bids) that imports DICOM (or'
    'BIDS) data into an internal database, tests the MRI protocols, anonymizes and preprocesses the scans, rates their'
    'quality, and exports them to BIDS. The readme files of the output are written again by each run (see'
    'cat_io_dcm2bids_readme), i.e. please do not add own notes to them.'
    ''
    sprintf('Last run: %s, %s',char(datetime('now','Format','yyyy-MM-dd HH:mm')),catver)
    sprintf('Settings: %s.',set)
    ''
    sprintf('%-28s %-62s %-14s %-9s %s','Directory','Content','Personal data','Sharing','Recreated')
    repmat('-',1,126)
    };
  D = dirs(Poutdir,job);
  for di = 1:size(D,1)
    txt{end+1} = deblank(sprintf('%-28s %-62s %-14s %-9s %s',D{di,1},D{di,4},D{di,5},D{di,6},D{di,7})); %#ok<AGROW>
  end
  txt = [txt; {
    ''
    'DO NOT CHANGE THE DATABASE (catDCM2BIDSdb): it is the basis of all other directories, which can be removed and'
    'are created again by the batch (with the database as input). The database refers to the original DICOM'
    'directories (for the image conversion and the DICOM export), i.e. keep the original data.'
    ''
    'Anonymization levels: 0 - none (BIDS, personal data), 1 - pseudonymized (BIDS-pseudonymized, labels with the'
    'site, PatientID, and scan date), 2 - anonymized (BIDS-anonymized, labels with the site, a number, and the'
    'time since the first session, the key stays in the database). Please check the data protection rules of your study'
    'before sharing.'
    ''}; {disclaimer(max(1,job.opts.anonymize))}];
end
