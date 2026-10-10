function def = cat_io_dcm2bids_defaults
%cat_io_dcm2bids_defaults. Default settings of cat_io_dcm2bids.
%  The defaults define the job structure of cat_io_dcm2bids, i.e., the
%  settings of the GUI (cat_conf_dcm2bids, that uses these values) and
%  further (expert/internal) settings that are not set by the GUI but can
%  be changed by the job in scripts (see cat_io_checkinopt). Functions
%  without access to the job use this function directly (e.g. the render
%  function or the output structure of the batch dependencies).
%  The quality ratings of the scans are defined by the protocols (qc*.json
%  files of the protocol directories).
%
%  def = cat_io_dcm2bids_defaults
%
%  See also cat_io_dcm2bids, cat_conf_dcm2bids.

  def.data                  = {};    % input DCM directories (add JSON/NII input later)
  def.outdir                = {pwd}; % main output directory
  def.subdir                = 'study'; % default
  def.BIDSdir               = 'BIDS';

  % dataset description (dataset_description.json and README.md)
  def.dataset.Name               = '';         % empty: study subdirectory
  def.dataset.Authors            = {};         % one author per cell, e.g. 'Jane Doe <jane.doe@uni.de>'
  def.dataset.HowToAcknowledge   = '';
  def.dataset.License            = 'CC0-1.0';  % SPDX identifier, '' for none
  def.dataset.README             = {''};       % text file with the head of the README.md
  def.dataset.Acknowledgements   = '';
  def.dataset.Funding            = {};
  def.dataset.EthicsApprovals    = {};
  def.dataset.ReferencesAndLinks = {};
  def.dataset.DatasetDOI         = '';

  def.dicts.Pprotocoldirs   = {};    % input DCM protocol directories
  def.dicts.Pcenterdict     = {};    % dictionary for center names (otherwise scanner ID)
  def.dicts.Pstudydict      = {};    % not implemented yet
  def.dicts.Psubjdict       = {};    % not implemented yet

  def.opts.ProtocolFileName = 1;     % replace protocol name by the filenames
                                     % of the evaluation protocols under protocol directories
  def.opts.gzipi            = 1;     % internal use of nii.gz (save disk space but maybe a bit slower)
  def.opts.gzipe            = 1;     % external use of nii.gz (save disk space but nonoptimal for SPM processing)
  def.opts.tolerance        = 5;     % tolerance in percent for MR parameters (does not help for ordinal variables)
  def.opts.Pdcm2nii         = '';    % dcm2niix executable ('' - detected if required, see cat_io_dcm2bids_importDCM)
  def.opts.output           = 2;     % 0-no, 1-only json, 2-json+nii (input protocols), 3-full
  %def.opts.studies          = '';    % study filter >> file
  def.opts.subIDform        = 2;     % (1) use only PatientID, i.e. sub-PID
                                     % (2) add the ScannerID, i.e. sub-SITE-PID
                                     % (3) add also StudyName, i.e. sub-STUDY-SITE-PID
                                     % (4) i.e. sub-SITE-STUDY-PID
  def.opts.anonymize        = 1;     % imaging data: 0-no, 1-light-anat-only,  2-strong-all, 3-skull-strip-light?
                                     % meta data:    0-no, 1-basic-reduced,    2-strong-only-coded

  def.opts.hrrealign        = 0;     % highres-realignment for further use? difficult to keep this persistent
  def.opts.rerun            = 0;     % rerun - overwrite existing
  def.opts.tableformat      = 'csv'; % csv/tsv
  def.opts.checkgzipi       = 1;     % look for unzipped data and zip it
  def.opts.runqc            = 1;     % run QC (all data)
  def.opts.preprocessing    = 1;     % run segmentation (anat)
  def.opts.denoise          = 0;     % do denoising
  def.opts.ignoreScouts     = 1;     % remove localizer/scouts ASAP (also derived and 2D images by dcm2niix -i y)
  def.opts.protocolsubdirs  = 0;     % use additional sub-directories to separate between protocols
  def.opts.BIDSsep          = '';    % Use extra BIDS-incompatible separator such as - to separate the site
                                     % sub-SITE-SUBJECT

  % orientation rating of each session (RMSE of the offset of the image 
  % origin and of the rotation to MNI space, see writeSessionAffines in 
  % cat_io_dcm2bids_pp) and coverage rating of each scan (missing brain 
  % coverage, see cat_io_dcm2bids_overlap) with the [best worst] measure 
  % for the marks 1 and 6 (similar to the protocol qc*.json files)
  def.opts.orientation.OTM  = [ 6 30];   % orientation translation measure: RMSE of the origin offset in mm 
                                         % (i.e. an offset of 10 and 50 mm along one axis)
  def.opts.orientation.ORM  = [ 3 12];   % orientation rotation measure: RMSE of the rotation in degree
                                         % (i.e. a rotation of 5 and 20 degree around one axis)
  def.opts.coverage.BBM     = [ 0 20];   % coverage measure: missing coverage in percent of the expected 
                                         % brain (whole brain or protocol mask <protocol>_msk.nii)
  def.opts.coverage.fwhm    = 4;         % smoothing of the coverage maps in mm

  % render slices of each scan into BIDS-report/render/... (see cat_io_dcm2bids_render)
  def.opts.render.run       = 1;       % render images
  def.opts.render.source    = 3;       % 1-raw, 2-derivatives, 3-both
  def.opts.render.sessions  = 0;       % 1-only first session per subject, 0-all sessions
  def.opts.render.slicemode = 'mni';   % 'mni' - affine, 'mnirigid' - rigid registration to MNI space
  def.opts.render.sort      = 'name';  % order of scans: 'name' or 'SQR' (subjects by worst rating)
  def.opts.render.slices    = [3 10; 2 0; 1 0]; % R1 overview slices [orientation MNI-position]
                                       % with 3-axial, 2-coronal, 1-sagittal
  def.opts.render.tiles     = [4 5];   % R1 tiles per page [x y]: [3 4], [4 5], [5 7], or [6 8]
  def.opts.render.rows      = 1;       % R2 scan-row pages per image type (one row of slices per scan)
  def.opts.render.subjects  = 1;       % R3 subject reports (PDF with one row of slices per scan)
  % internal render settings (not in the GUI)
  def.opts.render.infoside  = 'left';  % side of the information panel in R2/R3: 'left' or 'right'
  def.opts.render.rowimg    = 1;       % write the slice row of each scan as image (e.g. for HTML)
  def.opts.render.cache     = 1;       % use and write the slice cache
  def.opts.render.fov       = 200;     % in-plane field of view in mm
  def.opts.render.res       = [];      % in-plane resolution in mm ([] - 300/dpi, i.e. 1 mm for 300 dpi)
  def.opts.render.interp    = [];      % interpolation: 0-nearest, 1-trilinear, 2-7 B-spline degree
                                       % ([] - trilinear for dpi<=300, otherwise 4th degree B-spline)
  def.opts.render.dpi       = 300;     % print resolution

  % datatypes and (lower case) suffixes of the output structure that is
  % used by the batch dependencies (see BIDSoutputs in cat_io_dcm2bids_bids),
  % for now, we keep this a bit shorter
  def.BIDSoutputs = {
    'anat', {'t1w','t2w','flair','pdw'}; %,'t2star','t1map','t2map','t2starmap','pdt2','mt','mtr'};
    'func', {'bold'}; %,'sbref'};
    'dwi',  {'dwi'}; %,'sbref'};
    'fmap', {'epi'}; %,'fieldmap','phase'};
    };
end
