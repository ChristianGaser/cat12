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
  def.dataset.phenotypic         = {};         % tables (CSV/TSV) with phenotypic data of the sessions (see 
                                               % cat_io_dcm2bids_phenotypic), exported to the sessions files
  def.dataset.phenotypickeys     = [];         % key columns of the tables (numbers, empty: detected by their names)

  def.dicts.Pprotocoldirs   = {};    % input DCM protocol directories
  def.dicts.Pcenterdict     = {};    % dictionary for center names (otherwise scanner ID)
  def.dicts.Pstudydict      = {};    % not implemented yet
  def.dicts.Psubjdict       = {};    % not implemented yet

  def.opts.ProtocolFileName = 1;     % replace protocol name by the filenames
                                     % of the evaluation protocols under protocol directories
  def.opts.gzipi            = 1;     % internal use of nii.gz (save disk space but maybe a bit slower)
  def.opts.gzipe            = 1;     % external use of nii.gz (save disk space but nonoptimal for SPM processing)
  def.opts.tolerance        = 0.5;   % tolerance in percent for numeric MR parameters (only for rounding, as the same
                                     % protocol results in the same values, see testProtocols in cat_io_dcm2bids_protocols)
  def.opts.Pdcm2nii         = '';    % dcm2niix executable ('' - detected if required, see cat_io_dcm2bids_importDCM)
  def.opts.output           = 2;     % 0-no, 1-only json, 2-json+nii (input protocols), 3-full
  %def.opts.studies          = '';    % study filter >> file
  def.opts.sitelength       = 3;     % letters of the city in the site abbreviation without site dictionary (city, 
                                     % manufacturer, and field strength, e.g. HeiSI30T, see siteRegistry in 
                                     % cat_io_dcm2bids_protocols)
  def.opts.sitevendor       = 2;     % letters of the manufacturer in the site abbreviation (0 - none, 1 - e.g. S, 
                                     % 2 - e.g. SI for Siemens, GE, PH for Philips, see vendorCode in cat_io_dcm2bids_protocols)
  def.opts.subIDform        = 2;     % (1) use only PatientID, i.e. sub-PID
                                     % (2) add the ScannerID, i.e. sub-SITE-PID
                                     % (3) add also StudyName, i.e. sub-STUDY-SITE-PID
                                     % (4) i.e. sub-SITE-STUDY-PID
  def.opts.anonymize        = 1;     % anonymization (BIDS directory BIDS-<level> and DICOM export): 0 - no, 1 - pseudonymized 
                                     % (site and PatientID, scan date, defaced anatomical scans), 2 - anonymized (site and
                                     % number, time since the first session, shifted dates, defaced scans, see 
                                     % cat_io_dcm2bids_anonymize)
  def.opts.sesunit          = 'w';   % unit of the session labels and dates of the anonymization level 2, i.e. the time 
                                     % since the first session of the subject (rounded): 'd' - days, 'w' - weeks, 
                                     % 'm' - months, 'y' - years, e.g. ses-w0000 and ses-w0052 one year later

  def.opts.dicomexport      = 0;     % anonymized DICOM export of the imported scans (cleaned headers, defaced images,
                                     % birth date by opts.anonymize, see cat_io_dcm2bids_export): 0 - no, 1 - scans with
                                     % conform or semiconform protocol (all scans without protocols), 2 - all scans
  def.opts.hrrealign        = 0;     % highres-realignment for further use? difficult to keep this persistent
  def.opts.rerun            = 0;     % rerun - overwrite existing
  def.opts.tableformat      = 'csv'; % csv/tsv
  def.opts.scantypes        = {'anat/T1w','anat/T2w','dwi/dwi','func/bold'}; % typical scan types (datatype/
                                     % suffix) for the counts of the site table and the columns of the session
                                     % table (see cat_io_dcm2bids_tables and cat_io_dcm2bids_report)
  def.opts.checkgzipi       = 1;     % look for unzipped data and zip it
  def.opts.runqc            = 1;     % run QC (all data)
  def.opts.preprocessing    = 1;     % basic preprocessing of all scans: 0 - none, 1 - processing (segmentation 
                                     % of anat, DTI of dwi, rs-fMRI measures of func, see cat_io_dcm2bids_pp_anat*), 
                                     % 2 - processing with derivatives (maps)
  def.opts.denoise          = 0;     % denoising (sanlm) of anat before the segmentation and of the DWI volumes 
                                     % before the tensor fit (see cat_io_dcm2bids_pp_dwi)
  def.opts.pplevel          = 4;     % preprocessing level, i.e. accuracy and time of the DTI and rs-fMRI measures 
                                     % (see levelParameters in cat_io_dcm2bids_pp_anat): 1 - full quality (SPM default 
                                     % parameters, reference), 2 - optimized (fast parameters, slight resolution 
                                     % limits), 3 - quick review (lower resolution, only the first DWI and BOLD
                                     % scan), 4 - preview; continuous levels with about 26, 9, 4, and 2.5 min for a
                                     % session with T1w, 2 DWI, and rs-fMRI scans (4 for the development, later 2);
                                     % one segmentation per session for all levels (T1w, otherwise UNIT1, T2w, or
                                     % FLAIR, see anatScan in cat_io_dcm2bids_pp_anat)
  def.opts.dtimaps          = 1;     % DTI maps for opts.preprocessing = 2 (see cat_io_dcm2bids_pp_dwi): 1 - FA,
                                     % DEC (direction-encoded colors), and the normalized FA (wFA, e.g. for voxel-
                                     % based analyses), 2 - also MD, AD, and RD (their white matter values are
                                     % always estimated)
  def.opts.ignoreScouts     = 1;     % remove localizer/scouts ASAP (also derived and 2D images by dcm2niix -i y)
  def.opts.protocolsubdirs  = 0;     % use additional sub-directories to separate between protocols
  def.opts.BIDSsep          = '';    % Use extra BIDS-incompatible separator such as - to separate the site
                                     % sub-SITE-SUBJECT

  % orientation rating of each session (RMSE of the offset of the image 
  % origin and of the rotation to MNI space, see writeSessionAffines in 
  % cat_io_dcm2bids_pp_anat), coverage rating of each scan (missing brain
  % coverage, see cat_io_dcm2bids_overlap), and phase-encoding reversal
  % rating of scans with opposite phase-encoding direction with the
  % [best worst] measure for the marks 1 and 6 (similar to the protocol
  % qc*.json files)
  def.opts.orientation.OTM  = [ 6 30];   % orientation translation measure: RMSE of the origin offset in mm 
                                         % (i.e. an offset of 10 and 50 mm along one axis)
  def.opts.orientation.ORM  = [ 3 12];   % orientation rotation measure: RMSE of the rotation in degree
                                         % (i.e. a rotation of 5 and 20 degree around one axis)
  def.opts.coverage.BBM     = [ 0 20];   % coverage measure: missing coverage in percent of the expected 
                                         % brain (whole brain or protocol mask <protocol>_msk.nii)
  def.opts.coverage.fwhm    = 4;         % smoothing of the coverage maps in mm
  def.opts.reversal.PRM     = [.80 .95]; % phase-encoding reversal measure: correlation of the first volumes of scans
                                         % with opposite phase-encoding direction (e.g. AP/PA field maps), where high
                                         % values indicate a missing reversal (see checkPEreversal in cat_io_dcm2bids_pp_anat)

  % protocol fit rating of each scan (PFR) of the protocol fit measure PFM,
  % i.e., the weighted number of deviations from the closest protocol (0.5
  % for each numeric field outside the tolerance, 1 for each text field or
  % list entry), and the protocol status by the PFR (see testProtocols in
  % cat_io_dcm2bids_protocols), where the PFM scale is defined by the marks
  % 0.5 and 4.5 (the best mark and the limit of semiconform scans)
  def.opts.protocol.PFM     = [ 0 3.5];  % PFM for the marks 0.5 and 4.5 (i.e. about 1.1 marks per deviation point)
  def.opts.protocol.status  = [1.5 4.5]; % PFR limits of the protocol status: conform (PFR < 1.5, scans are named by
                                         % the protocol), semiconform, and nonconform (PFR > 4.5, unknown protocol)
  % critical parameters, which are never the same if they differ (10 points,
  % i.e., nonconform), as well as the image component of the ImageType
  % (magnitude/phase) and the direction in the series name (e.g. PA of
  % "SpinEchoFieldMap_PA") that has to fit to the direction in the protocol
  % name, and key parameters of the contrast and resolution (at least 1
  % point, whereas other numeric deviations count 0.5)
  def.opts.protocol.critical = {'PhaseEncodingDirection','InPlanePhaseEncodingDirectionDICOM', ...
    'PulseSequenceDetails','ScanningSequence','MRAcquisitionType','MagneticFieldStrength', ...
    'EchoNumber','BValueShells'};
  def.opts.protocol.key      = {'EchoTime','EchoTime1','EchoTime2','RepetitionTime','InversionTime', ...
    'FlipAngle','SliceThickness','VoxelSize','MultibandAccelerationFactor','SliceOrder'};

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

  % report figure of each BIDS directory (see cat_io_dcm2bids_reportfig): A4
  % overview page and one page per site (for several sites) with diagrams 
  % of the report tables in <BIDSdir>-report/report_<BIDSsubdir>_<design>.pdf
  def.opts.reportfig.run       = 1;            % create the report figure
  def.opts.reportfig.design    = {'CAT_bright','CAT_dark','DZPG_bright'}; % designs (one PDF each, 
                                              % DCM2BIDS/colormaps/<design>.json with the logos, rating icons, colors, 
                                              % and colormaps): 'CAT_bright', 'CAT_dark', and/or 'DZPG_bright'
  def.opts.reportfig.fontsize  = 8;            % font size of the text in pt
  def.opts.reportfig.fontscale = [0.96 0.8];   % relative font size of header-like elements (titles of the diagrams) 
                                              % and of axes, table entries, and legends
  def.opts.reportfig.agerange  = [0 100];      % age range of the age pyramid in years
  def.opts.reportfig.agestep   = 5;            % bin size of the age pyramid in years

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
