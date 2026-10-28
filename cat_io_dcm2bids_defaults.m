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
  def.nproc                 = 0;     % parallel background MATLAB processes for the subject packages (0 - none, i.e.
                                     % this process; each process needs about 4 GB RAM, see cat_io_dcm2bids_parallel)

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
  def.opts.output           = 2;     % 0-tables only, 1-only json, 2-json+nii of scans with (semi)conform protocol, 3-all
  %def.opts.studies          = '';    % study filter >> file
  def.opts.sitelength       = 3;     % letters of the city in the site abbreviation without site dictionary (city, 
                                     % manufacturer, and field strength, e.g. HeiSI30T, see siteRegistry in 
                                     % cat_io_dcm2bids_protocols)
  def.opts.sitevendor       = 2;     % letters of the manufacturer in the site abbreviation (0 - none, 1 - e.g. S, 
                                     % 2 - e.g. SI for Siemens, GE, PH for Philips, see vendorCode in cat_io_dcm2bids_protocols)
  def.opts.subIDform        = 2;     % (1) use only PatientID, i.e. sub-PID (links the scanners of a subject in 
                                     %     the BIDS output and for the T1w of other sessions, only for unique IDs)
                                     % (2) add the site, i.e. sub-SITE-PID
                                     % (3) add the StudyName, i.e. sub-STUDY-PID
                                     % (4) sub-SITE-STUDY-PID, (5) sub-STUDY-SITE-PID
  def.opts.anonymize        = 1;     % anonymization (BIDS directory BIDS-<level> and DICOM export): 0 - no, 1 - pseudonymized 
                                     % (site and PatientID, scan date, defaced anatomical scans), 2 - anonymized (site and
                                     % number, time since the first session, shifted dates, defaced scans, see 
                                     % cat_io_dcm2bids_anonymize)
  def.opts.sesunit          = 'w';   % unit of the session labels and dates of the anonymization level 2, i.e. the time 
                                     % since the first session of the subject (rounded): 'd' - days, 'w' - weeks, 
                                     % 'm' - months, 'y' - years, e.g. ses-w0000 and ses-w0052 one year later

  def.opts.agelimit         = [1 inf]; % age limit of the sessions in years ([min max], see ageLimit in cat_io_dcm2bids):
                                     % sessions outside the limit or without age (e.g. age 0, the undefined age of
                                     % many DICOMs) are not processed and are listed in the table excluded_sessions of
                                     % catDCM2BIDSdb-tables, [-inf inf] keeps all sessions (also without age)

  def.opts.dicomexport      = 0;    % anonymized DICOM export of the imported scans (cleaned headers, defaced images,
                                     % birth date by opts.anonymize, see cat_io_dcm2bids_export): 0 - no, 1 - scans with
                                     % conform or semiconform protocol (all scans without protocols), 2 - all scans
  def.opts.exportmode       = 1;     % BIDS and DICOM export: 1 - extend (the exports of all runs, i.e. the scans of this
                                     % run are added), 2 - only selected (the exports contain only the scans of this run,
                                     % i.e. the directories <BIDS>-selected and catDCM2BIDSdb-dicom-selected are written
                                     % again by each run, whereas the database is always extended)
  def.opts.hrrealign       = 0;     % highres-realignment for further use? difficult to keep this persistent
  def.opts.rerun            = 0;     % rerun - overwrite existing
  def.opts.tableformat      = 'csv'; % csv/tsv
  def.opts.scantypes        = {'anat/T1w','anat/T2w','dwi/dwi','func/bold'}; % typical scan types (datatype/
                                     % suffix) for the counts of the site table and the columns of the session
                                     % table (see cat_io_dcm2bids_tables and cat_io_dcm2bids_report)
  def.opts.checkgzipi       = 1;     % look for unzipped data and zip it
  def.opts.workdir          = 'auto';% work directory of the unzipped images of the processing (see workSetup in 
                                     % cat_io_dcm2bids_helper): 'auto' - the temporary directory of the system if
                                     % its disk is at least 2 times faster than the database disk (e.g. internal SSD
                                     % vs. external hard disk), 'native' - in the database, or a directory
  def.opts.workcache        = 8;     % maximum size of the unzipped images of the work directory per process in GB
  def.opts.workroot         = '';    % (internal: work directory of this run, see workSetup)
  def.opts.runqc            = 1;     % run QC (all data)
  def.opts.preprocessing    = 1;     % basic preprocessing of all scans: 0 - none, 1 - processing (segmentation 
                                     % of anat, DTI of dwi, rs-fMRI measures of func, see cat_io_dcm2bids_pp_anat*), 
                                     % 2 - processing with derivatives (maps), whereas the maps of 1 stay in the 
                                     % database and are used for the renderings and report (derivatives index)
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
  def.opts.dtimix           = 1;     % mixed DTI of several DWI scans of a session (e.g. the standard and advanced
                                     % protocol of the DZPG, see processMix in cat_io_dcm2bids_pp_dwi)
  def.opts.dtimaps          = 1;     % DTI maps for opts.preprocessing = 2 (see cat_io_dcm2bids_pp_dwi): 1 - FA,
                                     % DEC (direction-encoded colors), and the normalized FA (wFA, e.g. for voxel-
                                     % based analyses), 2 - also MD, AD, and RD (their white matter values are
                                     % always estimated)
  def.opts.ignoreScouts     = 1;     % remove localizer/scouts ASAP (also derived images except the MP2RAGE uniform 
                                     % image and T1 map, and 2D images, see ignoredDerived in cat_io_dcm2bids_db)
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
  def.opts.orientation.OTM  = [10 60];   % orientation translation measure: RMSE of the origin offset in mm 
                                         % (i.e. an offset of about 17 and 100 mm along one axis; typical 
                                         % sessions of the test data have 10-27 mm)
  def.opts.orientation.ORM  = [ 5 30];   % orientation rotation measure: RMSE of the rotation in degree
                                         % (i.e. a rotation of about 9 and 50 degree around one axis; typical
                                         % sessions of the test data have 5-15 degree)
  def.opts.coverage.fwhm    = 4;         % smoothing of the coverage maps in mm
  % expected coverage of each datatype (otherwise func, see cat_io_dcm2bids_overlap): brain with CSF 
  % (csf 1) or GM and WM (0), without the outer erode mm, and the weight of the cerebellum and 
  % brainstem (of the expected map), i.e. anatomical scans should cover the whole brain, fMRI not 
  % the CSF and less the cerebellum, and DWI also not the outer 2 mm of the cortex (low signal), 
  % with the coverage measure BBM (missing coverage in percent of the expected brain, limited by 
  % the protocol mask <protocol>_msk.nii if available) for the marks 1 and 6, i.e. anatomical 
  % scans should miss nothing (a partially missing cerebellum of about 1.5% gives about 3), 
  % whereas fMRI and DWI are more tolerant (median of 0.1-1.6% and 1.2% in NKI-TRT, i.e. about 
  % 1.2-2.3 and 1.8), and protocols with a smaller coverage need a mask for good ratings
  def.opts.coverage.anat    = struct('csf',1,'erode',0,'cerebellum',1,  'BBM',[0 4]);
  def.opts.coverage.func    = struct('csf',0,'erode',0,'cerebellum',0.5,'BBM',[0 6]);
  def.opts.coverage.dwi     = struct('csf',0,'erode',2,'cerebellum',0.5,'BBM',[0 8]);
  def.opts.reversal.PRM     = [.80 .95]; % phase-encoding reversal measure: correlation of the first volumes of scans
                                         % with opposite phase-encoding direction (e.g. AP/PA field maps), where high
                                         % values indicate a missing reversal (see checkPEreversal in cat_io_dcm2bids_pp_anat)

  % default QC definitions of each BIDS datatype for scans without the QC
  % definition of a protocol (qc<protocol>.json), i.e. without protocol 
  % definitions or with a nonconform protocol: [best worst] measures for the
  % marks 1 and 6 (NaN to skip a measure, see qualityRating in
  % cat_io_dcm2bids_qc), broader than the DZPG definitions, as they cover
  % different protocols: 0.75 times the median of the measure (mark 1) and
  % 2.5 times the median or its 90th percentile (mark 6) of the scans of the
  % NKI-RS, NKI-TRT, and Buchert test data with the default definition
  % (2026-10, n = 470-1290 per datatype; the DZPG definitions use the same
  % factors of the DZPG pilot data without percentile, see DCM2BIDS/DZPG3T)
  def.opts.qcdefault.anat   = struct('BSM',nan,          'WSM',nan,           'ISR',[0.0667 0.222],'NSR',[0.00127 0.039], 'RES',[0.8 2.0]);
  def.opts.qcdefault.func   = struct('BSM',[0.065 0.237],  'WSM',nan,           'ISR',[0.0509 0.170],'NSR',[0.0246 0.0819],'RES',[2.0 4.0]);
  def.opts.qcdefault.dwi    = struct('BSM',[0.0622 0.213],'WSM',[0.0299 0.0995],'ISR',[0.0347 0.116],'NSR',[0.0503 0.168], 'RES',[1.5 3.0]);
  def.opts.qcdefault.fmap   = struct('BSM',[0.00 0.15],    'WSM',nan,           'ISR',[0.0688 0.229],'NSR',[0.0199 0.0663],'RES',[2.0 4.0]);
  def.opts.qcdefault.other  = struct('BSM',[0.10 0.30],'WSM',nan,       'ISR',[0.10 0.25],'NSR',[0.03 0.10], 'RES',[2.0 4.0]);

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
  def.opts.render.subfilter = {};      % only the subject and session reports (R3) of these subjects (BIDS labels, 
                                       % e.g. of a subject package, see cat_io_dcm2bids_parallel), with the slice
                                       % cache for R1/R2 (empty - all subjects and outputs)
  def.opts.render.dbdir     = '';      % database directory with the index of the derivatives that were not exported
                                       % (opts.preprocessing = 1, see registerDerivative in cat_io_dcm2bids_helper)
                                       % ('' - catDCM2BIDSdb of the output directory)

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
  def.opts.reportfig.agerange  = [];           % age range of the age pyramid and the brain markers in years (empty: 
                                              % the ages of the sample rounded to 10 years for a range of at least 20 
                                              % years, to 5 years for 10 to 20 years, and to 1 year below)
  def.opts.reportfig.agestep   = [];           % bin size of the age pyramid in years (empty: 1, 2, 5, or 10 years
                                              % for at most 20 bins)
  def.opts.reportfig.maxsites  = 20;           % maximum number of sites shown separately (site table columns, 
                                              % site labels, and site pages), more sites only together 

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
