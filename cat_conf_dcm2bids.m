function dcm2bids = cat_conf_dcm2bids(expert)
%cat_conf_dcm2bids. Batch definition to convert DICOM to BIDS. 

  if ~exist('expert','var')
    expert = cat_get_defaults('extopts.expertgui'); 
  end
  def = cat_io_dcm2bids_defaults; % default values (see cat_io_dcm2bids_defaults)

  % define input
  datadir          = cfg_files;
  datadir.tag      = 'data';
  datadir.name     = 'Input Directories';
  datadir.filter   = 'dir';
  datadir.ufilter  = '.*';
  datadir.num      = [1 Inf];
  datadir.help     = {
   ['Select directories with DICOM data. The directories and all their sub-directories are imported ' ...
    'into the database of the output directory (outputdir/studyDBdir/catDCM2BIDSdb) and then exported ' ...
    'to BIDS. The import first reads only the DICOM headers (fast) to store the JSON sidecars of new ' ...
    'scans together with their DICOM directory, whereas the images are converted only if they are ' ...
    'required (not for the JSON-only output levels). Directories that were imported before are not ' ...
    'read again (only with the rerun option, opts.rerun, in scripts). Series that are split over several ' ...
    'directories are incomplete and not imported (error message). ']
    ''
   ['To export already imported data again (e.g. with other protocols or if the DICOM data is no ' ...
    'longer available), select directories of this database: ']
    '  catDCM2BIDSdb                  .. all scans of the database'
    '  catDCM2BIDSdb/sub-*[/ses-*[/snr-*]] .. all scans of a subject, session, or a single scan'
    ''
   ['Directories of the database of another output directory are not supported. ' ...
    'DICOM and database directories can be combined, where each scan is only processed once. ']
    };

  % output directory
  outdir            = cfg_files;
  outdir.tag        = 'outdir';
  outdir.name       = 'Output Directory';
  outdir.filter     = 'dir';
  outdir.ufilter    = '.*';
  outdir.num        = [1 1];
  % box-drawing characters for the directory tree (├─ entry, └─ last entry, │ continuation)
  b = [char(9500) char(9472) ' ']; e = [char(9492) char(9472) ' ']; v = [char(9474) '  ']; s = '   ';
  outdir.help       = {[ ...
    'Select a directory where files are written to. ' ...
    'The batch will create a subdirectory "catDCM2BIDS" with converted NIFTI/JSON data. ' ...
    'In case of protocolsubdirs, the batch will create a protocol-conform and a non-conform BIDS structure.' ...
    'It will also create a subdirectory with used MR protocols and final reports. ']
    ''
    'outputdir/studyDBdir/'
   ['  ' b 'BIDS']
   ['  ' v e 'BIDS[-protocoldirs]']
   ['  ' v s e '[BIDS-protocoldir]']
   ['  ' v s s b '[derivatives']
   ['  ' v s s v e 'catDCM2BIDS]']
   ['  ' v s s v s e 'sub-[siteID-]patientID']
   ['  ' v s s v s s e 'ses-[date|seriesID]']
   ['  ' v s s v s s s e '{anat|dwi|func|fmaps|...}']
   ['  ' v s s v s s s s e '*sub-*_ses-*_*.{nii|json} .. processed data']
   ['  ' v s s e 'sub-[siteID-]patientID']
   ['  ' v s s s e 'ses-[date|seriesID]']
   ['  ' v s s s s e '{anat|dwi|func|fmaps|...}']
   ['  ' v s s s s s e 'sub-*_ses-*_*.{nii|json} .. raw images']
   ['  ' b 'BIDS-private          .. critical patient information (name, day of birth)']
   ['  ' b 'BIDS-report           .. report figure (report_*.pdf), tables, and renderings']
   ['  ' b 'catDCM2BIDSdb         .. main "database"']
   ['  ' e 'catDCM2BIDSdb-tables  .. main database overview tables']
   ''
    };
  clear b e v s;
  
  subdir            = cfg_entry;
  subdir.tag        = 'subdir';
  subdir.name       = 'Study Directory';
  subdir.strtype    = 's';
  subdir.num        = [0 Inf];
  subdir.val        = {def.subdir};
  subdir.help       = {
   ['The directory is created within the chosen output directory. ' ...
    'To use it also in BIDS use only letters and digits. ' ...
    'If no name is given no subdirectory is created. ']
    ''
    };




  % Dictionaries 
  % =======================================================================
  % define directory with protocol filter
  protocoldir          = cfg_files;
  protocoldir.tag      = 'Pprotocoldirs';
  protocoldir.name     = 'Protocol Directories';
  protocoldir.filter   = 'dir';
  protocoldir.ufilter  = '.*';
  protocoldir.val      = {{''}};
  protocoldir.num      = [0 Inf];
  protocoldir.help     = {
   ['Select directories with JSON files with specified DICOM entries to filter for relevant protocols. ' ...
    'If no directory is specified then all data will be exported to a BIDS directory. ' ...
    'The file name of the JSON filter file together with the series number ' ...
    'will be used to specify the BIDS acquisition entity "acq-###FILENAME", where scans of one acquisition ' ...
    '(e.g. SBRef and bold, field-map magnitude and phase, AP and PA) get the series number of the first one ' ...
    'and the file name without the information of the suffix and further entities (e.g. AP for dir-AP).']
    ''
    'E.g., a "myt1w.json" with:'
    '  {'
    '    "SeriesDescription":                   "mprage_sag_0p8mm",'
  	'    "SliceThickness":                      0.8,'
  	'    "EchoTime":                            0.00222,'
  	'    "RepetitionTime":                      2.4,'
  	'    "InversionTime":                       1.03,'
  	'    "FlipAngle":                           8'
    '  }'
    ''
   ['It is possible to define a quality file with additional prefix "qc" that ' ...
    'defines the range of the following quality measures:']
   ... ['To convert the quality measures into a standardized rating, each protocol JSON file ' ...
   ... 'requires an addition file with prefix "qc" that includes the scaling range for each ' ...
   ... 'quality measure ranging from perfect to unacceptable quality: ']
    ''
    'E.g., a "qcmyt1w.json" with:'
    '  {'
    '    "BSM": NaN,'
  	'  	 "WSM": NaN,'
  	'  	 "ISR": [0.10, 0.30],'
  	'    "NSR": [0.03, 0.09],'
  	'  	 "RES": [0.80, 1.00]'
    '  }'
    ''
   ['with Between Scan Movement (BSM), Within Scan Movement (WSM) for high-dimensional data ' ...
    '(e.g., functional and diffusion data but also anatomical rescans); ' ...
    'Inhomogeneity Signal Ratio (ISR), Noise Signal Ratio (NSR), and ' ...
    'RMS resolution value RES of the voxel dimension. ']
    }; 
  

  % define directory with protocol filter
  %%%%%% NOT FULLY IMPLEMENTED YET  
  studydict          = cfg_files;
  studydict.tag      = 'Pstudydict';
  studydict.name     = 'Study Dictionary File (expert)';
  studydict.filter   = 'any';
  studydict.ufilter  = '.*\.json$';
  studydict.val      = {{''}};
  studydict.num      = [0 Inf];
  studydict.hidden   = expert<1;
  studydict.help     = {
    ['In case of multiple centers it is possible to add this to the BIDS subject code to improve readability. ' ...
    'Define and link a json file that specifies your "StudySerialNumber" ' ...
    'and defines the desired "StudyAbbreviation". '] 
    ''
    'E.g., a "mystudies.json" with:'
    '  ['
    '    {'
    '      "StudySerialNumber":      "000815",'
    '      "StudyAbbreviation":      "Motion"'
    '    }'
    '    {'
    '      "StudyNumber":           "000007",'
    '      "StudyAbbreviation":      "ADNI"'
    '    }'
    '  ]'
    ''
    };   

  % Subject table
  % csv file to redefine IDs and add further data 
  %
  subjectdict          = cfg_files;
  subjectdict.tag      = 'Psubjectdict';
  subjectdict.name     = 'Subject Dictionary Table (expert)';
  subjectdict.filter   = 'any';
  subjectdict.ufilter  = '.*\.csv$';
  subjectdict.val      = {{''}};
  subjectdict.num      = [0 Inf];
  subjectdict.hidden   = expert<1;
  subjectdict.help     = {
    'Integration of phenotypical data by CSV tables with subject-specific "PatientID" or session-specific "StudyNumber".'
    ''
    'E.g., a "subjects.csv" with:'
    'PatientID, MMSE, GROUP'
    '000000001, 30,   0'
    '000000002, 12,   1'
    ''
    }; 

  % select files with center information 
  % automatic naming based on the site name and scanner field strength.
  centerdict          = cfg_files;
  centerdict.tag      = 'Pcenterdict';
  centerdict.name     = 'Center Dictionary File (expert)';
  centerdict.filter   = 'any';
  centerdict.ufilter  = 'site.*\.json$';
  centerdict.val      = {{''}};
  centerdict.num      = [0 Inf];
  centerdict.hidden   = true;
  centerdict.help     = {
    [ ...
    'In case of multiple centers it is possible to include a centerID into ' ...
    'the subject code to avoid overlap of center specific subject IDs. ' ...
    'Define and link a json file that specifies your "DeviceSerialNumber" ' ...
    'and defines the desired "InstitutionAbbreviation". ' ...
    ] 
    ''
    'E.g., a "mysites.json" with:'
    '  ['
    '    {'
    '      "DeviceSerialNumber":           "000815",'
    '      "InstitutionAbbreviation":      "JE",'
    '      "Color":                        "#A82117"'
    '    },'
    '    {'
    '      "DeviceSerialNumber":           "000007",'
    '      "InstitutionAbbreviation":      "NA"'
    '    },'
    '    {'
    '      "Name":                         "partner site Halle-Jena-Magdeburg",'
    '      "Cities":                       ["Halle","Jena","Magdeburg"],'
    '      "Color":                        "#A82117"'
    '    }'
    '  ]'
    ''
    ['The optional "Color" (RGB hex) is the color of the site in the report figure. Entries with "Cities" and ' ...
     '"Color" (without "DeviceSerialNumber") define the color of all scanners of these cities (city of the ' ...
     '"InstitutionAddress", also of scanners that are not in the dictionary, "*" for all other cities) without ' ...
     'own color, e.g. of all scanners of a partner site, where scanners with the same color differ slightly in ' ...
     'brightness. See also DCM2BIDS/sites.json. ']
    ''
    };
    


  % Options:
  % =======================================================================

  % subjectID definition
  subIDform         = cfg_menu;
  subIDform.tag     = 'subIDform';
  subIDform.name    = 'Subject ID form (expert)';
  subIDform.help    = {
    ['Definition of the BIDS subject ID with PatientID (PID) only or as combination ' ...
    'with the SITE, as DeviceSerialNumber or recoded by a fitting entry in the ' ...
    '"Center Dictionary File". The SITE entry is used by default as the PID is given ' ...
    'by a center and might be not unique in multicenter studies. ']
    };
  subIDform.labels  = {
    'sub-PID', ...
    'sub-SITE-PID', ...
    };
  subIDform.values  = {1,2};
  if expert > 1 % extended version 
    % to run multiple studies a study dictionary json/tsv file could be used
    subIDform.name    = 'Subject ID form (developer)';
    subIDform.labels  = {
      'sub-PID', ...
      'sub-SITE-PID', ...
      'sub-STUDY-PID', ...
      'sub-SITE-STUDY-PID', ...
      'sub-STUDY-SITE-PID'
      };
    subIDform.values  = {1,2,3,4,5};
    subIDform.help    = [subIDform.help; {
      'The STUDY is defined by the GUI entry here. '; 
      }];
  end
  subIDform.val     = {def.opts.subIDform};
  subIDform.hidden  = expert<1;



  % === not implemented yet ===
  ProtocolFileName         = cfg_menu;
  ProtocolFileName.tag     = 'ProtocolFileName';
  ProtocolFileName.name    = 'Use Fitting Protocol Filter File Name';
  ProtocolFileName.labels  = {'Yes','No'};
  ProtocolFileName.values  = {1,0};
  ProtocolFileName.val     = {def.opts.ProtocolFileName};
  ProtocolFileName.hidden  = true; %expert<1;
  ProtocolFileName.help    = {
    'Redefine the name of a protocol by the filename of the fitting protocol filter.'
    };

  % study selector/filter - NOT WORKING YET
  studies         = cfg_entry;
  studies.tag     = 'studies';
  studies.name    = 'Study Selector/Filter (expert)';
  studies.strtype = 's';
  studies.num     = [0 Inf];
  studies.val     = {''};
  studies.hidden  = true; %expert<1;
  studies.help    = {
    'Specify the export of specific studies by studyID or the defined study abbreviations.' ''};

  % zipping
  gzipi            = cfg_menu;
  gzipi.tag        = 'gzipi';
  gzipi.name       = 'GZIP Internal Images (expert)';
  gzipi.labels     = {'Yes','No'};
  gzipi.values     = {1,0};
  gzipi.val        = {def.opts.gzipi};
  gzipi.hidden     = expert<1;
  gzipi.help       = {'GZIP the NIFTI images of the internal database (catDCM2BIDSdb) to save space.'};

  gzipe            = cfg_menu;
  gzipe.tag        = 'gzipe';
  gzipe.name       = 'GZIP Output Images';
  gzipe.labels     = {'Yes','No'};
  gzipe.values     = {1,0};
  gzipe.val        = {def.opts.gzipe};
  gzipe.help       = {'GZIP output NIFTI files in the BIDS directories. This might help in case of further SPM preprocesing.'};

  % limit output
  output           = cfg_menu;
  output.tag       = 'output';
  output.name      = 'Output level';
  if expert
    output.labels    = { ...
      'Overview Protocols/Studies (0)', ...
      'Overview + BIDS but only JSON (1)', ...
      'Overview + BIDS of Known Protocols/Studies (2)', ...
      'Overview + BIDS of All Protocols/Studies (3)'};
    output.values    = {0,1,2,3};
    output.val       = {def.opts.output};
    output.help      = {
      ['Use option 0 to import a subset of the data without running further processing and BIDS output ' ...
      'to get an overview of the used protocols and prepare your own protocol filter sets. '] 
      'Option 1 and 2 allows then to prepare the output of data that fits the defined protocols. ' 
      'Option 3 exports all data to BIDS. '
      ['The options 0 and 1 only use the DICOM header information (JSON), i.e., the images are not ' ...
      'converted (they are converted later if required by option 2 or 3). ']
      ''
      };
  else
    output.labels    = { ...
      'Only JSON (1)', ...
      'Only Fitting Protocols/Studies (2)', ...
      'All Protocols/Studies (2)'};
    output.values    = {1,2,3};
    output.val       = {def.opts.output};
    output.help      = {
      ['Use option 1 to fast import data without converting and processing of the image ' ...
       'to get an overview of the used protocols and prepare protocol filter sets. ' ...
       'Option 2 imports and process all data that fits the defined protocols, whereas option 3 consider all protocols. ']};
  end
  

  % anonymized DICOM export (see cat_io_dcm2bids_export)
  dicomexport        = cfg_menu;
  dicomexport.tag    = 'dicomexport';
  dicomexport.name   = 'Anonymized DICOM export';
  dicomexport.labels = {'No','Scans fitting the protocols','All imported scans'};
  dicomexport.values = {0,1,2};
  dicomexport.val    = {def.opts.dicomexport};
  dicomexport.help   = {
   ['Export of the DICOM files of the imported scans with cleaned headers and defaced images into the directory ' ...
    'catDCM2BIDSdb-dicom, e.g. to transfer the data to a central site that runs this batch again. ' ...
    'Option 1 exports the scans with conform or semiconform protocol (all scans without protocols), option 2 all ' ...
    'imported scans (localizers, scouts, derived and 2D images are not imported). ']
    ''
   ['Each run with new scans writes a directory <date-time> with one zip file per session (DICOM files in the ' ...
    'layout of the internal database) and a table of its content. Each scan is exported only once, i.e. the run ' ...
    'directories can be removed after the upload. The tables scans and sessions give an overview of all exports ' ...
    '(without personal data). ']
    ''
   ['Header: the names of the patient (replaced by the PatientID), physicians, and operators, addresses, phone ' ...
    'numbers, other IDs, and comments are cleaned, and the birth date is set to the 15th of the month (anonymization), ' ...
    'whereas the fields used by the batch are kept (e.g. PatientID, sex, age, dates, UIDs, scanner, protocol, and ' ...
    'private fields). A warning is given if parts of the patient name are found in other fields. The PatientID is ' ...
    'not changed and has to be a pseudonym! ']
    'Images: defaced as the BIDS export (anatomical scans, all scans for the extensive anonymization). '
   ['Limitations: only DICOM files with explicit VR little endian are exported, and UIH mosaic images (e.g. DTI ' ...
    'saved by slice) cannot be defaced yet, i.e. they are not exported if they have to be defaced. ']
    ''
    'THE EXPORT IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED. KEEP THE ORIGINAL DATA! '
    ''
    };

  % anonymize .. always required !
  anonymize         = cfg_menu;
  anonymize.tag     = 'anonymize';
  if expert
    anonymize.name    = 'Anonymization level (expert)';
    anonymize.labels  = {'No','Pseudonymized','Anonymized'}; 
    anonymize.values  = {0,1,2};
  else
    anonymize.name    = 'Anonymization level';
    anonymize.labels  = {'Pseudonymized','Anonymized'};
    anonymize.values  = {1,2};
  end
  anonymize.val     = {def.opts.anonymize};
  anonymize.help    = {
    'Anonymization of the BIDS export (directory BIDS-<level>) and the DICOM export: '
    ''
   ['No (0): personal data (subject labels with the PatientID, session labels with date and time, images with faces), ' ...
    'i.e. not for sharing (BIDS directory without level name). ']
   ['Pseudonymized (1): subject labels with the site and the PatientID (that has to be a pseudonym), session labels ' ...
    'with the scan date, age with one decimal, defaced anatomical images, no acquisition times in the NIfTI headers. ' ...
    'The DICOM export keeps the PatientID, dates, and UIDs, but cleans names, addresses, IDs, and comments (birth date ' ...
    'at the 15th of the month). ']
   ['Anonymized (2): subject labels with the site and a number (e.g. sub-BerSI30T00001), session labels with the time ' ...
    'since the first session of the subject (e.g. ses-w0000 and ses-w0052 one year later, unit in weeks by default, ' ...
    'see opts.sesunit in cat_io_dcm2bids_defaults), age in years, defaced images of all scans. The DICOM export uses ' ...
    'the subject and session labels as PatientID, PatientName, and StudyID, shifts the dates of a session to the ' ...
    'date of its label (first session at 1900-01-01) and its times to the time since its start (00:00:00, i.e. the ' ...
    'timing of scans and slices is kept), replaces the UIDs (also in private fields) and the scanner ' ...
    'and institution by the site abbreviation and city, removes the station, department, and study and procedure ' ...
    'descriptions, and replaces parts of the patient name, the PatientID, and the serial number in other fields by ' ...
    'X. The key (numbers and date ' ...
    'offsets) stays in the database (catDCM2BIDSdb/catDCM2BIDSanonymization.json) and the private tables, i.e. the ' ...
    'creator can still identify the subjects (e.g. in case of findings). The import of an anonymized DICOM export ' ...
    'gives the labels of the anonymized BIDS export (dates from 1850 to 1949 are treated as anonymized). ']
    ''
   ['Limitations: the anonymization is a best effort of the software without any warranty. Defaced (and also ' ...
    'skull-stripped) images, rare phenotypes, small sites, the phenotypic data of the user tables, and other ' ...
    'information can allow a re-identification. THE USER IS RESPONSIBLE FOR CHECKING THE DATA AND FOR THE DATA ' ...
    'PROTECTION RULES OF THE STUDY BEFORE SHARING. ']
    ''
    };
 
  % preprocessing
  preprocessing         = cfg_menu;
  preprocessing.tag     = 'preprocessing';
  preprocessing.name    = 'Run basic preprocessing';
  preprocessing.labels  = {'No','Yes (only report tables)','Yes (also images)'}; 
  preprocessing.values  = {0,1,2};
  preprocessing.val     = {def.opts.preprocessing};
  preprocessing.help    = { ...
   ['Run SPM preprocessing to estimate brain tissue volumes (anat), ' ...
    'diffusivity (FA/AD) and functional connectivity using SPM. ' ...
    'Export tables or also the maps: the segmentation of one anatomical scan per session (T1w, otherwise ' ...
    'UNIT1, T2w, or FLAIR, of several scans the one with the best quality rating: label map c0, bias-corrected ' ...
    'image m, and the modulated normalized GM map mwc1), the DTI maps (FA, the direction-encoded color map DEC ' ...
    'in head orientation, and the normalized FA wFA, e.g. for voxel-based analyses, optionally also MD, AD, and ' ...
    'RD, and the DEC map in the space of the T1w), and the distortion-corrected mean EPI of the rs-fMRI with the ' ...
    'connectivity matrix. The DWI and rs-fMRI scans are coregistered rigidly to the segmented scan of their ' ...
    'session and non-linearly only if it clearly improves the overlap (as the non-linear co-registration of CAT, ' ...
    'e.g. for remaining distortions of an EPI without distortion correction; the Dice overlap of the brain masks ' ...
    'is saved in the JSON and shown in the renderings, with a warning below 0.8). ']};

  % preprocessing level, i.e. the accuracy and processing time of the DTI 
  % and rs-fMRI measures (see levelParameters in cat_io_dcm2bids_pp_anat): the 
  % default user selects the optimized level or the preview, whereas the 
  % expert enters a (continuous) level with the evaluation table as help
  if expert
    pplevel             = cfg_entry;
    pplevel.tag         = 'pplevel';
    pplevel.name        = 'Preprocessing accuracy level (expert)';
    pplevel.strtype     = 'r';
    pplevel.num         = [1 1];
    pplevel.val         = {def.opts.pplevel};
    pplevel.help        = { ...
      ['Accuracy and processing time of the DTI and rs-fMRI measures as continuous level from 1 to 4 (the ' ...
       'segmentation of the anatomical scans is the same for all levels, see levelParameters in cat_io_dcm2bids_pp_anat): ']
      '  1 - full quality: SPM default parameters of the distortion correction (SCOPE) at the native resolution (reference)'
      '  2 - optimized: fast settings of SCOPE at 2.5 mm'
      '  3 - quick review: SCOPE at 4 mm and only the first DWI and rs-fMRI scan per session (one segmentation per session for all levels)'
      '  4 - preview: SCOPE at 6 mm (rough maps)'
      ''
      ['Evaluation with a pilot session (T1w, 2 DWI, and rs-fMRI scans): processing time, differences to level 1 ' ...
       '(WM FA, tSNR, and segregation), and correlation with level 1 (FA map and connectivity matrix): ']
      '  level           time    WM FA   r(FA map)  tSNR   segregation  r(conn)'
      '  1 full quality  26 min  0.372     1        37.8     0.731        1'
      '  1.5             17 min  +0.000    0.98     -0.1    +0.007        0.999'
      '  2 optimized      9 min  +0.001    0.98     +0.0    +0.019        0.999'
      '  2.5              7 min  +0.001    0.96     -0.1    +0.014        0.998'
      '  3 quick review   4 min  +0.001    0.95     +0.1    +0.018        0.997'
      '  3.5            2.5 min  +0.001    0.95     -0.2    +0.007        0.996'
      '  4 preview      2.5 min  +0.001    0.93     +0.2    +0.008        0.996'
      ''
      ['I.e., the global measures are the same for all levels, whereas the maps differ slightly. Nevertheless, ' ...
       'the data of a study should be processed with the same level. Results of a lower level (higher accuracy) ' ...
       'are kept. ']
      };
  else
    pplevel             = cfg_menu;
    pplevel.tag         = 'pplevel';
    pplevel.name        = 'Preprocessing accuracy';
    pplevel.labels      = {'Optimized','Preview'};
    pplevel.values      = {2,4};
    [~,li]              = min(abs([2 4] - def.opts.pplevel)); 
    pplevel.val         = pplevel.values(li);
    pplevel.help        = { ...
      ['Accuracy and processing time of the DTI and rs-fMRI measures: "Optimized" for the final results ' ...
       '(about 9 min for a session with T1w, 2 DWI, and rs-fMRI scans) and "Preview" for a fast review ' ...
       '(about 2.5 min, with the same global measures, e.g. the white matter FA, but rougher maps). The data of ' ...
       'a study should be processed with the same accuracy. Results of the higher accuracy are kept. ']};
  end

  % DTI maps of the preprocessing (see cat_io_dcm2bids_pp_dwi)
  dtimaps         = cfg_menu;
  dtimaps.tag     = 'dtimaps';
  dtimaps.name    = 'DTI maps (expert)';
  dtimaps.labels  = {'FA, DEC, and wFA','Also MD, AD, and RD'};
  dtimaps.values  = {1,2};
  dtimaps.val     = {def.opts.dtimaps};
  dtimaps.hidden  = expert < 1; 
  dtimaps.help    = { ...
    ['DTI maps of the basic preprocessing with images: the fractional anisotropy FA, the direction-encoded ' ...
     'color map DEC, and the normalized FA wFA (e.g. for voxel-based analyses), optionally also the mean, axial, ' ...
     'and radial diffusivity MD, AD, and RD (their white matter values are always estimated for the tables). ']
    };

  protocolsubdirs         = cfg_menu;
  protocolsubdirs.tag     = 'protocolsubdirs';
  protocolsubdirs.name    = 'Use subdirectories to separate protocols (expert)';
  protocolsubdirs.labels  = {'No','Yes'}; 
  protocolsubdirs.values  = {0,1};
  protocolsubdirs.val     = {def.opts.protocolsubdirs};
  protocolsubdirs.hidden  = expert < 1; 
  protocolsubdirs.help    = { ...
    'Use subdirectories to separate the protocols of the main protocol directories. '
    };

  tolerance          = cfg_entry;   
  tolerance.tag      = 'tolerance';
  tolerance.name     = 'Tolerance';
  tolerance.help     = {
      'Tolerance of numeric protocol parameters in percent. The same protocol results in the same values (also for different scanner software versions), i.e., the tolerance is only required for the rounding of values (e.g. a pixel bandwidth of 2289 or 2290 Hz/px, i.e., 0.04%). Larger values would hide real changes of the protocol, e.g., an inversion time of 1000 instead of 1030 ms (3%). '
    };
  tolerance.strtype  = 'e';
  tolerance.num      = [1 1];
  tolerance.hidden   = expert<1; 
  tolerance.val      = {def.opts.tolerance};

  % rendering of slices of each scan (BIDS-report/render/..., see 
  % cat_io_dcm2bids_render): the aim is the visual identification of outliers 
  % in large sets of similar images, i.e., all scans of one image type are 
  % rendered at the same MNI positions. The default user only decides about 
  % the rendering (the rendering is relatively fast and prepares all outputs 
  % at once), the subfields are expert options to adapt the organization. 
  rsource         = cfg_menu;
  rsource.tag     = 'source';
  rsource.name    = 'Data (expert)';
  rsource.labels  = {'Raw BIDS data','Derivatives (catDCM2BIDS)','Both'};
  rsource.values  = {1,2,3};
  rsource.val     = {def.opts.render.source};
  rsource.hidden  = expert<1;
  rsource.help    = {
    'Render the raw BIDS images, the derivatives (e.g. the segmentation), or both. '};

  rsessions         = cfg_menu;
  rsessions.tag     = 'sessions';
  rsessions.name    = 'Sessions (expert)';
  rsessions.labels  = {'Only first session per subject','All sessions'};
  rsessions.values  = {1,0};
  rsessions.val     = {def.opts.render.sessions};
  rsessions.hidden  = expert<1;
  rsessions.help    = {
    'Render only the first session of each subject or all sessions. '};

  rslicemode        = cfg_menu;
  rslicemode.tag    = 'slicemode';
  rslicemode.name   = 'MNI registration (expert)';
  rslicemode.labels = {'Affine','Rigid'};
  rslicemode.values = {'mni','mnirigid'};
  rslicemode.val    = {def.opts.render.slicemode};
  rslicemode.help   = {
   ['The slices are defined in MNI space using the affine registration of the session ' ...
    '(from the SPM segmentation or an additional affine registration of an anatomical scan), ' ...
    'either with the full affine transformation that also scales the brain to a similar size ' ...
    '(default) or only its rigid part that keeps the original size. Scans without registration ' ...
    'are not rendered (black tile marked by "no MNI"). ' ...
    'Colored lines show the coordinate planes of the original (scanner) space: x=0 (red), y=0 (green), and z=0 (blue). ']
    };

  rsort           = cfg_menu;
  rsort.tag       = 'sort';
  rsort.name      = 'Order of scans (expert)';
  rsort.labels    = {'Subject and session name','Quality rating (worst first)'};
  rsort.values    = {'name','SQR'};
  rsort.val       = {def.opts.render.sort};
  rsort.hidden    = expert<1;
  rsort.help      = {
   ['Order of the scans on the overview and scan-row pages, either by subject and session name, ' ...
    'or the subjects by their worst overall quality rating (SQR) with unrated subjects at the end. ' ...
    'The sessions of a subject stay together in both cases. Note that the order by rating can ' ...
    'differ between image types and that the rating cannot be fully trusted. ']};

  % R1: overview pages with one slice per scan
  rtiles          = cfg_menu;
  rtiles.tag      = 'tiles';
  rtiles.name     = 'Overview tiles per page (expert)';
  rtiles.labels   = {'3x4','4x5','5x7','6x8'};
  rtiles.values   = {[3 4],[4 5],[5 7],[6 8]};
  rtiles.val      = {def.opts.render.tiles};
  rtiles.hidden   = expert<1;
  rtiles.help     = {
    'Number of tiles (columns x rows) per overview page in A4 portrait format. '};

  rslorient        = cfg_menu;
  rslorient.tag    = 'orient';
  rslorient.name   = 'Orientation';
  rslorient.labels = {'Axial','Coronal','Sagittal'};
  rslorient.values = {3,2,1};
  rslorient.val    = {3};
  rslorient.help   = {'Orientation of the slice. '};

  rslpos           = cfg_entry;
  rslpos.tag       = 'pos';
  rslpos.name      = 'Position (mm)';
  rslpos.strtype   = 'r';
  rslpos.num       = [1 1];
  rslpos.val       = {10};
  rslpos.help      = {'MNI coordinate of the slice in mm, i.e., x for sagittal, y for coronal, and z for axial slices. '};

  rslice           = cfg_branch;
  rslice.tag       = 'slices'; % the repeat is harvested with this tag as structure array
  rslice.name      = 'Slice';
  rslice.val       = {rslorient, rslpos};
  rslice.help      = {'Orientation and MNI position of a slice. '};

  % default slices (axial z=10, coronal y=0, and sagittal x=0)
  defsl  = def.opts.render.slices; 
  rslval = cell(1,size(defsl,1)); 
  for si = 1:size(defsl,1)
    rslo = rslorient; rslo.val = {defsl(si,1)}; 
    rslp = rslpos;    rslp.val = {defsl(si,2)}; 
    rslval{si} = rslice; rslval{si}.val = {rslo, rslp}; 
  end
  rslices          = cfg_repeat;
  rslices.tag      = 'slices';
  rslices.name     = 'Overview slices (expert)';
  rslices.values   = {rslice};
  rslices.val      = rslval;
  rslices.num      = [0 Inf];
  rslices.hidden   = expert<1;
  rslices.help     = {
   ['Slices of the overview pages, where each slice gives its own pages with one tile per scan. ' ...
    'The overview pages give a very brief overview of many images but are biased by the selected ' ...
    'slice, i.e., they are useful for a rough global review. The default slices are ' ...
    'axial z=10 (basal ganglia and ventricles), coronal y=0 (subcortical structures), and ' ...
    'sagittal x=0 (between the hemispheres, shows the offset of the origin). Further interesting ' ...
    'slices are coronal y=-60 (symmetric cut through the cerebellum) and sagittal x=-30 and x=30 ' ...
    '(frontal, parietal, temporal, and cerebellar areas and especially the hippocampus of both ' ...
    'hemispheres). No slice means no overview pages. Each tile shows the subject, the session with the ' ...
    'quality rating, and the scan (e.g. the acquisition label with series number and protocol). ']};
  clear defsl rslval rslo rslp si

  % R2: scan-row pages with one row of slices per scan
  rrows           = cfg_menu;
  rrows.tag       = 'rows';
  rrows.name      = 'Scan-row pages (expert)';
  rrows.labels    = {'No','Yes'};
  rrows.values    = {0,1};
  rrows.val       = {def.opts.render.rows};
  rrows.hidden    = expert<1;
  rrows.help      = {
   ['Pages of each image type with one row per scan. Each scan has a header with its BIDS filename ' ...
    'followed by an information panel (quality ratings, image size, acquisition parameters, and ' ...
    'registration) and the slices axial z=10, sagittal x=0, and coronal y=0. ' ...
    'There are 6 scans per A4 portrait page. ']};

  % R3: subject reports with all scans of a subject
  rsubjects        = cfg_menu;
  rsubjects.tag    = 'subjects';
  rsubjects.name   = 'Subject reports (expert)';
  rsubjects.labels = {'No','Yes'};
  rsubjects.values = {0,1};
  rsubjects.val    = {def.opts.render.subjects};
  rsubjects.hidden = expert<1;
  rsubjects.help   = {
   ['One PDF file per subject (BIDS-report/render/[protocol/]subjects/render_sub-*.pdf) with the scan ' ...
    'rows (as the scan-row pages) of selected images of the subject: the original anatomical scans of the ' ...
    'segmentation (e.g. T1w), the label map of the segmentation, the normalized GM map, the ' ...
    'direction-encoded color map of the DWI on the T1w, and the mean EPI with the network of the rs-fMRI ' ...
    '(similar images of several scans together), otherwise all scans ordered by session, datatype (anat, ' ...
    'dwi, func, perf, others), and name. Field maps are not rendered. ']};

  norender        = cfg_const;
  norender.tag    = 'norender';
  norender.name   = 'No';
  norender.val    = {0};
  norender.help   = {'No rendering. '};

  dorender        = cfg_branch;
  dorender.tag    = 'render';
  dorender.name   = 'Yes';
  dorender.val    = {rsource, rsessions, rslicemode, rsort, rtiles, rslices, rrows, rsubjects};
  dorender.help   = {'Render slices of each scan. '};
  
  if expert
    render          = cfg_choice;
    render.tag      = 'render';
    render.name     = 'Render slices';
    render.values   = {norender, dorender};
    if def.opts.render.run, render.val = {dorender}; else, render.val = {norender}; end
  else
    render        = cfg_menu;
    render.tag    = 'render';
    render.name   = 'Render reports';
    render.labels = {'No','Yes'};
    render.values = {0,1};
    render.val    = {1};
    render.hidden = true; % it is nearly always required/helpful and relative cheap
    render.help   = {'Render slices of each scan. '};
  end
  render.help     = {
   ['Render slices of each scan in MNI space into PNG pages in the BIDS-report directory ' ...
    '(BIDS-report/render/[protocol/]datatype/subtype/) for the visual identification of outliers ' ...
    'in large sets of similar images. The scans of each image type are shown (1) on overview pages ' ...
    'with one slice per scan (one set of pages per slice) and (2) on scan-row pages with three ' ...
    'slices and further information per scan. In addition, (3) a PDF report per subject shows ' ...
    'the scan rows of all scans of the subject. ' ...
    'The labels are colored by the overall quality rating. 4D data is represented by its first volume. ']};


  % further possible parameter:
  %%%%%%%%  
  % - avoid/add BIDS field in the JSON files to further specify anonymizing settings?
  %   as string with +SubjectID or -SubjectID 
  % - use the Database as Project name also in the SubjectID
  %   > handling of multiple studies by dictionary
  % - flat to recode subjectIDs 
  %%%%%%%%  
 


  % dataset description (dataset_description.json, README.md)
  % =======================================================================
  dsname           = cfg_entry;
  dsname.tag       = 'Name';
  dsname.name      = 'Name';
  dsname.strtype   = 's';
  dsname.num       = [0 Inf];
  dsname.val       = {def.dataset.Name};
  dsname.help      = {'Name of the dataset (BIDS "Name"). If empty, the study subdirectory name is used. '};

  dsauthors        = cfg_entry;
  dsauthors.tag    = 'Authors';
  dsauthors.name   = 'Authors';
  dsauthors.strtype = 's+';
  dsauthors.num    = [0 Inf];
  dsauthors.val    = {{''}};
  dsauthors.help   = {
   ['List of the authors of the dataset (BIDS "Authors"), one author per line, ideally with e-mail ' ...
    'address, e.g. "Jane Doe <jane.doe@uni-example.de>". ']};

  dshowto          = cfg_entry;
  dshowto.tag      = 'HowToAcknowledge';
  dshowto.name     = 'How to acknowledge';
  dshowto.strtype  = 's';
  dshowto.num      = [0 Inf];
  dshowto.val      = {def.dataset.HowToAcknowledge};
  dshowto.help     = {
   ['Text that describes how to acknowledge the dataset in publications (BIDS "HowToAcknowledge"), ' ...
    'e.g. a reference or a sentence for the acknowledgements. ']};

  dslicense        = cfg_menu;
  dslicense.tag    = 'License';
  dslicense.name   = 'License';
  dslicense.labels = {
    'CC0-1.0 (public domain, no restriction, most common)', ...
    'PDDL-1.0 (public domain for databases)', ...
    'CC-BY-4.0 (free use with attribution of the authors)', ...
    'CC-BY-SA-4.0 (as CC-BY, derived data under the same license', ...
    'CC-BY-NC-4.0 (as CC-BY, but no commercial use)', ...
    'ODC-BY-1.0 (free use of the database with attribution)', ...
    'other / not specified'};
  dslicense.values = {'CC0-1.0','PDDL-1.0','CC-BY-4.0','CC-BY-SA-4.0','CC-BY-NC-4.0','ODC-BY-1.0',''};
  dslicense.val    = {def.dataset.License};
  dslicense.help   = {
   ['License of the dataset (BIDS "License") as SPDX identifier, as recommended by BIDS. For the import ' ...
    'of BIDS datasets, the most restrictive license of this setting and the source datasets is used ' ...
    '(with a warning). ']
    ''
    '  CC0-1.0 .. Creative Commons Zero: public domain, no restrictions (most common for open data)'
    '  PDDL-1.0 .. Open Data Commons Public Domain Dedication: public domain for databases'
    '  CC-BY-4.0 .. Creative Commons Attribution: free use with attribution of the authors'
    '  CC-BY-SA-4.0 .. Creative Commons Attribution-ShareAlike: as CC-BY, derived data under the same license'
    '  CC-BY-NC-4.0 .. Creative Commons Attribution-NonCommercial: as CC-BY, but no commercial use'
    '  ODC-BY-1.0 .. Open Data Commons Attribution: free use of the database with attribution'
    '  other / not specified .. no license entry (e.g. data under a data use agreement)'
    };

  dsreadme         = cfg_files;
  dsreadme.tag     = 'README';
  dsreadme.name    = 'README file';
  dsreadme.filter  = 'any';
  dsreadme.ufilter = '.*\.(txt|md|TXT|MD)$';
  dsreadme.num     = [0 1];
  dsreadme.val     = {{''}};
  dsreadme.help    = {
   ['Text file (txt/md) that defines the head of the README.md file of the BIDS dataset, e.g. with a ' ...
    'description of the study. For the import of BIDS datasets, it is extended by their README files ' ...
    'with the dataset name as header. ']};

  dsack            = cfg_entry;
  dsack.tag        = 'Acknowledgements';
  dsack.name       = 'Acknowledgements (expert)';
  dsack.strtype    = 's';
  dsack.num        = [0 Inf];
  dsack.val        = {def.dataset.Acknowledgements};
  dsack.hidden     = expert<1;
  dsack.help       = {'Text acknowledging contributions of individuals or institutions beyond the authors (BIDS "Acknowledgements"). '};

  dsfunding        = cfg_entry;
  dsfunding.tag    = 'Funding';
  dsfunding.name   = 'Funding (expert)';
  dsfunding.strtype = 's+';
  dsfunding.num    = [0 Inf];
  dsfunding.val    = {{''}};
  dsfunding.hidden = expert<1;
  dsfunding.help   = {'List of the sources of funding, e.g. grant numbers, one per line (BIDS "Funding"). '};

  dsethics         = cfg_entry;
  dsethics.tag     = 'EthicsApprovals';
  dsethics.name    = 'Ethics approvals (expert)';
  dsethics.strtype = 's+';
  dsethics.num     = [0 Inf];
  dsethics.val     = {{''}};
  dsethics.hidden  = expert<1;
  dsethics.help    = {'List of the ethics committee approvals of the research protocols, one per line, e.g. committee and approval number (BIDS "EthicsApprovals"). '};

  dsrefs           = cfg_entry;
  dsrefs.tag       = 'ReferencesAndLinks';
  dsrefs.name      = 'References and links (expert)';
  dsrefs.strtype   = 's+';
  dsrefs.num       = [0 Inf];
  dsrefs.val       = {{''}};
  dsrefs.hidden    = expert<1;
  dsrefs.help      = {'List of references to publications that contain information on the dataset, or links, one per line (BIDS "ReferencesAndLinks"). '};

  dsdoi            = cfg_entry;
  dsdoi.tag        = 'DatasetDOI';
  dsdoi.name       = 'Dataset DOI (expert)';
  dsdoi.strtype    = 's';
  dsdoi.num        = [0 Inf];
  dsdoi.val        = {def.dataset.DatasetDOI};
  dsdoi.hidden     = expert<1;
  dsdoi.help       = {'The Digital Object Identifier of the dataset (not the corresponding paper), e.g. "doi:10.18112/openneuro.ds000001.v1.0.0" (BIDS "DatasetDOI"). '};

  % phenotypic data of the sessions (see cat_io_dcm2bids_phenotypic)
  dspheno          = cfg_files;
  dspheno.tag      = 'phenotypic';
  dspheno.name     = 'Phenotypic data tables';
  dspheno.filter   = 'any';
  dspheno.ufilter  = '.*\.(csv|tsv|CSV|TSV)$';
  dspheno.num      = [0 Inf];
  dspheno.val      = {{''}};
  dspheno.help     = {
   ['CSV or TSV tables with data of the scan sessions (e.g. clinical or behavioral data) with one row per session ' ...
    '(or per subject). The rows are matched to the sessions of the internal database by their key columns, and ' ...
    'all other columns are added to the BIDS sessions files (sub-<label>_sessions.tsv). The key columns are not ' ...
    'exported, but please check the other columns for personal data. ']
    ''
   ['Key columns (if not given by their numbers): PatientID (patientid, pid), name (name, patientname, firstname, ' ...
    'lastname, surname, vorname, nachname; the words of all name columns and of the DICOM name in any order), birth ' ...
    'date (birthdate, dateofbirth, dob, geburtsdatum; year and month), scan date (scandate, date, datetime, ' ...
    'acquisitiondate, datum; day and the closest time), and visit (visit, session, timepoint; number of the session ' ...
    'of the subject by date). Dates as yyyy-mm-dd, yyyymmdd, dd.mm.yyyy, or dd/mm/yyyy (with optional time). ']
    ''
   ['Rows with a scan date or visit have to match exactly one session, rows without all sessions of one subject. ' ...
    'Sessions with several rows of a table get no data of this table. The data of a table are replaced by each ' ...
    'import of the table (by its file name). The match of each row is given in the private table ' ...
    'catDCM2BIDSdb-tables/phenotypic_<table>.csv. ']
    ''
    'E.g. a "visits.csv" with:'
    '  name, surname, scandate, MMSE, group'
    '  Max, Mustermann, 2025-11-20, 29, control'
    ''
    };

  dsphenokeys         = cfg_entry;
  dsphenokeys.tag     = 'phenotypickeys';
  dsphenokeys.name    = 'Key columns of the phenotypic tables';
  dsphenokeys.strtype = 'n';
  dsphenokeys.num     = [0 Inf];
  dsphenokeys.val     = {def.dataset.phenotypickeys};
  dsphenokeys.help    = {
   ['Numbers of the key columns of the phenotypic data tables for the matching to the sessions, e.g. 1:3 for name, ' ...
    'surname, and scan date. Empty: detection of the key columns by their names (see the help of the tables). ' ...
    'Key columns with other names are dates (scan or birth date), visits (numbers up to 99), or names. ']
    };

  dataset          = cfg_branch;
  dataset.tag      = 'dataset';
  dataset.name     = 'Dataset Descriptor';
  dataset.val      = {dsauthors, dsname, dshowto, dslicense, dsreadme, dsack, dsfunding, dsethics, dsrefs, dsdoi, ...
                      dspheno, dsphenokeys};
  dataset.help     = {
   ['Description of the BIDS dataset that is written to the dataset_description.json and README.md ' ...
    'files of the BIDS directory. BIDSVersion, DatasetType, and GeneratedBy (catDCM2BIDS and dcm2niix) are ' ...
    'set automatically. ']};


  % main fields
  % =======================================================================
  dicts            = cfg_branch;
  dicts.tag        = 'dicts';
  dicts.name       = 'Dictionary Files';
  dicts.val        = {protocoldir, centerdict, studydict, subjectdict, };
  dicts.help       = {'Optional files to replace standard variables by more handy names or codes. '}; 

  % report figure (see cat_io_dcm2bids_reportfig): designs of DCM2BIDS/colormaps
  % (all or one, other combinations by the defaults or a batch script)
  Pdes            = dir(fullfile(fileparts(mfilename('fullpath')),'DCM2BIDS','colormaps','*.json'));
  dnames          = regexprep({Pdes.name},'\.json$','');
  rdesign         = cfg_menu;
  rdesign.tag     = 'design';
  rdesign.name    = 'Designs';
  rdesign.labels  = [{'All designs'} strrep(dnames,'_',' ')];
  rdesign.values  = [{dnames} cellfun(@(d) {d},dnames,'UniformOutput',false)];
  ddef            = cellstr(def.opts.reportfig.design);
  if ~any(cellfun(@(v) isequal(v,ddef),rdesign.values))
    rdesign.labels = [{'Default designs'} rdesign.labels]; rdesign.values = [{ddef} rdesign.values];
  end
  rdesign.val     = {ddef};
  rdesign.help    = {
   ['Designs of the report figure, i.e. one PDF file per design (BIDS-report/report_<BIDSsubdir>_<design>.pdf): ' ...
    'the bright and dark CAT designs and the bright DZPG design with the logos, rating icons, colors, and ' ...
    'colormaps of DCM2BIDS/colormaps/<design>.json (further designs can be added there). ']};

  reportfig       = cfg_branch;
  reportfig.tag   = 'reportfig';
  reportfig.name  = 'Report figure';
  reportfig.val   = {rdesign};
  reportfig.help  = {
   ['A4 pages with the diagrams of the report tables of each BIDS directory (overview and one page per site): ' ...
    'number of sessions per month and per site, age and sex of the subjects, brain measures by age, sample ' ...
    'homogeneity, and quality ratings (see cat_io_dcm2bids_reportfig). ']};

  opts            = cfg_branch;
  opts.tag        = 'opts';
  opts.name       = 'Options';
  opts.val        = {ProtocolFileName, studies, anonymize, tolerance, preprocessing, pplevel, dtimaps, ...
                      gzipi, gzipe, output, dicomexport, protocolsubdirs, subIDform, render, reportfig};
  opts.help       = {['Parameters to control the processing and output of ' ...
                      'the converted data and its export to BIDS. ']}; 


  % batch
  % =======================================================================
  dcm2bids        = cfg_exbranch;
  dcm2bids.tag    = 'dcm2bids';
  dcm2bids.name   = 'DICOM2BIDS';
  dcm2bids.prog   = @cat_io_dcm2bids;
  dcm2bids.vout   = @vout_io_dcm2bids;
  dcm2bids.val    = {datadir, outdir, subdir, dicts, dataset, opts}; 
  dcm2bids.help   = { ...
   ['This batch uses DCM2NIIX to convert DICOM into NIFTI images with JSON sidecars. ' ...
    'It stores the converted data and reorganizes the output in BIDS. ' ...
    'It allows filtering for specific MRI protocols within a directory that outline relevant MR parameters within a JSON file. '] 
    ''
   ['Please check out the CAT subdirectory DCM2BIDS/DZPG3T for an example that includes the definition for ' ...
    'structural, functional, and diffusion scans used by the DZPG (German Center for Mental Health, https://www.dzpg.org/). ']
    ''
    'The batch also applies the SPM anonymization routine and runs a basic image quality control. '
    ''
   ['The anonymized DICOM export (option) writes the DICOM files of the imported scans with cleaned headers and ' ...
    'defaced images as zip files of the sessions, e.g. to transfer the data to a central site that runs this batch again. ']
    'THE DICOM EXPORT IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED! '
    };
end
function cdep = vout_io_dcm2bids(job)
%vout_io_dcm2bids. Dependencies for the converted raw BIDS images.
%  SPM defines the dependencies before the batch runs, i.e. without knowing
%  the data. Therefore, a fixed list of BIDS datatypes and suffixes is used
%  (see cat_io_dcm2bids('outputs')) that may also be empty after processing.
%  4D images are given as files (no volume expansion). 

  % gzipped output files (CAT can read them, but most SPM functions not)
  gz = '';
  if job.opts.gzipe, gz = ' (.gz)'; end

  out  = cat_io_dcm2bids('outputs'); 
  cdep = cfg_dep; ci = 0; 
  dts  = fieldnames(out); 
  for di = 1:numel(dts)
    sxs = fieldnames(out.(dts{di})); 
    for si = 1:numel(sxs)
      ci = ci + 1; 
      if ci > 1, cdep(ci) = cfg_dep; end
      cdep(ci).sname      = sprintf('%s %s%s', dts{di}, sxs{si}, gz);
      cdep(ci).src_output = substruct('.',dts{di},'.',sxs{si});
      cdep(ci).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end
end
