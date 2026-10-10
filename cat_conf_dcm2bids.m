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
   ['  ' b 'BIDS-report           .. overview tables']
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
    'The file name of the JSON filter file together with the sequence number ' ...
    'will be used to specify the BIDS acquisition field "acq-###-FILENAME"']
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
  centerdict          = cfg_files;
  centerdict.tag      = 'Pcenterdict';
  centerdict.name     = 'Center Dictionary File';
  centerdict.filter   = 'any';
  centerdict.ufilter  = 'site.*\.json$';
  centerdict.val      = {{''}};
  centerdict.num      = [0 Inf];
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
    '      "InstitutionAbbreviation":      "JE"'
    '    }'
    '    {'
    '      "DeviceSerialNumber":           "000007",'
    '      "InstitutionAbbreviation":      "NA"'
    '    }'
    '  ]'
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
  

  % anonymize .. always required !
  anonymize         = cfg_menu;
  anonymize.tag     = 'anonymize';
  if expert
    anonymize.name    = 'Anonymization level (expert)';
    anonymize.labels  = {'No','Yes - basic','Yes - extensive','Yes - extreme'}; 
    anonymize.values  = {0,1,2,3};
  else
    anonymize.name    = 'Anonymization level';
    anonymize.labels  = {'Basic','Extensive'};
    anonymize.values  = {1,2};
  end
  anonymize.val     = {def.opts.anonymize};
  anonymize.help    = {'Strength of the anonymization of DICOM header and image information. '};
 
  % preprocessing
  preprocessing         = cfg_menu;
  preprocessing.tag     = 'preprocessing';
  preprocessing.name    = 'Run basic SPM preprocessing (expert)';
  preprocessing.labels  = {'No','Yes','Yes (export to derivatives)'}; 
  preprocessing.values  = {0,1,2};
  preprocessing.val     = {def.opts.preprocessing};
  preprocessing.help    = { ...
   ['Run SPM preprocessing to estimate brain tissue volumes (anat), ' ...
    'diffusivity (FA/AD) and functional connectivity using SPM. ' ...
    'Export realigned maps. ']};

  % optimizing: denoising, bias-correction?, resampling to MNI with specific BB
  %{
  optimizing         = cfg_menu;
  optimizing.tag     = 'optimizing';
  optimizing.name    = 'Run optimization (expert)';
  optimizing.labels  = {'No','Yes'}; 
  optimizing.values  = {0,1};
  optimizing.val     = {0};
  optimizing.hidden  = expert<1; 
  optimizing.help    = { ...
    'Run optimization with denoising, slice-motion/bias-correction, and resampling in MNI.'};
  %}

  % preprocessing
  protocolsubdirs         = cfg_menu;
  protocolsubdirs.tag     = 'protocolsubdirs';
  protocolsubdirs.name    = 'Use subdirectories to separate protocols';
  protocolsubdirs.labels  = {'No','Yes'}; 
  protocolsubdirs.values  = {0,1};
  protocolsubdirs.val     = {def.opts.protocolsubdirs};
  protocolsubdirs.help    = { ...
    'Use subdirectories to separate the protocols of the main protocol directories. '
    };

  tolerance          = cfg_entry;   
  tolerance.tag      = 'tolerance';
  tolerance.name     = 'Tolerance';
  tolerance.help     = {
      'Tolerance level of input parameters in percent. E.g., if the TR time should be 2.1 but is 2.2.'
    };
  tolerance.strtype  = 'e';
  tolerance.num      = [1 1];
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
    'hemispheres). No slice means no overview pages. ']};
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
    'rows (as the scan-row pages) of all scans of the subject ordered by session, datatype (anat, dwi, ' ...
    'func, fmap, perf, others), and name. ']};

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

  render          = cfg_choice;
  render.tag      = 'render';
  render.name     = 'Render slices';
  render.values   = {norender, dorender};
  if def.opts.render.run, render.val = {dorender}; else, render.val = {norender}; end
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
  dslicense.labels = {'CC0-1.0','PDDL-1.0','CC-BY-4.0','CC-BY-SA-4.0','CC-BY-NC-4.0','ODC-BY-1.0','other / not specified'};
  dslicense.values = {'CC0-1.0','PDDL-1.0','CC-BY-4.0','CC-BY-SA-4.0','CC-BY-NC-4.0','ODC-BY-1.0',''};
  dslicense.val    = {def.dataset.License};
  dslicense.help   = {
   ['License of the dataset (BIDS "License") as SPDX identifier, as recommended by BIDS. For the import ' ...
    'of BIDS datasets, the most restrictive license of this setting and the source datasets is used ' ...
    '(with a warning). ']
    ''
    '  CC0-1.0      .. Creative Commons Zero: public domain, no restrictions (most common for open data)'
    '  PDDL-1.0     .. Open Data Commons Public Domain Dedication: public domain for databases'
    '  CC-BY-4.0    .. Creative Commons Attribution: free use with attribution of the authors'
    '  CC-BY-SA-4.0 .. Creative Commons Attribution-ShareAlike: as CC-BY, derived data under the same license'
    '  CC-BY-NC-4.0 .. Creative Commons Attribution-NonCommercial: as CC-BY, but no commercial use'
    '  ODC-BY-1.0   .. Open Data Commons Attribution: free use of the database with attribution'
    '  other / not specified .. no license entry (e.g. data under a data use agreement)'
    };

  dsreadme         = cfg_files;
  dsreadme.tag     = 'README';
  dsreadme.name    = 'README head';
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

  dataset          = cfg_branch;
  dataset.tag      = 'dataset';
  dataset.name     = 'Dataset';
  dataset.val      = {dsname, dsauthors, dshowto, dslicense, dsreadme, dsack, dsfunding, dsethics, dsrefs, dsdoi};
  dataset.help     = {
   ['Description of the BIDS dataset that is written to the dataset_description.json and README.md ' ...
    'files of the BIDS directory. BIDSVersion, DatasetType, and GeneratedBy (catDCM2BIDS and dcm2niix) are ' ...
    'set automatically. ']};


  % main fields
  % =======================================================================
  dicts            = cfg_branch;
  dicts.tag        = 'dicts';
  dicts.name       = 'Dictionary files';
  dicts.val        = {protocoldir, centerdict, studydict, subjectdict, };
  dicts.help       = {'Optional files to replace standard variables by more handy names or codes. '}; 

  opts            = cfg_branch;
  opts.tag        = 'opts';
  opts.name       = 'Options';
  opts.val        = {ProtocolFileName, studies, anonymize, tolerance, preprocessing, ...
                      gzipi, gzipe, output, protocolsubdirs, subIDform, render};
  opts.help       = {['Parameters to control the processing and output of ' ...
                      'the converted data and its export to BIDS. ']}; 


  % batch
  % =======================================================================
  dcm2bids        = cfg_exbranch;
  dcm2bids.tag    = 'dcm2bids';
  dcm2bids.name   = 'DICOM2BIDS';
  dcm2bids.prog   = @cat_io_dcm2bids;
  dcm2bids.vout   = @vout_io_dcm2bids;
  dcm2bids.val    = {datadir, outdir, subdir, dataset, dicts, opts}; 
  dcm2bids.help   = { ...
   ['This batch uses DCM2NIIX to convert DICOM into NIFTI images with JSON sidecars. ' ...
    'It stores the converted data and reorganizes the output in BIDS. ' ...
    'It allows filtering for specific MRI protocols within a directory that outline relevant MR parameters within a JSON file. '] 
    ''
   ['Please check out the CAT subdirectory DCM2BIDS/DZPG3T for an example that includes the definition for ' ...
    'structural, functional, and diffusion scans used by the DZPG (German Center for Mental Health, https://www.dzpg.org/). ']
    ''
    'The batch also applies the SPM anonymization routine and runs a basic image quality control. ' ...
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
