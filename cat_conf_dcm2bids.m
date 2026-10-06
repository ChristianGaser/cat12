function dcm2bids = cat_conf_dcm2bids(expert)
%cat_conf_dcm2bids. Batch definition to convert DICOM to BIDS. 

  if ~exist('expert','var')
    expert = cat_get_defaults('extopts.expertgui'); 
  end

  % define input
  datadir          = cfg_files;
  datadir.tag      = 'data';
  datadir.name     = 'Input Directories';
  datadir.filter   = 'dir';
  datadir.ufilter  = '.*';
  datadir.num      = [1 Inf];
  datadir.help     = {'Select directory with DICOM data directories.'}; 
% what do I do in case of already imported data those raw files are not available any longer?
% >> empty input to output just all internal files?
% >> selection of internal DCM2NIIX directories with extra handling

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
  subdir.val        = {'study'};
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
  subIDform.val     = {2};
  subIDform.hidden  = expert<1;



  % === not implemented yet ===
  ProtocolFileName         = cfg_menu;
  ProtocolFileName.tag     = 'ProtocolFileName';
  ProtocolFileName.name    = 'Use Fitting Protocol Filter File Name';
  ProtocolFileName.labels  = {'Yes','No'};
  ProtocolFileName.values  = {1,0};
  ProtocolFileName.val     = {1};
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
  gzipi.val        = {1};
  gzipi.hidden     = expert<1;
  gzipi.help       = {'GZIP internal NIFTI files in DCM2NIIX import directory to save space.'};

  gzipe            = cfg_menu;
  gzipe.tag        = 'gzipe';
  gzipe.name       = 'GZIP Output Images';
  gzipe.labels     = {'Yes','No'};
  gzipe.values     = {1,0};
  gzipe.val        = {1};
  gzipe.hidden     = expert<0;
  gzipe.help       = {'GZIP output NIFTI files in the BIDS directories. This might help in case of further SPM preprocesing.'};

  % limit output
  output           = cfg_menu;
  output.tag       = 'output';
  output.name      = 'Output level';
  output.labels    = { ...
    'Overview Protocols/Studies (0)', ...
    'Overview + BIDS but only JSON (1)', ...
    'Overview + BIDS of Known Protocols/Studies (2)', ...
    'Overview + BIDS of All Protocols/Studies (3)'};
  output.values    = {0,1,2,3};
  output.val       = {2};
  output.help      = {
   ['Use option 0 to import a subset of the data without running further processing and BIDS output ' ...
    'to get an overview of the used protocols and prepare your own protocol filter sets. '] 
    'Option 1 and 2 allows then prepare the output the only data that fits the defined protocols. ' 
    'Option 3 exports all Data to BIDS. '
    ''
    };

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
  anonymize.val     = {1};
  anonymize.help    = {'Strength of the anonymization of DICOM header and image information. '};
 
  % preprocessing
  preprocessing         = cfg_menu;
  preprocessing.tag     = 'preprocessing';
  preprocessing.name    = 'Run basic SPM preprocessing (expert)';
  preprocessing.labels  = {'No','Yes','Yes (export to derivatives)'}; 
  preprocessing.values  = {0,1,2};
  preprocessing.val     = {1};
  preprocessing.hidden  = expert<1; 
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
  protocolsubdirs.val     = {0};
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
  tolerance.val      = {5};

  % rendering of one slice per scan (BIDS-report/render/...)
  % only the on/off choice is visible for default users, the subfields
  % are expert options
  rsource         = cfg_menu;
  rsource.tag     = 'source';
  rsource.name    = 'Data (expert)';
  rsource.labels  = {'Raw BIDS data','Derivatives (catDCM2BIDS)','Both'};
  rsource.values  = {1,2,3};
  rsource.val     = {3};
  rsource.hidden  = expert<1;
  rsource.help    = {
    'Render the raw BIDS images, the derivatives (e.g. the segmentation), or both. '};

  rsessions         = cfg_menu;
  rsessions.tag     = 'sessions';
  rsessions.name    = 'Sessions (expert)';
  rsessions.labels  = {'Only first session per subject','All sessions'};
  rsessions.values  = {1,0};
  rsessions.val     = {0};
  rsessions.hidden  = expert<1;
  rsessions.help    = {
    'Render only the first session of each subject (one tile per subject) or all sessions. '};

  rtiles          = cfg_menu;
  rtiles.tag      = 'tiles';
  rtiles.name     = 'Tiles per page (expert)';
  rtiles.labels   = {'3x4','4x5','5x7','6x8'};
  rtiles.values   = {[3 4],[4 5],[5 7],[6 8]};
  rtiles.val      = {[4 5]};
  rtiles.hidden   = expert<1;
  rtiles.help     = {
    'Number of tiles (columns x rows) per page in A4 portrait format. '};

  rorient         = cfg_menu;
  rorient.tag     = 'orient';
  rorient.name    = 'Orientation (expert)';
  rorient.labels  = {'Axial','Coronal','Sagittal'};
  rorient.values  = {3,2,1};
  rorient.val     = {3};
  rorient.hidden  = expert<1;
  rorient.help    = {
    'Orientation of the rendered slice in world space. '};

  rslicemode        = cfg_menu;
  rslicemode.tag    = 'slicemode';
  rslicemode.name   = 'Slice position (expert)';
  rslicemode.labels = {'World space','Image center'};
  rslicemode.values = {'world','center'};
  rslicemode.val    = {'world'};
  rslicemode.hidden = expert<1;
  rslicemode.help   = {
   ['The slice is defined in world space (by the image orientation matrix ' ...
    'without registration) and therefore shows positioning differences. ' ...
    'Alternatively, the slice is placed through the center of each image, ' ...
    'which is better to compare anatomy and image quality. ']};

  rslice          = cfg_entry;
  rslice.tag      = 'slice';
  rslice.name     = 'Slice [mm] (expert)';
  rslice.strtype  = 'r';
  rslice.num      = [1 1];
  rslice.val      = {0};
  rslice.hidden   = expert<1;
  rslice.help     = {
    'Position of the slice in mm (in world space or relative to the image center). '};

  norender        = cfg_const;
  norender.tag    = 'norender';
  norender.name   = 'No';
  norender.val    = {0};
  norender.help   = {'No rendering. '};

  dorender        = cfg_branch;
  dorender.tag    = 'render';
  dorender.name   = 'Yes';
  dorender.val    = {rsource, rsessions, rtiles, rorient, rslicemode, rslice};
  dorender.help   = {'Render one slice per scan. '};

  render          = cfg_choice;
  render.tag      = 'render';
  render.name     = 'Render slices';
  render.values   = {norender, dorender};
  render.val      = {dorender};
  render.help     = {
   ['Render one slice of each scan into tiled PNG pages in the BIDS-report directory ' ...
    '(BIDS-report/render/[protocol/]datatype/subtype/) as quick visual check. ' ...
    'Each tile shows the subject, session and overall quality rating. ' ...
    '4D data is represented by its first volume. ']};


  % further possible parameter:
  %%%%%%%%  
  % - avoid/add BIDS field in the JSON files to further specify anonymizing settings?
  %   as string with +SubjectID or -SubjectID 
  % - use the Database as Project name also in the SubjectID
  %   > handling of multiple studies by dictionary
  % - flat to recode subjectIDs 
  %%%%%%%%  
 


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
  dcm2bids.val    = {datadir, outdir, subdir, dicts, opts}; 
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
  try
    if job.opts.gzipe, gz = ' (.gz)'; end
  end

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
