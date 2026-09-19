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
  datadir.help     = {'Select directory with DICOM or NII/JSON data.'}; 
% what do I do in case of already imported data those raw files are not available any longer?
% >> selection of internal DCM2NIIX dir >> need special handling

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
  protocoldir.help     = {[ ...
    'Select directory with JSON files with specified DICOM entries to filter for relevant protocols. ' ...
    'If no directory is specified then all data will be exported to a BIDS directory. ' ...
    'The file name of the JSON filter file together with the sequence number ' ...
    'will be used to specify the BIDS acquisition filed "acq-###-FILENAME"']
    ''
    'Eg. a "myt1w.json" with:'
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
    'define the range of the following quality measures:']
    ''
    '  {'
    '    "SeriesDescription":                   "mprage_sag_0p8mm",'
    '    "SliceThickness":                      0.8,'
    '    "EchoTime":                            0.00222,'
    '    "RepetitionTime":                      2.4,'
    '    "InversionTime":                       1.03,'
    '    "FlipAngle":                           8'
    '  }'
    }; 
  

  % define directory with protocol filter
  studydict          = cfg_files;
  studydict.tag      = 'Pstudydict';
  studydict.name     = 'Study Dictonary File';
  studydict.filter   = 'any';
  studydict.ufilter  = '.*\.json$';
  studydict.val      = {{''}};
  studydict.num      = [0 Inf];
  studydict.hidden   = expert<1;
  studydict.help     = {
    [ ...
      'In case of multiple studies it is posible to include this into the subject code. ' ...
      'Define and link a json file that specify your "DeviceSerialNumber" and defines the wished "InstitutionAbbreviation". ' ...
    ] 
    ''
    'Eg. a "mysites.json" with:'
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

  % define directory with protocol filter
  subjectdict          = cfg_files;
  subjectdict.tag      = 'Psubjectdict';
  subjectdict.name     = 'Subject Dictonary Table';
  subjectdict.filter   = 'any';
  subjectdict.ufilter  = '.*\.csv$';
  subjectdict.val      = {{''}};
  subjectdict.num      = [0 Inf];
  subjectdict.hidden   = expert<1;
  subjectdict.help     = {
    'Integration of phenotypical data by CSV tables with subject-specific "PatientID" or session-specific "StudyNumber".'
    ''
    'Eg. a "subjects.csv" with:'
    'PatientID, MMSE, GROUP'
    '000000001, 30,   0'
    '000000002, 12,   1'
    ''
    }; 

  % select files with study information 
  centerdict          = cfg_files;
  centerdict.tag      = 'Pcenterdict';
  centerdict.name     = 'Center Dictonary File';
  centerdict.filter   = 'any';
  centerdict.ufilter  = '.*\.json$';
  centerdict.val      = {{''}};
  centerdict.num      = [0 Inf];
  centerdict.help     = {
    ['In case of multiple studies it is posible to add this to improve ' ...
    'readability by short names that could be used a variable. ' ...
    'Define and link a json file that specify your "StudySerialNumber" ' ...
    'and defines the wished "StudyAbbreviation". '] 
    ''
    'Eg. a "mystudies.json" with:'
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


  % =======================================================================

  % subjecttable ?
  % csv file to redefine IDs and add other phenotypical data 

  % subjectIDsetup?
  %  sub-ID
  %  sub-CITE-SID
  %  sub-STUDY-SID
  %  sub-CITE-STUDY-SID % omni-center
  %  sub-STUDY-CITE-SID % multi-center 
  % >>>> add both to participants
  subIDform         = cfg_menu;
  subIDform.tag     = 'subIDform';
  subIDform.name    = 'Subject ID form (expert)';
  if expert>1
    subIDform.name    = 'Subject ID form (developer)';
    subIDform.labels  = {
      'sub-PID', ...
      'sub-SITE-PID', ...
      'sub-STUDY-PID', ...
      'sub-SITE-STUDY-PID', ...
      'sub-STUDY-SITE-PID'
      };
    subIDform.values  = {1,2,3,4,5};
    subIDform.help    = {
      'PID=PatientID, SITE=SannerID, STUDY=StudyID .'
      };
  else
    subIDform.labels  = {
      'sub-PID', ...
      'sub-SITE-PID', ...
      };
    subIDform.values  = {1,2};
  end
  subIDform.val     = {2};
  subIDform.hidden  = expert<1;
  subIDform.help    = {
   ['Defintion of the BIDS subject ID with PatientID (PID) only or as combination ' ...
    'with the SITE, as DeviceSerialNumber or recoded by a fitting entry in the ' ...
    '"Center Dictonary File". The SITE entry is used by default as the PID is given ' ...
    'by a center and might be not unique in multicenter studies. ']
    };


  


  % output directory
  outdir            = cfg_files;
  outdir.tag        = 'outdir';
  outdir.name       = 'Output directory';
  outdir.filter     = 'dir';
  outdir.ufilter    = '.*';
  outdir.num        = [1 1];
  %outdir.val{1}     = {''};
  outdir.help       = {[ ...
    'Select a directory where files are written to. ' ...
    'The batch will create a protocoll-conform and a non-conform BIDS structure, and a (temporary) subdirectory with converted NIFTI/JSON data. ' ...
    'It will also create a subdirectory with used MR protocols and final reports. ']};

  subdir            = cfg_entry;
  subdir.tag        = 'subdir';
  subdir.name       = 'Project Name';
  subdir.strtype    = 's';
  subdir.num        = [0 Inf];
  subdir.val        = {'CATBIDS'};
  subdir.help       = {
    'The directory is created within the choosen output directory. If no name is given no subdirecty is created. ' ''};


  % === not implemented yet ===
  ProtocolFileName         = cfg_menu;
  ProtocolFileName.tag     = 'ProtocolFileName';
  ProtocolFileName.name    = 'Use Fitting Protocol Filter File Name';
  ProtocolFileName.labels  = {'Yes','No'};
  ProtocolFileName.values  = {1,0};
  ProtocolFileName.val     = {1};
  ProtocolFileName.hidden  = expert<1;
  ProtocolFileName.help    = {
    'Redefine the name of a protocol by the filename of the fitting protocol filter.'
    };

  % not really needed as the temp-dir is in principle the internal data-base
  deltemp         = cfg_menu;
  deltemp.tag     = 'deltemp';
  deltemp.name    = 'Detele Temporary Files';
  deltemp.labels  = {'Yes','No'};
  deltemp.values  = {1,0};
  deltemp.val     = {0};
  deltemp.hidden  = expert<1;
  deltemp.help    = {''};

  % study selector/filter - NOT WORKING YET
  studies         = cfg_entry;
  studies.tag     = 'studies';
  studies.name    = 'Study Selector/Filter (expert)';
  studies.strtype = 's';
  studies.num     = [0 Inf];
  studies.val     = {''};
  studies.hidden  = expert<1;
  studies.help    = {
    'Specify the export of specify studies by studyID or the defined study appreviations.' ''};

  % zipping
  gzipi            = cfg_menu;
  gzipi.tag        = 'gzipi';
  gzipi.name       = 'GZIP Internal Images (expert)';
  gzipi.labels     = {'Yes','No'};
  gzipi.values     = {1,0};
  gzipi.val        = {1};
  gzipi.hidden     = expert<1;
  gzipi.help       = {'GZIP internal NIFTI files in DCM2NIIX improt directory to save space.'};

  gzipe            = cfg_menu;
  gzipe.tag        = 'gzipe';
  gzipe.name       = 'GZIP Output Images';
  gzipe.labels     = {'Yes','No'};
  gzipe.values     = {1,0};
  gzipe.val        = {1};
  gzipe.hidden     = expert<0;
  gzipe.help       = {'GZIP output NIFTI files in the BIDS directories.'};

  verbose          = cfg_menu;
  verbose.tag      = 'verbose';
  verbose.name     = 'Be verbose (expert)';
  verbose.labels   = {'No','Yes - basic','Yes - extensive'};
  verbose.values   = {0,1,2};
  verbose.val      = {1};
  verbose.hidden   = expert<1;
  verbose.help     = {
   ['Comand line output level as one row per scan (basic) or with addition ' ...
    'processing information (extensive) for debugging. ']};
 
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
  output.val       = {1};
  output.help      = {
   ['Use option 0 to get an overview of (a subset) of your data with DICOM ' ...
    'convertion but without further processing and BIDS output. The reported ' ...
    'protocols can then be used to define ones own protocol filter sets. '] 
    'Option 1 and 2 allows then prepare the BIDS data for the given protocol sets. '
    ''
    };

  % avoid BIDS field
  

  % anonymize .. always required !
  anonymize         = cfg_menu;
  anonymize.tag     = 'anonymize';
  if expert
    anonymize.name    = 'Anonymization level (expert)';
    anonymize.labels  = {'No','Yes - basic','Yes - extensive'};
    anonymize.values  = {0,1,2};
  else
    anonymize.name    = 'Anonymization level';
    anonymize.labels  = {'Basic','Extensive'};
    anonymize.values  = {1,2};
  end
  anonymize.val     = {1};
  anonymize.hidden  = expert<0;
  anonymize.help    = {'Strength of the anonymization of DICOM header and image information. '};

  % optimize .. denoise, bias-corrected, realign?-BB? ...
  % ... we are not doing this yet

  % runpreprocessing ... 
  % ... not yet, but for instance full packages "SPM-CAT-diff-..."

  % QC .. create derivatives directory with BIDS structure and write their the QC values as JSON
  
  dicts            = cfg_branch;
  dicts.tag        = 'dicts';
  dicts.name       = 'Dictonary files';
  dicts.val        = {protocoldir, centerdict, studydict, subjectdict, };
  dicts.help       = {'Voluntary files to replace standard variables by more handy names or codes. '}; 

  opts            = cfg_branch;
  opts.tag        = 'opts';
  opts.name       = 'Options';
  opts.val        = {ProtocolFileName, studies, deltemp, anonymize, gzipi, gzipe, verbose, output, subIDform};
  opts.help       = {'Parameters to control the selection of input files. '}; 

  % main field
  dcm2bids        = cfg_exbranch;
  dcm2bids.tag    = 'dcm2bids';
  dcm2bids.name   = 'DICOM2BIDS';
  dcm2bids.prog   = @cat_io_dcm2bids;
  %dcm2bids.vout   = @vout_io_dcm2nii;
  dcm2bids.val    = {datadir, outdir, subdir, dicts, opts}; 
  dcm2bids.help   = {[ ...
   ['Convert DICOM via dcm2niix into NIFTI with JSON sidecar. ' ...
    'Reorganize the converted data into BIDS and filter for given protocols specified by JSON files. ' ...
    'Select your consensus protocols (eg. DZPG) or extract available protocols from your DICOMs and ' ...
    'define your own default protocol set(s) by using the output option 0 that stopts after DICOM import. '] ...
    'For data quality control basic image quality measures are extracted that can be used to define protocol-' ...
    'specific ranges to create ratings that can simply be averaged.' 
    ]}; 
end
function cdep = vout_io_dcm2nii
% connect to BIDS2PP batch
  cdep = cfg_dep;

  cdep(end).sname      = 'conform BIDS';
  cdep(end).src_output = substruct('.','avg','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  cdep(end+1) = cfg_dep;
  cdep(end).sname      = 'unconform BIDS';
  cdep(end).src_output = substruct('.','avg','()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

end