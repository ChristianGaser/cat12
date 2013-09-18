function vbm_db(action,FT,DBhdir,DBname,opt,varargin)
% ______________________________________________________________________
% WARNING: THIS IS A PRIVATE FUNCTION THAT IS STILL IN DEVELOPMENT! 
%          IT ONLY WORKS ON OUR SERVERS, AND WILL MAYBE BECOME PART OF
%          THE VBM TOOLBOX IN THE FUTURE.
% ______________________________________________________________________
%
% vbm_db(action,FT,DBhdir,DBname,opt,varargin)
% 
%   action = {'import','update'} 
%   DBhdir = vbm database home directory
%   DBname = vbm databse name (default = 'vbmDB')
%   FT     = input files 
%   opt    = inport options
%
% vbm_db('import','/Volumes/MyBook/NISALS/home/graz','/Volumes/MyBook/MRData','vbm_DB')
%
% Unkown fields are replaced by a similar number of the letter X.
%
% hier muss rein, wie ich von einer fremnden struktur die nötigen infos extrahieren kann
%
% CSV:            SORTIEREN???
% ALLGEMEIN:      VERSIONIEREN???
% 
% SUBVERZEICHNISSE FÜR 
%  RAW
%  PRE
%  VIEW
%
% Subpath in the resultdirectory an filename:
%   Project, GroupID, Center,
%   ScannerManufacturer, MagneticFieldStrength, Sequence, Contrast, 
%   SubjectID, YOB, DOB, Sex, Age, Handness,
%   Scandate, Scantime
%
% Example = {{'SubjectID'} {'Contrast'} {'ScanDate' '-' 'ScanTime'} for
% MaxMueller_T1_20120809-145322
%
% TODO: - write help
%       - scantime problem (n different fields for the time and non 
%         gives the real time of the image) 
%       - IDs
%       - anonymize center function
%       - MR sequence table
%       * subfields 
%         - further dicom data
%         - einbinden von [DOB,DOS,age,txt]=CompleteDates(DOB,DOS,age)
%       * deface => FreeSurfer | after preprocessing
%       - csv-export > csv-import
% __________________________________________________________________________________________________
% $Id$
 

 
  % action... can set this latter based on DBhdir and the 
  % -> if these files are within DBdir than update else import
  if ~exist('action','var'), error('MATLAB:vbm_db','Need action!\n'); end
  if ~exist('opt','var'), opt = struct(); end
  
  
  % vbm databse name/directory and homedirectory
  switch action
    case 'import' % import to RAW directory
      if ~exist('DBhdir','var') || ~exist(DBhdir,'dir') 
        DBhdir = spm_select(1,'dir','select database home directory'); 
      end
      if ~exist('DBname','var')
        DBname = 'vbmDB';
      end
    case 'update' % updates of the RAW directory
      if ~exist('DBhdir','var') || ~exist(DBhdir,'dir') || ...
         ~exist('DBname','var') || ~exist(fullfile(DBhdir,DBname),'dir') 
        DBdir = spm_select(1,'dir','select database directory'); 
        [DBhdir,DBname] = fileparts(DBdir);  
      end
    %case 'RAW2PRE'   % creation of preprocessing directories
    %case 'PRE2VIEW'  % creation of views on the preprocessing directories
  end
  DBdir = fullfile(DBhdir,DBname);
  
  
  
  % get files (or dirs?)
  if ~exist('FT','var') || isempty(FT),
    switch action
      case 'import' % import2RAW
        FT = spm_select([1 Inf],'dir','select import directories');
      case 'update'
        FT = spm_select([1 Inf],'image','select image files for update','',DBdir);
      case 'delete'
    end
  end
  FT = cellstr(FT); for fti=1:numel(FT), if FT{fti}(end-1)==',', FT{fti}(1:end-2); end; end
  
    
  
  % DEFAULT OPTIONS:
  mdef=vbm_db_defaults; 
  if strcmp(action,'import'), 
  % There are special import options, to definie special input ways like XML (ADNI) or a special csv
  % file (OASIS). This is not necessary, if we update something in the vbm database.
    if ~isfield(opt,'importdb');
      if exist('varargin','var') && numel(varargin)>0 && ~isempty(varargin{1})
        opt.importdb = varargin{1};
      else
        opt.importdb = '';
      end
    end
    switch opt.importdb
      case 'ADNI',      def=vbm_db_defaults_ADNI; 
      case 'INDI',      def=vbm_db_defaults_INDI;
      case 'NISALS',    def=vbm_db_defaults_NISALS;
      case 'DICOM',     def=vbm_db_defaults_DICOM; % sicher? ==own??
      case 'NIFTI',     def=vbm_db_defautls_NIFTI; % sicher? ==own??
      case 'own?';      if exist('varargin','var') && numel(varargin)>1 && ~isempty(varargin{2}), def=eval(varargin{2}); end
      otherwise,        def=vbm_db_defaults; % sicher? ==own??
    end
  end
  def = checkinopt(mdef,def);

  
  % das hier fliegt wohl raus... 
  % ____________________________________________________________________
  % meta data options
  def.meta.prefix      = 'data_';
  def.meta.ScanIDtype  = 'uint32';
  def.meta.usecsv      = 1; % [0 = no, 1 = yes]  % file for updates
  def.meta.usemysql    = 1; % [0 = no, 1 = yes]  % file for updates
  def.meta.importorder = {'csv','fname','dcm'}; 

  % Anwendungsverzeichnisse von MRIcron, Freesurfer, FSL und MySQL
  def.MRIcron     = '/Applications/MRIcron/';
  def.FreeSurfer  = '/Applications/Freesurfer5.1';
  def.FSL         = '/Applications/FSL4.1';
  def.MySQL       = '/usr/local/mysql-5.5.28-osx10.6-x86_64/';
  
  % mysql options
  def.mysql.path    = def.MySQL;
  def.mysql.user    = 'root';
  def.mysql.pwd     = 'sql4NISALS';
  def.mysql.db      = 'vbmDB';
  def.mysql.tname   = 'vbmDB';
  def.mysql.tid     = 'ScanID';
  
  % anonymize
  def.anomize     = 1;                         % Anonymisieren der Bilder
  


  
  % Subpath in the resultdirectory an filename:
  %   Project, GroupID, Center,
  %   ScannerManufacturer, MagneticFieldStrength, Sequence, Contrast, 
  %   SubjectID, YOB, DOB, Sex, Age, Handness,
  %   Scandate, Scantime
  %
  %
  % Example:
  %   input: {{'SubjectID'} {'Contrast'} {'ScanDate' '-' 'ScanTime'} 
  %   ouput: MaxMueller_T1_20120809-145322
  %
  % Although, filenames will make it easier to identify images, they
  % have may be changed, if the meta informations were corrected...
  def.pathname = {
   % {'Center'}
    {'Contrast'}
  };
  def.filename = {
    {'Project'}
   % {'GroupID'} 
    {'Center'} 
    {'ScannerManufacturer' 'MagneticFieldStrength'}
    {'SubjectID'}
   % {'YOB' 'SubjectID'}
   % {'Sex' 'Age' 'Handness'}
    {'Contrast' '-' 'Contrast2'} % '-' 'Sequence'}
   % {'ScanDate' '-' 'ScanTime'}
    {'ScanID'};  % this field is necessary for unique idententification!
  }; 
  def.Contrast = {'T1w','T2w','PDw'}; % contrast for import {'T1w','T2w','PDw','DTI',...}

  if ~exist('opt','var'), opt = struct(); end
  opt = checkinopt(opt,def);

  opt.initPATH = sprintf( [...
    'export MRIcron=%s;' ...
    'export FSLDIR=%s; sh ${FSLDIR}/etc/fslconf/fsl.sh;' ...
    'export FREESURFER_HOME=%s; ${FREESURFER_HOME}/SetUpFreeSurfer.sh;' ...
    'MNI_DIR=${FREESURFER_HOME}/mni; ' ...
    'PATH=$PATH:$FSLDIR/bin:$MNI_DIR/bin:$MRIcron; ' ...
    'source $FREESURFER_HOME/FreeSurferEnv.sh; ' ...
    'export FSLOUTPUTTYPE=NIFTI NIFTIexport=NIFTI; ' ...
    ],opt.MRIcron,opt.FSL,opt.FreeSurfer);
  
  
  
  % variables given by input variables (if these fields are specified in 
  % opt, than we overwrite them)
  opt.DBhdir  = DBhdir;
  opt.DBname  = DBname;
  opt.DBdir   = fullfile(opt.DBhdir,opt.DBname);
  opt.DBcsv   = fullfile(opt.DBdir,[opt.DBname '.csv']); % main csv-file ... das hier wird wohl meine export variable werden
  opt.DBRAW   = fullfile(opt.DBdir,'RAW');
  opt.DBTMP   = fullfile(tempdir,'vbmDBimport');
  opt.vbmDBfile  = fullfile(opt.DBdir,'vbmDBfile.mat');  
  
  
  %opt.DBPRE   = fullfile(opt.DBdir,'PRE');
  %opt.DBVIEW  = fullfile(opt.DBdir,'VIEW');
  
  
  
  % create main directories
  if ~exist(opt.DBdir,'dir'), mkdir(opt.DBdir); end
  if ~exist(opt.DBRAW,'dir'), mkdir(opt.DBRAW); end
  %if ~exist(opt.DBPRE,'dir'), mkdir(opt.DBPRE); end
  %if ~exist(opt.DBVIEW,'dir'), mkdir(opt.DBVIEW); end
  
  
  switch action
    case 'import'
      dicom_import(FT,opt)
      nifti_import(FT,opt)    
    case 'update'
      % hier wählich ich aus dem db-verzeichnis und muss nach source 
      % (mat,cvs,mysql) daten auslesen, aktualisieren und files und dirs 
      % umsortieren
        
      % set update source
      if exist('varargin','var') && numel(varargin)>0 && ~isempty(varargin{1})
        opt.source = varargin{1};
      else 
        opt.source = 'csv';
      end

      % get/check input MR-files
      vbm_db_update(FT,opt);
    case 'delete'
    otherwise
      error('MATLAB:vbm_db','Unkown action ''%s''',action);
  end
end

function vbm_db_update(FT,opt)
  % set mat-filenames
  Fmeta = cell(size(FT));
  if opt.meta.usemat  
    for f=1:numel(FT), [pp,ff]=fileparts(FT{f}); Fmeta{f} = fullfile(opt.DBdir,[opt.meta.prefix 'mat'],[opt.meta.prefix ff '.mat']); end
  else
    for f=1:numel(FT), [pp,ff]=fileparts(FT{f}); Fmeta{f} = fullfile(pp,[opt.meta.prefix ff '.mat']); end
  end
      
  % get new data
  meta=struct();
  switch opt.source
    case 'mat'
      % load mat file
      for f=1:numel(FT)
        tmp = load(Fmeta{f},'meta');
        fn  = fieldnames(tmp.meta);
        for fni=1:numel(fn), meta(f).(fn{fni}) = tmp.meta.(fn{fni}); end 
      end
    case 'csv'
      % load csv file > convert to struct
      % hier wäre wohl das hauptfile gut zur editierung geeignet...
      meta = vbm_io_struct(opt.DBcsv);
    case 'mysql'
      % get entries
      
      meta = vbm_io_mysql({FT(:).fname}',opt.mysql);
  end
  
  % rename file and update meta
  for f=1:numel(FT)
    % path and filename
    fpath = metaname(opt,meta(f),'pathname',opt.DBdir);
    fname = metaname(opt,meta(f),'filename'); meta(f).fname = fname;

    % rename nifti
    movefile(FT{f},[fullfile(fpath,fname),'.nii']); 
  
    % update mat-file (delete old, write new)
    if strcmp(opt.source,{'csv','mysql'})
      save(Fmeta{f},'meta');
      % delete old meta-file and write a new one
      delete(Fmeta{f});
      if opt.separateMetaFile
        if ~exist(fullfile(opt.DBdir,'meta'),'dir'); mkdir(fullfile(nifti_dir,'meta')); end
        save(fullfile(opt.DBdir,'meta',['meta_' fname '.mat']),'metadicom','meta');
      else
        save(fullfile(fpath,['meta_' fname '.mat']),'metadicom','meta');
      end
    end
    
    % update csv-file
    if strcmp(opt.source,{'mat','mysql'})
      % hier musst du die alten einträge durch neue ersetzen!!!
        vbm_io_csv(opt.csv,meta);
    end
    
    % update mysql
    if strcmp(opt.source,{'mat','csv'})
      % hier muss du nur ein mysql-update ausführen 
      vbm_io_mysql(opt,meta); 
    end      
  end
end
function def=vbm_db_defaults
  def.dcm2nii.reorienate    = 1; % [0=nothing| 1*reorientat images | 2=reorientate & crop images]
 
  def.maxflength.SubjectID  = 10;
  def.maxflength.Sequence   = 20;
  def.maxflength.SubjectID  = 10;
  def.maxflength.Sequence   = 20;

  def.use_vbmDBfile = 1; 
  
  def.Contrasts = {
    'T1w' {'T1','GR','GR\IR','SE\IR','MPRAGE','MPRage','mpr','*tfl3d1'}
    'T2w' {'T2','EP','ep_b0'}
    'PDw' {'PD'}
    'DTI' {'DTI','DWI'}
    'LOC' {'scout','localizer','Survey'} 
  };
  def.badfields = {'PatientName' 'PatientsName'}; 
  def.metafields = {
  % default_chararcter for equal_field lenght 'X' '=' '.'
  % 
  % MF().Category     = '';                     % Category of S.(MF.Category).(MF.Name)
  % MF().Name         = '';                     % Name of the field
  % MF().Datatype     = '';                     % Datatype of the field (for the database)
  % MF().Fieldlength  = 8;                      % 
  % MF().InputNames   = {''};                   % Name of the dicom fiels
  % MF().FilePos      = Position???             % Position of the filename
  % MF().InputPref    = {'XML','CSV','File','Dicom'}; % Input preference
  % MF().Content().New   = '';
  % MF().Content().Txt   = {''};
  % MF().Content().Value = value;
  % MF().Content().Range = [low high]; 
  %
  % Category   MetaName   DicomName                   MetaContent < DicomContent     
  % Subject    Contrast   {Protocoll,Series,...}      {'T1w'  {'T1','GR','GR\IR','SE\IR','MPRAGE','MPRage','mpraxial15','*tfl3d1'}}
  %
  % Only Dicom Data described here will be storred, because other field may contain private data!
  % The subgroup field allows grouping of variables and can be used to create further tables of a 
  % final relational database.
  %
  % 'OLDNAME'                         'NEWNAME'                 'SUBGROUP'
  %
  %File Information:
    'Filename'                        'original_fname'          '' 
    'StudyDate'                       ''                        'times'
    'SeriesDate'                      ''                        'times'
    'AcquisitionDate'                 ''                        'times'
    'StudyTime'                       ''                        'times'
    'SeriesTime'                      ''                        'times'
    'AcquisitionTime'                 ''                        'times'
    'InstanceCreationDate'            ''                        'times'
    'InstanceCreationTime'            ''                        'times'
    'PerformedProcedureStepStartDate' ''                        'times'
    'PerformedProcedureStepEndDate'   ''                        'times'
    'PerformedProcedureStepStartTime' ''                        'times'
    'PerformedProcedureStepEndTime'   ''                        'times'
  % Scanner Information:
    'Manufacturer'                    ''                        'scanner'
    'ManufacturerModelName'           ''                        'scanner'
    'MagneticFieldStrength'           ''                        'scanner'
    'TimeOfLastCalibration'           ''                        'scanner'
    'StationName'                     ''                        'scanner'
    'SoftwareVersion'                 ''                        'scanner'
  % Site Information
    'InstitutionName'                 ''                        'site'
    'InstitutionAdress'               ''                        'site'
    'ReferringPhysiciansName'         ''                        'site'
    'StationName'                     ''                        'site'
  % Subject Information     
    'PatientBirthDate'                'DOB'                     'sub'
    'PatientAge'                      'Age'                     'sub'
    'PatientSex'                      'Sex'                     'sub'
    'PatientWeight'                   'weight'                  'sub'
    'StudyID'                         ''                        'sub'
  % Sequence Information (Scan Protocoll); 
    'ProtocolName'                    ''                        'image'
    'SequenceName'                    ''                        'image'
    'SequenceVariant'                 ''                        'image'
    'ScanningSequence'                ''                        'image'
    'AngioFlag'                       ''                        'image'
    'ScanOptions'                     ''                        'image'
    'MRAcquisitionType'               ''                        'image'
    'RepetitionTime'                  ''                        'image'
    'EchoTime'                        ''                        'image'
    'InversionTime'                   ''                        'image'
    'FlipAngle'                       ''                        'image'
    'SAR'                             ''                        'image'
  % image properties
    'SeriesDescription'               ''                        'image'
    'Width'                           ''                        'image'
    'High'                            ''                        'image'
    'SliceThickness'                  ''                        'image'
    'Rows'                            ''                        'image'
    'Columns'                         ''                        'image'
    'BitDepth'                        ''                        'image'
    'PatientPosition'                 ''                        'image'
    'AcquisitionMatrix'               ''                        'image'
    'PercentSampling'                 ''                        'image'
    'PercentPhaseFieldOfView'         ''                        'image'
    'ImagePositionPatient'            ''                        'image'
    'ImageOrientationPatient'         ''                        'image'
    'BodyPartExamined'                ''                        'image'
  };
  def.manufacturer = {
    'GE'        'GE'
    'Siemens'   'SI'
    'Philips'   'SP'
    'Hitachi'   'HI'
    'Toshiba'   'TO'
  };

  % http://www.imaios.com/en/e-Courses/e-MRI/MRI-Sequences/Sequences-acronyms
  % Philips   Siemens      GE           Hitachi     Toshiba
  %
  def.sequences = {
    {'SE'           'SE'         'SE'         'SE'        'SE'        } 'SE'        % Spin Echo (SE)
    {'Multi SE'     'Multi écho MS'         'SE'         'SE'        'SE'        } 'MSE'       % Multi Echo SE
    {'Turbo SE'     'SE'         'SE'         'SE'        'SE'        } 'FSE'       % Fast SE
    {'SSH-TSE'      'SE'         'SE'         'SE'        'SE'        } 'UFSE'      % Ultra fast SE
    {'UFSE'         'SE'         'SE'         'SE'        'SE'        } 'UFSE'      % Ultra fast SE
    {'IR'           'SE'         'SE'         'SE'        'SE'        } 'IR'        % IR
    {'IR TSE'       'SE'         'SE'         'SE'        'SE'        } 'IR'        % IR
    {'STIR'         'SE'         'SE'         'SE'        'SE'        } 'STIR'      % STIR
    {'STIR TSE'     'SE'         'SE'         'SE'        'SE'        } 'STIR'      % STIR     
    {'FLAIR'        'SE'         'SE'         'SE'        'SE'        } 'FLAIR'     % FLAIR
    {'FLAIR TSE'    'SE'         'SE'         'SE'        'SE'        } 'FLAIR'     % FLAIR
    {'FFE'          'SE'         'SE'         'SE'        'SE'        } 'GE'        % Gradient Echo (GE)
    {'T1-FFE'       'SE'         'SE'         'SE'        'SE'        } 'UFGE'      % Ultra fast GE
    {'T2-FFE'       'SE'         'SE'         'SE'        'SE'        } 'UFGE'      % Ultra fast GE
    {'THRIVE'       'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'IR-TFE'       'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'FFE'          'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'T2-FFE'       'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'T2'           'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'Balanced FFE' 'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'FFE'          'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'FFE'          'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
    {'FFE'          'SE'         'SE'         'SE'        'SE'        } 'UFGE'     
  };
end
function def=vbm_db_defaults_NISALS
  def.project               = 'NISALS';
  
  % IDs for groups, center, Contrasts, bad-field, ...
  def.groups = {
    'controls'      'HC';
    'patients-als'  'ALS';
  };
  def.center = { % directory name, centerID, anonymize center
    'Unkown'        'UKN'   0; % special group - if we don't know it 
    'Other'         'OTH'   0; % special group - if we do not like to know it :D
    'Bad'           'BAD'   0; % special group - bad datasets
    'Unkown'        'XXX'   0; % special group - internal if we don't know it
  };
 

end

function dicom_import(input_dir,opt) 

  % 1.  Find all dicoms files and sort them by there home directories.
  % 1.1 Find all files within the subdirectories and remove all non-dicom files from this list.
  %     input_dir = char(input_dir); 
  files = findfiles(input_dir,'*')';
  for fi=numel(files):-1:1
    [p,f,e]=fileparts(files{fi});
    fdata=dir(files{fi}); % remove to small files and check access
    if isempty(f) || f(1)=='.' || isempty(fdata) || fdata.bytes < 256 || ...
      sum(strcmpi(e,{'.ima','.dcm','.dic',''}))==0  || ...
      sum(strcmpi(e,{'.nii','.img','.hdr','.mnc','jpg','txt','png','xls','doc','docx'}))
      files(fi)='';
    end
  end
  if isempty(files), fprintf(1,'No DICOMs in %s!\n',char(input_dir)); return; end

  
  % 1.2 GroupID the files based on there home directories.
  dfiles{1}{1} = files{1}; pdir=fileparts(files{1}); di=1;
  for fi=2:numel(files);
    p=fileparts(files{fi});
    if ~strcmp(p,pdir), di=di+1; pdir=p; dfiles{di}{1}=files{fi}; else dfiles{di}{end+1}=files{fi}; end %#ok<AGROW>
  end
  clear files;
  

  % 2.  Run convertation for each subject %tempname
  tmpdir = opt.DBTMP; if ~exist(tmpdir,'dir'); mkdir(tmpdir); end; olddir=pwd; cd(tmpdir); 
  for j=numel(dfiles):-1:1
    % 2.1 read all dicom header
    warning off; dicoms = spm_dicom_headers(char(dfiles{j})); warning on; %#ok<WNON,WNOFF>
    if isempty(dicoms), continue; end

    
    % 2.2 call mricron dcm convert
    % -a anonymize, -b begin 4d volume, -d, -d date, -e events (s###a###), -f source-filename, -g gzip
    % -i ID, -l end 4d volume, -n nifti, -o output dir, -p protocol, -s spm2/5 (analyse only), v ?
    % -r reorientate, -x reorientate and crop
    reorienate = containers.Map({0,1,2},{'-r n -x n','-r y -x n','-r y -x y'});
    %fprintf('\n\n\n>>\n');
    [SR1,SR2] = system(sprintf([sprintf('PATH=$(echo %s:$PATH);',opt.MRIcron)  ...
                                'dcm2nii -a n -c n -d y -e y -f Y -g n -i n ' ...
                                '-n y -p y -v y %s -o %s %s;'], ...
                                reorienate(opt.dcm2nii.reorienate),tmpdir, ...
                                ([char(dfiles{j}') repmat(' ',size(dfiles{j},2),1)])'));  %#ok<NASGU>
    %fprintf('\n<<\n\n\n');
    if opt.dcm2nii.reorienate>0
      nii_files=findfiles(tmpdir,'o*.nii')';
      for fi=1:numel(nii_files), 
        [pp,ff,ee]=fileparts(nii_files{fi}); 
        movefile(nii_files{fi},fullfile(pp,[ff(2:end),ee]));
      end
    else
      nii_files=findfiles(tmpdir,'co*.nii')';  
      for fi=1:numel(nii_files), delele(nii_files{fi}); end
    end
    if opt.dcm2nii.reorienate>1
      nii_files=findfiles(tmpdir,'co*.nii')';  
      for fi=1:numel(nii_files), 
        [pp,ff,ee]=fileparts(nii_files{fi}); 
        movefile(nii_files{fi},fullfile(pp,[ff(3:end),ee]));
      end
    else
      nii_files=findfiles(tmpdir,'co*.nii')';  
      for fi=1:numel(nii_files), delete(nii_files{fi}); end
    end
    nii_files = findfiles(tmpdir,'*.nii')';

    if exist(opt.vbmDBfile,'file') && opt.use_vbmDBfile
      load(opt.vbmDBfile);
    else
      DB = {};
    end

    % 2.3 process each scan separately 
    %     extract diocm-meta-data, analyse filestructure and insert them into the database
    for d=1:numel(nii_files)
      if ~isempty(nii_files{d})
        meta = get_dicom_meta(dicoms,nii_files{d},opt); %mx{j,d}=meta;
   
        % deidentification
        if opt.anomize && isfield(meta,'Contrast') && strcmp(meta.Contrast,'T1w');
          try
            [SR1,SR2] = system(sprintf('%s  mri_deface %s %s %s %s;', ... 
              opt.initPATH, nii_files{d}, ...
              '$FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca', ...
              '$FREESURFER_HOME/average/face.gca', ...
              nii_files{d}));
          end
          %RS(SR1,SR2); % error message
        end
        
        DB=insert(nii_files{d},meta,opt,DB);
        clear meta;
      end
    end
    
    % save DB variable
    save(opt.vbmDBfile,'DB');
    
    
    % delete old niftis
    delete *;
  end
  rmdir(tmpdir,'s');
  cd(olddir);
end  
function nifti_import(input_dir,opt)
  % find analyse files and convert to nifti
  
  imgfiles = findfiles(input_dir,'*.img')';
  if ~isempty(imgfiles)
    tmpdir = fullfile(tempdir,'vbmDBimport'); 
    if ~exist(tmpdir,'dir'); mkdir(tmpdir); end;
    olddir=pwd; cd(tmpdir); 
  end
   
  for i=1:min(numel(imgfiles),inf) % ##################################
    [pp,ff]=fileparts(imgfiles{i});
    if ~exist(fullfile(pp,[ff '.nii']),'file') 
      if exist(fullfile(pp,[ff '.hdr']),'file')
        Vhdr = spm_vol(imgfiles{i});
        V    = spm_read_vols(Vhdr);
        Vhdr.fname = fullfile(tmpdir,[ff '.nii']);
        spm_write_vol(Vhdr,V);
        niifiles{i} = Vhdr.fname;
      else
        imgfiles{i} = ''; 
      end 
    end
  end
  
  if exist(opt.vbmDBfile,'file') && opt.use_vbmDBfile
    load(opt.vbmDBfile);
  else
    DB = {};
  end
  
  % find nifties, extract nifti-meta-data (filestructure) and insert them into the database
  if ~isempty(imgfiles) && ~exist(fullfile(pp,[ff '.nii']),'file') 
    if isempty(niifiles), fprintf(1,'No NIFTIs in %s!\n',char(input_dir)); return; end
    for d=1:numel(niifiles)
      %try
      if ~isempty(imgfiles{d})
        
        meta = get_nifti_meta(imgfiles{d},opt); % need the original data
        
        % deidentification
        if opt.anomize && strcmp(meta.Contrast,'T1w');
          [pp,ff,ee]=fileparts(niifiles{d});
          aniifiles{d}=fullfile(pp,['a' ff ee]);
          try
            [SR1,SR2] = system(sprintf('%s  mri_deface %s %s %s %s;', ... 
              opt.initPATH, niifiles{d}, ...
              '$FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca', ...
              '$FREESURFER_HOME/average/face.gca', ...
              aniifiles{d}));
          end
          %RS(SR1,SR2); % error message
        else 
          aniifiles{d} = niifiles{d};
        end

        DB = insert(aniifiles{d},meta,opt,DB);
      end  
      
    % save DB variable
    save(opt.vbmDBfile,'DB');
   
     % catch
     %   fprintf('Error: %s\n',imgfiles{d});
     % end
    end
    rmdir(tmpdir,'s');
  else
    niifiles = findfiles(input_dir,'*.nii')';
    if isempty(niifiles), fprintf(1,'No NIFTIs in %s!\n',char(input_dir)); return; end
    for d=1:numel(niifiles)
      meta = get_nifti_meta(niifiles{d},opt); % need the original data
      
      % deidentification
      if opt.anomize && strcmp(meta.Contrast,'T1w');
        [pp,ff,ee]=fileparts(niifiles{d});
        aniifiles{d}=fullfile(pp,['a' ff ee]);
        try
          [SR1,SR2] = system(sprintf('%s  mri_deface %s %s %s %s;', ... 
            opt.initPATH, niifiles{d}, ...
            '$FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca', ...
            '$FREESURFER_HOME/average/face.gca', ...
            aniifiles{d}));
        end
        %RS(SR1,SR2); % error message
      else 
        aniifiles{d} = niifiles{d};
      end
      
      DB = insert(aniifiles{d},meta,opt,DB);
    end
  end
end

function meta = create_meta()
% fixed fieldlength ...
% fieldname, oldfieldnames, dicomfieldname, fieldtype, fieldlenght / default / filenamestruct
  metastruct = { ...
    ... % filename and database varialbes
    'oldfname'     , 'char'    , ...
    'fpath'        , 'char'    , ...
    'fname'        , 'char'    , ...
    'ScanID'       , 'uint32'  , ...
    ... % project, center and investigator 
    'Project'      , 'char'    , ...
    'Center'       , 'char'    , ...
    ... %'Investigator'        , 'char'    , ...
    ... % Group varaibles
    ... % name, disease, 
    ... % subject varialbes
    'Subject'      , 'char'    , ...
    'GroupID'      , 'char'    , ...
    'Sex'          , 'char'    , ...
    ... % examination variables
    'Date'         , 'char'    , ...
    'Time'         , 'char'    , ...
    'DateTime'     , 'char'    , ...
    'Age'          , 'double'  , ...
    'DOB'          , 'double'  , ...
    'YOB'          , 'double'  , ...
    'Tests'        , 'char'    , ...
    ... % scanner variables
    'Manufacturer' , 'char'    , ...
    'Model'        , 'char'    , ...
    'fieldstrength', 'double'  , ...
    ... % protocoll variables
  };
  metastructinit = [metastruct(1,:)'; metastruct(2,:)'];
  meta = eval(sprintf('struct(%s%s)',metastructinit(1),sprintf(',%s',metastructinit(2:end))));
  
  %{
   'Filename'                        'original_fname'          '' 
    'StudyDate'                       ''                        'times'
    'SeriesDate'                      ''                        'times'
    'AcquisitionDate'                 ''                        'times'
    'StudyTime'                       ''                        'times'
    'SeriesTime'                      ''                        'times'
    'AcquisitionTime'                 ''                        'times'
    'InstanceCreationDate'            ''                        'times'
    'InstanceCreationTime'            ''                        'times'
    'PerformedProcedureStepStartDate' ''                        'times'
    'PerformedProcedureStepEndDate'   ''                        'times'
    'PerformedProcedureStepStartTime' ''                        'times'
    'PerformedProcedureStepEndTime'   ''                        'times'
  % Scanner Information:
    'Manufacturer'                    ''                        'scanner'
    'ManufacturerModelName'           ''                        'scanner'
    'MagneticFieldStrength'           ''                        'scanner'
    'TimeOfLastCalibration'           ''                        'scanner'
    'StationName'                     ''                        'scanner'
    'SoftwareVersion'                 ''                        'scanner'
  % Site Information
    'InstitutionName'                 ''                        'site'
    'InstitutionAdress'               ''                        'site'
    'ReferringPhysiciansName'         ''                        'site'
    'StationName'                     ''                        'site'
  % Subject Information     
    'PatientBirthDate'                'DOB'                     'sub'
    'PatientAge'                      'Age'                     'sub'
    'PatientSex'                      'Sex'                     'sub'
    'PatientWeight'                   'weight'                  'sub'
    'StudyID'                         ''                        'sub'
  % Sequence Information (Scan Protocoll); 
    'ProtocolName'                    ''                        'image'
    'SequenceName'                    ''                        'image'
    'SequenceVariant'                 ''                        'image'
    'ScanningSequence'                ''                        'image'
    'AngioFlag'                       ''                        'image'
    'ScanOptions'                     ''                        'image'
    'MRAcquisitionType'               ''                        'image'
    'RepetitionTime'                  ''                        'image'
    'EchoTime'                        ''                        'image'
    'InversionTime'                   ''                        'image'
    'FlipAngle'                       ''                        'image'
    'SAR'                             ''                        'image'
  % image properties
    'SeriesDescription'               ''                        'image'
    'Width'                           ''                        'image'
    'High'                            ''                        'image'
    'SliceThickness'                  ''                        'image'
    'Rows'                            ''                        'image'
    'Columns'                         ''                        'image'
    'BitDepth'                        ''                        'image'
    'PatientPosition'                 ''                        'image'
    'AcquisitionMatrix'               ''                        'image'
    'PercentSampling'                 ''                        'image'
    'PercentPhaseFieldOfView'         ''                        'image'
    'ImagePositionPatient'            ''                        'image'
    'ImageOrientationPatient'         ''                        'image'
    'BodyPartExamined'                ''                        'image'
  %}
end
function meta = get_dicom_meta(dicoms,nii_file,opt)

  % 2.3.1 find the correct header based on the filename
  %       if you can't find one, print this and continue with the next image
  %nii_file
  %dicoms{1}
  %fprintf('\n%s\n\n',nii_file);
   % [pp1,ff1]=fileparts(nii_file);
  
  for dc=1:numel(dicoms)
    [pp,ff]=fileparts(dicoms{dc}.Filename);
    if ~isempty(strfind(nii_file,ff)) || ~isempty(strfind(nii_file,strrep(ff,'.',''))) 
      dicom_meta=dicoms{dc}; break;
    end
  end
  if ~exist('dicom_meta','var') || isempty(dicom_meta)
    fprintf('No DICOM HDR: %s\n',nii_file); meta=struct(); return 
  end
  
  % Dicom Information
  Vnii=spm_vol(nii_file);
  
  % RESOLUTION
  % SpacingBetweenSlices, SliceThickness
  meta.DimX = Vnii(1).dim(1);
  meta.DimY = Vnii(1).dim(2);
  meta.DimZ = Vnii(1).dim(3);
  if any(Vnii(1).dim(1:3) == 1), meta.Contrast = 'LOC'; return;  end
  
  
  % to analays path and filename
  fnameparts = eval(['{'' ' lower(strrep(dicom_meta.Filename,'/',''',''')) '''}'])'; fnameparts(1)=[]; 
 
  % dicom and nifti filename
  if isfield(dicom_meta,'Filename'), meta.dicomname = dicom_meta.Filename; end
  [pp,ff] = fileparts(nii_file); meta.dcm2niiname = ff; 
  
  
  % 2.3.2 set meta fields
  % PROJECT
  meta.Project = opt.project;


  % GROUP
  groupn=0; for i=1:size(opt.groups,1), if sum(strcmpi(opt.groups{i,1},fnameparts))>0, groupn=i; end; end 
  if groupn==0, meta.GroupID = 'XX'; else meta.GroupID = opt.groups{groupn,2};  end


  % CENTER
  centern=0; for i=1:size(opt.center,1), if strfind(lower(strrep(dicom_meta.Filename,'/',' ')),lower(opt.center{i,1}))>0, centern=i; end; end 
  if centern==0, meta.center = 'XXX'; else  meta.Center = opt.center{centern,2}; end


  % DATE
  if isfield(dicom_meta,'PerformedProcedureStepStartTime') && isfield(dicom_meta,'AcquisitionDate')
    meta.ScanDayTime = sprintf('%s-%s%s%s',...
                       datestr(dicom_meta.AcquisitionDate,'yyyymmdd'), ...
                       sprintf('%02d', floor(dicom_meta.PerformedProcedureStepStartTime/3600)), ...
                       sprintf('%02d', floor(rem(dicom_meta.PerformedProcedureStepStartTime/60,60))), ...
                       sprintf('%02d', floor(rem(dicom_meta.PerformedProcedureStepStartTime,60))));
  else
    [pp,ff]    = fileparts(dicom_meta.Filename);
    [pp,ffnii] = fileparts(nii_file);
    if      ffnii(1)=='o', meta.ScanDayTime = ffnii(2:strfind(ffnii,ff)-1);
    elseif  ffnii(1)=='c', meta.ScanDayTime = ffnii(3:strfind(ffnii,ff)-1);  
    %elseif  meta.ScanDayTime = ffnii(3:strfind(ffnii,ff)-1); 
    else    meta.ScanDayTime = 'XXXXXXXX-XXXXXX';
    end
  end 
  meta.DOS=meta.ScanDayTime(1:8);
  meta.TOS=meta.ScanDayTime(10:end);
 
  % SubjectID (directory name)
  % - PatientID
  % - PatientsName
  try
    if centern>0, SubjectID  = fnameparts{find(strcmpi(opt.center{centern},fnameparts)==1,1,'first')+3};
    else          SubjectID  = repmat('X',[1,opt.maxflength.SubjectID]);
    end
  catch 
    error('MATLAB:vbm_db:unknown_center', ...
         ['There is no center that match for the path ''%s''.' ...
          'Pleace add the center to the ''opt.center'' varialbe!'], ...
          dicom_meta.Filename);
  end
  meta.SubjectID = SubjectID(1:min(opt.maxflength.SubjectID,numel(SubjectID)));
  % check if this exist?!
  % you have to read and compare both volumes!!!


  % SEX
  if      isfield(dicom_meta,'PatientsSex'), meta.Sex = dicom_meta.PatientsSex(1);
  elseif  isfield(dicom_meta,'PatientSex'),  meta.Sex = dicom_meta.PatientSex;
  else    meta.Sex = 'X';
  end
  % DOB, YOB
  if      isfield(dicom_meta,'PatientsBirthDate'), meta.DOB = datestr(dicom_meta.PatientsBirthDate,'yyyymmdd');
  elseif  isfield(dicom_meta,'PatientBirthDate'),  meta.DOB = dicom_meta.PatientBirthDate;
  else    meta.DOB = 'XXXXXXXX';
  end
  try
    meta.YOB = meta.DOB(1:4);
  catch
    fprintf('%s\n',meta.DOB);
    meta.DOB = 'XXXXXXXX';
    meta.YOB = 'XXXX';
  end
    
  % AGE
  if      isfield(dicom_meta,'PatientsAge'), meta.Age = dicom_meta.PatientsAge(2:3);
  elseif  isfield(dicom_meta,'PatientAge'),  meta.Age = dicom_meta.PatientAge(2:3);
  else    meta.Age = 'XX';
  end
  % AGE (floatingpoint years)
  [tmp,tmp,meta.AgeY] = CompleteDates(meta.DOB,meta.DOS,'');
  % AGE (days)
  % Weight
  


  % HANDNESS
  % here we need some addition meta data, because there are dicom informations about handness.
  meta.Handness = 'X';


  % WEIGTING
  Contrast = lower(strrep(dicom_meta.Filename,'/',' ')); 
  if isfield(dicom_meta,'ScanningSequence'), Contrast = [Contrast ' ' dicom_meta.ScanningSequence]; end 
  if isfield(dicom_meta,'Protocol'),         Contrast = [Contrast ' ' dicom_meta.Protocol]; end 
  if isfield(dicom_meta,'ProtocolName'),     Contrast = [Contrast ' ' dicom_meta.ProtocolName]; end
  if isfield(dicom_meta,'SequenceName'),     Contrast = [Contrast ' ' dicom_meta.SequenceName]; end
  contrastn=0; 
  for filesi=1:size(opt.Contrasts,1)
    for jj=1:numel(opt.Contrasts{filesi,2})
      if sum(~isempty(strfind(lower(Contrast),lower(opt.Contrasts{filesi,2}{jj}))))>0, contrastn=filesi; end;
    end
  end 
  if contrastn==0, meta.Contrast = 'XXX'; else  meta.Contrast = opt.Contrasts{contrastn,1}; end
  if strcmp('LOC',meta.Contrast), return; end
  if contrastn==0, fprintf('Unknow Contrast %s\n',Contrast); 
    return;
  end
  % scan pareamter:  TE TR ... 
  %   Fletcher LM et al, Magnetic Resonance in Medicine 1993; 29: 623-630 (for 1.5 Tesla)
  % tissue parameter:           CSF       GM        WM       Skull    Mussel    Fat
  %   - proton density (pd)     
  %   - t1-relaxiontime (s):    0.8-20    0.76-1.08 1.09-2.15 0.5-2.2 0.95-1.82 0.2-075
  %   - t2-relaxion-time (ms):  110-2000  61-100    61-109    50-165  20-67     53-94
  if isfield(dicom_meta,'EchoTime'),       meta.TE = dicom_meta.EchoTime;       else meta.TE = nan; end
  if isfield(dicom_meta,'RepetitionTime'), meta.TR = dicom_meta.RepetitionTime; else meta.TR = nan; end
  if isfield(dicom_meta,'InversionTime'),  meta.TI = dicom_meta.InversionTime;  else meta.TI = nan; end
  % ImagedNucleus???
  % NumberofPhaseEncodingSteps???
  % TransmittingCoil
  % ImagingFrequency
  % PercentPhaseFieldofView
  % Fettunterdrückung???
  % FOV???
  if isfield(dicom_meta,'FlipAngle'),      meta.FlipAngle = dicom_meta.FlipAngle; else meta.FlipAngle = nan; end
  % Echo numbers
  % Echozugwinkel???
  if     meta.TE<20 && meta.TR<2000, meta.Contrast2 = 'T1w';
  elseif meta.TE>20 && meta.TR>2000, meta.Contrast2 = 'T2w';
  elseif meta.TE<20 && meta.TR>2000, meta.Contrast2 = 'PDw'; 
  elseif any(isnan([meta.TE,meta.TR])), meta.Contrast2 = 'LOC'; meta.Contrast = 'LOC'; return;
  else                                  meta.Contrast2 = 'LOC';
  end 
%   if ~strcmp(meta.Contrast,meta.Contrast2)
%     meta;
%   end
  
  
%  meta.ScanTime = TE * Nph * Nac (*Npart)
  % PD:   TR >> T1, TE << T2
  % T1:   TR ~  T1, TE << T2
  % T2:   TR >> T1, TE ~  T2
  
  % QA?
  if isfield(dicom_meta,'SAR'),  meta.SAR = dicom_meta.SAR; end 
  if isfield(dicom_meta,'NumberofAverages'), meta.TR = dicom_meta.NumberofAverages; end  
  
    
  % SEQUENCE, for equal length use fill up with zeros => '%010s'
  meta.Sequence = repmat('X',[1,opt.maxflength.Sequence]);
  if isfield(dicom_meta,'SequenceName')
    Sequence = dicom_meta.SequenceName;
  elseif isfield(dicom_meta,'ScanningSequence')
    Sequence = dicom_meta.ScanningSequence;
    if isfield(dicom_meta,'SequenceVariant') && ...
      ~strcmp(dicom_meta.ScanningSequence,dicom_meta.SequenceVariant)
        Sequence = [Sequence '-' dicom_meta.SequenceVariant];
    end
    if isfield(dicom_meta,'ProtocolName') && ...
      ~strcmp(dicom_meta.ScanningSequence,dicom_meta.ProtocolName)
        Sequence = [Sequence '-' dicom_meta.ProtocolName]; 
    end

  elseif isfield(dicom_meta,'ProtocolName')
    Sequence = dicom_meta.ProtocolName;
  end
  meta.Sequence = fliplr(sprintf(sprintf('%%%ds',opt.maxflength.Sequence), ...
    fliplr(Sequence(1:min(opt.maxflength.Sequence,numel(Sequence))))));

  % you may need a check if the image still exist (read and compare)


  % SCANNER: 
  % MANUFACTURER:
  manu=0; for i=1:size(opt.manufacturer,1), if strfind(lower(strrep(dicom_meta.Manufacturer,'/',' ')),lower(opt.manufacturer{i,1}))>0, manu=i; end; end 
  if manu==0, meta.ScannerManufacturer = 'XX'; else  meta.ScannerManufacturer = opt.manufacturer{manu,2}; end
  % MODELL_
  if isfield(dicom_meta,'ManufacturersModelName'),   meta.ScannerModel    = dicom_meta.ManufacturersModelName; 
  else                                               meta.ScannerModel    = ''; 
  end
  % SORTWARE:
  if isfield(dicom_meta,'SoftwareVersions'),         meta.SoftwareVersions = dicom_meta.SoftwareVersions;
  else                                               meta.SoftwareVersions = '';
  end
  % MAGNETICFIELDSTRENG
  if isfield(dicom_meta,'MagneticFieldStrength')
    meta.MagneticFieldStrength = sprintf('%02.0fT',10*dicom_meta.MagneticFieldStrength);
  else 
    meta.MagneticFieldStrength = 'XXT'; 
  end

  

  
  
  % INSTITUTION
  if isfield(dicom_meta,'InstitutionName'),     meta.InstitutionName   = dicom_meta.InstitutionName;
  else                                          meta.InstitutionName   = '';
  end
  if isfield(dicom_meta,'InstitutionAddress'),  meta.InstitutionAddress = dicom_meta.InstitutionAddress;
  else                                          meta.InstitutionAddress = '';
  end


end
function meta = get_nifti_meta(nii_file,opt)
 fnameparts = eval(['{'' ' lower(strrep(nii_file,'/',''',''')) '''}'])'; fnameparts(1)=[]; 

  % Dicom Information
  Vnii=spm_vol(nii_file);
  
  % RESOLUTION
  % SpacingBetweenSlices, SliceThickness
  meta.DimX = Vnii(1).dim(1);
  meta.DimY = Vnii(1).dim(2);
  meta.DimZ = Vnii(1).dim(3);
  if any(Vnii(1).dim(1:3) == 1), meta.Contrast = 'LOC'; return;  end
  
  % STUDY
  meta.Project = opt.project;

  % GROUP
  groupn=0; for i=1:size(opt.groups,1), if sum(strcmpi(opt.groups{i,1},fnameparts))>0, groupn=i; end; end 
  if groupn==0, meta.GroupID = 'unkown'; else meta.GroupID = opt.groups{groupn,2};  end

  % CENTER
  centern=0; for i=1:size(opt.center,1),  if sum(strcmpi(strrep(fnameparts,'/',' '),opt.center{i,1}))>0, centern=i; end; end 
  if centern==0, meta.Center = 'XXX'; else  meta.Center = opt.center{centern,2}; end


  % DATE
  meta.StudyTime = 'XXXXXXXX-XXXXXX';
  meta.ScanDate  = 'XXXXXXXX';
  meta.ScanTime  = 'XXXXXX'; 

  % SUBJECTID
  if centern>0, meta.Subjectdir = fnameparts{find(strcmpi(opt.center{centern},fnameparts)==1,1,'first')+3};
  else          meta.Subjectdir = 'unknown';
  end
  if strcmp('Leuven',meta.Center); [ppp,fff]=fileparts(fnameparts{end}); meta.Subjectdir = fff; end
  meta.SubjectID = meta.Subjectdir(1:min(opt.maxflength.SubjectID,numel(meta.Subjectdir)));
  meta.Sex = 'X';
  meta.YOB = 'XXXX';
  meta.Age = 'XX';
  meta.Handness = 'X';


  % CONTRAST
  contrastn = 0; 
  Contrast  = lower(strrep(nii_file,'/',' ')); 
  for filesi=1:size(opt.Contrasts,1)
    for jj=1:numel(opt.Contrasts{filesi,2})
      if sum(~isempty(strfind(lower(Contrast),lower(opt.Contrasts{filesi,2}{jj}))))>0, contrastn=filesi; end;
    end
  end 
  if contrastn==0, meta.Contrast = 'XXX'; else  meta.Contrast = opt.Contrasts{contrastn,1}; end
  meta.Contrast2 = meta.Contrast;
%  weightn=0; for i=1:size(opt.Contrasts,2), if ~isempty(cell2mat(strfind(fnameparts,lower(opt.Contrasts{i})))); weightn=i; end; end 
%   if weightn==0, meta.Contrast = 'XX'; else meta.Contrast = opt.Contrasts{weightn}; end
%   meta.Contrast = fnameparts{find(cellfun('isempty',strfind(fnameparts,lower(Contrasts{weightn})))==0,1,'last')};
%   if any(strcmp(meta.Contrast,{'DTI','DWI'})), meta.Contrast = [meta.Contrast '-' fnameparts{end}]; continue;  end

  meta.Sequence = meta.Contrast;

  % MANUFACTURER
  manu=0; for i=1:size(opt.manufacturer,2), if ~isempty(cell2mat(strfind(fnameparts,lower(opt.manufacturer{i}))))>0, manu=i; end; end 
  if manu==0, meta.ScannerManufacturer = 'XX'; else  meta.ScannerManufacturer = opt.manufacturer{manu,2}; end


  % MAGNETICFIELDSTRENG
  meta.MagneticFieldStrength = 'XXT'; 
end

function [ScanID,DB] = create_ScanID(meta,RAWdir,opt,DB)
% this function has to check a lot of things...
% - images with the same name, but different volume size are no rescans - maybe 
%   * scouts (less than 10 pixel in one direction)      => irgnore
%   * differnt resolution (clear)                       => ????
%   
% - images have the same resolution, weigthing and date => rescans
% - images have the same resolution, weigthing and date => followups

% first, we need all existing files and there metafiles in our database 
% (because, it should be a unique, decentral concept we have to find them!)

  if isfield(meta,'Contrast') && any(strcmp(opt.Contrast,meta.Contrast))    % import only if the contrast match
                                                                            % further fields: multiscan-handling
                                                    
    if ~exist('DB','var') || isempty(DB)  
      stime = clock;                                               
      RAWs = findfiles(RAWdir,'*.nii')'; 
      FmetaRAWs = cell(size(RAWs)); ScanIDs = zeros(size(RAWs)); ImportIDs = cell(size(RAWs));
  
      if numel(RAWs)>0
        for RAWsi=1:numel(RAWs)
          [pp,ff]          = fileparts(RAWs{RAWsi});
          FmetaRAWs{RAWsi} = fullfile(pp,[opt.meta.prefix ff '.xml']);
          if exist(FmetaRAWs{RAWsi},'file')
            tmp              = vbm_io_xml(FmetaRAWs{RAWsi}); 
            ScanIDs(RAWsi)   = double(tmp.ScanID);
            ImportIDs{RAWsi} = tmp.ImportID;
          end
        end
        DB = [num2cell(ScanIDs),ImportIDs];
        save(opt.vbmDBfile,'DB');

        fprintf('Found %d scans in %d seconds and save the enties in ''%s''!\n',...
          size(DB,1),etime(stime,clock),opt.vbmDBfile);
      else
        fprintf('##########');
        DB = {};  
      end
      
    end
    
    %if size(DB,1)==10,  DB(1:min(10,size(DB,1)),:), end
    if ~isempty(DB)
      % now, we have get the meta data and check, if they fit to our new case
      metaequal = find(~cellfun('isempty',strfind(DB(:,2),meta.ImportID)),1,'first');
          
      % check if there is a file with the same id
      if ~isempty(metaequal)
        ScanID = DB{metaequal,1};    
      else
        ScanID = max(cell2mat(DB(:,1)))+1; %str2double(char(
        DB = [DB;{ScanID,meta.ImportID}];
      end
    else
      ScanID = 1;
      DB = {ScanID,meta.ImportID};
    end
    ScanID = cast(ScanID,opt.meta.ScanIDtype);
  else
    ScanID = cast(0,opt.meta.ScanIDtype);
  end
end
function DB = insert(nifti,meta,opt,DB)
  % check database for a old/new ScanID depending on the import options like Contrast... 
  if isfield(meta,'dcm2niiname'), meta.ImportID = meta.dcm2niiname; 
  else                            meta.ImportID = nifti;
  end
  [meta.ScanID,DB] = create_ScanID(meta,opt.DBRAW,opt,DB);
  
  % import only if we get a valid id
  if isfield(meta,'ScanID') && meta.ScanID>0  
 
    % pathname and filename
    meta.fpath = metaname(opt,meta,'pathname',opt.DBRAW);
    meta.fname = metaname(opt,meta,'filename'); 
    vbmDBnii = fullfile(meta.fpath,[meta.fname '.nii']);
    vbmDBxml = fullfile(meta.fpath,[opt.meta.prefix meta.fname '.xml']);
    
    
    % move niftis from tmp-dir, write the meta-file and export to mysql
    if exist(nifti,'file')
      if ~exist(meta.fpath,'dir'); mkdir(meta.fpath); end %; end
      
      % move new
      copyfile(nifti,vbmDBnii ); 

      % write meta data to the nii-file
      vbm_io_struct(vbmDBxml,'ScanID',meta);

      % write meta data to mySQL database
      vbm_io_mysql(opt.mysql,meta);
    
      fprintf('%s\n',fullfile(meta.fpath,meta.fname));
    
    end
    
    save(opt.vbmDBfile,'DB');
  end
end

function s=FieldName(s) % fieldlength
  switch class(s),
    case {'double','single'}, s=num2str(s,'%08.2f'); 
    case {'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}
      s=num2str(s,sprintf('%%0%dd',length(num2str(intmax(class(s))))));
    case 'char'
      ts = isstrprop(s,'alphanum'); ts(strfind(s,'.'))=1; s(ts); %???
      for i='#_*/\ ', s = strrep(s,i,''); end 
    case 'cell'
      s=FieldName(char(s));
 end
 
end
function [DOB,DOS,age]=CompleteDates(DOB,DOS,age)
% calculates one missing variable
%   DOB = DayOfBirth
%   DOS = DayOfScan
%   age in years 
%   txt ... add string that explain estimated values...
%
% [value type
% type: 0=
% age_d  = age in days
% age_y  = age in years
% age_ry = 
%
% DOB = DayOfBirth 
% YOB = YearOfBirth
% DOS = DayOfScan
% YOS = YearOfScan

  missed=[isempty(DOB),isempty(DOS),isempty(age)];

  switch char(missed*'m')
    case 'm  '
      DOB = datestr(datenum(DOS,'yyyy-mm-dd') - 365*age,'yyyy-mm-dd');
      DOB = DOB + sum(mod(str2double(DOB(1:4)):1:str2double(DOS(1:4)),4)==3)/365;
    case ' m '
      DOS = datestr(datenum(DOB,'yyyy-mm-dd') + 365*age,'yyyy-mm-dd');
      DOS = DOS + sum(mod(str2double(DOB(1:4)):1:str2double(DOS(1:4)),4)==3)/365;
    case '  m'
      age = (datenum(DOS,'yyyy-mm-dd') - datenum(DOB,'yyyy-mm-dd') - ...
            sum(mod(str2double(DOB(1:4)):1:str2double(DOS(1:4)),4)==3))/365; % schaltjahre mit einem tag weniger
    case '   '
      % test it
  end
end
function fname=metaname(opt,meta,cfield,fnameI)
    % correct for special characters in the metafiles to avoid problems with the filename
    FN=fieldnames(meta); for fn=1:numel(FN), meta.(FN{fn})=FieldName(meta.(FN{fn})); end

    if ~exist('fnameI','var'), fname = ''; else fname = fnameI; end
    for pi=1:numel(opt.(cfield))
      subdir=''; 
      for pj=1:numel(opt.(cfield){pi})
        if any(strcmp({'-','_','.'},opt.(cfield){pi}{pj}))
          subdir=[subdir fieldname(opt.(cfield){pi}{pj})];  %#ok<AGROW>
        else
          if isnumeric(meta.(opt.(cfield){pi}{pj}))
            subdir=[subdir num2str(meta.(opt.(cfield){pi}{pj}),'%d')];  %#ok<AGROW>
          else
            subdir=[subdir meta.(opt.(cfield){pi}{pj})];  %#ok<AGROW>
          end
        end
      end
      if exist('fnameI','var'),
        fname = fullfile(fname,subdir); 
      else
        fname = [fname '_' subdir];%#ok<AGROW>
      end
    end
    if ~exist('fnameI','var'), fname(1)=[]; end
end

 