function vbm_db_batch(DCMdir,MAINdir)
% ______________________________________________________________________
% WARNING: THIS IS A PRIVATE FUNCTION THAT IS STILL IN DEVELOPMENT! 
%          IT ONLY WORKS ON OUR SERVERS, AND WILL MAYBE BECOME PART OF
%          THE VBM TOOLBOX IN THE FUTURE.
% ______________________________________________________________________
% Aufruf von vbm_db zur Erstellung der blabla Datenbank:
% ______________________________________________________________________
% Bei vbm_db handelt es um eine dezentral Datenbank, d.h. jede Datei hat
% ein XML-File mit den zusätzlichen Informationen und es gibt kein
% zentrales Datenbankfile. Die XML-Metadaten besitzen einen Prefix, so 
% das als Dateiblock besser von den Bildern abgrenzbar sind. Diese 
% Designentscheidung beruht primär auf der dezentral orientierten 
% Struktur von SPM!
% vbm_db durchsucht das DCMdir nach DICOM und NIFTI Files, konvertiert
% diese, schreibt sie als NIFTI in eine definierbare Verzeichnisstruktur 
% und fügt die Zusatzdaten als XML-File dazu. 
% ACHTUNG vbm_db ist in der Entwicklung und sorgt aktuell nicht für
% Datenkonsistenz - also nicht Vergleichbar mit einem DBMS!
% ______________________________________________________________________
% $Id$


% Verzeichnisse:
% ______________________________________________________________________
% Die Daten werden aus dem DCMdir in das Datenbankverzeichnis
% DBdir/DBname importiert. Das DCMdir darf ein Subverzeichnis sein, so
% das z.B. nur die Daten aus diesem Unterverzeichnis importiert werden.
  if ~exist('MAINdir','var') || ~exist(MAINdir,'dir')
    host = system('hostname -s');
    switch host
      case 'nedigs10'
        MAINdir = '../MRI/MRI_blabla';
      otherwise
        MAINdir = '../MRI/MRI_blabla';
    end     
  end

  % Importerzeichnis mit den RAW-Daten 
  if ~exist('DCMdir','var') || ~exist(DCMdir,'dir')
    DCMdir = fullfile(MAINdir,'repo_freezes/20130917/'); % DICOM/ANALYSE RAW Data 
  end
  DBdir  = fullfile(MAINdir,'jenana');                                % Hauptverzeichnis der Datenbank
  DBname = fullfile(MAINdir,'nisalsDB');                              % Name der Datenbank (und nicht unbedingt des Projektes!)
  
  opt.DBTMP   = fullfile(MAINdir,'jenana/vbmDBimport');               % temporäres Importverzeichnis
  
  % Anwendungsverzeichnisse von MRIcron, Freesurfer, FSL und MySQL
  opt.MATLAB      = fullfile(MAINdir,'../matlab/');
  opt.VBM         = fullfile(MAINdir,'../matlab/SPM8R4290_VBM12+');
  opt.MRIcron     = fullfile(MAINdir,'../mricron');
  opt.FreeSurfer  = '/usr/local/freesurfer';
  opt.FSL         = fullfile(MAINdir,'../fsl5');
  opt.MySQL       = '/usr/bin'; % /usr/bin/...
  opt.initPATH    = sprintf( [...
      'export MRIcron=%s;' ...
      'export FSLDIR=%s; sh ${FSLDIR}/etc/fslconf/fsl.sh;' ...
      'export FREESURFER_HOME=%s; ${FREESURFER_HOME}/SetUpFreeSurfer.sh;' ...
      'MNI_DIR=${FREESURFER_HOME}/mni; ' ...
      'PATH=$PATH:$FSLDIR/bin:$MNI_DIR/bin:$MRIcron; ' ...
      'source $FREESURFER_HOME/FreeSurferEnv.sh; ' ...
      'export FSLOUTPUTTYPE=NIFTI NIFTIexport=NIFTI; ' ...
      ],opt.MRIcron,opt.FSL,opt.FreeSurfer);
  
    
  % MY-SQL options
  % ______________________________________________________________________
  opt.mysql.path    = opt.MySQL ;          % Verzeichnis der MYSQL DB
  opt.mysql.user    = '..';                % DB Nutzer
  opt.mysql.pwd     = '..';                % DB Nutzer Passwort 
  opt.mysql.db      = '..';                % Name der SQL Datenbank
  opt.mysql.tname   = DBname;              % Name der SQL Tabelle
  opt.mysql.tid     = 'ScanID';            % ID Feld eines Bildes
    
  
  % other options
  % ______________________________________________________________________
  opt.importdb    = 'blabla';                  % Datenbankimporttyp 
  opt.meta.prefix = 'data_';                   % Metadaten-Prefix 
  opt.anomize     = 1;                         % Anonymisieren der Bilder
  opt.Contrast    = {'T1w','T2w','PDw','DTI'}; % Import nur von ... {'T1w','T2w','PDw','DTI',...}
  
    
  % Pfad und Dateinamen durch Kombination mit Dateidaten:
  % ______________________________________________________________________
  %   Project, GroupID, Center,
  %   ScannerManufacturer, MagneticFieldStrength, Sequence, Contrast, 
  %   SubjectID, YOB, DOB, Sex, Age, Handness,
  %   Scandate, Scantime
  %
  % Beispiel
  %   input: {{'SubjectID'} {'Contrast'} {'ScanDate' '-' 'ScanTime'} 
  %   ouput: MaxMueller_T1_20120809-145322
  %
  % Da sich alle Felder sich ändern könnten sollte man an dieser Stelle mit der 
  % Benennung ruhig geizig sein.
  % ______________________________________________________________________
  opt.pathname = {
    {'Contrast'}
    {'Center'}
  };
  opt.filename = {
    {'Project'}
    {'GroupID'}  % subject groups like HC, AD, ALS, ...
    {'Center'} 
    {'ScannerManufacturer' 'MagneticFieldStrength'}
    {'SubjectID'}
  % {'YOB' 'SubjectID'}
  % {'Sex' 'Age'} % 'Handness'}
    {'Contrast'} % '-' 'Contrast2'} % '-' 'Sequence'}
  % {'ScanDate' '-' 'ScanTime'}
    {'ScanID'};  % this field is necessary for unique idententification!
  }; 
  
  
  % Center
  % Zuordnung der Zentren zu spezifischen IDs
  % ______________________________________________________________________
    opt.center = { % directory name, centerID
      'blablabla'      'bla'
      'blublublu'      'blu'
      ... die Sachen hier unten sind unwichtig...
      'Unkown'        'UKN' % special group - if we don't know it 
      'Other'         'OTH' % special group - if we do not like to know it :D
      'Bad'           'BAD' % special group - bad datasets
      'Unkown'        'XXX' % special group - internal if we don't know it
    };
  
  % vbm zum Pfad dazufügen
  addpath(genpath(opt.VBM));
  addpath(genpath(opt.MATLAB));
  
  % Hauptaufruf:
  %fprintf('vbm_db(''import'',''%s'',''%s'',''%s'');\n',DCMdir,DBdir,DBname);
  vbm_db('import',DCMdir,DBdir,DBname,opt)
  
end