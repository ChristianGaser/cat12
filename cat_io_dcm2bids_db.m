function varargout = cat_io_dcm2bids_db(action,varargin)
%cat_io_dcm2bids_db. Internal database of cat_io_dcm2bids (catDCM2BIDSdb).
%  Each imported scan is stored in a scan directory of the database
%    sub-<DeviceSerialNumber>-<PatientID>/ses-<StudyID>[-<ScanDate>]/snr-<SeriesNumber>-<ProtocolName>/
%  with the JSON sidecar "folder=protocol=time=series.json", the import 
%  information "catDCM2BIDSimport_folder=protocol=time=series.json" (DICOM
%  directory or source file of a BIDS dataset), and the image (only if it 
%  was converted). These functions define the database layout, the identity
%  of the scans, and the access to the import information, i.e., only this 
%  file has to be adapted for another kind of database.
%
%  varargout = cat_io_dcm2bids_db(action,varargin)
%
%  Actions (see the help of the local functions):
%    Pdbdir = cat_io_dcm2bids_db('getDBdir',V)
%    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Pdir)
%    V = cat_io_dcm2bids_db('assurePatientDCMfields',V)
%    id = cat_io_dcm2bids_db('getScanID',V,fnameparts)
%    key = cat_io_dcm2bids_db('scanKey',V)
%    [Pdcm, Psnr, Popt, Psrc, Psrcsnr] = cat_io_dcm2bids_db('getImportIndex',Pmdbdir)
%    S = cat_io_dcm2bids_db('readImportJSON',P)
%    cat_io_dcm2bids_db('writeImportJSON',P,S)
%    cat_io_dcm2bids_db('checkGzipi',Pdbdir,gzipi,run)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'getDBdir','getDBscanJSONs','assurePatientDCMfields','getScanID', ...
          'scanKey','getImportIndex','readImportJSON','writeImportJSON', ...
          'checkGzipi'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_db:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Pdbdir = getDBdir(V)
%getDBdir. Scan directory in the database (relative to the database), with 
%  the scan date in the session directory if available.
  ses = sprintf('ses-%s', V.StudyID); 
  if ~isnat(datetime(V.ScanDate)), ses = sprintf('%s-%s', ses, datetime( V.ScanDate , 'Format','yyyyMMdd')); end
  Pdbdir = fullfile( ...
    sprintf('sub-%s-%s',   V.DeviceSerialNumber, V.PatientID), ses, ...
    sprintf('snr-%04d-%s', V.SeriesNumber, V.ProtocolName)); 
end
% =========================================================================
function Pjson = getDBscanJSONs(Pdir)
%getDBscanJSONs. dcm2niix sidecars of a scan directory of the database. 
%  These are the JSON files with the dcm2niix name pattern 
%  "folder=protocol=time=series" without the files of the import and the
%  further processing (catDCM2BIDSimport_*, catDCM2BIDSqc_*, ...). 
  Pjson = cat_vol_findfiles(Pdir,'*=*=*=*.json',struct('depth',1)); 
  Pjson( strncmp(spm_file(Pjson,'basename'),'catDCM2BIDS',11) ) = []; 
end
% =========================================================================
function V = assurePatientDCMfields(V) 

  % correction of special characters
  Lold = {'ä','ü','ö'};    L2old = {'ß','Ã¤','Ã¼','Ã¶','ÃŸ','Ã„','Ã–','Ãœ'}; 
  Lnew = {'ae','ue','oe'}; L2new = {'ss','ae','ue','oe','ss','Ae','Oe','Ue'}; 
  FN   = {'PatientName','InstitutionName','InstitutionalDepartmentName','InstitutionAddress', ...
          'ProcedureStepDescription','ReferringPhysicianName','ImageComments'};
  FN   = intersect(FN,fieldnames(V)); 
  for fni = 1:numel(FN)
    V.(FN{fni}) = cat_io_strrep(V.(FN{fni}),L2old,L2new);
    V.(FN{fni}) = cat_io_strrep(V.(FN{fni}),Lold,Lnew);
    V.(FN{fni}) = cat_io_strrep(V.(FN{fni}),upper(Lold),upper(Lnew));
  end

  
  %% fields that are required later but may be missing (e.g. in case of 
  %  already anonymized DICOMs)
  % subject ID: use the name or the study UID if no PatientID is available,
  % to avoid that different subjects get the same ID
  if ~isfield(V,'PatientID') || isempty(V.PatientID)
    if isfield(V,'PatientName') && ~isempty(V.PatientName)
      V.PatientID = regexprep(char(V.PatientName),'[^a-zA-Z0-9]','');
    elseif isfield(V,'StudyInstanceUID') && ~isempty(V.StudyInstanceUID)
      uid = regexprep(char(V.StudyInstanceUID),'[^0-9]',''); 
      V.PatientID = ['UID' uid(max(1,end-9):end)];
    else
      V.PatientID = 'NA'; 
    end
  end
  if isnumeric(V.PatientID), V.PatientID = num2str(V.PatientID); end
  
  % numeric fields
  FN = {'SeriesNumber'}; 
  for fni = 1:numel(FN)
    if ~isfield(V,FN{fni}) || isempty(V.(FN{fni})), V.(FN{fni}) = 0; end
  end

  % need this for the session
  if ~isfield(V,'ScanDate') 
    try
      V.ScanDate = datetime(V.AcquisitionDateTime,'Format','uuuu-MM-dd');
    catch
      V.ScanDate = NaT('Format','uuuu-MM-dd'); % missing or unreadable date
    end
  end

  % (re)estimate age
  if isfield(V,'PatientBirthDate') && ~isempty(V.PatientBirthDate)
    try
      ScanDate     = datetime(V.ScanDate,'Format','uuuu-MM-dd');
      BirthDate    = datetime(V.PatientBirthDate,'Format','uuuu-MM-dd');
      V.PatientAge = char(duration( ScanDate - BirthDate, 'Format','y'));
      V.PatientAge = str2double(V.PatientAge(1:end-4)); 
    catch
      V.PatientAge = NaN; 
    end
  elseif ~isfield(V,'PatientAge') 
    V.PatientAge = NaN;
  end

  % numeric fields
  FN = {'PatientAge', 'PatientWeight'}; %, 'PatientHeight'};
  for fni = 1:numel(FN)
    if ~isfield(V,FN{fni}) || isempty(V.(FN{fni})), V.(FN{fni}) = NaN; end
  end

  % text fields (after the age estimation that uses the PatientBirthDate)
  FN = {'PatientSex','PatientName','PatientBirthDate','AcquisitionDateTime', ...
        'StudyID','DeviceSerialNumber','SeriesDescription'};
  for fni = 1:numel(FN)
    if ~isfield(V,FN{fni}) || isempty(V.(FN{fni})), V.(FN{fni}) = 'NA'; end
    if isnumeric(V.(FN{fni})), V.(FN{fni}) = num2str(V.(FN{fni})); end
  end
end
% =========================================================================
function id = getScanID(V,fnameparts)
%getScanID. Identity of a converted scan, i.e., the DICOM series and the 
%  dcm2niix image suffix of the series number (last part of the filename 
%  "folder=protocol=time=series", e.g. 7, 7a, or 17_e2). The folder name is
%  ignored as the same series can be converted from different directories. 
  if isfield(V,'SeriesInstanceUID') && ~isempty(V.SeriesInstanceUID)
    id = sprintf('%s|%s', V.SeriesInstanceUID, fnameparts{end}); 
  else
    % without series UID (e.g. BIDS input): site/dataset, subject, and name
    id = sprintf('%s|%s|%s', V.DeviceSerialNumber, V.PatientID, strjoin(fnameparts(2:end),'=')); 
  end
end
% =========================================================================
function key = scanKey(V)
%scanKey. Subject/session/series identifier to count images of one series.
  if ~isfield(V,'SeriesNumber') || isempty(V.SeriesNumber), key = ''; return; end
  key = sprintf('%d',V.SeriesNumber);
  FN  = {'DeviceSerialNumber','PatientID','StudyID','StudyInstanceUID'}; % site/dataset and session also for BIDS input
  for fni = 1:numel(FN)
    if isfield(V,FN{fni}), key = [key '|' char(string(V.(FN{fni})))]; end %#ok<AGROW>
  end
end
% =========================================================================
function [Pdcm, Psnr, Popt, Psrc, Psrcsnr] = getImportIndex(Pmdbdir)
%getImportIndex. DICOM directories (Pdcm), their scan directories (Psnr), 
%  and dcm2niix options (Popt) of the DICOM imports, and source files (Psrc)
%  and scan directories (Psrcsnr) of BIDS imports of the import information 
%  (catDCM2BIDSimport_*.json) in the database. 
  Pimp = cat_vol_findfiles(Pmdbdir,'catDCM2BIDSimport_*.json',struct('depth',4)); 
  Pdcm = {}; Psnr = {}; Popt = {}; Psrc = {}; Psrcsnr = {}; 
  for ii = 1:numel(Pimp)
    try
      I = readImportJSON(Pimp{ii}); 
      if isfield(I,'DICOMdirectory')
        Pdcm{end+1,1} = I.DICOMdirectory; Popt{end+1,1} = I.dcm2niixOptions; %#ok<AGROW>
        Psnr{end+1,1} = fileparts(Pimp{ii}); %#ok<AGROW>
      elseif isfield(I,'SourceFile')
        Psrc{end+1,1} = I.SourceFile; Psrcsnr{end+1,1} = fileparts(Pimp{ii}); %#ok<AGROW>
      end
    end
  end
end
% =========================================================================
function S = readImportJSON(P)
%readImportJSON. Read the import information (UTF-8, e.g. for paths with
%  special characters). 
  fid = fopen(P,'r','n','UTF-8'); txt = fread(fid,'*char')'; fclose(fid); 
  S   = jsondecode(txt); 
end
% =========================================================================
function writeImportJSON(P,S)
%writeImportJSON. Write the import information (UTF-8). 
  cat_io_dcm2bids_helper('writeJSON',P,S); 
end
% =========================================================================
function checkGzipi(Pdbdir, gzipi, run)
%checkGzipi. Assure the internal gzip status of the NIfTIs of the database 
%  defined by gzipi (e.g. after interruptions of previous imports).
  pverb = 1; 
  if exist(Pdbdir,'dir') && run
    if gzipi
      Punpacked = cat_vol_findfiles(Pdbdir,'*.nii'); 
      for fi = 1:numel(Punpacked)
        if ~exist([Punpacked{fi} '.gz'],'file')
          if pverb, fprintf('Updating database storing gzipped NIFTIs... '); pverb = 0; end
          gzip(Punpacked{fi});
        end
        delete(Punpacked{fi}); 
      end
    else
      Ppacked = cat_vol_findfiles(Pdbdir,'*.nii.gz'); 
      for fi = 1:numel(Ppacked)
        if ~exist([Ppacked{fi}(1:end-3)],'file')
          if pverb, fprintf('Updating database storing NIFTIs... '); pverb = 0; end
          gunzip(Ppacked{fi});
        end
        delete(Ppacked{fi}); 
      end
    end
  end
  if pverb == 0, fprintf(' done.\n'); end % if something was printed before 
end
