function varargout = cat_io_dcm2bids_importDCM(action,varargin)
%cat_io_dcm2bids_importDCM. DICOM import of cat_io_dcm2bids (dcm2niix).
%  The DICOM headers of new scans are imported into the database 
%  (importDICOMheaders, dcm2niix header only, i.e. fast), whereas the images
%  are converted only if they are required (convertImages, e.g. not for the
%  JSON-only outputs). convertImages also copies the images of imported BIDS
%  datasets (see cat_io_dcm2bids_importBIDS). The dcm2niix executable 
%  (opts.Pdcm2nii) is detected only if DICOM data has to be read or converted
%  (see getDCM2NIIX), i.e., dcm2niix is not required for BIDS inputs or the 
%  export of already converted data.
%
%  varargout = cat_io_dcm2bids_importDCM(action,varargin)
%
%  Actions (see the help of the local functions):
%    Psnr = cat_io_dcm2bids_importDCM('importDICOMheaders',Pdcmdirs,Pmdbdir,Ptmp,opts)
%    msg = cat_io_dcm2bids_importDCM('convertImages',Psnr,Pmdbdir,Ptmp,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'importDICOMheaders','convertImages'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_importDCM:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Psnr = importDICOMheaders(Pdcmdirs, Pmdbdir, Ptmp, opts)
%importDICOMheaders. Import the DICOM headers of new scans into the database.
%  The DICOM header of each directory (without its sub-directories) is 
%  converted to JSON sidecars by dcm2niix (header only, i.e. fast, see 
%  runDCM2NIIX). Localizers and scouts are ignored. New scans get the scan 
%  directory of the database 
%    sub-DeviceSerialNumber-PatientID/ses-StudyID-ScanDate/snr-SeriesNumber-ProtocolName 
%  with the JSON sidecar "folder=protocol=time=series.json" and the import 
%  information "catDCM2BIDSimport_folder=protocol=time=series.json" with the
%  DICOM directory and the dcm2niix options (used to convert the image only
%  if required, see convertImages). The import information of a directory 
%  is written after all its scans, i.e., an interrupted directory is read 
%  again. DICOM directories that were imported before (with the same 
%  options) are not read again (only with opts.rerun). 
%  The same scan (series and image) from another DICOM directory is a copy
%  if both directories contain the same files (not imported again) or 
%  otherwise a series that is split over several directories, which is 
%  incomplete and not imported (error message), as dcm2niix converts each 
%  directory separately and does not warn about incomplete series. 
%  Psnr is the sorted list of the scan directories of all scans of the 
%  DICOM directories (new and previously imported). 

  Psnr = {}; nnew = 0; nold = 0; 
  if isempty(Pdcmdirs), return; end
  fprintf('DCM2BIDS - import DICOM headers of %d directories\n', numel(Pdcmdirs)); 
  
  % DICOM directories imported before (with the same dcm2niix options)
  [Pimpdcm, Pimpsnr, Pimpopt] = cat_io_dcm2bids_db('getImportIndex',Pmdbdir); 
  dopt  = dcm2niixOptions(opts); 
  ver   = ''; % dcm2niix version (only if required)
  seen  = containers.Map('KeyType','char','ValueType','any'); % scans of this run
  flist = containers.Map('KeyType','char','ValueType','any'); % files of directories
  split = {}; % IDs of split series 

  for di = 1:numel(Pdcmdirs)
    % previously imported directory: scans of the database
    imported = strcmp(Pimpdcm, Pdcmdirs{di}) & strcmp(Pimpopt, dopt); 
    if any(imported) && ~opts.rerun
      Psnr = [Psnr; Pimpsnr(imported)]; %#ok<AGROW>
      nold = nold + 1; 
      continue
    end

    % directories without files (e.g. only with sub-directories)
    if isempty(dirFiles(Pdcmdirs{di},flist)), continue; end

    % header only conversion of the directory
    if isempty(ver)
      if isempty(opts.Pdcm2nii), opts.Pdcm2nii = getDCM2NIIX; end % only if required
      [~,ver] = system(sprintf('"%s" --version',opts.Pdcm2nii)); 
      ver = regexp(ver,'v[0-9]+\.[0-9]+\.[0-9]+','match','once'); 
    end
    Ptmpdi = fullfile(Ptmp,sprintf('hdr%06d',di)); 
    runDCM2NIIX(Pdcmdirs{di}, Ptmpdi, opts, 1); 
    Pjson  = cat_vol_findfiles(Ptmpdi,'*.json',struct('depth',1)); 
    rec    = struct('Pimp',{},'S',{}); % import information of this directory
    for ji = 1:numel(Pjson)
      fnameparts = strsplit(spm_file(Pjson{ji},'basename'),'='); 
      if numel(fnameparts) < 4 || any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        continue
      end
      V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pjson{ji})); 
      V.ProtocolName = fnameparts{2}; 
      Pdbdirpath = fullfile(Pmdbdir,cat_io_dcm2bids_db('getDBdir',V)); 
      Pjsondb    = fullfile(Pdbdirpath,spm_file(Pjson{ji},'filename')); 
      scanid     = cat_io_dcm2bids_db('getScanID',V,fnameparts); 
      if any(strcmp(split,scanid)), continue; end

      % the same scan from this run or in the database (also with another 
      % name, i.e. from another DICOM directory) 
      Pexist = ''; Pprev = ''; 
      if isKey(seen,scanid)
        Pexist = seen(scanid).Pjsondb; Pprev = seen(scanid).dir; 
      elseif exist(Pdbdirpath,'dir')
        Pother = cat_io_dcm2bids_db('getDBscanJSONs',Pdbdirpath); 
        same   = cellfun(@(P) strcmp(cat_io_dcm2bids_db('getScanID',cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(P)), ...
          strsplit(spm_file(P,'basename'),'=')), scanid), Pother); 
        if any(same)
          Pexist = Pother{find(same,1)}; 
          Pimp   = spm_file(Pexist,'prefix','catDCM2BIDSimport_'); 
          if exist(Pimp,'file'), I = cat_io_dcm2bids_db('readImportJSON',Pimp); Pprev = I.DICOMdirectory; end
        end
      end
      if ~isempty(Pexist) && ~isempty(Pprev) && ~strcmp(Pprev,Pdcmdirs{di})
        if exist(Pprev,'dir') && ~isequal(dirFiles(Pprev,flist),dirFiles(Pdcmdirs{di},flist))
          % split series: incomplete, not imported 
          cat_io_cprintf('err',['  Series %d "%s" (%s) is split over several directories and therefore ' ...
            'incomplete - not imported:\n    %s\n    %s\n'], V.SeriesNumber, V.ProtocolName, ...
            fileparts(cat_io_dcm2bids_db('getDBdir',V)), Pprev, Pdcmdirs{di}); 
          split{end+1} = scanid; %#ok<AGROW>
          if isKey(seen,scanid) && seen(scanid).new
            revokeScan(seen(scanid).Pjsondb); nnew = nnew - 1; 
          end
          continue
        end
        if exist(Pprev,'dir')
          % copy of the data: not imported again
          Psnr{end+1,1} = fileparts(Pexist); %#ok<AGROW>
          continue
        end
        % the DICOM directory of the existing scan is no longer available: 
        % use this directory for the image conversion 
      end

      % new scan: JSON sidecar and import information (written after the 
      % directory is complete)
      if isempty(Pexist), Pexist = Pjsondb; end
      new = ~exist(Pexist,'file'); 
      if ~exist(fileparts(Pexist),'dir'), mkdir(fileparts(Pexist)); end
      if new || opts.rerun
        copyfile(Pjson{ji}, Pexist); 
        nnew = nnew + new; 
      end
      rec(end+1).Pimp = spm_file(Pexist,'prefix','catDCM2BIDSimport_'); %#ok<AGROW>
      rec(end).S      = struct( ...
        'DICOMdirectory',  Pdcmdirs{di}, ...
        'ImportDate',      char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
        'dcm2niix',        ver, ...
        'dcm2niixOptions', dopt, ...
        'ScanID',          scanid); 
      seen(scanid)  = struct('dir',Pdcmdirs{di},'Pjsondb',Pexist,'new',new); 
      Psnr{end+1,1} = fileparts(Pexist); %#ok<AGROW>
    end

    % import information of the complete directory
    for ri = 1:numel(rec)
      if exist(strrep(rec(ri).Pimp,'catDCM2BIDSimport_',''),'file') % not revoked
        cat_io_dcm2bids_db('writeImportJSON',rec(ri).Pimp, rec(ri).S); 
      end
    end
    rmdir(Ptmpdi,'s'); 
  end

  % existing scan directories (without revoked split series)
  Psnr = unique(Psnr); 
  Psnr = Psnr( cellfun(@(P) exist(P,'dir') && ~isempty(cat_io_dcm2bids_db('getDBscanJSONs',P)), Psnr) ); 
  fprintf('  %d new scans, %d directories imported before\n\n', nnew, nold); 
end
% =========================================================================
function F = dirFiles(P,flist)
%dirFiles. Sorted names of the (non-hidden) files of a directory (cached in
%  the containers.Map flist). 
  if isKey(flist,P), F = flist(P); return; end
  D = dir(P); 
  F = sort({D( ~[D.isdir] & ~strncmp({D.name},'.',1) ).name}); 
  flist(P) = F; %#ok<NASGU> handle object
end
% =========================================================================
function revokeScan(Pjsondb)
%revokeScan. Remove a scan that was imported in this run (sidecar and import
%  information, and the scan directory if it is empty). 
  Pimp = spm_file(Pjsondb,'prefix','catDCM2BIDSimport_'); 
  if exist(Pjsondb,'file'), delete(Pjsondb); end
  if exist(Pimp,'file'), delete(Pimp); end
  D = dir(fileparts(Pjsondb)); 
  if all([D.isdir]), rmdir(fileparts(Pjsondb)); end
end
% =========================================================================
function msg = convertImages(Psnr, Pmdbdir, Ptmp, opts)
%convertImages. Convert the missing images of the scans of the scan 
%  directories Psnr. The DICOM directory of the import information 
%  (catDCM2BIDSimport_*.json) of each scan is converted once with its 
%  dcm2niix options. The image (and bval/bvec) and its JSON sidecar of the 
%  full conversion replace the sidecar of the header-only import (they can 
%  differ, e.g. for name conflicts of dcm2niix within a directory, where 
%  the further images "...a" are added as new scans). 
%  msg is a containers.Map with the reason for scans without image. 
  msg = containers.Map('KeyType','char','ValueType','char'); 
  if opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end

  % scans without image grouped by DICOM directory and dcm2niix options
  grp = containers.Map('KeyType','char','ValueType','any'); 
  for si = 1:numel(Psnr)
    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{si}); 
    for ji = 1:numel(Pjson)
      if exist(spm_file(Pjson{ji},'ext',niiext),'file'), continue; end
      Pimp = spm_file(Pjson{ji},'prefix','catDCM2BIDSimport_'); 
      if ~exist(Pimp,'file'), msg(Pjson{ji}) = 'no import information'; continue; end
      I = cat_io_dcm2bids_db('readImportJSON',Pimp); 
      if isfield(I,'SourceFile') 
        % BIDS import: copy the image (and bval/bvec)
        if ~exist(I.SourceFile,'file'), msg(Pjson{ji}) = 'source image not available'; continue; end
        copyImage(I.SourceFile, spm_file(Pjson{ji},'ext',niiext)); 
        if isfield(I,'SourceBval') && ~isempty(I.SourceBval), copyfile(I.SourceBval, spm_file(Pjson{ji},'ext','.bval')); end
        if isfield(I,'SourceBvec') && ~isempty(I.SourceBvec), copyfile(I.SourceBvec, spm_file(Pjson{ji},'ext','.bvec')); end
        continue
      end
      if ~exist(I.DICOMdirectory,'dir'), msg(Pjson{ji}) = 'DICOM directory not available'; continue; end
      key = [I.DICOMdirectory '|' I.dcm2niixOptions]; 
      if isKey(grp,key), grp(key) = [grp(key); Pjson(ji)]; else, grp(key) = Pjson(ji); end
    end
  end
  if grp.Count == 0, return; end
  fprintf('DCM2BIDS - convert images of %d DICOM directories\n', grp.Count);
  if isempty(opts.Pdcm2nii), opts.Pdcm2nii = getDCM2NIIX; end % only if required

  keys = grp.keys; 
  for gi = 1:numel(keys)
    Pdcm   = keys{gi}(1:find(keys{gi}=='|',1,'last')-1); 
    gopt   = keys{gi}(find(keys{gi}=='|',1,'last')+1:end); 
    opts2  = opts; opts2.ignoreScouts = contains(gopt,'-i y'); % options of the import
    Pout   = fullfile(Ptmp,sprintf('img%06d',gi)); 
    runDCM2NIIX(Pdcm, Pout, opts2, 0); 

    % assign the converted images to the scans of this DICOM directory
    Pconv = cat_vol_findfiles(Pout,'*.json',struct('depth',1)); 
    for ci = 1:numel(Pconv)
      name       = spm_file(Pconv{ci},'basename'); 
      fnameparts = strsplit(name,'='); 
      if ~exist(fullfile(Pout,[name niiext]),'file') || numel(fnameparts) < 4 || ...
          any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        continue
      end
      V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pconv{ci})); 
      V.ProtocolName = fnameparts{2}; 
      Pdbdirpath = fullfile(Pmdbdir,cat_io_dcm2bids_db('getDBdir',V)); 
      Pjsondb    = fullfile(Pdbdirpath,[name '.json']); 
      if exist(Pjsondb,'file')
        if exist(fullfile(Pdbdirpath,[name niiext]),'file'), continue; end % image exists
        Pimp = spm_file(Pjsondb,'prefix','catDCM2BIDSimport_'); 
        if ~exist(Pimp,'file'), continue; end
        I = cat_io_dcm2bids_db('readImportJSON',Pimp); 
        if ~strcmp(I.DICOMdirectory,Pdcm), continue; end % scan of another directory
      else
        % further image of a dcm2niix name conflict (e.g. "7a" of "7") of a 
        % scan of this directory
        Pbase = fullfile(Pdbdirpath,[regexprep(name,'[a-z]$','') '.json']); 
        if strcmp(Pbase,Pjsondb) || ~exist(Pbase,'file'), continue; end
        I = cat_io_dcm2bids_db('readImportJSON',spm_file(Pbase,'prefix','catDCM2BIDSimport_')); 
        if ~strcmp(I.DICOMdirectory,Pdcm), continue; end
        I.ScanID = cat_io_dcm2bids_db('getScanID',V,fnameparts); 
        cat_io_dcm2bids_db('writeImportJSON',spm_file(Pjsondb,'prefix','catDCM2BIDSimport_'), I); 
      end
      % image and the sidecar of the full conversion 
      for ext = {niiext,'.bval','.bvec'}
        if exist(fullfile(Pout,[name ext{1}]),'file')
          movefile(fullfile(Pout,[name ext{1}]), fullfile(Pdbdirpath,[name ext{1}])); 
        end
      end
      copyfile(Pconv{ci}, Pjsondb); 
    end
    for ji = 1:numel(grp(keys{gi}))
      Pj = grp(keys{gi}); 
      if ~exist(spm_file(Pj{ji},'ext',niiext),'file'), msg(Pj{ji}) = 'no image after DICOM conversion'; end
    end
    rmdir(Pout,'s'); 
  end
  fprintf('\n'); 
end
% =========================================================================
function copyImage(Psrc, Pdst)
%copyImage. Copy a NIfTI image with the gzip format of the destination name.
  srcgz = endsWith(Psrc,'.gz'); dstgz = endsWith(Pdst,'.gz'); 
  if srcgz == dstgz
    copyfile(Psrc, Pdst); 
  elseif dstgz
    copyfile(Psrc, Pdst(1:end-3)); gzip(Pdst(1:end-3)); delete(Pdst(1:end-3)); 
  else
    copyfile(Psrc, [Pdst '.gz']); gunzip([Pdst '.gz']); delete([Pdst '.gz']); 
  end
end
% =========================================================================
function runDCM2NIIX(Pdcm, Pout, opts, headeronly)
%runDCM2NIIX. Convert the DICOM files of the directory Pdcm (without its 
%  sub-directories, as all sub-directories are processed separately) into 
%  the directory Pout with the file names "folder=protocol=time=series". 
%  With headeronly, only the JSON sidecars are written (dcm2niix -b o, fast), 
%  otherwise also the images (and bval/bvec) in the gzip format of 
%  opts.gzipi. With opts.ignoreScouts, dcm2niix also ignores derived, 
%  localizer, and 2D images (-i y). 
  if ~exist(Pout,'dir'), mkdir(Pout); end
  if headeronly, b = 'o'; else, b = 'y'; end
  if opts.gzipi, gz = 'y'; else, gz = 'n'; end
  cmd = sprintf('"%s" %s -z %s -b %s -o "%s" "%s"', ...
    opts.Pdcm2nii, dcm2niixOptions(opts), gz, b, Pout, Pdcm); 
  [status,cmdout] = system(cmd); %#ok<ASGLU>

  % assure the gzip status of the images
  if ~headeronly
    if opts.gzipi
      P = cat_vol_findfiles( Pout , '*.nii' ,struct('depth',1)); 
      for fi=1:numel(P), gzip(P{fi}); delete(P{fi}); end
    else
      P = cat_vol_findfiles( Pout , '*.nii.gz' ,struct('depth',1)); 
      for fi=1:numel(P), gunzip(P{fi}); delete(P{fi}); end
    end
  end
end
% =========================================================================
function str = dcm2niixOptions(opts)
%dcm2niixOptions. Common dcm2niix options of the header and image conversion
%  (file names, no Philips precise scaling, no BIDS anonymizing, no 
%  sub-directories, and with opts.ignoreScouts no derived/localizer/2D images).
  if opts.ignoreScouts, ign = 'y'; else, ign = 'n'; end
  str = sprintf('-f "%%f=%%p=%%t=%%s" -p n -ba n -d 0 -i %s', ign); 
end
% =========================================================================
function P = getDCM2NIIX
%getDCM2NIIX. Path of the dcm2niix executable at its default installation
%  (MRIcroGL on macOS, MRIcron on Windows, /usr or /opt on Linux), with an
%  error if it is not available. It is only called if DICOM data has to be
%  read or converted and no path is given by opts.Pdcm2nii.
  if ismac
    P = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix'; 
  elseif ispc
    P = 'C:\Program Files\MRIcron\dcm2niix.exe';
    if ~exist(P,'file')
      P = 'C:\Program Files (x86)\MRIcron\dcm2niix.exe';
    end
  else
    [~,P] = system('find /usr /opt -type f -name "dcm2niix" 2>/dev/null');
  end
  P = deblank(P); 
  if ~exist(P,'file')
    %%%%%%% maybe use input path or include it in CAT?
    error('cat_io_dcm2bids:noDcm2niix', ...
      'Cannot find dcm2niix. Please install from: \n  %s\n\n', ... 
      spm_file('https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage' , ...
        'link','https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage')); 
  end
end
