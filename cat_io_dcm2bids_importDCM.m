function varargout = cat_io_dcm2bids_importDCM(action,varargin)
%cat_io_dcm2bids_importDCM. DICOM import of cat_io_dcm2bids (dcm2niix).
%  The DICOM headers of new scans are imported into the database 
%  (importDICOMheaders, dcm2niix header only, i.e. fast), whereas the images
%  are converted only if they are required (convertImages, e.g. not for the
%  JSON-only outputs), i.e. for the outputs with images directly after the
%  header of each new DICOM directory (files read once from the disk), 
%  otherwise by the subject packages (e.g. images of a former JSON-only 
%  import). convertImages also copies the images of imported BIDS
%  datasets (see cat_io_dcm2bids_importBIDS). The dcm2niix executable 
%  (opts.Pdcm2nii) is detected only if DICOM data has to be read or converted
%  (see getDCM2NIIX), i.e., dcm2niix is not required for BIDS inputs or the 
%  export of already converted data.
%
%  varargout = cat_io_dcm2bids_importDCM(action,varargin)
%
%  Actions (see the help of the local functions):
%    Psnr = cat_io_dcm2bids_importDCM('importDICOMheaders',Pdcmdirs,Pmdbdir,Ptmp,opts)
%    msg = cat_io_dcm2bids_importDCM('convertImages',Psnr,Pmdbdir,Ptmp,opts[,Psub,quiet])
%    P = cat_io_dcm2bids_importDCM('getDCM2NIIX')
%
%  See also cat_io_dcm2bids.

  switch action
    case {'importDICOMheaders','convertImages','getDCM2NIIX'}
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
%  again. Each read directory is registered in the directory index of the
%  database (catDCM2BIDSimportdirs.json, also directories without scans,
%  e.g. with only localizers, see dirIndex) with its files (number, size,
%  and the last modification), i.e., DICOM directories that were imported
%  before (with the same options) are not read again (only with opts.rerun
%  or changed files, e.g. an incomplete transfer, with a warning, where the
%  images of existing scans are kept). For outputs with images (opts.output
%  > 1), the new scans of a directory are converted directly after its
%  header (see convertImages), i.e., the files are read once from the disk
%  (cached by the system, see prefetchFiles).
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
  
  % DICOM directories imported before (with the same dcm2niix options):
  % scans of the import information and the directory index
  [Pimpdcm, Pimpsnr, Pimpopt] = cat_io_dcm2bids_db('getImportIndex',Pmdbdir);
  X     = dirIndex('read',Pmdbdir);
  dopt  = dcm2niixOptions(opts);
  ver   = ''; % dcm2niix version (only if required)
  seen  = containers.Map('KeyType','char','ValueType','any'); % scans of this run
  flist = containers.Map('KeyType','char','ValueType','any'); % files of directories
  split = {}; % IDs of split series
  physio = cell(0,1); % physiological recordings (see attachPhysio)
  convert = opts.output > 1; % images of the new scans (see convertImages)
  nread = 0; nchanged = 0;

  spm_progress_bar('init',numel(Pdcmdirs),'DICOM header import','DICOM directories')
  for di = 1:numel(Pdcmdirs)
    spm_progress_bar('Set',di)

    % directories without files (e.g. only with sub-directories)
    F = dirFiles(Pdcmdirs{di},flist);
    if isempty(F.name), continue; end

    % previously imported directory with the same files (directory index or,
    % for databases of former versions, the import information): scans of
    % the database
    imported = strcmp(Pimpdcm, Pdcmdirs{di}) & strcmp(Pimpopt, dopt);
    xi = find(strcmp({X.DICOMdirectory}, Pdcmdirs{di}) & strcmp({X.dcm2niixOptions}, dopt), 1);
    changed = ~isempty(xi) && ~(X(xi).Files == F.Files && X(xi).Bytes == F.Bytes && strcmp(X(xi).Modified,F.Modified));
    if ~opts.rerun && ~changed && (~isempty(xi) || any(imported))
      if isempty(xi), X = dirIndex('add',X,Pdcmdirs{di},dopt,F,nnz(imported)); end
      Psnr = [Psnr; Pimpsnr(imported)]; %#ok<AGROW>
      nold = nold + 1;
      continue
    end
    if changed && ~opts.rerun
      cat_io_cprintf('warn',['  The files of the DICOM directory changed since its import (e.g. an incomplete ' ...
        'transfer), i.e. it is imported again (the images of existing scans are kept, use opts.rerun ' ...
        'to convert them again):\n    %s\n'],Pdcmdirs{di});
      nchanged = nchanged + 1;
    end
    tdi = tic; nnewdi = nnew;

    % header only conversion of the directory
    if isempty(ver)
      if isempty(opts.Pdcm2nii), opts.Pdcm2nii = getDCM2NIIX; end % only if required
      [~,ver] = system(sprintf('"%s" --version',opts.Pdcm2nii)); 
      ver = regexp(ver,'v[0-9]+\.[0-9]+\.[0-9]+','match','once'); 
    end
    Ptmpdi = fullfile(Ptmp,sprintf('hdr%06d',di)); 
    optsh = opts; if convert, optsh.prefetch = inf; end % (whole files for the conversion below)
    runDCM2NIIX(Pdcmdirs{di}, Ptmpdi, optsh, 1);
    Pjson  = cat_vol_findfiles(Ptmpdi,'*.json',struct('depth',1)); 
    rec    = struct('Pimp',{},'S',{}); % import information of this directory
    for ji = 1:numel(Pjson)
      % physiological recordings of Siemens PhysioLog series (dcm2niix
      % "_physio" files without patient and session fields, also of the
      % header-only conversion) are attached to their fMRI scan after all
      % directories (see attachPhysio)
      if endsWith(spm_file(Pjson{ji},'basename'),'_physio')
        Pph = fullfile(Ptmp,'physio',sprintf('%06d',di));
        if ~exist(Pph,'dir'), mkdir(Pph); end
        for ext = {'.json','.tsv.gz'}
          if exist(spm_file(Pjson{ji},'ext',ext{1}),'file'), movefile(spm_file(Pjson{ji},'ext',ext{1}),Pph); end
        end
        physio{end+1,1} = fullfile(Pph,spm_file(Pjson{ji},'filename')); %#ok<AGROW>
        continue
      end
      fnameparts = cat_io_dcm2bids_db('nameParts',spm_file(Pjson{ji},'basename'));
      if numel(fnameparts) < 4 || any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        continue
      end
      V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pjson{ji})); 
      if opts.ignoreScouts && cat_io_dcm2bids_db('ignoredDerived',V), continue; end % (derived images, see ignoredDerived)
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
          cat_io_dcm2bids_db('nameParts',spm_file(P,'basename'))), scanid), Pother); 
        if any(same)
          Pexist = Pother{find(same,1)}; 
          Pimp   = spm_file(Pexist,'prefix','catDCM2BIDSimport_'); 
          if exist(Pimp,'file'), I = cat_io_dcm2bids_db('readImportJSON',Pimp); Pprev = I.DICOMdirectory; end
        end
      end
      if ~isempty(Pexist) && ~isempty(Pprev) && ~strcmp(Pprev,Pdcmdirs{di})
        if exist(Pprev,'dir') && ~isequal(dirFiles(Pprev,flist).name,F.name)
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
    Psnrdi = {};
    for ri = 1:numel(rec)
      if exist(strrep(rec(ri).Pimp,'catDCM2BIDSimport_',''),'file') % not revoked
        cat_io_dcm2bids_db('writeImportJSON',rec(ri).Pimp, rec(ri).S);
        Psnrdi{end+1,1} = fileparts(rec(ri).Pimp); %#ok<AGROW>
      end
    end
    rmdir(Ptmpdi,'s');

    % images of the new scans while the files are cached (no prefetch)
    if convert && nnew > nnewdi
      opts0 = opts; opts0.prefetch = 0;
      convertImages(unique(Psnrdi), Pmdbdir, Ptmp, opts0, {}, 1);
    end
    X = dirIndex('add',X,Pdcmdirs{di},dopt,F,numel(Psnrdi));
    nread = nread + 1;
    if mod(nread,20) == 0, dirIndex('write',Pmdbdir,X); end % (interruptions)
    fprintf('  %4d/%d  %-60s %5d files %8.1f MB  %3d new scans  %6.1f s\n', di, numel(Pdcmdirs), ...
      spm_str_manip(Pdcmdirs{di},'a60'), F.Files, F.Bytes/2^20, nnew - nnewdi, toc(tdi));
  end
  spm_progress_bar('Clear')
  dirIndex('write',Pmdbdir,X);
  if ~isempty(physio), attachPhysio(physio, unique([Pimpsnr; Psnr])); end

  % existing scan directories (without revoked split series)
  Psnr = unique(Psnr);
  Psnr = Psnr( cellfun(@(P) exist(P,'dir') && ~isempty(cat_io_dcm2bids_db('getDBscanJSONs',P)), Psnr) );
  fprintf('  %d new scans, %d directories read (%d changed), %d directories imported before\n\n', ...
    nnew, nread, nchanged, nold);
end
% =========================================================================
function F = dirFiles(P,flist)
%dirFiles. Sorted names of the (non-hidden) files of a directory with their
%  number, size in bytes, and the last modification (cached in the
%  containers.Map flist).
  if isKey(flist,P), F = flist(P); return; end
  D = dir(P); D = D( ~[D.isdir] & ~strncmp({D.name},'.',1) );
  F = struct('name',{sort({D.name})},'Files',numel(D),'Bytes',sum([D.bytes]),'Modified','');
  if ~isempty(D), F.Modified = char(datetime(max([D.datenum]),'ConvertFrom','datenum','Format','yyyy-MM-dd''T''HH:mm:ss')); end
  flist(P) = F; %#ok<NASGU> handle object
end
% =========================================================================
function X = dirIndex(action,varargin)
%dirIndex. Directory index of the DICOM import (catDCM2BIDSimportdirs.json
%  of the database): the read DICOM directories with their dcm2niix options,
%  files (number, size in bytes, and the last modification, see dirFiles),
%  the number of their scans (0 for directories without scans, e.g. only
%  localizers), and the import date.
%    X = dirIndex('read',Pmdbdir)
%    X = dirIndex('add',X,Pdcm,dopt,F,nscans)   (replaces an entry of Pdcm)
%        dirIndex('write',Pmdbdir,X)
  X0 = struct('DICOMdirectory',{},'dcm2niixOptions',{},'Files',{},'Bytes',{},'Modified',{}, ...
    'Scans',{},'ImportDate',{});
  switch action
    case 'read'
      X = X0; P = fullfile(varargin{1},'catDCM2BIDSimportdirs.json');
      if exist(P,'file')
        try
          X = cat_io_dcm2bids_db('readImportJSON',P);
          if ~isempty(X), X = orderfields(reshape(X,1,[]),X0); else, X = X0; end
        catch e
          cat_io_cprintf('warn','  Cannot read the directory index "%s" (%s), i.e. directories are read again.\n',P,e.message);
          X = X0;
        end
      end
    case 'add'
      [X,Pdcm,dopt,F,nscans] = deal(varargin{:});
      X(strcmp({X.DICOMdirectory},Pdcm) & strcmp({X.dcm2niixOptions},dopt)) = [];
      X(end+1) = struct('DICOMdirectory',Pdcm,'dcm2niixOptions',dopt,'Files',F.Files,'Bytes',F.Bytes, ...
        'Modified',F.Modified,'Scans',nscans,'ImportDate',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));
    case 'write'
      [Pmdbdir,X] = deal(varargin{:});
      cat_io_dcm2bids_helper('writeJSON',fullfile(Pmdbdir,'catDCM2BIDSimportdirs.json'),X);
  end
end
% =========================================================================
function revokeScan(Pjsondb)
%revokeScan. Remove a scan that was imported in this run (sidecar, image of
%  the conversion of the import, and import information, and the scan
%  directory if it is empty).
  Pimp = spm_file(Pjsondb,'prefix','catDCM2BIDSimport_');
  for ext = {'.json','.nii','.nii.gz','.bval','.bvec'}
    if exist(spm_file(Pjsondb,'ext',ext{1}),'file'), delete(spm_file(Pjsondb,'ext',ext{1})); end
  end
  if exist(Pimp,'file'), delete(Pimp); end
  D = dir(fileparts(Pjsondb));
  if all([D.isdir]), rmdir(fileparts(Pjsondb)); end
end
% =========================================================================
function attachPhysio(Pphysio, Psnr)
%attachPhysio. Attach the physiological recordings (dcm2niix sidecars
%  "folder=protocol=time=series_recording-<type>_physio.json" with their
%  tsv.gz) to their fMRI scan of the scan directories Psnr, i.e., the scan
%  with the same protocol name and study time (time part of the name) and
%  the closest lower series number (e.g. series 13 of the PhysioLog of the
%  fMRI series 11 with its SBRef 10). The files are stored as associated
%  files of the scan (catDCM2BIDSassoc_<sidecar>_recording-<type>_physio.*)
%  that are exported with the BIDS name of the scan (see processScans in
%  cat_io_dcm2bids_process). Recordings without (unique) scan are skipped.

  % sidecars of the scan directories with their name parts
  Pjson = {};
  for si = 1:numel(Psnr)
    if exist(Psnr{si},'dir'), Pjson = [Pjson; cat_io_dcm2bids_db('getDBscanJSONs',Psnr{si})]; end %#ok<AGROW>
  end
  parts = cellfun(@(P) cat_io_dcm2bids_db('nameParts',spm_file(P,'basename')), Pjson, 'UniformOutput', false);
  ok    = cellfun(@numel,parts) == 4; Pjson = Pjson(ok); parts = parts(ok);
  sn    = cellfun(@(p) str2double(regexp(p{4},'^\d+','match','once')), parts);

  for phi = 1:numel(Pphysio)
    p   = cat_io_dcm2bids_db('nameParts',spm_file(Pphysio{phi},'basename'));
    psn = str2double(regexp(p{end},'^\d+','match','once'));
    rec = regexprep(p{end},'^\d+[a-z]*_',''); % e.g. recording-cardiac_physio
    cand = [];
    if numel(p) == 4 && ~isempty(p{3}) && ~isnan(psn)
      cand = find(cellfun(@(q) strcmp(q{2},p{2}) && strcmp(q{3},p{3}), parts) & sn < psn);
      cand = cand(sn(cand) == max(sn(cand)));
    end
    if isempty(cand) || numel(unique(cellfun(@fileparts,Pjson(cand),'UniformOutput',false))) > 1
      cat_io_cprintf('warn','  Physiological recording "%s" without (unique) fMRI scan - not imported.\n', ...
        spm_file(Pphysio{phi},'basename'));
      continue
    end
    for ci = 1:numel(cand) % (e.g. the echoes of multi-echo fMRI)
      for ext = {'.json','.tsv.gz'}
        if exist(spm_file(Pphysio{phi},'ext',ext{1}),'file')
          copyfile(spm_file(Pphysio{phi},'ext',ext{1}), fullfile(fileparts(Pjson{cand(ci)}), ...
            sprintf('catDCM2BIDSassoc_%s_%s%s',spm_file(Pjson{cand(ci)},'basename'),rec,ext{1})));
        end
      end
    end
  end
end
% =========================================================================
function msg = convertImages(Psnr, Pmdbdir, Ptmp, opts, Psub, quiet)
%convertImages. Convert the missing images of the scans of the scan 
%  directories Psnr. The DICOM directory of the import information 
%  (catDCM2BIDSimport_*.json) of each scan is converted once with its 
%  dcm2niix options. The image (and bval/bvec) and its JSON sidecar of the 
%  full conversion replace the sidecar of the header-only import (they can 
%  differ, e.g. for name conflicts of dcm2niix within a directory, where 
%  the further images "...a" are added as new scans). The warnings of 
%  dcm2niix are added to the sidecar of each image (ConversionWarnings, 
%  see conversionWarnings, an empty list without warnings), also for scans
%  of former versions (checked once by a header-only conversion of their 
%  DICOM directory and a full conversion only if dcm2niix gives warnings),
%  with the merged echoes of incomplete echo images (see mergeEchoes). 
%  Psub limits the assigned images to these subject directories of the 
%  database (e.g. the subject package of a parallel process, as a DICOM 
%  directory can contain scans of several subjects that are converted by 
%  their own process; empty for all). quiet suppresses the progress output
%  (e.g. for the conversion of each directory of the import). 
%  msg is a containers.Map with the reason for scans without image. 
  msg = containers.Map('KeyType','char','ValueType','char'); 
  if ~exist('Psub','var') || isempty(Psub), Psub = {}; else, Psub = cellfun(@(p) cat_io_dcm2bids_helper('normPath',p),cellstr(Psub),'UniformOutput',false); end
  if ~exist('quiet','var'), quiet = 0; end
  if opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end

  % scans without image grouped by DICOM directory and dcm2niix options
  % (and scans with image but without the dcm2niix warnings in chk)
  grp = containers.Map('KeyType','char','ValueType','any'); 
  chk = containers.Map('KeyType','char','ValueType','any'); 

  for si = 1:numel(Psnr)
    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{si}); 
    for ji = 1:numel(Pjson)
      if exist(spm_file(Pjson{ji},'ext',niiext),'file')
        if ~contains(cat_io_dcm2bids_helper('readText',Pjson{ji}),'"ConversionWarnings"')
          Pimp = spm_file(Pjson{ji},'prefix','catDCM2BIDSimport_'); 
          if exist(Pimp,'file')
            I = cat_io_dcm2bids_db('readImportJSON',Pimp); 
            if isfield(I,'DICOMdirectory') && exist(I.DICOMdirectory,'dir') % (not for BIDS imports)
              key = [I.DICOMdirectory '|' I.dcm2niixOptions]; 
              if isKey(chk,key), chk(key) = [chk(key); Pjson(ji)]; else, chk(key) = Pjson(ji); end
            end
          end
        end
        continue
      end
      Pimp = spm_file(Pjson{ji},'prefix','catDCM2BIDSimport_'); 
      if ~exist(Pimp,'file'), msg(Pjson{ji}) = 'no import information'; continue; end
      I = cat_io_dcm2bids_db('readImportJSON',Pimp); 
      if isfield(I,'SourceFile') 
        % BIDS import: copy the image (and bval/bvec)
        if ~exist(I.SourceFile,'file')
          msg(Pjson{ji}) = 'source image not available'; 
          [~,lnk] = system(sprintf('test -L "%s" && echo link',I.SourceFile)); 
          if contains(lnk,'link') % e.g. a DataLad/git-annex file without content
            msg(Pjson{ji}) = 'source image not available (link without file, e.g. DataLad: run "datalad get")'; 
          end
          continue
        end
        copyImage(I.SourceFile, spm_file(Pjson{ji},'ext',niiext)); 
        if isfield(I,'SourceBval') && ~isempty(I.SourceBval), copyWritable(I.SourceBval, spm_file(Pjson{ji},'ext','.bval')); end
        if isfield(I,'SourceBvec') && ~isempty(I.SourceBvec), copyWritable(I.SourceBvec, spm_file(Pjson{ji},'ext','.bvec')); end
        continue
      end
      if ~exist(I.DICOMdirectory,'dir'), msg(Pjson{ji}) = 'DICOM directory not available'; continue; end
      key = [I.DICOMdirectory '|' I.dcm2niixOptions]; 
      if isKey(grp,key), grp(key) = [grp(key); Pjson(ji)]; else, grp(key) = Pjson(ji); end
    end
  end
  if grp.Count == 0 && chk.Count == 0, return; end
  if isempty(opts.Pdcm2nii), opts.Pdcm2nii = getDCM2NIIX; end % only if required
  checkWarnings(chk, Ptmp, opts, niiext); 
  if grp.Count == 0, return; end
  if ~quiet, fprintf('DCM2BIDS - convert images of %d DICOM directories\n', grp.Count); end

  keys = grp.keys; 
  if ~quiet, spm_progress_bar('init',numel(keys),'DICOM image import','DICOM directories'); end
  for gi = 1:numel(keys)
    if ~quiet, spm_progress_bar('Set',gi); end

    Pdcm   = keys{gi}(1:find(keys{gi}=='|',1,'last')-1); 
    gopt   = keys{gi}(find(keys{gi}=='|',1,'last')+1:end); 
    opts2  = opts; opts2.dcm2niixIgnore = contains(gopt,'-i y'); % options of the import
    Pout   = fullfile(Ptmp,sprintf('img%06d',gi)); 
    out    = runDCM2NIIX(Pdcm, Pout, opts2, 0); 
    W      = conversionWarnings(out, Pdcm); 
    Pnew   = {}; % database sidecars of the converted images

    % assign the converted images to the scans of this DICOM directory
    Pconv = cat_vol_findfiles(Pout,'*.json',struct('depth',1)); 
    for ci = 1:numel(Pconv)
      name       = spm_file(Pconv{ci},'basename'); 
      fnameparts = cat_io_dcm2bids_db('nameParts',name); 
      if ~exist(fullfile(Pout,[name niiext]),'file') || numel(fnameparts) < 4 || ...
          any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        continue
      end
      V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pconv{ci})); 
      V.ProtocolName = fnameparts{2}; 
      Pdbdirpath = fullfile(Pmdbdir,cat_io_dcm2bids_db('getDBdir',V)); 
      Pjsondb    = fullfile(Pdbdirpath,[name '.json']); 
      if ~isempty(Psub) && ~any(strcmp(cat_io_dcm2bids_helper('normPath',fileparts(fileparts(Pdbdirpath))),Psub))
        continue % scan of another subject package
      end
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
      if isKey(W,name), w = W(name); else, w = {}; end
      setSidecarField(Pjsondb,'ConversionWarnings',w); 
      Pnew{end+1,1} = Pjsondb; %#ok<AGROW>
    end
    if ~isempty(Pnew), mergeEchoes(Pdcm, Pout, opts2, W, out, Pnew, niiext); end
    for ji = 1:numel(grp(keys{gi}))
      Pj = grp(keys{gi}); 
      if ~exist(spm_file(Pj{ji},'ext',niiext),'file'), msg(Pj{ji}) = 'no image after DICOM conversion'; end
    end
    rmdir(Pout,'s'); 
  end
  if ~quiet, spm_progress_bar('Clear'); fprintf('\n'); end
end
% =========================================================================
function W = conversionWarnings(cmdout, Pdcm)
%conversionWarnings. Warnings of dcm2niix of each converted image (map of 
%  the image name to a cell of warnings) from the console output cmdout of 
%  a full conversion (dcm2niix prints the warnings of an image before its
%  line "Convert ... as <image>", which is missing for header-only
%  conversions). Paths and DICOM file names are replaced (no personal data,
%  e.g. for the BIDS sidecar). 
  W = containers.Map('KeyType','char','ValueType','any'); 
  L = regexp(cmdout,'\r?\n','split'); w = {}; 
  for li = 1:numel(L)
    l = strtrim(L{li}); 
    if startsWith(l,'Warning:')
      l = strrep(l,Pdcm,'<DICOMdirectory>'); 
      l = regexprep(l,'(?<![\w.])(/|[A-Za-z]:[\\/])[^\s''"]*','<path>'); 
      l = regexprep(l,'[^\s''"]+\.(dcm|DCM|ima|IMA)\>','<file>'); 
      if ~any(strcmp(w,l)), w{end+1} = l; end %#ok<AGROW>
    else
      t = regexp(l,'^Convert \d+ DICOM as (.*) \(\d+x','tokens','once'); 
      if ~isempty(t)
        [~,name] = fileparts(t{1}); W(name) = w; w = {}; 
      end
    end
  end
end
% =========================================================================
function geom = geometryWarning(w)
%geometryWarning. True if one of the dcm2niix warnings w describes a wrong
%  geometry of the image, i.e. missing slices or varying slice distances 
%  (e.g. the NKI-RS field maps with one of the two echoes of each slice, 
%  see mergeEchoes).
  geom = any(cat_io_contains(w,{'Interslice distance varies','Missing images'})); 
end
% =========================================================================
function setSidecarField(Pjson, fn, val)
%setSidecarField. Add (or replace) the field fn with the value val in the 
%  JSON sidecar Pjson. A new field is inserted into the text, i.e. the 
%  other fields of the dcm2niix sidecar are kept as they are. 
  txt = cat_io_dcm2bids_helper('readText',Pjson); 
  if contains(txt,['"' fn '"'])
    V = jsondecode(txt); V.(fn) = val; 
    cat_io_dcm2bids_helper('writeJSON',Pjson,V); 
    return
  end
  if iscell(val) && isempty(val), enc = '[]'; else, enc = jsonencode(val); end
  i = find(txt == '}',1,'last'); 
  pre = deblank(txt(1:i-1)); 
  if ~endsWith(pre,'{'), pre = [pre ',']; end
  cat_io_dcm2bids_helper('writeText',Pjson,sprintf('%s\n\t"%s": %s\n}\n',pre,fn,enc)); 
end
% =========================================================================
function mergeEchoes(Pdcm, Pout, opts, W, cmdout, Pjson, niiext)
%mergeEchoes. Image of the merged echoes (dcm2niix -m y) of series whose
%  echo images have missing slices with the complementary slices in the 
%  other echo, e.g. the NKI-RS field maps with one echo of each slice (the
%  echoes with the same instance number got the same DICOM file name at 
%  the export, i.e. one of them is lost). The merged image is complete but
%  with the mixed echo times of the slices, i.e. only for the QC and the 
%  alignment (see cat_io_dcm2bids_qc and scannerAlignment). It is saved as
%    catDCM2BIDSmerged_<scan without the echo suffix>.nii(.gz)
%  in the database scan directory and named in the sidecars of the echo 
%  images (ConversionMergedEchoes), if dcm2niix gives no geometry warning
%  for it (e.g. for duplicated slices of complete echoes). W are the 
%  warnings (see conversionWarnings) and cmdout the console output of the
%  conversion of the DICOM directory Pdcm, and Pjson the database sidecars 
%  of the scans of this directory.
  if ~contains(cmdout,'echo varies') || ~any(cellfun(@(k) geometryWarning(W(k)),W.keys)), return; end
  opts.dcm2niixMerge = true; Pm = [Pout '_merge']; 
  Wm = conversionWarnings(runDCM2NIIX(Pdcm, Pm, opts, 0), Pdcm); 
  names = spm_file(Pjson,'basename'); 
  for k = Wm.keys
    if geometryWarning(Wm(k{1})) || ~exist(fullfile(Pm,[k{1} niiext]),'file'), continue; end
    base = regexprep(k{1},'_e\d+$',''); 
    tgt  = Pjson( strcmp(regexprep(names,'_e\d+$',''),base) & ...
      cellfun(@(n) isKey(W,n) && geometryWarning(W(n)),names) ); 
    if isempty(tgt), continue; end
    Pmerged = fullfile(fileparts(tgt{1}),['catDCM2BIDSmerged_' base niiext]); 
    movefile(fullfile(Pm,[k{1} niiext]),Pmerged); 
    for ti = 1:numel(tgt), setSidecarField(tgt{ti},'ConversionMergedEchoes',spm_file(Pmerged,'filename')); end
    cat_io_cprintf('warn','  %s: echoes with missing slices merged for the QC and alignment (%s)\n', ...
      base,spm_file(Pmerged,'filename')); 
  end
  if exist(Pm,'dir'), rmdir(Pm,'s'); end
end
% =========================================================================
function checkWarnings(chk, Ptmp, opts, niiext)
%checkWarnings. Warnings of dcm2niix for scans with image but without them
%  in their sidecar (ConversionWarnings, e.g. of former versions) grouped
%  in chk by DICOM directory and dcm2niix options (see convertImages): a
%  header-only conversion of the directory (fast), and only if dcm2niix 
%  gives warnings a full conversion into a temporary directory to assign 
%  them to the images (see conversionWarnings) and to merge incomplete 
%  echoes (see mergeEchoes). The images of the database are kept.
  if chk.Count == 0, return; end
  fprintf('DCM2BIDS - check dcm2niix warnings of %d DICOM directories\n', chk.Count);
  keys = chk.keys; nw = 0; 
  for gi = 1:numel(keys)
    Pdcm  = keys{gi}(1:find(keys{gi}=='|',1,'last')-1); 
    gopt  = keys{gi}(find(keys{gi}=='|',1,'last')+1:end); 
    opts2 = opts; opts2.dcm2niixIgnore = contains(gopt,'-i y'); 
    Pout  = fullfile(Ptmp,sprintf('chk%06d',gi)); 
    out   = runDCM2NIIX(Pdcm, Pout, opts2, 1); 
    W     = containers.Map('KeyType','char','ValueType','any'); 
    if contains(out,'Warning:')
      rmdir(Pout,'s'); out = runDCM2NIIX(Pdcm, Pout, opts2, 0); 
      W = conversionWarnings(out, Pdcm); 
    end
    Pjson = chk(keys{gi}); 
    for ji = 1:numel(Pjson)
      name = spm_file(Pjson{ji},'basename'); 
      if isKey(W,name), w = W(name); else, w = {}; end
      setSidecarField(Pjson{ji},'ConversionWarnings',w); nw = nw + ~isempty(w); 
    end
    mergeEchoes(Pdcm, Pout, opts2, W, out, Pjson, niiext); 
    if exist(Pout,'dir'), rmdir(Pout,'s'); end
  end
  fprintf('  %d scan(s) with dcm2niix warnings\n', nw);
end
% =========================================================================
function copyImage(Psrc, Pdst)
%copyImage. Copy a NIfTI image with the gzip format of the destination name
%  (writable, see copyWritable).
  srcgz = endsWith(Psrc,'.gz'); dstgz = endsWith(Pdst,'.gz'); 
  if srcgz == dstgz
    copyWritable(Psrc, Pdst); 
  elseif dstgz
    copyWritable(Psrc, Pdst(1:end-3)); gzip(Pdst(1:end-3)); delete(Pdst(1:end-3)); 
  else
    copyWritable(Psrc, [Pdst '.gz']); gunzip([Pdst '.gz']); delete([Pdst '.gz']); 
  end
end
% =========================================================================
function copyWritable(Psrc, Pdst)
%copyWritable. Copy of a file (the content of a symbolic link) that is
%  writable, also for read-only sources (e.g. the locked files of DataLad
%  datasets), as the database files are updated.
  copyfile(Psrc, Pdst); 
  fileattrib(Pdst,'+w','u'); 
end
% =========================================================================
function cmdout = runDCM2NIIX(Pdcm, Pout, opts, headeronly)
%runDCM2NIIX. Convert the DICOM files of the directory Pdcm (without its 
%  sub-directories, as all sub-directories are processed separately) into 
%  the directory Pout with the file names "folder=protocol=time=series". 
%  With headeronly, only the JSON sidecars are written (dcm2niix -b o, fast), 
%  otherwise also the images (and bval/bvec) in the gzip format of 
%  opts.gzipi. Localizers, scouts, derived images (but not the uniform 
%  image and T1 map of MP2RAGE, see ignoredDerived in cat_io_dcm2bids_db), 
%  and 2D images are ignored by the import (opts.ignoreScouts), i.e. not 
%  by dcm2niix (-i n; imports of former versions with -i y are converted
%  with their options). cmdout is the console output of dcm2niix (e.g. for
%  its warnings, see conversionWarnings). The files are read sequentially
%  before (see prefetchFiles): the whole files for the images, the first
%  MB for the header (about the part that dcm2niix reads), or the bytes of
%  opts.prefetch (internal, e.g. inf for the header import with the image
%  conversion afterwards and 0 for files that are already cached).
  if ~exist(Pout,'dir'), mkdir(Pout); end
  if isfield(opts,'prefetch'), nb = opts.prefetch; elseif headeronly, nb = 2^20; else, nb = inf; end
  prefetchFiles(Pdcm,nb);
  if headeronly, b = 'o'; else, b = 'y'; end
  if opts.gzipi, gz = 'y'; else, gz = 'n'; end
  cmd = sprintf('"%s" %s -z %s -b %s -o "%s" "%s"', ...
    opts.Pdcm2nii, dcm2niixOptions(opts), gz, b, Pout, Pdcm); 
  [~,cmdout] = system(cmd);

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
function prefetchFiles(Pdcm,nb)
%prefetchFiles. Read the files of the DICOM directory Pdcm (without its
%  sub-directories, the first nb bytes of each file, all for inf, none for
%  0) once sequentially, i.e. into the file cache of the
%  operating system, before dcm2niix reads them. The import is limited by
%  the disk and not by the CPU (dcm2niix needs about 2 s for 4000 cached 
%  files), and dcm2niix reads the files in a seek-intensive order, i.e. the 
%  sequential reading before is faster on hard disks (e.g. external USB 
%  disks: header import about 1.8 times and image conversion about 1.3 
%  times faster), whereas it costs little on SSDs. Parallel processes 
%  are slower on hard disks (concurrent seeks).
  if nb <= 0, return; end
  D = dir(Pdcm); D = D(~[D.isdir]);
  for fi = 1:numel(D)
    fid = fopen(fullfile(Pdcm,D(fi).name),'r');
    if fid < 0, continue; end
    n = 0;
    while ~feof(fid) && n < nb, n = n + numel(fread(fid,min(2^24,nb - n),'*uint8')); end
    fclose(fid);
  end
end
% =========================================================================
function str = dcm2niixOptions(opts)
%dcm2niixOptions. Common dcm2niix options of the header and image conversion
%  (file names, no Philips precise scaling, no BIDS anonymizing, no 
%  sub-directories, and all images, as the import selects the scans, e.g.
%  the MP2RAGE uniform image and T1 map that are derived images; -i y only
%  for the conversion of imports of former versions with these options).
  if isfield(opts,'dcm2niixIgnore') && opts.dcm2niixIgnore, ign = 'y'; else, ign = 'n'; end
  str = sprintf('-f "%%f=%%p=%%t=%%s" -p n -ba n -d 0 -i %s', ign); 
  if isfield(opts,'dcm2niixMerge') && opts.dcm2niixMerge, str = [str ' -m y']; end % (only for mergeEchoes)
end
% =========================================================================
function P = getDCM2NIIX
%getDCM2NIIX. Path of the dcm2niix executable, with an error if it is not
%  available. It is only called if DICOM data has to be read or converted
%  and no path is given by opts.Pdcm2nii: 
%   (1) the release of dcm2niix tested with DCM2BIDS (DCM2BIDS/dcm2niix.json,
%       updated by the GitHub workflow .github/workflows/dcm2niix_update.yml)
%       in the cache DCM2BIDS/dcm2niix/<version>/<os> (otherwise in the 
%       MATLAB prefdir, if CAT is not writable), downloaded from GitHub on
%       first use and verified by its SHA-256 checksum (see downloadDCM2NIIX)
%   (2) a manually copied binary in DCM2BIDS/dcm2niix/dcm2niix.<os> (e.g.
%       for systems without internet access)
%   (3) the default installations (MRIcroGL on macOS, MRIcron on Windows, 
%       /usr or /opt on Linux), with a warning if their version is not the
%       tested one
  persistent Pcache % (one download test per MATLAB session)
  if ~isempty(Pcache) && exist(Pcache,'file'), P = Pcache; return; end

  Pcat = fileparts(mfilename('fullpath')); 
  if ismac,     os = 'mac'; exe = 'dcm2niix'; 
  elseif ispc,  os = 'win'; exe = 'dcm2niix.exe'; 
  else,         os = 'lnx'; exe = 'dcm2niix';     % (x86-64 Linux)
  end
  R = struct('version',''); 
  try R = jsondecode(fileread(fullfile(Pcat,'DCM2BIDS','dcm2niix.json'))); end %#ok<TRYNC>

  % (1) tested release (cache or download)
  if ~isempty(R.version)
    Pdirs = {fullfile(Pcat,'DCM2BIDS','dcm2niix',R.version,os), fullfile(prefdir,'CAT_dcm2niix',R.version,os)}; 
    for di = 1:numel(Pdirs)
      P = fullfile(Pdirs{di},exe); 
      if exist(P,'file'), Pcache = P; return; end
    end
    for di = 1:numel(Pdirs)
      try
        P = downloadDCM2NIIX(R,os,exe,Pdirs{di}); Pcache = P; return
      catch e
        msg = e.message; 
        if ~contains(e.identifier,'notWritable'), break; end % (otherwise try the next directory)
      end
    end
    cat_io_cprintf('warn',['  dcm2niix %s not available (%s), use another installation \n' ...
      '  (or copy the dcm2niix of %s into DCM2BIDS/dcm2niix/%s/%s).\n'],R.version,msg, ...
      sprintf('https://github.com/%s/releases/tag/%s',R.repository,R.version),R.version,os); 
  end

  % (2) manually copied binary
  ext = struct('mac','maca64','win','w32','lnx','a64'); 
  P = fullfile(Pcat,'DCM2BIDS','dcm2niix',['dcm2niix.' ext.(os)],exe); 

  % (3) default installations
  if ~exist(P,'file')
    if ismac
      P = '/Applications/MRIcroGL.app/Contents/Resources/dcm2niix'; 
    elseif ispc
      P = 'C:\Program Files\MRIcron\dcm2niix.exe';
      if ~exist(P,'file')
        P = 'C:\Program Files (x86)\MRIcron\dcm2niix.exe';
      end
      if ~exist(P,'file')
        P = 'C:\Program Files (x86)\RIcroGL\Resources\dcm2niix.exe';
      end
    else
      [~,P] = system('find /usr /opt -type f -name "dcm2niix" 2>/dev/null');
      P = strtok(P,newline); % (first match)
    end
  end

  P = deblank(P); 
  if ~exist(P,'file')
    error('cat_io_dcm2bids:noDcm2niix', ...
     ['Cannot find dcm2niix (%s). \n ' ...
      'Please install from: \n  %s\n\n'], ... 
      P, spm_file('https://github.com/rordenlab/dcm2niix/releases' , ...
           'link','https://github.com/rordenlab/dcm2niix/releases')); 
  end
  if ~isempty(R.version)
    [~,ver] = system(sprintf('"%s" --version',P)); ver = regexp(ver,'v[0-9]+\.[0-9]+\.[0-9]+','match','once'); 
    if ~strcmp(ver,R.version)
      cat_io_cprintf('warn','  Use dcm2niix %s (%s), not the version %s tested with DCM2BIDS.\n',ver,P,R.version); 
    end
  end
  Pcache = P; 
end
% =========================================================================
function P = downloadDCM2NIIX(R,os,exe,Pdir)
%downloadDCM2NIIX. Download the dcm2niix release R (see DCM2BIDS/dcm2niix.json)
%  of the system os (mac, lnx, win) from GitHub into the directory Pdir 
%  (see webDownload), with the test of the SHA-256 checksum of the zip file
%  (see sha256File) and of the version of the executable exe. 
  A   = R.assets.(os); 
  url = sprintf('https://github.com/%s/releases/download/%s/%s',R.repository,R.version,A.name); 
  try
    if ~exist(Pdir,'dir'), mkdir(Pdir); end
    Ptest = fullfile(Pdir,'.writetest'); fid = fopen(Ptest,'w'); fclose(fid); delete(Ptest); 
  catch
    error('cat_io_dcm2bids:notWritable','directory "%s" is not writable',Pdir); 
  end
  fprintf('DCM2BIDS - download dcm2niix %s from GitHub (%s)\n',R.version,url); 
  Pzip = fullfile(Pdir,A.name); 
  try
    webDownload(Pzip,url); 
    h = sha256File(Pzip); 
    if ~strcmpi(h,A.sha256)
      error('cat_io_dcm2bids:checksum','wrong SHA-256 checksum of %s (%s instead of %s)',A.name,h,A.sha256); 
    end
    unzip(Pzip,Pdir); delete(Pzip); 
    P = fullfile(Pdir,exe); 
    if ~exist(P,'file'), error('cat_io_dcm2bids:noExe','no %s in %s',exe,A.name); end
    if ~ispc, system(sprintf('chmod +x "%s"',P)); end
    [~,ver] = system(sprintf('"%s" --version',P)); 
    if ~contains(ver,R.version)
      error('cat_io_dcm2bids:version','the downloaded dcm2niix reports another version (%s)',strtrim(ver)); 
    end
  catch e
    if exist(Pzip,'file'), delete(Pzip); end
    if exist(Pdir,'dir'), rmdir(Pdir,'s'); end % (no incomplete installation)
    D = dir(fileparts(Pdir)); if all(ismember({D.name},{'.','..'})), rmdir(fileparts(Pdir)); end % (empty version directory)
    rethrow(e); 
  end
end
% =========================================================================
function webDownload(P,url)
%webDownload. Download of url into the file P by MATLAB, with the system 
%  certificates or by curl of the system, if the certificates of MATLAB do
%  not verify the server (e.g. new root certificates of the GitHub release
%  server or a proxy/VPN with its own certificates in the system store). 
%  The download is tested by the SHA-256 checksum (see downloadDCM2NIIX).
  if ismac
    Pcrt = {'/etc/ssl/cert.pem'}; 
  elseif ispc
    Pcrt = {}; 
  else
    Pcrt = {'/etc/ssl/certs/ca-certificates.crt','/etc/pki/tls/certs/ca-bundle.crt','/etc/ssl/cert.pem'}; 
  end
  Pcrt = Pcrt(cellfun(@(f) exist(f,'file') > 0,Pcrt)); 
  msg  = {}; 

  % MATLAB with its own and with the system certificates
  for ci = 0:numel(Pcrt)
    try
      if ci == 0, wo = weboptions('Timeout',60); else, wo = weboptions('Timeout',60,'CertificateFilename',Pcrt{ci}); end
      websave(P,url,wo); return
    catch e
      msg{end+1} = sprintf('MATLAB: %s',strtok(e.message,newline)); %#ok<AGROW>
      if exist(P,'file'), delete(P); end
    end
  end

  % curl of the system (macOS, Linux, Windows 10+ with the Windows certificates)
  [st,out] = system(sprintf('curl -fsSL --retry 2 --connect-timeout 60 -o "%s" "%s"',P,url)); 
  if st == 0 && exist(P,'file'), return; end
  if exist(P,'file'), delete(P); end
  msg{end+1} = sprintf('curl: %s',strtrim(out)); 
  error('cat_io_dcm2bids:download','download failed (%s)', ...
    strjoin(regexprep(msg,'"https://[^"]*"','(URL)'),'; ')); 
end
% =========================================================================
function h = sha256File(P)
%sha256File. SHA-256 checksum (lowercase hex) of the file P (Java, otherwise
%  shasum on macOS/Linux or certutil on Windows).
  try
    fid = fopen(P,'r'); b = fread(fid,inf,'*uint8'); fclose(fid); 
    md = java.security.MessageDigest.getInstance('SHA-256'); 
    md.update(typecast(b,'int8')); 
    h = lower(reshape(dec2hex(typecast(md.digest,'uint8'),2)',1,[])); 
  catch
    if ispc
      [~,out] = system(sprintf('certutil -hashfile "%s" SHA256',P)); 
    else
      [~,out] = system(sprintf('shasum -a 256 "%s"',P)); 
    end
    h = lower(regexp(out,'[0-9a-fA-F]{64}','match','once')); 
  end
end
