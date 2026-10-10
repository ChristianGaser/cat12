function varargout = cat_io_dcm2bids_import(action,varargin)
%cat_io_dcm2bids_import. Input sources and import stage of cat_io_dcm2bids.
%  The input directories are classified as DICOM directories, BIDS datasets,
%  or directories of the database of the output before any output is 
%  written (getInputSources). Then, the DICOM headers of new scans and the 
%  scans of BIDS datasets are imported into the database and the images are
%  converted (DICOM) or copied (BIDS) if required (importSources, see 
%  cat_io_dcm2bids_importDCM and cat_io_dcm2bids_importBIDS).
%
%  varargout = cat_io_dcm2bids_import(action,varargin)
%
%  Actions (see the help of the local functions):
%    Psrc = cat_io_dcm2bids_import('getInputSources',Pdata,Pmdbdir,ignoreScouts)
%    [Psnr,convmsg] = cat_io_dcm2bids_import('importSources',Psrc,Pmdbdir,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'getInputSources','importSources'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_import:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Psrc = getInputSources(Pdata, Pmdbdir, ignoreScouts)
%getInputSources. Classify the input directories. 
%  (1) DICOM directories: the directory and all its sub-directories (each 
%      without its sub-directories, see runDCM2NIIX in 
%      cat_io_dcm2bids_importDCM) are imported into the database 
%      (Psrc.dicom), where directories of a database are ignored.
%  (2) Directories of the database of this output (Pmdbdir) to export 
%      already imported data, e.g. with other protocols or if the DICOM 
%      data is no longer available: the database itself, subject (sub-*), 
%      session (ses-*), or scan (snr-*) directories give all scan 
%      directories below (Psrc.database). 
%  (3) BIDS datasets or directories within a BIDS dataset (Psrc.bids with 
%      the dataset root and the selection, see cat_io_dcm2bids_importBIDS),
%      where derivatives are ignored and the sourcedata directory is 
%      imported as DICOM directory. BIDS datasets within DICOM directories 
%      are imported as BIDS datasets. 
%  Directories of another database give an error, as merging databases 
%  is not supported. The lists are sorted and without repetitions. 

  fs   = regexptranslate('escape',filesep); 
  Psrc = struct('dicom',{{}},'database',{{}},'bids',struct('root',{},'sel',{})); 
  if ignoreScouts
    cat_io_cprintf('blue','\n  Skip all localizer and scouts!\n\n'); 
  end
  for di = 1:numel(Pdata)
    Pdi = regexprep(Pdata{di},['[' fs ']+$'],''); % no final filesep
    
    % database directory: last path component "catDCM2BIDSdb" and the rest
    tok = regexp(cat_io_dcm2bids_helper('absPath',Pdi), ['^(.*' fs 'catDCM2BIDSdb)(' fs '.*)?$'], 'tokens', 'once'); 

    % (3) BIDS dataset or a directory within a BIDS dataset (with the 
    % dataset_description.json), where the derivatives are ignored and the
    % sourcedata directory is imported as DICOM directory
    if isempty(tok)
      Pdi = cat_io_dcm2bids_helper('normPath',Pdi); 
      Pbr = bidsRoot(Pdi); 
      if ~isempty(Pbr)
        rel = regexp(Pdi(numel(Pbr)+1:end),'[^/\\]+','match','once'); 
        if strcmp(rel,'derivatives')
          cat_io_cprintf('warn','  Derivatives of BIDS datasets are not imported: %s\n',Pdi); 
          continue
        elseif ~strcmp(rel,'sourcedata')
          Psrc.bids(end+1) = struct('root',Pbr,'sel',Pdi); 
          if isempty(rel) && exist(fullfile(Pbr,'sourcedata'),'dir')
            cat_io_cprintf('blue',['  The BIDS dataset "%s" has a sourcedata directory that may contain the DICOM ' ...
              'data, which can be selected as separate input (but gives other subject IDs).\n'], spm_file(Pbr,'basename')); 
          end
          continue
        end
      end
    end

    % (1) DICOM directory: the directory and all its sub-directories (with 
    % absolute paths as they are stored in the database)
    if isempty(tok)
      Pdirsdi = [{Pdi}; cat_vol_findfiles(Pdi,'*',struct('dirs',1))]; 
      Pdirsdi(cat_io_contains(Pdirsdi,'catDCM2BIDSdb')) = []; 
      % BIDS datasets within the directory (imported as BIDS, without their
      % sub-directories, i.e. also without the sourcedata)
      % (datasets within other datasets, e.g. derivatives, are not imported)
      Pdirsdi = sort(Pdirsdi); % parent directories first
      isbids = cellfun(@(P) exist(fullfile(P,'dataset_description.json'),'file') > 0, Pdirsdi); 
      for bri = find(isbids)'
        if isempty(Pdirsdi{bri}), continue; end % within another dataset
        Psrc.bids(end+1) = struct('root',Pdirsdi{bri},'sel',Pdirsdi{bri}); 
        Pdirsdi( strncmp(Pdirsdi,[Pdirsdi{bri} filesep],numel(Pdirsdi{bri})+1) ) = {''}; 
      end
      Pdirsdi(isbids | cellfun('isempty',Pdirsdi)) = []; 
      if ignoreScouts
        Pdirsdi(cat_io_contains(lower(Pdirsdi),{'localizer','scout'})) = []; 
      end
      Psrc.dicom = [Psrc.dicom; Pdirsdi]; 
      continue
    end

    % (2) database directory (only of this output)
    if ~strcmp( cat_io_dcm2bids_helper('canonPath',tok{1}) , cat_io_dcm2bids_helper('canonPath',Pmdbdir) )
      error('cat_io_dcm2bids:foreignDatabase', ...
        ['The selected directory "%s" is not part of the database of the output directory ' ...
         '"%s". Select directories of this database or change the output directory.'], Pdi, Pmdbdir); 
    end
    if isempty(tok{2}), rel = cell(1,0); else, rel = strsplit(tok{2}(2:end),filesep); end
    levels = {'sub-','ses-','snr-'}; 
    if numel(rel) > 3 || ~all(cellfun(@(r,l) strncmp(r,l,4), rel, levels(1:numel(rel))))
      error('cat_io_dcm2bids:noDatabaseDirectory', ...
        ['The selected directory "%s" is not a subject (sub-*), session (ses-*), or scan ' ...
         '(snr-*) directory of the database.'], Pdi); 
    end
    Pdir = fullfile(Pmdbdir,rel{:}); 
    if numel(rel) == 3
      Pdirsdi = {Pdir}; 
    else
      Pdirsdi = cat_vol_findfiles(Pdir,'snr-*',struct('dirs',1,'depth',3 - numel(rel))); 
    end
    Psrc.database = [Psrc.database; Pdirsdi]; 
  end
  Psrc.dicom    = unique(Psrc.dicom); 
  Psrc.database = unique(Psrc.database); 
  if ~isempty(Psrc.bids)
    [~,ui] = unique(strcat({Psrc.bids.root},'|',{Psrc.bids.sel})); Psrc.bids = Psrc.bids(sort(ui)); 
  end
  if ~isempty(Psrc.database)
    cat_io_cprintf('blue','  Export %d database scan directories.\n\n', numel(Psrc.database)); 
  end
end
% =========================================================================
function R = bidsRoot(P)
%bidsRoot. Root directory of the BIDS dataset (with dataset_description.json) 
%  that contains the directory P (or empty). 
  R = ''; 
  while ~isempty(P)
    if exist(fullfile(P,'dataset_description.json'),'file'), R = P; return; end
    Pp = fileparts(P); 
    if strcmp(Pp,P), return; end
    P = Pp; 
  end
end
% =========================================================================
function [Psnr,convmsg] = importSources(Psrc,Pmdbdir,opts)
%importSources. Import the input sources Psrc (see getInputSources) into 
%  the database Pmdbdir and convert (or copy) the images if required. 
%  (1) The DICOM headers of new scans (fast, dcm2niix header only) and the 
%      scans of BIDS datasets are imported into the database, see 
%      importDICOMheaders in cat_io_dcm2bids_importDCM and 
%      cat_io_dcm2bids_importBIDS. 
%  (2) The images of these scans are converted (DICOM) or copied (BIDS) if 
%      required (not for the JSON-only outputs 0 and 1), see convertImages 
%      in cat_io_dcm2bids_importDCM. 
%  Psnr are the scan directories of all selected scans (DICOM, BIDS, and 
%  database input) and convmsg is a containers.Map with the reason of scans
%  without image. 

  % temporary directory of this run for the dcm2niix conversions (removed 
  % at the end, also after errors or interruptions; directories of other 
  % runs are only removed if they are older than one day)
  [~,tmpname] = fileparts(tempname); 
  Ptmp     = fullfile(Pmdbdir,['catDCM2BIDStmp_' tmpname]); 
  Dtmp     = dir(fullfile(Pmdbdir,'catDCM2BIDStmp*')); 
  for ti = 1:numel(Dtmp)
    if Dtmp(ti).isdir && datetime(Dtmp(ti).datenum,'ConvertFrom','datenum') < datetime('now') - days(1)
      rmTmpDir(fullfile(Pmdbdir,Dtmp(ti).name)); 
    end
  end
  tmpclean = onCleanup(@() rmTmpDir(Ptmp)); % removes Ptmp at the end

  % (1) import the DICOM headers of new scans and the scans of BIDS datasets
  % into the database and get the scan directories of all selected scans 
  % (DICOM, BIDS, and database input)
  Psnr = unique([cat_io_dcm2bids_importDCM('importDICOMheaders',Psrc.dicom, Pmdbdir, Ptmp, opts); ...
                 cat_io_dcm2bids_importBIDS(Psrc.bids, Pmdbdir, opts); Psrc.database]); 

  % (2) convert the images of these scans if required (not for the JSON-only
  % outputs 0 and 1)
  if opts.output > 1
    convmsg = cat_io_dcm2bids_importDCM('convertImages',Psnr, Pmdbdir, Ptmp, opts); 
  else
    convmsg = containers.Map('KeyType','char','ValueType','char'); 
  end
end
% =========================================================================
function rmTmpDir(Ptmp)
%rmTmpDir. Remove the temporary conversion directory of the database. 
  if exist(Ptmp,'dir'), rmdir(Ptmp,'s'); end
end
