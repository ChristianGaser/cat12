function varargout = cat_io_dcm2bids_process(action,varargin)
%cat_io_dcm2bids_process. Processing of the scans of cat_io_dcm2bids.
%  The scans of the scan directories of the database (sorted, i.e. the 
%  sessions one after another) are tested by the protocols, named and 
%  exported in BIDS, rated by the QC, and each session is preprocessed at
%  its end (segmentation of one anatomical scan, registration, coverage, 
%  DTI and rs-fMRI measures), and the sessions without anatomical scan at
%  the end with the T1w of another session of the subject (processScans).
%  The scans are processed in subject packages, i.e. all scans of a 
%  subject (PatientID) in one call, also by parallel processes (see 
%  cat_io_dcm2bids_parallel). Outputs for all subjects (tables of the 
%  database, participants, DICOM export, reports) are written by 
%  cat_io_dcm2bids with the results R of all packages. 
%
%  varargout = cat_io_dcm2bids_process(action,varargin)
%
%  Actions (see the help of the local functions):
%    R = cat_io_dcm2bids_process('processScans',C,Psnr,convmsg)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_parallel.

  switch action
    case {'processScans'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_process:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function R = processScans(C,Psnr,convmsg)
%processScans. Process the scans of the scan directories Psnr (sorted, all
%  scans of the subjects, e.g. a subject package of cat_io_dcm2bids_parallel)
%  with the context C of the run (see cat_io_dcm2bids: job, directories, 
%  sites, protocols, and the anonymization registry) and the conversion 
%  messages convmsg (see convertSources in cat_io_dcm2bids_import). 
%  R contains the records of the scans for the tables of the database (with
%  the deferred protocol variants varg and subject table rows dbV), the 
%  fields of the DICOM export, the participant data, and the references of
%  BIDS inputs (see cat_io_dcm2bids). 
  job = C.job; Poutdir = C.Poutdir; Pmdbdir = C.Pmdbdir; 
  Pprodictdirs = C.Pprodictdirs; sites = C.sites; protocols = C.protocols; 
  anon = C.anon; tol = job.opts.tolerance; 

  %% basic initialization that might have to be extended
  Vjson = cell(1,numel(Psnr)); sci = 0;
  sub = Vjson; ses = Vjson; datatype = Vjson; pro = Vjson; acq = Vjson; acqtime = Vjson; 
  sn = Vjson; task = Vjson; suffix = Vjson; site = Vjson; scankey = Vjson; scanid = Vjson;
  ent = Vjson; lab = Vjson; acqsn = Vjson; acqnr = Vjson; rec = Vjson; % records of the tables of the database
  Panon = Vjson; Pdbdir = Vjson; Pdbjson = Vjson; pstatus = Vjson;
  Pdbdirpath = Vjson; BIDSpathd = Vjson; Pdbnii = Vjson;
  BIDSpath = Vjson; BIDSdir = Vjson; BIDSfile = Vjson; BIDSsub = Vjson; Pprot = Vjson;
  PID = ''; sni = 0; stime = datetime('now'); bidsroots = {}; convver = {}; 
  srcmap   = containers.Map('KeyType','char','ValueType','any'); % new names of BIDS inputs
  intended = {}; % exported JSONs with IntendedFor references
  partdata = containers.Map('KeyType','char','ValueType','any'); % participant data
  anatrun  = containers.Map('KeyType','char','ValueType','any'); % anatomical scans of the sessions of each subject (this run)
  sesnoanat = {}; % sessions without anatomical scan (DTI and rs-fMRI at the end of the run)
  QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]); 
  varg = Vjson; dbV = {}; % deferred protocol variants and subject table rows (see scanRecord)
  if job.opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end
  for fdiri = 1:numel(Psnr)
    
    % (3) process the scans of each scan directory of the database
    % =====================================================================
    Pdirfi = Psnr{fdiri}; 
    Pjson  = cat_io_dcm2bids_db('getDBscanJSONs',Pdirfi); 
    Pnii   = spm_file(Pjson,'ext',niiext); 
    for fscni = 1:numel(Pnii) % for every scan
      sci = sci + 1;

      % get nii2dcm information 
      Vjson{sci} = cat_io_json(Pjson{fscni});
      if isfield( Vjson{sci} , 'ImageTypeText')
        if isfield( Vjson{sci} , 'ImageType')
          Vjson{sci}.ImageType = unique([reshape(cellstr(Vjson{sci}.ImageType),[],1); ...
                                         reshape(cellstr(Vjson{sci}.ImageTypeText),[],1)]); 
        else
          Vjson{sci}.ImageType = Vjson{sci}.ImageTypeText; 
        end
      end
      Vjson{sci} = cat_io_dcm2bids_db('assurePatientDCMfields',Vjson{sci});
      fnameparts = cat_io_dcm2bids_db('nameParts',spm_file(Pjson{fscni},'basename')); 
      Vjson{sci}.ProtocolName = fnameparts{2};

      % import information, e.g. the datatype, suffix, and participant data 
      % of imported BIDS datasets 
      Isrc = struct(); 
      if exist(spm_file(Pjson{fscni},'prefix','catDCM2BIDSimport_'),'file')
        try %#ok<TRYNC>
          Isrc = cat_io_dcm2bids_db('readImportJSON',spm_file(Pjson{fscni},'prefix','catDCM2BIDSimport_')); 
        end
      end
      isbids = isfield(Isrc,'SourceType') && strcmp(Isrc.SourceType,'BIDS'); 
      if isbids, bidsroots{end+1} = Isrc.SourceDataset; end %#ok<AGROW>
      if isfield(Isrc,'dcm2niix') && ~isempty(Isrc.dcm2niix), convver{end+1} = Isrc.dcm2niix; end %#ok<AGROW>


      % create table header (command line output only)
      % ===================================================================
      if ~strcmp(PID, Vjson{sci}.PatientID)
        sni = sni + 1;

        if sci>1
          fprintf('%s\n',repmat('-',1,154)); 
          fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
          stime = datetime('now'); 
        end

        sname = sprintf('sub-%s-%s',Vjson{sci}.DeviceSerialNumber,strrep(Vjson{sci}.PatientID,'_',''));
        if numel(sname) > 39, sname = ['..' sname(end-36:end)]; end % keep the end of long names
        cat_io_cprintf([0 0.2 .8],'\n%-45s',sprintf('%4d) %s', sni, sname));
        if job.opts.output
          if isempty(Pprodictdirs) || isempty(Pprodictdirs{1})
            fprintf('%15s : %-50s%8s |%5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              ' ',' ','BSM','WSM','ISR','NSR','RES','SQR'); 
          else
            fprintf('%15s : %-50s%8s |%5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              'matching-protocol','PFR','BSM','WSM','ISR','NSR','RES','SQR'); 
          end
        else
          fprintf('%15s : %-50s%8s |\n','dicom-protocol','matching-protocol','PFR'); 
        end
        fprintf('%s\n',repmat('-',1,154)); 
        PID = Vjson{sci}.PatientID; 
      end



      % ignore localizer and scout scans
      % ===================================================================
      fname1 = sprintf('%3d) %s', Vjson{sci}.SeriesNumber, spm_str_manip(strrep(sprintf('%s_%s_%s', ...
        strrep(Vjson{sci}.PatientID, strrep(char(Vjson{sci}.ScanDate),'-',''),''), ...
        strrep(char(Vjson{sci}.ScanDate),'-',''), Vjson{sci}.ProtocolName),'__','_'),'l55'));
      if any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        cat_io_cprintf([.5 .5 .5],'%60s : ignore localizer/scout\n',fname1); 
        sci = sci - 1; % reset import counter
        continue; 
      end
      % image of the scan (required for the outputs 2 and 3, see convertImages in cat_io_dcm2bids_importDCM)
      if job.opts.output > 1 && ~exist(Pnii{fscni},'file')
        if isKey(convmsg,Pjson{fscni}), msg = convmsg(Pjson{fscni}); else, msg = 'not converted'; end
        cat_io_cprintf('err','%60s : no image (%s)\n',fname1,msg); 
        sci = sci - 1; % reset import counter
        continue; 
      end
      % derived images (but not the uniform image and T1 map of MP2RAGE, 
      % see ignoredDerived in cat_io_dcm2bids_db), e.g. of imports of 
      % former versions or of BIDS inputs
      if job.opts.ignoreScouts && cat_io_dcm2bids_db('ignoredDerived',Vjson{sci})
        cat_io_cprintf([.5 .5 .5],'%60s : ignore derived image\n',fname1); 
        sci = sci - 1; % reset import counter
        continue; 
      end
      % ignore 2D data
      if exist(Pnii{fscni},'file')
        Vsz = dir(Pnii{fscni});
        if isempty(Vsz) || ~isfield(Vsz,'bytes')
          cat_io_cprintf([.5 .5 .5],'%60s : ignore 2D data\n',fname1); 
          sci = sci - 1; % reset import counter
          continue; 
        end
        if Vsz.bytes/1024 < 1000 
          try
            evalc('V = spm_vol(Pnii{fscni});'); 
          catch
            cat_io_cprintf('err','%60s : unreadable image\n',fname1); 
            sci = sci - 1; % reset import counter
            continue; 
          end
          if numel(V(1).dim)>2 && any(V(1).dim < 5) 
            cat_io_cprintf([.5 .5 .5],'%60s : ignore 2D data\n',fname1); 
            sci = sci - 1; % reset import counter
            continue; 
          end
        end
        % image dimensions [x y z volumes] and voxel size [x y z] (header only)
        % for the protocol tables and test and the BIDS sidecar
        [vx,dim] = cat_io_dcm2bids_helper('niiHeaderInfo',Pnii{fscni});
        Vjson{sci}.Dimensions = dim;
        Vjson{sci}.VoxelSize  = round(vx,4);
        % b-values of DWI for the b-value shells of the protocol test (see
        % normalizeFields in cat_io_dcm2bids_protocols)
        Pbval = spm_file(Pjson{fscni},'ext','.bval');
        if exist(Pbval,'file')
          Vjson{sci}.BValues = sscanf(fileread(Pbval),'%f')';
        end
      end

      % process each scan only once, i.e., the same series (and dcm2niix 
      % image suffix, e.g. 7a or 17_e2 for multiple images of a series), 
      % e.g. if a scan directory contains copies from different DICOM 
      % directories (older imports)
      scanid{sci} = cat_io_dcm2bids_db('getScanID',Vjson{sci},fnameparts); 
      if any(strcmp(scanid(1:sci-1), scanid{sci}))
        cat_io_cprintf([.5 .5 .5],'%60s : ignore repeated scan\n',fname1); 
        sci = sci - 1; % reset import counter
        continue; 
      end
    


      %% create DB structure
      % to get an overview of imported data to define coding/filter files
      % ===================================================================
      % * overview directories/files with csv/tsv-tables and/or json files:
      %    - scans     (internal list to avoid double imports)
      %    - sessions  
      %    - subjects  (maybe relevant later to add basic phenotype data, 
      %                 e.g., subject- or session-specific)   
      %    - studies   (to organize different projects)
      %    - protocols (to create own protocol filter lists)
      %
      % * Data is internally stored by in the dbdir with subject, session
      %   and series-number. It might be possible to extend the IDs by 
      %   name, date, protocol for human readability but let's start simple.
      % ===================================================================
      Pdbdirpath{sci} = Pdirfi; 
      Pdbdir{sci}     = Pdirfi(numel(Pmdbdir)+2:end); 
      Pdbnii{sci}     = Pnii{fscni}; % image in the database (e.g. for the session registration)
      Pdbjson{sci}    = Pjson{fscni}; % sidecar in the database (e.g. for the DICOM export)
      dbV{end+1} = Vjson{sci}; %#ok<AGROW> % subject table of the database (written after all packages, see cat_io_dcm2bids)

      % test protocols
      [match,pname,mismatchstr,Pnfailedid,Pfname,PT] = ...
        cat_io_dcm2bids_protocols('testProtocols',Pjson{fscni},Vjson{sci},protocols,tol,job,sites);
      Pprot{sci}   = Pfname;    % (closest) protocol file, e.g. for the coverage mask
      pstatus{sci} = PT.status; % protocol status, e.g. for the DICOM export

      % site definition: the label of the BIDS export (subject labels, 
      % sessions files, DICOM export; without location for the anonymization
      % level 2, see setupSites) and of the tables of the database (siteDB)
      site{sci} = cat_io_dcm2bids_protocols('setupSites',Vjson{sci},sites,job.opts);
      siteDB    = cat_io_dcm2bids_protocols('setupSites',Vjson{sci},sites);
  
      % study definition %%% need later refinement
      studynam = cat_io_dcm2bids_bids('bidsLabel',job.subdir);


      % main BIDS fields (subject, session, datatype, weighting, ...)
      % ===================================================================
      % subject: prefix (site and/or study) and PatientID, or the number of
      % the subject (anonymized), without the prefix for an ID that starts
      % with it (e.g. of an anonymized DICOM export, i.e. the same label)
      bs = job.opts.BIDSsep; dayoffset = 0; 
      if job.opts.anonymize > 1 % number and date offset of the subject (anonymized)
        [~,dbsub] = fileparts(fileparts(fileparts(Pdirfi))); 
        [NPID,dayoffset] = cat_io_dcm2bids_anonymize('subject',anon,dbsub); 
      else
        NPID = strrep(cat_io_dcm2bids_bids('bidsLabel',Vjson{sci}.PatientID),'_',''); 
      end
      sitelab = cat_io_dcm2bids_bids('bidsLabel',site{sci}); 
      switch job.opts.subIDform
        case 1, pre = {};                  % sub-PID
        case 2, pre = {sitelab};           % sub-SITE-PID
        case 3, pre = {studynam};          % sub-STUDY-PID
        case 4, pre = {sitelab,studynam};  % sub-SITE-STUDY-PID
        case 5, pre = {studynam,sitelab};  % sub-STUDY-SITE-PID
      end
      if ~isempty(pre) && numel(NPID) > numel([pre{:}]) && startsWith(NPID,[pre{:}])
        NPID = NPID(numel([pre{:}]) + 1:end); 
      end
      sub{sci} = ['sub-' strjoin([pre {NPID}],bs)]; %#ok<*SAGROW>
      
      % session and its acquisition time of the sessions file (by the 
      % anonymization level: date and time, date and full hour, or date of
      % the session label and the 3-hour interval of the time of the day, 
      % see acqTime in cat_io_dcm2bids_anonymize)
      ses{sci}     = cat_io_dcm2bids_bids('sessionLabel',fileparts(Pdirfi),fnameparts{3},job.opts,dayoffset); 
      acqtime{sci} = cat_io_dcm2bids_anonymize('acqTime',fnameparts{3},dayoffset,job.opts.sesunit,job.opts.anonymize); 

      if job.opts.output == 0 % only the tables of the database
        [datatype{sci},suffix{sci},ent{sci}] = scanType(pname,Vjson{sci},Isrc,isbids); 
        [rec{sci},varg{sci}] = scanRecord(Vjson{sci},scanid{sci},siteDB,fileparts(Pdbdir{sci}), ...
          datatype{sci},suffix{sci},ent{sci},PT,''); 
        fprintf('\n'); continue
      end
      if job.opts.output == 2 && strcmp(PT.status,'nonconform') 
        % output 2: only scans with a conform or semiconform protocol (PFR 
        % <= opts.protocol.status(2), all scans without protocols), i.e. 
        % scans with a nonconform protocol are only part of the tables of 
        % the database (output 3 for all scans)
        [datatype{sci},suffix{sci},ent{sci}] = scanType(pname,Vjson{sci},Isrc,isbids); 
        [rec{sci},varg{sci}] = scanRecord(Vjson{sci},scanid{sci},siteDB,fileparts(Pdbdir{sci}), ...
          datatype{sci},suffix{sci},ent{sci},PT,''); 
        cat_io_cprintf([.5 .5 .5],' nonconform protocol (not exported, see output 3)\n'); continue
      end

      % protocol directory 
      if match > 0
        if job.opts.protocolsubdirs
          BIDSsubdir = fullfile(job.BIDSdir, job.BIDSsubdir, spm_file(protocols{match,1},'basename'));
        else
          BIDSsubdir = fullfile(job.BIDSdir, job.BIDSsubdir);
        end
        pro{sci}   = pname;
      else
        if isempty(Pprodictdirs) || isempty( Pprodictdirs{1} ) || ~job.opts.protocolsubdirs
          BIDSsubdir  = fullfile(job.BIDSdir, job.BIDSsubdir);
        else
          BIDSsubdir  = fullfile(job.BIDSdir, job.BIDSsubdir, 'mismatch');
        end
        pro{sci} = fnameparts{2}; % protocol name (with separators for setupSuffix)
      end
      BIDSsub{sci} = BIDSsubdir; 

      % datatype, suffix, further entities (e.g. dir, echo, part), and the 
      % label of the acq entity, i.e., the protocol name without the 
      % information of the suffix and entities (e.g. SBRef, phase, echo, 
      % direction), and the task 
      [datatype{sci},suffix{sci},ent{sci},label] = scanType(pro{sci},Vjson{sci},Isrc,isbids); 
      lab{sci} = cat_io_dcm2bids_bids('bidsLabel',label,1); 
      series = pro{sci}; 
      FN = {'ProtocolName', 'SeriesDescription', 'SequenceName'};
      for fni = 1:numel(FN)
        if isfield(Vjson{sci},FN{fni}), series = [series ' ' Vjson{sci}.(FN{fni})]; end %#ok<AGROW>
      end
      task{sci} = cat_io_dcm2bids_bids('setupTask',datatype{sci}, series );
      if isbids % BIDS input: task of the source dataset
        if isfield(Isrc.Entities,'task'), task{sci} = ['_task-' Isrc.Entities.task]; else, task{sci} = ''; end
      end
      

      % define BIDS naming
      % ===================================================================
      BIDSdir{sci}   = fullfile(sub{sci}, ses{sci}, datatype{sci}); 
      BIDSpath{sci}  = fullfile(Poutdir, BIDSsubdir, BIDSdir{sci}); 
      BIDSpathd{sci} = fullfile(cat_io_dcm2bids_helper('derivativesDir',fullfile(Poutdir, BIDSsubdir)), BIDSdir{sci}); 
      % %%%%%%%%%%%%%%%%%%%%%%%%%% refine run definition 
      % acq entity with the series number and the label (a run entity would 
      % require all scans of the session, as only the second of two equal 
      % scans shows that the first one needs a run number): 
      % - images of one series with the same suffix and entities (e.g. 
      %   further images of dcm2niix) get also the image number 
      % - images of one acquisition (consecutive scans of the session with 
      %   the same datatype and label, e.g. SBRef and bold, field-map 
      %   magnitude and phase, or AP and PA field maps) get the series 
      %   number of the first one, as BIDS associates such files by their
      %   equal entities (e.g. phasediff and magnitude1, AP and PA)
      sn{sci}      = sprintf('%03.0f',Vjson{sci}.SeriesNumber); 
      scankey{sci} = cat_io_dcm2bids_db('scanKey',Vjson{sci}); % stored, as cleanupVjson removes the Patient fields
      prev = find( strcmp( scankey(1:sci-1) , scankey{sci} ) ); 
      numbered = any( strcmp( suffix(prev) , suffix{sci} ) & strcmp( ent(prev) , ent{sci} ) ); 
      if ~numbered && ~isbids && numel(Pjson) > fscni
        % look ahead to the next file of the series (but no repeated scan, 
        % e.g. the same series converted from another directory)
        Pjsonnext = spm_file(Pjson{fscni+1},'path',Pdirfi); 
        if exist(Pjsonnext,'file')
          Vjsonnext = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pjsonnext)); 
          fnamenext = cat_io_dcm2bids_db('nameParts',spm_file(Pjsonnext,'basename')); 
          if strcmp( cat_io_dcm2bids_db('scanKey',Vjsonnext) , scankey{sci} ) && ...
            ~any(strcmp( scanid(1:sci) , cat_io_dcm2bids_db('getScanID',Vjsonnext,fnamenext) ))
            [suffixnext,entnext] = cat_io_dcm2bids_bids('setupSuffix',datatype{sci},pro{sci},Vjsonnext); 
            numbered = strcmp(suffixnext,suffix{sci}) && strcmp(entnext,ent{sci}); 
          end
        end
      end
      if numbered, acqnr{sci} = sprintf('%s%d', bs, numel(prev) + 1); else, acqnr{sci} = ''; end
      acqsn{sci} = sn{sci}; 
      if sci > 1 && isempty(acqnr{sci}) && isempty(acqnr{sci-1}) && ...
          strcmp(BIDSpath{sci-1},BIDSpath{sci}) && strcmp(lab{sci-1},lab{sci}) 
        % name with the series number of the previous scan, if not yet used
        BIDSfile1 = sprintf('%s_%s%s_acq-%s%s%s%s_%s.json', sub{sci}, ses{sci}, task{sci}, ...
          acqsn{sci-1}, bs, lab{sci}, ent{sci}, suffix{sci}); 
        if ~any(strcmp( BIDSfile(strcmp(BIDSpath(1:sci-1),BIDSpath{sci})) , BIDSfile1 ))
          acqsn{sci} = acqsn{sci-1}; 
        end
      end
      acq{sci}      = sprintf('acq-%s%s%s%s', acqsn{sci}, acqnr{sci}, bs, lab{sci}); 
      BIDSfile{sci} = sprintf('%s_%s%s_%s%s_%s.json', ...
        sub{sci}, ses{sci}, task{sci}, acq{sci}, ent{sci}, suffix{sci});
      if ~match && ~isempty(protocols) && ~isempty(protocols{1})
        for pfi = 1:numel(Pnfailedid)
          pdir      = fullfile( BIDSpath{sci} , 'matchfiles' );
          if ~exist(pdir,'dir'), mkdir(pdir); end
          BIDSfile2 = spm_file(BIDSfile{sci},'ext','.csv'); 

          try
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','current','required'}; ...
              mismatchstr{ Pnfailedid(pfi) }] ); 
          catch
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','current','required'}; ...
              {mismatchstr(Pnfailedid(pfi),1),'mlt. values','mlt. values'}] ); 
          end
        end
      end

      % record of the scan for the tables of the database (see 
      % cat_io_dcm2bids_tables, with the QC rating added below)
      [rec{sci},varg{sci}] = scanRecord(Vjson{sci},scanid{sci},siteDB,fileparts(Pdbdir{sci}), ...
        datatype{sci},suffix{sci},ent{sci},PT,BIDSfile{sci}); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data preparation 
% - anonymizing, basic-preprocessing, QC .. all outputs are  and would have to be packed 

      %% (anonymized) input file for internal processing (not zipped!) 
      if job.opts.output > 1
        Panon{fscni} = cat_io_dcm2bids_pp_anat('anonymize',Pnii{fscni}, job.opts, datatype{sci}); 
        

        % Basic QC with 
        %  - BSM (Between Scan Motion) in 4D data (average correction of the realignment)
        %  - WSM (Within Scan Motion) in 4D data (image variance between scans)
        %  - ISR (Inhomogeneity to Signal Ratio)
        %  - NSR (Noise to Signal Ratio)
        %  - RES (RMS of voxel resolution)
        if job.opts.runqc
          [Prnii{fscni},QM(fscni),Pqc{fscni}] = cat_io_dcm2bids_qc(spm_file(Panon{fscni}), datatype{sci}, job.opts, Pfname); %#ok<AGROW>
          if exist(spm_file(Pqc{fscni},'ext','json'),'file')
            Jqc = cat_io_json(spm_file(Pqc{fscni},'ext','json')); 
            if isfield(Jqc,'qualityrating') && isfield(Jqc.qualityrating,'SQR') && isnumeric(Jqc.qualityrating.SQR) && ~isempty(Jqc.qualityrating.SQR)
              rec{sci}.SQR = double(Jqc.qualityrating.SQR); 
            end
          end
        else
          fprintf('\n'); % end of the scan line (otherwise ended by the QC values)
        end
        % (the segmentation of the anatomical scan of the session follows at
        % the end of the session, see anatScan in cat_io_dcm2bids_pp_anat)
      else
        fprintf('\n');
      end

      
      
      

      %% create result dir and copy files
      % ===================================================================
      if job.opts.output > 0
        if ~exist(BIDSpath{sci},'dir'), mkdir(BIDSpath{sci}); end
        if job.opts.output > 1
          % remove older derivatives of the scan, as all its derivatives are
          % written again by this run (e.g. maps of former versions such as
          % the noise maps of the QC or segmentations of further scans), 
          % but not the results of the preprocessing of former runs without
          % preprocessing in this run (prefix before the BIDS name, e.g. c0,
          % m, mwc1, or catDCM2BIDSsegus_, see exportSegmentation and the
          % DTI and rs-fMRI measures; realigned 4D data r* of former 
          % versions are removed)
          Dold = dir(fullfile(BIDSpathd{sci},['*' spm_file(BIDSfile{sci},'basename') '*']));
          Dold = Dold(~[Dold.isdir]);
          if ~job.opts.preprocessing
            pre  = regexprep({Dold.name},'sub-.*$','');
            Dold = Dold(cellfun('isempty',regexp(pre,'^(c0|m|mwc1|neuromorphometrics_|catDCM2BIDS(segus|dwi|dwimix|rsfmri)_.*)$','once')));
          end
          for di = 1:numel(Dold), delete(fullfile(Dold(di).folder,Dold(di).name)); end
          % scan, b-values and directions, QC, and the realignment of 4D data
          % (rp_*.txt of the QC, without resliced images, see 
          % cat_io_dcm2bids_qc; the segmentation of the session is exported 
          % at the end of the session, see exportSegmentation)
          der = [0 0 0 1 1];
          pre = {'','','','catDCM2BIDSqc_','rp_'};
          ext = {niiext,'.bval','.bvec','.json','.txt'};
          for ei = 1:numel(ext)
            if exist( spm_file(Pjson{fscni},'ext',ext{ei},'prefix',pre{ei}) , 'file' ) || ...
               exist( spm_file(Pjson{fscni},'ext',ext{ei},'prefix',[pre{ei} 'anon_']) , 'file' )
              if der(ei), BIDSpath0 = BIDSpathd{sci}; else, BIDSpath0 = BIDSpath{sci}; end
              if ~exist(BIDSpath0,'dir'), mkdir(BIDSpath0); end
              
              file = spm_file(Pjson{fscni},'ext',ext{ei},'prefix',pre{ei}); 
              if exist( spm_file(file,'prefix','anon_'),'file')
                file = spm_file(file,'prefix','anon_'); 
              elseif exist( spm_file(Pjson{fscni},'ext',ext{ei},'prefix',[pre{ei} 'anon_']),'file')
                file = spm_file(Pjson{fscni},'ext',ext{ei},'prefix',[pre{ei} 'anon_']);
              end
    
              if strcmp(ext{ei},niiext)
                cat_io_dcm2bids_helper('exportNii',file, ...
                  fullfile(BIDSpath0,[pre{ei} spm_file(BIDSfile{sci},'basename')]), job.opts); 
              else
                copyfile( file, spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext',ext{ei},'prefix',pre{ei})); 
              end
            end
          end

          % denoised average of 4D data of the QC (see cat_io_dcm2bids_qc and
          % cat_io_dcm2bids_render), but not the noise maps of older versions
          Pavg = [regexprep(Panon{fscni},'\.nii(\.gz)?$','') '_denoisedavg' niiext];
          if exist(Pavg,'file')
            if ~exist(BIDSpathd{sci},'dir'), mkdir(BIDSpathd{sci}); end
            cat_io_dcm2bids_helper('exportNii',Pavg, fullfile(BIDSpathd{sci}, ...
              ['catDCM2BIDSqc_denoisedavg_' spm_file(BIDSfile{sci},'basename')]), job.opts);
          end

          % protocol fit measure and rating (not part of the SQR) and the 
          % protocol test (DICOM protocol name, closest protocol, status, and 
          % deviations) in the QC JSON of the scan (see testProtocols in 
          % cat_io_dcm2bids_protocols and cat_io_dcm2bids_report) 
          if ~isnan(PT.PFM)
            Pqcj = fullfile(BIDSpathd{sci},['catDCM2BIDSqc_' spm_file(BIDSfile{sci},'basename') '.json']); 
            if exist(Pqcj,'file'), J = cat_io_json(Pqcj); else, J = struct(); end
            J.qualitymeasures.PFM = PT.PFM; J.qualityrating.PFR = PT.PFR; 
            J.protocol = struct('ProtocolName',Vjson{sci}.ProtocolName,'SeriesNumber',Vjson{sci}.SeriesNumber, ...
              'protocol',PT.protocol, ...
              'status',PT.status,'deviations',{cellfun(@(f,v,p) sprintf('%s: %s (%s)',f,v,p), ...
              PT.deviations(:,1),PT.deviations(:,2),PT.deviations(:,3),'UniformOutput',false)}); 
            if ~exist(BIDSpathd{sci},'dir'), mkdir(BIDSpathd{sci}); end
            cat_io_json(Pqcj,J); 
          end
        end

        % anonymize and save json file
        Vjsonp = cat_io_dcm2bids_protocols('cleanupVjson',Vjson{sci},'bids'); 
        if ~isempty(task{sci}) && ~isfield(Vjsonp,'TaskName') % required by BIDS
          Vjsonp.TaskName = regexprep(task{sci},'^_task-',''); 
        end
        if contains(ent{sci},'_part-phase') && ~isfield(Vjsonp,'Units') % required by BIDS
          Vjsonp.Units = 'arbitrary'; % e.g. -4096 to 4095 of Siemens phase images 
        end
        % image dimensions and voxel size (if the image is available)
        if isfield(Vjson{sci},'Dimensions'), Vjsonp.Dimensions = Vjson{sci}.Dimensions; end
        if isfield(Vjson{sci},'VoxelSize'),  Vjsonp.VoxelSize  = Vjson{sci}.VoxelSize;  end
        % scanner fields of the anonymization level 2 (no model and coil, 
        % only the type of the software, see cat_io_dcm2bids_anonymize)
        Vjsonp = cat_io_dcm2bids_anonymize('sidecar',Vjsonp,job.opts.anonymize); 
        cat_io_json( spm_file( fullfile(BIDSpath{sci},BIDSfile{sci}),'ext','.json'), Vjsonp); 

        % files that belong to the scan (BIDS input, e.g. events, physio, 
        % aslcontext) with the new name (see bidsAssociated in cat_io_dcm2bids_importBIDS)
        newbase = regexprep(spm_file(BIDSfile{sci},'basename'),'_[^_]+$',''); 
        Passoc  = dir(fullfile(Pdirfi,['catDCM2BIDSassoc_' spm_file(Pjson{fscni},'basename') '_*'])); 
        for ai = 1:numel(Passoc)
          rest = Passoc(ai).name(numel(['catDCM2BIDSassoc_' spm_file(Pjson{fscni},'basename') '_'])+1:end); 
          copyfile(fullfile(Pdirfi,Passoc(ai).name), fullfile(BIDSpath{sci},[newbase '_' rest])); 
        end

        % new names of BIDS inputs and IntendedFor references (remapped after 
        % the export of all scans)
        if isbids
          if job.opts.gzipe, ext0 = '.nii.gz'; else, ext0 = '.nii'; end
          srcmap(cat_io_dcm2bids_helper('normPath',Isrc.SourceFile)) = struct('root',BIDSsubdir,'sub',sub{sci}, ...
            'rel',fullfile(ses{sci},datatype{sci},[spm_file(BIDSfile{sci},'basename') ext0])); 
          if isfield(Vjson{sci},'IntendedFor') && isfield(Vjsonp,'IntendedFor')
            intended{end+1} = struct('json',spm_file(fullfile(BIDSpath{sci},BIDSfile{sci}),'ext','.json'), ...
              'root',BIDSsubdir,'sub',sub{sci},'srcroot',Isrc.SourceDataset, ...
              'srcsub',['sub-' Isrc.Entities.sub],'entries',{cellstr(Vjson{sci}.IntendedFor)}); %#ok<AGROW>
          end
        end


        %%%% all protocols

        %%%% all subjects


      end


      %% write private and participant data
      % ===================================================================
      % The private.tsv should contain fields that are removed in the BIDS
      % processing such as the real Patient name and his birth data etc. 
      % It might be saved in another directory to avoid unwanted uploading?
      % ===================================================================
%spm_file(Pprodictdirs,'basename')      
      % participant data of the participants.tsv and private table, which 
      % are written once at the end (see cat_io_dcm2bids_bids), with the data 
      % of the first scan (minimum age) and further columns of BIDS imports 
      if isbids && isfield(Isrc,'Participant'), extra = Isrc.Participant; else, extra = []; end
      pkey = [BIDSsubdir '|' sub{sci}]; 
      if ~isKey(partdata,pkey) || Vjson{sci}.PatientAge < partdata(pkey).V.PatientAge || ...
          (isnan(partdata(pkey).V.PatientAge) && ~isnan(Vjson{sci}.PatientAge))
        partdata(pkey) = struct('BIDSsubdir',BIDSsubdir,'sub',sub{sci},'V',Vjson{sci},'extra',extra); 
      end

    end

    % end of a session (the scan directories are sorted, i.e., the scans of a
    % session are processed one after another): segmentation of one
    % anatomical scan of the session, affine registration of the
    % session to MNI space with the orientation rating of the session and
    % the coverage rating of its scans (see writeSessionAffines in 
    % cat_io_dcm2bids_pp_anat), the phase-encoding reversal test of its scans 
    % (checkPEreversal), and the session data of the subject (sessions.tsv) 
    if job.opts.output > 0 && (fdiri == numel(Psnr) || ~strcmp(fileparts(Psnr{fdiri+1}),fileparts(Pdirfi)))
      sesids = find(strcmp(cellfun(@fileparts,Pdbdirpath(1:sci),'UniformOutput',false),fileparts(Pdirfi)) & ...
        ~cellfun('isempty',BIDSsub(1:sci))); 
      SOR = nan; 
      if job.opts.output > 1 && ~isempty(sesids)
        % segmentation of the anatomical scan of the session (one per session,
        % e.g. the T1w with the best quality rating, see anatScan in 
        % cat_io_dcm2bids_pp_anat) and its derivatives 
        ia = sesids(cat_io_dcm2bids_pp_anat('anatScan',Pdbnii(sesids),datatype(sesids),suffix(sesids))); 
        if ~isempty(ia) && job.opts.preprocessing
          try
            cat_io_dcm2bids_pp_anat('segmentanat',cat_io_dcm2bids_pp_anat('anonymize',Pdbnii{ia},job.opts,datatype{ia}), ...
              datatype{ia},job.opts); 
            exportSegmentation(Pdbnii{ia},BIDSpathd{ia},BIDSfile{ia},job.opts); 
          catch e
            cat_io_cprintf('err','%60s : segmentation failed: %s\n',spm_file(BIDSfile{ia},'basename'),e.message); 
          end
        end
        % anatomical scans of each subject (by the BIDS label, i.e. also of
        % other scanners for subject labels without site, see opts.subIDform)
        % for the sessions without anatomical scan, whose DTI and rs-fMRI 
        % are processed at the end of the run (see below)
        if ~isempty(ia)
          if isKey(anatrun,sub{ia}), anatrun(sub{ia}) = [anatrun(sub{ia}) ia]; else, anatrun(sub{ia}) = ia; end
        end
        noanat = isempty(ia) && job.opts.preprocessing && any(ismember(datatype(sesids),{'dwi','func'})); 
        if noanat
          sesnoanat{end+1} = sesids; %#ok<AGROW>
        else
          S = cat_io_dcm2bids_pp_anat('writeSessionAffines',Pdbdirpath(sesids), Pdbnii(sesids), datatype(sesids), ...
            suffix(sesids), BIDSpathd(sesids), sub(sesids), ses(sesids), job.opts, BIDSfile(sesids), Pprot(sesids)); 
          if ~isempty(S), SOR = S(1).SOR; end
        end
        cat_io_dcm2bids_pp_anat('checkPEreversal',Pdbnii(sesids),Vjson(sesids),datatype(sesids), ...
          BIDSpathd(sesids),BIDSfile(sesids),job.opts); 
        if noanat
          cat_io_cprintf([.5 .5 .5],'%60s : no anatomical scan (DTI and rs-fMRI at the end of the run)\n', ...
            sprintf('session %s',spm_str_manip(spm_file(fileparts(Pdirfi),'basename'),'l50'))); 
        elseif job.opts.preprocessing % DTI and rs-fMRI measures of the session (see cat_io_dcm2bids_pp_dwi/_rsfmri)
          cat_io_dcm2bids_pp_dwi('processSession',Pdbnii(sesids),Vjson(sesids),datatype(sesids), ...
            suffix(sesids),BIDSpathd(sesids),BIDSfile(sesids),job.opts); 
          cat_io_dcm2bids_pp_rsfmri('processSession',Pdbnii(sesids),Vjson(sesids),datatype(sesids), ...
            suffix(sesids),BIDSpathd(sesids),BIDSfile(sesids),job.opts); 
        end
      end
      for bsd = reshape(unique(BIDSsub(sesids)),1,[]) % each BIDS (protocol) directory of the session
        id = sesids(find(strcmp(BIDSsub(sesids),bsd{1}),1,'last')); 
        cat_io_dcm2bids_bids('writeSubjectTSV',Vjson{id}, sub{id}, ses{id}, Poutdir, bsd{1}, job.opts.anonymize, SOR, site{id}, ...
          cat_io_dcm2bids_phenotypic('sessionData',fileparts(Pdbdirpath{id})), acqtime{id}); 
      end
    end
  end 
  if sci>1
    fprintf('%s\n',repmat('-',1,154)); 
    fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
  end

  % sessions without anatomical scan (e.g. longitudinal studies with one 
  % T1w): DTI and rs-fMRI with the segmented anatomical scan of the closest
  % session (scan date, also of a later session) of the subject in this 
  % run, with a warning (see method.T1w_other_session of the DTI and 
  % rs-fMRI, the tables, and the renderings), and the registration of the
  % session by the coregistration of its DTI or rs-fMRI to this scan (see 
  % coregSessionAffine in cat_io_dcm2bids_pp_anat) 
  if ~isempty(sesnoanat)
    fprintf('\nDCM2BIDS - sessions without anatomical scan\n'); 
    for ni = 1:numel(sesnoanat)
      sesids = sesnoanat{ni}; pses = sesids; R = []; SOR = nan; 
      sesline = [sub{sesids(1)} '_' ses{sesids(1)}]; 
      if isKey(anatrun,sub{sesids(1)})
        ix = anatrun(sub{sesids(1)}); 
        dd = arrayfun(@(k) abs(days(Vjson{k}.ScanDate - Vjson{sesids(1)}.ScanDate)),ix); dd(isnan(dd)) = inf;
        [~,m] = min(dd); ix = ix(m); pses = [sesids ix]; 
        cat_io_cprintf('warn','%60s : DTI and rs-fMRI with %s\n',sesline,spm_file(BIDSfile{ix},'basename')); 
      else
        cat_io_cprintf([.5 .5 .5],'%60s : no anatomical scan of the subject\n',sesline); 
      end
      cat_io_dcm2bids_pp_dwi('processSession',Pdbnii(pses),Vjson(pses),datatype(pses), ...
        suffix(pses),BIDSpathd(pses),BIDSfile(pses),job.opts); 
      cat_io_dcm2bids_pp_rsfmri('processSession',Pdbnii(pses),Vjson(pses),datatype(pses), ...
        suffix(pses),BIDSpathd(pses),BIDSfile(pses),job.opts); 
      if numel(pses) > numel(sesids)
        R = cat_io_dcm2bids_pp_anat('coregSessionAffine',Pdbnii(sesids),fileparts(Pdbdirpath{ix}), ...
          Pdbnii{ix},datatype{ix},suffix{ix},ses{ix},job.opts); 
      end
      S = cat_io_dcm2bids_pp_anat('writeSessionAffines',Pdbdirpath(sesids), Pdbnii(sesids), datatype(sesids), ...
        suffix(sesids), BIDSpathd(sesids), sub(sesids), ses(sesids), job.opts, BIDSfile(sesids), Pprot(sesids), R); 
      if ~isempty(S), SOR = S(1).SOR; end
      for bsd = reshape(unique(BIDSsub(sesids)),1,[]) % (session orientation rating in the sessions file)
        id = sesids(find(strcmp(BIDSsub(sesids),bsd{1}),1,'last')); 
        cat_io_dcm2bids_bids('writeSubjectTSV',Vjson{id}, sub{id}, ses{id}, Poutdir, bsd{1}, job.opts.anonymize, SOR, site{id}, ...
          cat_io_dcm2bids_phenotypic('sessionData',fileparts(Pdbdirpath{id})), acqtime{id}); 
      end
    end
  end

  % results of the scans for the outputs of all packages (see cat_io_dcm2bids)
  % (fields of the scans with sci elements, e.g. without BIDS name for output 0)
  fs = @(c) [c(1:min(end,sci)) cell(1,sci - min(numel(c),sci))]; 
  R  = struct('rec',{fs(rec)},'varg',{fs(varg)},'dbV',{dbV},'Pdbjson',{fs(Pdbjson)}, ...
    'datatype',{fs(datatype)},'pstatus',{fs(pstatus)},'site',{fs(site)},'sub',{fs(sub)}, ...
    'ses',{fs(ses)},'BIDSfile',{fs(BIDSfile)},'BIDSsub',{fs(BIDSsub)},'intended',{intended}, ...
    'srcmap',srcmap,'partdata',partdata,'bidsroots',{bidsroots},'convver',{convver}); 

  if 0 
      job.opts.reg2MNI = 2; % 0-no, 1-yes,mat-sidecar, 2-yes,permanent
      if job.opts.reg2MNI
      %% apply registration to MNI space
      %  We can do this only after full import as we need to use the best data 
      %  to estimate the best way to get a good rigid registration setup. 
      %  Besides the anatomical scans also the localizer/scouts would be handy :D 
        fprintf('\nDCM2BIDS - adjust orientation\n')
        Psesdb  = spm_file(Pdbdirpath(1:sci),'path');  
        [Psesdbu,u1,u2] = unique(Psesdb); 
        Psesdbu0 = Pdbdir(1:sci); Psesdbu0 = Psesdbu0(u1); 
        for sesi = 1:numel(Psesdbu)
          fprintf('  Adjust session "%s"\n',Psesdbu0{sesi})
          %% get Affine transformation from the segmentation's
          Pseg8 = cat_vol_findfiles(Psesdbu{sesi},'*seg8.mat'); 
          for s8i = 1:numel(Pseg8)
            seg(s8i) = load(Pseg8{s8i}); %#ok<AGROW>
          end
          Affine = mean(cat(3, seg(:).Affine),3);
          imat   = spm_imatrix(Affine);                 
          Rigid  = spm_matrix(imat(1:6));   
          % apply the closest rigid transformation to bring the data closer to MNI
          Pc = spm_file( fullfile( BIDSpath(find(u2==sesi)) , BIDSfile(find(u2==sesi)) ),'ext',niiext); %#ok<FNDSB>
          for pci = 1:numel(Pc)
            fprintf('    %s\n',spm_fileparts( Pc{pci} ,'basename'))
    % maybe store the mat independently        
            Pci = cat_io_dcm2bids_helper('prepNii',Pc{pci},job.opts);
            Pmat = spm_file(Pci,'ext','mat');
            if exist(Pmat,'file'), delete(Pmat); end
            try
              M = spm_get_space(Pci);         % read current affine (header only) %%%%%% 4D?
            catch
              cat_io_cprintf('err','error reading "%s"\n',Pci);
              delete(Pci)
              continue
            end
            if job.opts.reg2MNI==1
              % soft application with NIFTI sidecar
              mat = Rigid * M; 
              save(spm_file(Pci,'ext','mat'),'mat');
              delete(Pci)
            else
              % hard application by updating the NFITI
              spm_get_space(Pci, Rigid * M);       % write new affine into the header
              cat_io_dcm2bids_helper('prepNiigz',Pci,job.opts);
            end
          end
    
        end
        fprintf('DCM2BIDS - adjust orientation done.\n')
      end
  end
end
% =========================================================================
% =========================================================================
function [datatype,suffix,ent,label] = scanType(pro,V,Isrc,isbids)
%scanType. BIDS datatype, suffix, further entities (e.g. dir, echo, part),
%  and acq label of a scan by its protocol name pro and sidecar V (with the
%  image fields, e.g. for protocol names without datatype, see 
%  setupDatatype and setupSuffix in cat_io_dcm2bids_bids), or the datatype, 
%  suffix, and entities of the source dataset of BIDS inputs. 
  if isbids
    datatype = Isrc.Datatype; 
  else
    datatype = cat_io_dcm2bids_bids('setupDatatype',pro,V); 
  end
  [suffix,ent,label] = cat_io_dcm2bids_bids('setupSuffix',datatype,pro,V);
  if isbids
    suffix = Isrc.Suffix; 
    ent    = cat_io_dcm2bids_bids('entityString',Isrc.Entities); 
  end
end
% =========================================================================
function [R,varg] = scanRecord(V,scanid,site,session,datatype,suffix,ent,PT,BIDSfile)
%scanRecord. Record of a scan for the tables of the database (see 
%  scanRecord in cat_io_dcm2bids_tables) without its protocol variant, 
%  which is set after all packages in the order of the scans (varg, see 
%  protocolVariant in cat_io_dcm2bids_protocols and cat_io_dcm2bids), as 
%  the variants are files of the database that parallel processes would 
%  write at the same time. 
  R    = cat_io_dcm2bids_tables('scanRecord',V,scanid,site,session,datatype,suffix,ent, ...
    '',PT,BIDSfile);
  varg = struct('PT',PT,'ProtocolName',V.ProtocolName,'datatype',datatype,'type',[suffix ent]); 
end
% =========================================================================
function exportSegmentation(Pdbnii,BIDSpathd,BIDSfile,opts)
%exportSegmentation. BIDS derivatives of the segmentation of the database
%  image Pdbnii (the anatomical scan of a session, see anatScan and
%  segmentanat in cat_io_dcm2bids_pp_anat): the volumes (JSON) and for
%  opts.preprocessing > 1 the label map c0, the bias-corrected image m, 
%  the modulated normalized GM map mwc1, and the atlas label map 
%  neuromorphometrics_ (of the anonymized or denoised image), which are 
%  otherwise only registered in the derivatives index of the database (for
%  the renderings and the report tables, see registerDerivative in 
%  cat_io_dcm2bids_helper).
  [pth,nam] = fileparts(regexprep(Pdbnii,'\.gz$',''));
  bnam = spm_file(BIDSfile,'basename');
  Pj   = fullfile(pth,['catDCM2BIDSsegus_' nam '.json']);
  if ~exist(Pj,'file'), return; end
  if ~exist(BIDSpathd,'dir'), mkdir(BIDSpathd); end
  copyfile(Pj,fullfile(BIDSpathd,['catDCM2BIDSsegus_' bnam '.json']));
  for pre = {'c0','m','mwc1','neuromorphometrics_'}
    for p2 = {'anon_','sanlm_anon_','sanlm_',''}
      P = fullfile(pth,[pre{1} p2{1} nam '.nii.gz']); if ~exist(P,'file'), P = P(1:end-3); end
      if exist(P,'file')
        cat_io_dcm2bids_helper('registerDerivative',P,fullfile(BIDSpathd,[pre{1} bnam]));
        if opts.preprocessing > 1
          cat_io_dcm2bids_helper('exportNii',P,fullfile(BIDSpathd,[pre{1} bnam]),opts); 
        end
        break
      end
    end
  end
end
