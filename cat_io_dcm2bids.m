function cat_io_dcm2bids(job)
%
%

% QUESTIONS: 
% * Use own subjectIDs ? 
%   No  - to avoid misalignment
%   Yes - to avoid human error - as number using private.tsv >> function 
%       - we will need this at least for the strong anonymizing setting
%       - a subject.tsv would also allow to integrate essential phenotypical  
%         data for groups or timepoints 
%
% * Flexibility
%   - tolerance parameter? but what about ordinary variables
%   - fallback option?
%
%
% TODO: 
% *** flag to write non-conform data or only protocols? 
%     eg. to extract only T1w data (see also recursive calls)
%     ... GUI done
%
% *** flag to extract only protocols (1) [and jsons (2)?] for protocol setup 
%     >> extract protocols!!! >> noQC, noPP, noAnon 
%
% *** QC (realignment) is slow ... maybe a flag ... 
%     if you have different levels you would need separation of files
%     0-no QC, 1-basic QC (just a subset or fast things?), 2-extend QC (full), 3-extend QC (for optimization/pp)
%
% * Addition protocol set directories?
%      protocols        anat1, ..., fMRI1...
%   >> study-defined    
%   >> subject-defined
%   >> session-defined
%
%
% *** overview dirs: 
%   * Protocols ... subdirs
%   * Studies (one-file with subjects and protocol-names)
%      JE2: subjects, images, StudyID
%   *** create a quality file based on the first scans (n>5)?  
%
% * Recursive calls / already converted/sorted data (i.e., nii/json import)
%
% * Non-conform/unknown protocol handling
%
% * Report Files (sub/scan) ****
%
% * write job-style into the outdir to control/monitor assure same structure 
%   and avoid weird combinations
%   - if outdir exist then check if the jog-style fits otherwise halt with 
%     message: (1) go on with old setting, (2) stop and change ouput dir
%
% * para-fields final DZPG protocol: 
%   - slice timing
%   - # slices
%   - head position ?
% * QC-test: head-orientation ?
%
%
% TODO-BONUS:
% * Anonymizing 
%   * own sub-ids
%   * avoid scan-dates on level 2.
%   * json-checks ...
%   * 4D call ******
%   - setting ... probably combine to keep it simple
%      - naming:  site: difficult as federated 
%                 sub:  0-sub-id,    1-sub-id,  2-new-id
%                 ses:  0-date-time, 1-date,    2-TP (order issues?)
%      - data:    0-no, 1-standard, 2-strong "json-cleanup"
%      - imgs:    0-no, 1-anat,     2-all
%
%
% * Further data: 
% =========================================================================
%   - spectroscopy
%   - MPM processing
%
%
% * Basic preprocessing, QC and anonymizing of (un)organized data?
% =========================================================================
%   - Preprocessing is a very complex issue that need general consensus! 
%     Yes, but we can make a basic suggestion that can be extended later. 
%     However, this might be irrelevant if things are done by the DZNE!
%     However, it should be an external batch that could be flagged here.
%       - T1w/T2w: SPM/CAT
%       - MPMs:    hMRI
%       - dMRI:    diffusion TB
%       - fMRI:    SPM
%   - Basic processing to detect neurological outliers?
%     This requires a normative model. 
%     It presents a logical step when the preprocessing is established!
%   - Basic optimization and inter-modality optimization/harmonization. 
%     This can be seen as additional part to further improve preprocessing, 
%     i.e., after preprocessing is established.
%       1) bias
%       2) denoising
%       3) contrast ("global mean")
%      (4) reorient & BB (rigid registration to MNI-space)
%   - How to add the data?
%     >> derivatives/TOOL/...
%     >> derivatives/catDCM2NII/../anat/[m|c0|mc1|y]*.nii(.gz) ... qc-files?
%                                  func/r*..
%
% * write report tables
% =========================================================================
%   - Overview for DZPG/Imaging, i.e., on line per study:
%     (Name, ID, number of (in)complete subjects, long/cross-design, modalities, avg. quality score)
%   - Overview for PI, i.e., one file per study/protocolset with on line per subject:
%     (SubID, ...)
%   - Detailed overview, i.e., one file for all protocolset
%
%    
%
% * protocol/QC check: 
% =========================================================================
%   - number of images/slices
%   - dcm convert error (missing dcm-files sliced etc.)
%
%
% * Features
% =========================================================================
%   - parallelization 
%   - import timer (subject-wise)
%
%
% * ISSUES & BUGS:
% =========================================================================
%  * BUG: handling of multiple runs ...
%  * avoid some read_vol for speed? 
%  * no qc-file message
%  * copy also qc-file for fits
%
%
% * TOTEST:
% =========================================================================
%   - get BIDS test cases
%   - QC Tests (MR-ART?, DZPG samples)

  def.data                  = {};    % input DCM directories (add JSON/NII input later)
  def.outdir                = {pwd}; % main output directory
  def.subdir                = 'CATDCM2NII'; % default
  
  def.dicts.Pprotocoldirs   = {};    % input DCM protocol directories
  def.dicts.Pcenterdict     = {};    % dictionary for center names (otherwise scanner ID)
  def.dicts.Pstudydict      = {};    % not implemented yet
  def.dicts.Psubjdict       = {};    % not implemented yet

  def.opts.ProtocolFileName = 1;     % replace protocol name by the filenames of the evaluation protocols under Pprotocoldirs 
  def.opts.gzipi            = 1;     % internal use of nii.gz (save disk space but maybe a bit slower) 
  def.opts.gzipe            = 1;     % external use of nii.gz (save disk space but nonoptimal for SPM processing)
  def.opts.tolerance        = 0;     % tolerance in percent for MR parameters (does not help for ordinal variables)
  def.opts.Pdcm2nii         = getDCM2NIIX; 
  def.opts.verb             = 1;     % 0-no, 1-minimal, 2-extensive
  def.opts.output           = 1; 
  def.opts.studies          = ''; 
  def.opts.subIDform        = 2;     % (1) use only PatientID, i.e. sub-PID 
                                     % (2) add the ScannerID, i.e. sub-SITE-PID 
  def.opts.anonymize        = 1;     % imaging data: 0-no, 1-light-anat-only,  2-strong-all, 3-skull-strip-light?
                                     % meta data:    0-no, 1-basic-reduced,    2-strong-onlycodeds
  def.opts.removeScanTime   = 1;
  def.opts.deletetmp        = 0;     % remove temporary directory (better not)
  
  def.opts.hrrealign        = 0;     % highres-realignment for further use
  def.opts.rerun            = 0;     % rerun - overwrite existing
  def.opts.tableformat      = 'csv'; % csv/tsv 
  
  if ~exist('job','var'), job = struct(); end
  job = cat_io_checkinopt(job,def); 

  job.opts.gzipi             = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

  %% ======================================================================
  if 0
    % my quick non-GUI test
    Pdcmdir      = {'/Users/robertdahnke/MRData/20260819 - CIRC-CAT - DCM2BIDS/DCM'};
    Pprodictdirs = {'/Users/robertdahnke/MRData/20260819 - CIRC-CAT - DCM2BIDS/protocols/DZPG'};
    Pcenterdict  = {'/Users/robertdahnke/MRData/20260819 - CIRC-CAT - DCM2BIDS/protocols/cites.json'}; 
    Pstudydict   = {'/Users/robertdahnke/MRData/20260819 - CIRC-CAT - DCM2BIDS/protocols/studies.json'}; 
    %Psubjdict    = {''}; 
    Poutdir      = '/Users/robertdahnke/MRData/20260819 - CIRC-CAT - DCM2BIDS/TMP5'; 
    tol          = 1; 
  else
    Pdcmdir      = job.data;
    Pprodictdirs = job.dicts.Pprotocoldirs;
    Pcenterdict  = job.dicts.Pcenterdict;
    Pprodictdirs = job.dicts.Pprotocoldirs;
    Pstudydict   = job.dicts.Pstudydict;
    %Psubjdict    = job.dicts.Psubjdict;
    Poutdir      = fullfile(job.outdir{1},job.subdir); 
    tol          = job.opts.tolerance;
  end
  if ~exist(Poutdir,'dir'), mkdir(Poutdir); end

  Pdbdirnam = 'catDCM2NIIdb';
  Popts = fullfile(Poutdir,Pdbdirnam,'catdcm2bids.mat');
  if ~checkPreviousSetting(Popts,job), return; end
  clear Popts; 

  % read site-names and protocol definitions
  sites     = getSites(Pcenterdict);
  %studies   = getStudies(Pstudydict);  %%%%%%%%%%%%
  %subj      = getSubjects(Psubjdict);  %%%%%%%%%%%%
  protocols = getProtocols(Pprodictdirs);  

  % get all sub-directories
  Pdcmdirs = {}; Pdcmdirs0 = {};
  for di = 1:numel(Pdcmdir)
    Pdcmdirsdi = cat_vol_findfiles(Pdcmdir{di},'*',struct('dirs',1));
    Pdcmdirsdi(cellfun(@(x) numel(x)>1,strfind(Pdcmdirsdi,Pdbdirnam))) = []; 
    Pdcmdirsdi(cat_io_contains(lower(Pdcmdirsdi),{'localizer','scout'})) = []; 
    Pdcmdirs   = [Pdcmdirs;  Pdcmdirsdi]; %#ok<AGROW>
    Pdcmdirs0  = [Pdcmdirs0; repmat(Pdcmdir(di),size(Pdcmdirsdi))]; %#ok<AGROW>
  end
  
  %%%%% special case of the internal database directories: DCM2NIIX case
  % DB subdirs
  Pdatadir     = fullfile(Poutdir,Pdbdirnam,'+data'); 
  Pprotocoldir = fullfile(Poutdir,Pdbdirnam,'+protocols'); 
  
  if ~exist(Pdatadir,'dir'),     mkdir(Pdatadir); end
  if ~exist(Pprotocoldir,'dir'), mkdir(Pprotocoldir); end
  
  % DB lists 
  Pscantable     = spm_file(fullfile(Pdatadir,'scans'),'ext',job.opts.tableformat); 
  Pstudytable    = spm_file(fullfile(Pdatadir,'studies'),'ext',job.opts.tableformat); 
  Psessiontable  = spm_file(fullfile(Pdatadir,'sessions'),'ext',job.opts.tableformat); 
  Psubjecttable  = spm_file(fullfile(Pdatadir,'subjects'),'ext',job.opts.tableformat); 
  Pprotocoltable = spm_file(fullfile(Pprotocoldir,'protocols'),'ext',job.opts.tableformat); 


  %% basic initialization that might have to be extended
  Vjson = cell(1,numel(Pdcmdirs)); jsoni = 0;
  sub = Vjson; ses = Vjson; datatype = Vjson; pro = Vjson; acq = Vjson; 
  run = Vjson; task = Vjson; suffix = Vjson; site = Vjson; 
  Panon = Vjson; Pp0 = Vjson; Pwc1 = Vjson; Pdbdir = Vjson;
  Pdbdirpath = Vjson;
  BIDSpath = Vjson; BIDSdir = Vjson; BIDSfile = Vjson;
  PID = ''; sni = 0; stime = datetime('now'); 
  QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]); 
  for fi = 1:numel(Pdcmdirs)
    
    % Convert DCM: 
    % =====================================================================
    Pdirfi = dcm2niix( Pdcmdirs{fi} , Pdcmdirs0{fi}, Poutdir, job.opts.Pdcm2nii, job.opts); 
%%%%% else gzip or gunzip depending all files
%%%%% special case of the internal database directories: DCM2NIIX case

  
    % get json files
    % =====================================================================
    if job.opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end
    Pjson  = cat_vol_findfiles(Pdirfi,'*.json',struct('depth',1));
    Pnii   = spm_file(Pjson,'ext',niiext); 
    for fj = 1:numel(Pnii) % for every scan
      jsoni = jsoni + 1;
  

      % get nii2dcm information 
      Vjson{jsoni} = cat_io_json(Pjson{fj});
      if isfield( Vjson{jsoni} , 'ImageTypeText')
        Vjson{jsoni}.ImageType = unique([Vjson{jsoni}.ImageType; Vjson{jsoni}.ImageTypeText ]); 
      end
      Vjson{jsoni} = assurePatientDCMfields(Vjson{jsoni});
      fnameparts   = strsplit(spm_file(Pjson{fj},'basename'),'='); 
      Vjson{jsoni}.ProtocolName = fnameparts{2};

      % ignore localizer and scout scans
      if any(cat_io_contains(lower(fnameparts),{'localizer','scout'}))
        continue; 
      end
      % ignore 2D data
      if exist(Pnii{fj},'file')
        Vsz = dir(Pnii{fj});
        if isempty(Vsz) || ~isfield(Vsz,'bytes'), continue; end
        if Vsz.bytes/1024 < 1000 
          evalc('V = spm_vol(Pnii{fj});'); 
          if numel(V.dim)>2 && any(V.dim < 2), continue; end
        end
      end

      % create table header
      if ~strcmp(PID, Vjson{jsoni}.PatientID)
        sni = sni + 1;

        if jsoni>1
          fprintf('%s\n',repmat('-',1,154)); 
          fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
          stime = datetime('now'); 
        end

        cat_io_cprintf([0 0.2 .8],'\n%-45s',sprintf('%4d) %s', sni, ...
          spm_str_manip(sprintf('sub-%s-%s',Vjson{jsoni}.DeviceSerialNumber, ...
            strrep(Vjson{jsoni}.PatientID,'_','')),'a38')));
        if job.opts.output
          if isempty(Pprodictdirs) || isempty(Pprodictdirs{1})
            fprintf('%-15s : %-50s%10s %5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              ' ',' ','BSM','WSM','ISR','NSR','RES','SQR'); 
          else
            fprintf('%-15s : %-50s%10s %5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              'matching-protocol','dismatch','BSM','WSM','ISR','NSR','RES','SQR'); 
          end
        end
        fprintf('%s\n',repmat('-',1,154)); 
        PID = Vjson{jsoni}.PatientID; 
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
      
      Pdbdir{jsoni} = fullfile( ...
        sprintf('sub-%s',   Vjson{jsoni}.PatientID), ...
        sprintf('ses-%s',   Vjson{jsoni}.StudyID), ...
        sprintf('snr-%03d', Vjson{jsoni}.SeriesNumber)); 

      % update 
      Pconv = spm_file(Pjson{fj},'ext',niiext,'path', Pdbdirpath{jsoni}); 
      if exist(Pconv,'file'), Pjson{fj} = Pconv; Pjson{fj} = Pconv; Pnii{fj} = spm_file(Pjson{fj},'ext',niiext); end
      ProcedureStepDescription = getFileString(Vjson{jsoni}.ProcedureStepDescription);

      % SCAN-LIST:
      %  - to avoid double entries ... but what to do ...
      %  - a flag (keep old / take new)
%%%%%%%  * this grows too strong >> internal variable + mat     
      if 0
        Thdr    = {'SeriesInstanceUID','DBpath'}; %,'IMPORTpath'}; 
        Tnewrow = {Vjson{jsoni}.SeriesInstanceUID, Pdbdir{jsoni}}; %, Pjson{fj}}; 
        updateTable(Pscantable,Thdr,Tnewrow,1,job.opts.rerun);
      end

      % SESSION-LIST:
      %  - to avoid double entries ... but what to do ...
      %  - a flag (keep old / take new)
%%%%%%%  * this one is usesless      
      if 0
        Thdr    = {'StudyInstanceUID','DBpath','ProcedureStepDescription'}; %,'IMPORTpath'}; 
        Tnewrow = {Vjson{jsoni}.StudyInstanceUID, spm_fileparts(Pdbdir{jsoni}), Vjson{jsoni}.ProcedureStepDescription}; %, Pjson{fj}}; 
        updateTable(Psessiontable,Thdr,Tnewrow,1,job.opts.rerun);
      end

      % SUBJECT-LIST ("constant" data):
%%%%%%%  * this one is nice but a session counter would be nice 
%%%%%%%  * maybe the scan period 
%%%%%%%  * import path (not robust :/)
      Thdr    = {'PatientID','PatientName', 'PatientSex','PatientBirthDate'}; 
      Tnewrow = {Vjson{jsoni}.PatientID, Vjson{jsoni}.PatientName, ...
                 Vjson{jsoni}.PatientSex, Vjson{jsoni}.PatientBirthDate}; 
      updateTable(Psubjecttable,Thdr,Tnewrow,1,job.opts.rerun);
      
      % STUDY-LIST (defined by ProcedureStepDescription with additional list of subjects):
      % - ProcedureStepDescription with subjects
      Pstudysubtable = spm_file(fullfile(Pdatadir,'study_subject_lists',ProcedureStepDescription),'ext',job.opts.tableformat); 
      Thdr    = {'Subjects'}; 
      Tnewrow = {Vjson{jsoni}.PatientID}; 
      nsub    = updateTable(Pstudysubtable,Thdr,Tnewrow,1,job.opts.rerun);
      % - study-list with number of subjects
      Thdr    = {'ProcedureStepDescription','nSubjects'}; 
      Tnewrow = {ProcedureStepDescription,nsub}; 
      updateTable(Pstudytable,Thdr,Tnewrow,1,job.opts.rerun);
% with counters for anat, func, fmap, dwi ?      
      
      % PROTOCOLS-LIST:
      datatype0         = setupDatatype(Vjson{jsoni}.ProtocolName); 
      ProtocolName      = getFileString(strrep(fnameparts{2},'_','-'));
      %Pprotocolsubtable = spm_file(fullfile(Pprotocoldir,ProcedureStepDescription,datatype0,ProtocolName),'ext',job.opts.tableformat); 
      Pprotocolsubtable = spm_file(fullfile(Pprotocoldir,datatype0,ProtocolName),'ext',job.opts.tableformat); 

      % save shorted protocol as json to be used as filter
      Vjson0 = cleanupVjson(Vjson{jsoni},'protocol'); 
      Pjson0 = spm_file( Pprotocolsubtable,'ext','.json'); 
      ProtocolName0 = ProtocolName; 
      if ~exist(Pjson0,'file')
        cat_io_json( Pjson0, Vjson0); 
      else
        %Pjsons = cat_vol_findfiles( fullfile(Pprotocoldir,ProcedureStepDescription,datatype0)  , '*.json');
        Pjsons = cat_vol_findfiles( fullfile(Pprotocoldir,datatype0)  , '*.json');
        Pfit = false(size(Pjsons)); 
        for pi = 1:numel(Pjsons)
          Vjson1 = cat_io_json( Pjsons{pi} ); 
          Pfit(pi) = structEqual(Vjson0,Vjson1); 
        end
        if ~any( Pfit )
          ci = 0;
          while exist(Pjson0,'file')
            ci = ci+1;
            Pjson0 = spm_file(Pprotocolsubtable,'suffix',sprintf('_%d',ci),'ext','.json'); 
            ProtocolName0 = spm_file(ProtocolName,'suffix',sprintf('_%d',ci)); 
          end
          cat_io_json( Pjson0, Vjson0); 
        end
      end
      Thdr    = {'Subjects'}; 
      Tnewrow = {Vjson{jsoni}.PatientID}; 
      nsub    = updateTable(Pprotocolsubtable,Thdr,Tnewrow,1,job.opts.rerun);

      % - study-list with number of subjects
      Thdr    = {'ProtocolName','nSubjects','ProcedureStepDescription', ...
                 'MRAcquisitionType','ScanningSequence', 'SequenceVariant', ...
                 'RepetitionTime', 'EchoTime', ...
                 'FlipAngle', 'SliceThickness', 'ImageType'}; 
      Vjson0 = Vjson{jsoni}; 
      for fni = 1:numel(Thdr)
        if ~isfield(Vjson0,Thdr{fni}), Vjson0.(Thdr{fni}) = ''; end
        if iscell(Vjson0.(Thdr{fni}))
          Vjson0.(Thdr{fni}) = strjoin(Vjson0.(Thdr{fni})); 
        end
      end
      Tnewrow = {sprintf('%s_%s',datatype0,ProtocolName0),nsub, Vjson0.ProcedureStepDescription, ...
        Vjson0.MRAcquisitionType, Vjson0.ScanningSequence,  Vjson0.SequenceVariant, ...
        Vjson0.RepetitionTime, Vjson0.EchoTime, ...
        Vjson0.FlipAngle, Vjson0.SliceThickness, Vjson0.ImageType}; 
      updateTable(Pprotocoltable,Thdr,Tnewrow,1,job.opts.rerun);


      % import files
      Pdbdirpath{jsoni} = fullfile(Poutdir,Pdbdirnam, Pdbdir{jsoni}); 
      if ~exist(Pdbdirpath{jsoni},'dir')
        mkdir(Pdbdirpath{jsoni}); 
        ext = {'.nii','.nii.gz','.bval','.bvec'};  
        for ei = 1:numel(ext)
          if exist(spm_file(Pjson{fj},'ext',ext{ei}),'file') 
            if ~exist( spm_file(Pjson{fj},'ext',ext{ei},'path', Pdbdirpath{jsoni}),'file')
              movefile( spm_file(Pjson{fj},'ext',ext{ei}) , Pdbdirpath{jsoni} );
            else
              delete( spm_file(Pjson{fj},'ext',ext{ei}) );
            end
          end
        end
        copyfile( Pjson{fj} , Pdbdirpath{jsoni} );
      end


      if job.opts.output == 0
        continue
      end
      Pjson{fj} = spm_file(Pjson{fj},'path', Pdbdirpath{jsoni});
      Pnii{fj}  = spm_file(Pnii{fj}, 'path', Pdbdirpath{jsoni});

      % test protocols
      [match,pname,dismatchstr,Pnfailedid,Pfname] = testProtcols(Pjson{fj},Vjson{jsoni},protocols,tol,job,sites);
  
      % site definition 
      site{jsoni} = setupSites(Vjson{jsoni},sites);
  

      % main BIDS fields (subject, session, datatype, weighting, ...)
      % ===================================================================
      % subject
      switch job.opts.subIDform
        case 1 % sub-PID
          sub{jsoni} = sprintf('sub-%s', strrep(Vjson{jsoni}.PatientID,'_','')); %#ok<*SAGROW>
        case 2 % sub-SITE-PID
          sub{jsoni} = sprintf('sub-%s-%s', site{jsoni}, strrep(Vjson{jsoni}.PatientID,'_','')); %#ok<*SAGROW>
      end

      % session
      if job.opts.anonymize > 1 % ses-StudyID
        ses{jsoni}  = sprintf('ses-%s',Vjson{jsoni}.StudyID);
      else % ses-date
        ses{jsoni}  = sprintf('ses-%s',fnameparts{3}(1:min(8,numel(fnameparts{3}))));
      end
      if numel(fnameparts{3}) > 8  &&  job.opts.anonymize==0
        ses{jsoni}  = sprintf('%s-%s',ses{jsoni},fnameparts{3}(9:end)); 
      end

      % protocol directory 
      if match
        BIDSsubdir  = sprintf('BIDS-%s',spm_file(protocols{match,1},'basename'));
        pro{jsoni}  = pname;
      else
        if isempty(Pprodictdirs) || isempty( Pprodictdirs{1} )
          BIDSsubdir  = 'BIDS';
          pro{jsoni}  = strrep(fnameparts{2},'_','-'); 
        else
          BIDSsubdir  = 'BIDS-dismatch';
          pro{jsoni}  = strrep(fnameparts{2},'_','-'); 
        end
      end

      % evaluate protocols
      datatype{jsoni} = setupDatatype(pro{jsoni}); 
      task{jsoni}     = setupTask(pro{jsoni} );
      % %%%%%%%%%%%%%%%%%%%%%%%%%% refine run definition 
      run{jsoni}      = sprintf('%03.0f',str2double(Vjson{jsoni}.SeriesNumber)); % fnameparts{4})); 
      acq{jsoni}      = sprintf('acq-%s-%s', run{jsoni}, pro{jsoni}); 
      % get suffix 
      [suffix{jsoni},acq{jsoni}] = setupSuffix(datatype{jsoni}, acq{jsoni}, Vjson{jsoni}.SeriesDescription);
      

      % define BIDS naming
      % ===================================================================
      BIDSdir{jsoni}  = fullfile(sub{jsoni}, ses{jsoni}, datatype{jsoni}); 
      BIDSpath{jsoni} = fullfile(Poutdir, BIDSsubdir, BIDSdir{jsoni}); 
      BIDSfile{jsoni} = sprintf('%s_%s_%s%s_%s.json', ...
        sub{jsoni}, ses{jsoni}, acq{jsoni}, task{jsoni}, suffix{jsoni});
      if ~match && ~isempty(protocols) && ~isempty(protocols{1})
        for pfi = 1:numel(Pnfailedid)
          %spm_file(char(Pprodictdir{Pnfailedid(pfi)}),'basename');
          pdir      = fullfile( BIDSpath{jsoni} , 'matchfiles' );
          if ~exist(pdir,'dir'), mkdir(pdir); end
          acq2      = sprintf('acq-%s-%s', run{jsoni}, pro{jsoni}); 
          BIDSfile2 = sprintf('%s_%s_%s%s_%s.csv', ...
            sub{jsoni}, ses{jsoni}, acq2, task{jsoni}, suffix{jsoni});

          try
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','is','should'}; dismatchstr{ Pnfailedid(pfi) }] ); 
          catch
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','is','should'}; {dismatchstr(Pnfailedid(pfi),1),'mlt. values','mlt. values'}] ); 
          end
        end
      end

      % add further fields?
      %Vjson{jsoni}.site = site

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data preparation 
% - anonymizing, basic-preprocessing, QC .. all outputs are  and would have to be packed 
      
      
      
      


      %% (anonymized) input file for internal processing (not zippered!) 
      if job.opts.output > 1
        Panon{fj} = anomize(Pnii{fj}, job.opts, datatype{jsoni}); 
      end    


      %% Basic QC:
      %  - NSR: Noise to Signal Ratio
      %  - ISR: Inhomogeneity to Signal Ratio
      %  - BSM: Between Scan Motion
      %  - RES: Resolution rating
      if job.opts.output > 1
        [Prnii{fj},QM(fj)] = runQC(spm_file(Panon{fj}), datatype{jsoni}, job.opts, Pfname); 
      else
        fprintf('\n');
      end


      %% Basic preprocessing 
      % segment anatomical scan ... what to do otherwise? how to save/fill data?
      if job.opts.output > 1  &&  strcmp(datatype{jsoni},'anat')
        [Pm{fj},Pp0{fj},Pwc1{fj}] = segmentanat(Panon{fj}, datatype{jsoni}, job.opts);
      end


      
      % write private and participant data
      % ===================================================================
      % The private.tsv should contain fields that are removed in the BIDS
      % processing such as the real Patient name and his birth data etc. 
      % It might be saved in another directory to avoid unwanted uploading?
      % ===================================================================
      writeScanReportTSV(Vjson{jsoni},sub{jsoni},site{jsoni},Poutdir,BIDSsubdir,'report'); %%%%%%%%%%
      writePrivateTSV(Vjson{jsoni},sub{jsoni},Poutdir,BIDSsubdir,'private')
      writeParticipantTSV(Vjson{jsoni},sub{jsoni},Poutdir,BIDSsubdir);
    

  % TODO: Create time point data file
  % Ptimepoints = fullfile(Poutdir,'timepoints.tsv');
 



      % create result dir and copy files
      % ===================================================================
      if job.opts.output > 0
        if ~exist(BIDSpath{jsoni},'dir'), mkdir(BIDSpath{jsoni}); end
        if job.opts.output > 1
          ext = {niiext,'.bval','.bvec'};  
          for ei = 1:numel(ext)
            if exist( spm_file(Pjson{fj},'ext',ext{ei}) , 'file' )
              if ~exist(BIDSpath{jsoni},'dir'), mkdir(BIDSpath{jsoni}); end
  
              file = spm_file(Pjson{fj},'ext',ext{ei}); 
              if exist( spm_file(file,'prefix','anon_'),'file')
                file = spm_file(file,'prefix','anon_'); 
              end
    
              if strcmp(ext{ei},niiext)
              % in case of niftis we try to avoid to copy to save time  
                if job.opts.gzipe
                  if ~exist(spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii.gz'),'file')
                    if ~job.opts.gzipi
                      gzip( file )
                      movefile( [file '.gz'], spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii.gz')); 
                    else
                      copyfile( file, spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii.gz')); 
                    end
                  end
                else
                  if ~exist(spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii'),'file')
                    if job.opts.gzipi
                      gunzip( file )
                      movefile( [file(1:end-3)], spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii')); 
                    else
                      copyfile( file, spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.nii')); 
                    end
                  end
                end
              else
                copyfile( file, spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext',ext{ei})); 
              end
            end
          end
        end

        % anonymize and save json file
        Vjson{jsoni} = cleanupVjson(Vjson{jsoni}); %%%%%%%%%%%%%%%%%%%%%%%
        cat_io_json( spm_file( fullfile(BIDSpath{jsoni},BIDSfile{jsoni}),'ext','.json'), Vjson{jsoni}); 


        %%%% all protocols

        %%%% all subjects


      end
    end
  end 
  fprintf('DCM2NII - import done.\n')

  if job.opts.output > 0
    % create final report form result dir
    writeSubReportTSV(Poutdir,BIDSsubdir,'report');
  
    % create/extend BIDS csv files ???
  
  
  
    % handling GZIP in output directory
    gzipunzipOutputdata(Poutdir,job.opts)
  end

end
% =========================================================================
function T = structEqual(S1,S2)
  if ~isstruct(S1), T = false; return; end
  if ~isstruct(S2), T = false; return; end

  FN1 = sort(fieldnames(S1)); 
  FN2 = sort(fieldnames(S2));
 
  if numel(FN1) ~= numel(FN2), T = false; return; end
  if any(~cellfun(@(x,y) strcmp(x,y), FN1, FN2)), T = false; return; end

  T = true;
  for fni = 1:numel(FN1)
    if (isnumeric( S1.(FN1{fni}) ) && isnumeric( S2.(FN1{fni}) )) || ...
       (islogical( S1.(FN1{fni}) ) && islogical( S2.(FN1{fni}) ))
      T = S1.(FN1{fni}) == S1.(FN1{fni}); 
    elseif ischar( S1.(FN1{fni}) ) && ischar( S2.(FN1{fni}) ) 
      T = strcmp(S1.(FN1{fni}),S1.(FN1{fni}));
    elseif isstruct( S1.(FN1{fni}) ) && isstruct( S2.(FN1{fni}) ) 
      T = structEqual(S1,S2);
    elseif iscellstr( S1.(FN1{fni}) ) && iscellstr( S2.(FN1{fni}) )   %#ok<ISCLSTR>
      T = strcmp(char(join(S1.(FN1{fni}))),char(join(S2.(FN1{fni})))); 
    elseif iscell( S1.(FN1{fni}) ) && iscell( S2.(FN1{fni}) )  
      error('structEqual:not implemented');
    else
      T = false; 
      return
    end
    if ~T; return; end
  end 
end
% =========================================================================
function entries = updateTable(Ptable,Thdr,Tnewrow,id,reimport)
  if exist(Ptable,'file')
    Tfiles = cat_io_csv(Ptable);
    entries = size(Tfiles,1)-1;
    
    % check if first ID entry already exists
    if isnumeric(Tfiles{2,1}) && ~reimport
      RID = Tnewrow{id};
      if ischar(RID), RID = str2double(Tnewrow{1}); end
      if any( [Tfiles{2:end,1}] == RID ) && ~reimport
        return
      end
    else
      if any( cat_io_contains( Tfiles(2:end,1), Tnewrow{id} ) & ...
          (cellfun(@(x)numel(x),Tfiles(2:end,1)) == numel(Tnewrow{id})) ) && ~reimport
        return
      end
    end

    % convert fields numbers if required
    for ci = 1:size(Tnewrow,2)
      if ischar(Tnewrow{ci}) && isnumeric(Tfiles{2,ci})
        Tnewrow{ci} = str2double(Tnewrow{ci}); 
      end
    end
    
    % add and sort rows
    Tfiles(end+1,:) = Tnewrow; 
    Tfiles(2:end,:) = sortrows(Tfiles(2:end,:),1);
    cat_io_csv(Ptable,Tfiles);
  else
    Tfiles = [Thdr;Tnewrow];
    cat_io_csv(Ptable,Tfiles);
  end
  entries = size(Tfiles,1)-1;
end
% =========================================================================
function str = getFileString(str)
  str = cat_io_strrep(str, ...
    {'ä','ü','ö',' '}, {'ae','ue','oe','_'});
  strd = double(str);
  str(strd<48 & strd~=45) = '_';
  str(strd>57 & strd<65)  = '_';
  str(strd>90 & strd<97 & strd~=95)  = '_';
  str(strd>122) = 'X';
end
% =========================================================================
function same = checkPreviousSetting(Popts,job)
  same = 1; 
  if exist(Popts,'file')
    load(Popts,'opts','dicts');
    FNopts = fieldnames(job.opts);
    for fni = 1:numel(FNopts)
      if ~isfield(opts,FNopts{fni})
        if isscalar(job.opts.(FNopts{fni})) && isscalar(opts.(FNopts{fni}))
          if job.opts.(FNopts{fni}) ~= opts.(FNopts{fni}), same = 0; end
        elseif isstring(job.opts.(FNopts{fni})) && isstring(opts.(FNopts{fni}))
          if ~strcmp(job.opts.(FNopts{fni}), opts.(FNopts{fni})), same = 0; end
        else
          same = 0; 
        end
      end
    end
    FNopts = fieldnames(job.dicts);
    for fni = 1:numel(FNopts)
      if ~isfield(dicts,FNopts{fni})
        if isscalar(job.dicts.(FNopts{fni})) && isscalar(dicts.(FNopts{fni}))
          if job.dicts.(FNopts{fni}) ~= dicts.(FNopts{fni}), same = 0; end
        elseif isstring(job.dicts.(FNopts{fni})) && isstring(dicts.(FNopts{fni}))
          if ~strcmp(job.dicts.(FNopts{fni}), dicts.(FNopts{fni})), same = 0; end
        else
          same = 0; 
        end
      end
    end
    if ~same 
      cat_io_cprintf('err', ...
        ['Error the underlying structure of the directory does not fit to the current settings. \n', ...
        'Use previous parameters or or change the output directory. \nMove on with prevous paramters?']);
    end
    %%
  else
    if ~exist(spm_file(Popts,'path'),'dir'), mkdir(spm_file(Popts,'path')); end
    opts = job.opts; dicts = job.dicts; 
    save(Popts,'opts','dicts');
    clear opts dicts; 
  end
end
% =========================================================================
function gzipunzipOutputdata(Poutdir,opts)
  BIDSsubdirs = cat_vol_findfiles( Poutdir , 'BIDS*', struct('dirs',1)); 
  if opts.gzipi == 1  &&  opts.gzipe == 0
    % gunzip all files in the result directory
    fprintf('DCM2NII - gunzip BIDS output data')
    for bi = 1:numel(BIDSsubdirs)
      Pniigz = cat_vol_findfiles( BIDSsubdirs{bi} , '*.nii.gz' ); 
      for fi = 1:numel(Pniigz), gunzip(Pniigz{fi}); delete(Pniigz{fi}); end
    end
    fprintf('done.\n')
  elseif opts.gzipi == 0  &&  opts.gzipe == 1
    % gzip all files in the result directory
    fprintf('DCM2NII - gzip BIDS output data')
    for bi = 1:numel(BIDSsubdirs)
      Pniigz = cat_vol_findfiles( BIDSsubdirs{bi} , '*.nii' ); 
      for fi = 1:numel(Pniigz), gzip(Pniigz{fi}); delete(Pniigz{fi}); end
    end
    fprintf('done.\n')
  end
end
% =========================================================================
function P = getDCM2NIIX
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
    error('Cannot find dcm2niix. Please install it from: \n  "%s"', ... 
      spm_file('https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage' , ...
        'link', 'https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage')); 
  end
end
% =========================================================================
function tmpdir = dcm2niix( Pdcmdirs , Pdcmdir, Poutdir, Pdcm2niix, opts)
  %tmpdir = strrep( Pdcmdirs , Pdcmdir, fullfile(Poutdir,'DCM2NIIX') ); 
  tmpdir = fullfile(Poutdir,'catDCM2NIIx',Pdcmdirs); 
  if ~exist(tmpdir,'dir') || opts.rerun
    mkdir(tmpdir); 
  
    if opts.gzipi, gz = '-z'; else, gz = ''; end

    % convert DCM in directory
    % - use = as more unique separator 
    cmd = sprintf('%s -f "%%f=%%p=%%t=%%s" -p n %s y -ba n -o "%s" "%s"', ...
      Pdcm2niix, gz, tmpdir,  Pdcmdirs); 
    [status,cmdout] = system(cmd); %#ok<ASGLU>
  end

  if opts.gzipi
    P = cat_vol_findfiles( tmpdir , '*.nii' ,struct('depth',1)); 
    for fi=1:numel(P), gzip(P{fi}); delete(P{fi}); end
  else
    P = cat_vol_findfiles( tmpdir , '*.nii.gz' ,struct('depth',1)); 
    for fi=1:numel(P), gunzip(P{fi}); delete(P{fi}); end
  end
end
% =========================================================================
function [match,pname,dismatchstr,Pnfailedid,Pfname] = testProtcols(Pjson,Vjson,protocols,tol,job,sites)
% test protocols

  if isempty(protocols{1,1})
  % no given protocols  
    match           = false;
    pname           = ''; 
    Pnfailedid      = {};
    Pfname          = ''; 
    pmatchs         = inf; 
    pmatch          = 0; 
    dismatchstr{1}  = cell(0,3); 
    
  else
    pmatch      = ones(1,size(protocols,1)); 
    pmatchs     = zeros(1,size(protocols,1)); 
    dismatchstr = cell(1,size(protocols,1)); 
    for pri = 1:size(protocols,1)
      dismatchstr{pri} = cell(0,3); 
      FNpi = fieldnames(protocols{pri,3});
      FNpi(cat_io_contains(FNpi,{'ConsistencyInfo','PulseSequenceDetails', ...
        'ImageComments','SequenceName','ProtocolName','SeriesDescription'})) = []; 
      pmatchn = 1; 
      pmatchs(pri) = numel(FNpi); 
      for fni = 1:numel(FNpi)
        if isfield( Vjson, FNpi{fni} ) 
          if strcmp(FNpi{fni}(1:2),'x_'), continue; end % comment
          
          if ( islogical( Vjson.(FNpi{fni})) || isnumeric( Vjson.(FNpi{fni})) ) && ...
             ( islogical( protocols{pri,3}.(FNpi{fni})) || isnumeric( protocols{pri,3}.(FNpi{fni})) )
            if numel( Vjson.(FNpi{fni}) ) == numel( protocols{pri,3}.(FNpi{fni}))
              pmatchn =  all( Vjson.(FNpi{fni}) >= protocols{pri,3}.(FNpi{fni})*(1-tol/100) & ...
                              Vjson.(FNpi{fni}) <= protocols{pri,3}.(FNpi{fni})/(1-tol/100)); 
            else
              pmatchn = 0; 
            end
          elseif ischar( Vjson.(FNpi{fni}) ) && ischar( protocols{pri,3}.(FNpi{fni}) )
            pmatchn = strcmp( Vjson.(FNpi{fni}) , protocols{pri,3}.(FNpi{fni}) ); 
          elseif iscell( Vjson.(FNpi{fni}) ) && iscell( protocols{pri,3}.(FNpi{fni}) )
            pmatchn = strcmp( char(Vjson.(FNpi{fni}(:))) , char(protocols{pri,3}.(FNpi{fni}(:))) ); 
          else
            pmatchn = 0; 
            % need refinement ! ... image type             
          end
          pmatch(pri) = pmatch(pri) & pmatchn; 
  
          if ~pmatchn
            if iscell( Vjson.(FNpi{fni}) ) 
              tstr1 = ''; 
              for ci = 1:numel( Vjson.(FNpi{fni}) )
                tstr1 = [tstr1 '+' Vjson.(FNpi{fni}){ci}]; %#ok<AGROW>
              end
              tstr2 = ''; 
              for ci = 1:numel( protocols{pri,3}.(FNpi{fni}) )
                tstr2 = [tstr2 '+' protocols{pri,3}.(FNpi{fni}){ci}]; %#ok<AGROW>
              end
            else
              if isscalar(Vjson.(FNpi{fni}))
                tstr1 = Vjson.(FNpi{fni});
              else
                tstr1 = sprintf('%dx%d %s', size(Vjson.(FNpi{fni}),1), ...
                  size(Vjson.(FNpi{fni}),2), class(Vjson.(FNpi{fni}))); 
              end
              if isscalar(Vjson.(FNpi{fni}))
                tstr2 = protocols{pri,3}.(FNpi{fni});
              else
                tstr2 = sprintf('%dx%d %s', size(protocols{pri,3}.(FNpi{fni}),1), ...
                  size(protocols{pri,3}.(FNpi{fni}),2), class(protocols{pri,3}.(FNpi{fni}))); 
              end
            end
            dismatchstr{pri} = [ dismatchstr{pri}; {FNpi{fni} tstr1 tstr2 }]; 
            
          end
        end
      end
      pmatch(pri) = pmatch(pri) & pmatchn; 
    end
  end
%% always print session 

  % display progress
  % =======================================================================
  fname1     = sprintf('%3d) %s', Vjson.SeriesNumber, spm_str_manip(strrep(sprintf('%s_%s_%s', ...
                strrep(Vjson.PatientID, strrep(char(Vjson.ScanDate),'-',''),''), ...
                strrep(char(Vjson.ScanDate),'-',''), Vjson.ProtocolName),'__','_'),'l55'));
  if all( ~pmatch )
  % unknown protocol  
    Pnfailed   = cellfun(@(x) size(x,1),dismatchstr);
    Pnfailedid = find(Pnfailed == min(Pnfailed) & Pnfailed < 6);
    pname      = Vjson.ProtocolName; 
    Pfname     = ''; 
    if job.opts.verb
      fprintf('%60s : ',fname1); 
      if isempty(protocols{1,1})
        datatype = setupDatatype(pname); 
        cat_io_cprintf([0 0 0.5],sprintf('%-50s%10s ', [datatype filesep pname], ''));
      elseif min(Pnfailed) < 6 
      % close protocol
        pname0 = spm_file(protocols{Pnfailedid(1),2},'basename'); 
        Pfname     = protocols{Pnfailedid(1),2}; 
        cat_io_cprintf([1 .5 0],sprintf('%-50s%10s ', ...
          pname0, sprintf('%2d/%2d',min(Pnfailed), numel(Pnfailed))));
      else
        cat_io_cprintf([.7 0 0],sprintf('%-50s%10s ', ...
          'Unknown protocol', sprintf('%2d/%2d',min(Pnfailed), numel(Pnfailed))));
      end
    end
    match = 0;
  else
    Pnfailedid = {}; 
    pname  = spm_file(protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2},'basename'); 
    Pfname = protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2}; 
    if job.opts.verb 
      pname0 = spm_str_manip( pname, 'l50');
      fprintf('%60s : ',fname1); 
      cat_io_cprintf([0 .5 0],sprintf('%-50s%10s ',pname0,''));
    end
    match = 1; 
  end
  if isempty(protocols{1,1}), return; end

  % save non-fitting protocols
  % =======================================================================
  % This should support central integration of protocols. 
  % So we create a directory with the full and a shorted version of the protocol. 
  % The shorted version includes all fields used so far in existing protocols. 
  % In case of close protocols we copy and adapt these and create a difference table.
  % Moreover, we create a list of all cases to see how often the protocol is used.
  % =======================================================================
  
  %% Pjson,Vjson,protocols,tol,opts
  if isempty(protocols{1,1})
    proSubdir = 'undefined'; 
  elseif match
    proSubdir = 'conform';
  elseif min(Pnfailed) < 6
    proSubdir = 'semiconform';
  else    
    proSubdir = 'nonconform'; 
  end

  proname = strrep(strrep(Vjson.ProtocolName,'_','-'),' ','-'); 
  if match
    prot  = pname;
  else
    prot  = proname; 
  end
  % evaluate protocols
  datatype = setupDatatype(proname); 

  % study


  site     = setupSites(Vjson,sites);
  Pprodir  = fullfile(job.outdir{1},job.subdir,'MRprotocols',datatype,proSubdir,proname);
  if ~exist(Pprodir,'dir'), mkdir(Pprodir); end

  Vjsonc   = cleanupVjson(Vjson); 
  if match == 0 && ~isempty(protocols{1,1})
    if min(Pnfailed) < 6
    % semi-match 
      for fi = 1:numel(Pnfailedid)
        Pprosubdir = fullfile(Pprodir,spm_file(protocols{Pnfailedid(fi),2},'basename')); 
        if ~exist(Pprosubdir,'dir'), mkdir(Pprosubdir); end
        Pprofilel = fullfile(Pprosubdir, ...
          sprintf('%s-%s-long.json', spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        cat_io_json(Pprofilel,Vjsonc);

        Pprofile = fullfile(Pprosubdir, ...
          sprintf('%s-%s.json', spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        % In this case the most similar protocol is used as template
        FNpi = intersect(fieldnames(protocols{Pnfailedid(fi),3}),fieldnames(Vjson));
        for fni = 1:numel(FNpi)
          Vjson3.(FNpi{fni}) = Vjson.(FNpi{fni});
        end
        cat_io_json(Pprofile,Vjson3);

        Pprofiled = fullfile(Pprosubdir,sprintf('%s-%s-diff.csv', ...
          spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        cat_io_csv(Pprofiled,dismatchstr{pri});

        Pqc = spm_file(protocols{Pnfailedid(fi),2},'prefix','qc');
        if exist(Pqc,'file')
          Pprofileqc = fullfile(Pprosubdir,sprintf('pc%s-%s.json', ...
            spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
          copyfile(Pqc,Pprofileqc);
        end
      end
    end
  end


  % create a list of cases to see how often the protocol is used
  Pprofilec = fullfile(Pprodir,sprintf('%s-%s-list.csv',proname,site));
  if exist(Pprofilec,'file')
    Yc = cat_io_csv(Pprofilec);
    Yc{end+1,1} = Pjson;
    Yc = unique(Yc);
  else
    Yc{1,1} = Pjson; 
  end
  cat_io_csv(Pprofilec,Yc);
  
  % save the full unknown protocol
  Pprofilel = fullfile(Pprodir,sprintf('%s-%s-long.json',proname,site));
  cat_io_json(Pprofilel,Vjsonc);

  % save a shorted version as starting point   
  Pprofile  = fullfile(fileparts(Pjson),sprintf('%s-%s.json',proname,site));
  Pprofile2 = fullfile(Pprodir,sprintf('%s-%s.json',proname,site));
  FN = cellfun(@(x) fieldnames(x), protocols(:,3) , 'UniformOutput', false );
  FNM = {}; for fni=1:numel(FN), FNM = [FNM; FN{fni}]; end; FNM = unique(FNM);
  FNM = intersect(FNM,fieldnames(Vjson));
  for fni = 1:numel(FNM)
    Vjson4.(FNM{fni}) = Vjson.(FNM{fni});
  end
  Vjson4 = cleanupVjson(Vjson4); 
  cat_io_json(Pprofile, Vjson4);
  cat_io_json(Pprofile2,Vjson4);

end
% =========================================================================
function sites = getSites(Pcenterdict)
  sites = {}; 
  if ~isempty(Pcenterdict) 
    for di = 1 % to support more you would have to match the fields 
      if ~isempty(Pcenterdict{di}) 
        if ~exist(Pcenterdict{di},'file')
          error('cat_io_dcm2bids:Pcenterdict','Center dictonary file "%s" is not existing.', Pcenterdict{di});
        end
        sites = [sites; struct2cell(cat_io_json(Pcenterdict{1}))']; %#ok<AGROW>
      end
    end
  end
end
% =========================================================================
function protocols = getProtocols(Pprodictdirs)
  protocols = cell(numel(Pprodictdirs),3); pdi = 0; 
  for di = 1:numel(Pprodictdirs)
    if isempty(Pprodictdirs{di}), continue; end  
    if ~exist(Pprodictdirs{di},'dir')
      error('cat_io_dcm2bids:protocoldir','Protocol directory %d  "%s" does not exist.\n',di,Pprodictdirs{di}); 
    end
    Pprodictsubdirs = cat_vol_findfiles(Pprodictdirs{di},'*.json');
    Pprodictsubdirs(cat_io_contains(Pprodictsubdirs,[filesep 'qc'])) = []; 
    for fi = 1:numel(Pprodictsubdirs)
      pdi = pdi + 1; 
      protocols{pdi,1} = spm_file(Pprodictdirs{di},'basename'); 
      protocols{pdi,2} = Pprodictsubdirs{fi}; 
      protocols{pdi,3} = cat_io_json(Pprodictsubdirs{fi});
    end
  end
  if size(protocols,1)<1
    Pstr = ''; for di = 1:numel(Pprodictdirs), Pstr = sprintf('%s  %s\n',Pprodictdirs{di}); end 
    error('cat_io_dcm2bids:noProtocols','No Protocols found in:\n %s',Pstr);
  end
end
% =========================================================================
function site = setupSites(Vjson,sites)
  if ~isempty(sites)
    DevNam = cat_io_contains( sites(:,1) , Vjson.DeviceSerialNumber );
  else
    DevNam = 0; 
  end
  if any(DevNam)
    site = sites{find(DevNam,1,'first'),2};
  else
    site = Vjson.DeviceSerialNumber; 
  end
end
% =========================================================================
function Vjson = cleanupVjson(Vjson,type)
% remove Patient specific entries. type={'basic'*|'protocol'}.

  if ~exist('type','var'), type = 'basic'; end

  FN    = fieldnames(Vjson); 
  Vjson = rmfield(Vjson, FN(cat_io_contains(FN,'Patient')));  
  
  % maybe 2-3 levels of data reduction ?
  % (0) all DCM2NII, (1) light cleanup, (2) severe cleanup

  switch type
    case 'basic'
      % negative list, i.e. remove these entries
      RFn = {
        ... 'SeriesInstanceUID'; 'StudyInstanceUID'; 'StudyID'; 
        ... 'ProcedureStepDescription'; 
        'BodyPartExamined';
        }; 
    case 'protocol'
      % positive list, i.e. keep only these entries
      % this should be for an overview of various files
      RFn = {
        'Modality'; 
        'MagneticFieldStrength'; 
        'Manufacturer';
        'PatientPosition';
        'MRAcquisitionType';
        'ScanningSequence';
        'SequenceVariant';
        'ScanOptions';
        'ImageType';
        'NonlinearGradientCorrection';
        'SliceThickness';
        'EchoTime';
        'RepetitionTime';
        'SpoilingState';
        'InversionTime';
        'FlipAngle';
        'PartialFourier';
        'BaseResolution';
        'DiffusionScheme';
        'PhaseResolution';
        'CoilString';
        'RefLinesPE';
        'CoilCombinationMethod';
        'MatrixCoilMode';
        'PercentPhaseFOV';
        'PercentSampling';
        'PhaseEncodingSteps';
        'AcquisitionMatrixPE';
        'DwellTime';
        'ReconMatrixPE';
        'InPlanePhaseEncodingDirectionDICOM';
        'DerivedVendorReportedEchoSpacing';
        'DwellTime';
        ... 
        'PulseSequenceName'; 
        'ProcedureStepDescription';
        ...
        'PhaseEncodingDirection';
        ...'SliceTiming';
        'InPlanePhaseEncodingDirectionDICOM';
        'PixelBandwidth';
        'TotalReadoutTime';
        'EffectiveEchoSpacing';
        'DerivedVendorReportedEchoSpacing';
        'ParallelReductionFactorInPlane';
        'BandwidthPerPixelPhaseEncode';
        'EchoTrainLength';
        'MultibandAccelerationFactor';
        }; 
      RFn = setdiff( fieldnames(Vjson), RFn ); 
  end
  RFn   = intersect( RFn , fieldnames(Vjson) ); 
  Vjson = rmfield(Vjson,RFn); 
end
% =========================================================================
function V = assurePatientDCMfields(V) 

  % text fields
  FN = {'PatientSex'};
  for fni = 1:numel(FN)
    if ~isfield(V,FN{fni})
      V.(FN{fni}) = 'NA';
    end
  end

  % need this for the session
  if ~isfield(V,'ScanDate') 
    if isfield(V,'AcquisitionDateTime')
      V.ScanDate = datetime(V.AcquisitionDateTime,'Format','uuuu-MM-dd');
    else
      V.ScanDate = 'NaN';
    end
  end

  % estimate age if not given but possible 
  if ~isfield(V,'PatientAge') && isfield(V,'PatientBirthDate') && ...
      isfield(V,'ScanDate') 
    ScanDate     = datetime(V.ScanDate,'Format','uuuu-MM-dd');
    BirthDate    = datetime(V.PatientBirthDate,'Format','uuuu-MM-dd');
    V.PatientAge = char(duration( ScanDate - BirthDate, 'Format','y'));
    V.PatientAge = str2double(V.PatientAge(1:end-4)); 
  else
    V.PatientAge = NaN;
  end

  % numeric fields
  FN = {'PatientAge', 'PatientWeight'}; %, 'PatientHeight'};
  for fni = 1:numel(FN)
    if ~isfield(V,FN{fni})
      V.(FN{fni}) = NaN;
    end
  end
end
% =========================================================================
function datatype = setupDatatype(pro)
  if cat_io_contains( lower(pro) , {'mpr','mp2r','t1w','t2w','pdw','flair','inv','uni','t1','t2','tse'} )
    datatype  = 'anat'; 
  elseif cat_io_contains( lower(pro) , {'fmri','rest','bold','func','rmri','rs','tb'} )
    datatype  = 'func'; 
  elseif cat_io_contains( lower(pro) , {'dti','dwi','diff','dmri'} )
    datatype  = 'dwi';
  elseif cat_io_contains( lower(pro) , {'mrs'} )
    datatype  = 'mrs';
  elseif cat_io_contains( lower(pro) , {'field','fieldmap','mag','b0','b1'} )
    datatype  = 'fmap';
  else
    % unknown protocols
    cat_io_cprintf('red', sprintf('\n  Unkown BIDS datatype (anat/func/...) for protocol "%s"\n', lower(pro)) ) 
    datatype = 'other';
  end
end
% =========================================================================
function [suffix,pro] = setupSuffix(datatype,pro,ser)
  suffix = '';
  if strcmp( datatype , 'anat')
    if cat_io_contains( lower(pro) , {'t1map'})
      suffix = 'T1map'; 
    elseif cat_io_contains( lower(pro) , {'t2map'})
      suffix = 'T2map'; 
    elseif cat_io_contains( lower(pro) , {'t2starmap'})
      suffix = 'T2starmap'; 
    elseif cat_io_contains( lower(pro) , {'pdt2'})
      suffix = 'PDT2'; 
    elseif cat_io_contains( lower(pro) , {'mtr'})
      suffix = 'MTR'; 
    elseif cat_io_contains( lower(pro) , {'mt'})
      suffix = 'MT'; 
    elseif cat_io_contains( lower(pro) , {'mpr','mp2r','t1'})
      suffix = 'T1w'; 
    elseif cat_io_contains( lower(pro) , {'t2'}) && cat_io_contains( lower(pro) , {'star'})
      suffix = 'T2star';
    elseif cat_io_contains( lower(pro) , {'t2'})
      suffix = 'T2w';
    elseif cat_io_contains( lower(pro) , {'pd'})
      suffix = 'PDw';
    elseif cat_io_contains( lower(pro) , {'flair'})
      suffix = 'FLAIR';
    end
  
  elseif strcmp( datatype , 'dwi')
    if cat_io_contains( lower(pro) , {'sbref'}) || cat_io_contains( lower(ser) , {'sbref'})
      suffix = 'sbref';
    else
      suffix = 'dwi';
    end
    %pro = cat_io_strrep(pro,{'sbref'},{''});
  
  elseif strcmp( datatype , 'fmap')
    if cat_io_contains( lower(pro) , {'fieldmap'})
      suffix = 'fieldmap'; 
    else
      suffix = 'epi'; 
    end
  
  elseif strcmp( datatype , 'func')  
    if cat_io_contains( lower(pro) , {'sbref'}) || cat_io_contains( lower(ser) , {'sbref'})
      suffix = 'sbref'; 
    else
      suffix = 'bold'; 
    end
  end
end
% =========================================================================
function task = setupTask(datatype,pro)
  if strcmp( datatype , 'func')
    if cat_io_contains( lower(pro) , {'rs','rest'})
      task = sprintf('_task-rest');
    elseif cat_io_contains( lower(pro) , {'motor'})
      task = sprintf('_task-motor');
    elseif cat_io_contains( lower(pro) , {'lang'})
      task = sprintf('_task-language');
    elseif cat_io_contains( lower(pro) , {'stroop'})
      task = sprintf('_task-stroop');
    elseif cat_io_contains( lower(pro) , {'nback'})
      task = sprintf('_task-nback');
    elseif cat_io_contains( lower(pro) , {'memory'})
      task = sprintf('_task-memory');
    elseif cat_io_contains( lower(pro) , {'faces'})
      task = sprintf('_task-faces');
    elseif cat_io_contains( lower(pro) , {'reward'})
      task = sprintf('_task-reward');
    elseif cat_io_contains( lower(pro) , {'gambling'})
      task = sprintf('_task-gambling');
    elseif cat_io_contains( lower(pro) , {'oddball'})
      task = sprintf('_task-oddball');
    elseif cat_io_contains( lower(pro) , {'attention'})
      task = sprintf('_task-attention');
    elseif cat_io_contains( lower(pro) , {'inhibition'})
      task = sprintf('_task-inhibition');
    elseif cat_io_contains( lower(pro) , {'go'})
      task = sprintf('_task-gonogo');
    elseif cat_io_contains( lower(pro) , {'flanker'})
      task = sprintf('_task-flanker');
    elseif cat_io_contains( lower(pro) , {'social'})
      task = sprintf('_task-social');
    elseif cat_io_contains( lower(pro) , {'pain'})
      task = sprintf('_task-pain');
    elseif cat_io_contains( lower(pro) , {'auditory'})
      task = sprintf('_task-auditory');
    elseif cat_io_contains( lower(pro) , {'visual','picture'})
      task = sprintf('_task-visual');
    elseif cat_io_contains( lower(pro) , {'movie'})
      task = sprintf('_task-movie');
    elseif cat_io_contains( lower(pro) , {'audio'})
      task = sprintf('_task-audio');
    elseif cat_io_contains( lower(pro) , {'working'})
      task = sprintf('_task-workingmemory');
    else
      task = sprintf('_task-other');
    end
  else
    task = '';
  end
end
% =========================================================================
function writeScanReportTSV(Vjson,sub,site,Poutdir,BIDSsubdir,subdir)
% Create report file with scans per row for further evaluation!
% - for all files
% - for each protocol set (extraction of the main report)
return

  Preport = fullfile(Poutdir, subdir, ...
    sprintf('scanreport_site-%s_study-%s_protocoldir-%s.tsv', site, BIDSsubdir));
  if exist(Preport,'file')
    Treport = cat_io_csv(Preport, '','', struct('delimiter','\t')); 
    pidpa = find(matches( Treport(2:end,1) , sub )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Treport,1) + 1; end
    
  else
    Treport = {
      'participant_id','participant_sex','participant_age',... 
      ... 'center_id','study_id','study_date', ...
      ... scanner_model_field-strength , ...
      'mr_name','mr_para'}; 
    pidpa = 2;
  end
  Treport(pidpa,:) = {
    sub, Vjson.PatientSex, round(Vjson.PatientAge), ...
    ... site, Vjson.StudyID, Vjson.StudyID, ...     
    Vjson.ProtocolName, BIDSsubdir};
  Treport(2:end,:) = sortrows(Treport(2:end,:)); 
  
  cat_io_csv(Preport,Treport,'','',struct('delimiter','\t')); 

  % subreports

end
% =========================================================================
function writeSubReportTSV(sub,Poutdir,BIDSsubdir,subdir)
% Create report file with one subjects per row and simplified scan data
% this one could run on the scanreport-tsv files

return

  Preport = fullfile(Poutdir, subdir, ...
    sprintf('subjectreport_site-%s_study-%s_protocoldir-%s.tsv', Vjson.center, BIDSsubdir));
  if exist(Preport,'file')
    Treport = cat_io_csv(Preport, '','', struct('delimiter','\t')); 
    pidpa = find(matches( Treport(2:end,1) , sub )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Treport,1) + 1; end
  else
    Treport = {
      'participant_id','participant_sex','participant_age',... 
      'center_id','study_id','study_date','study_timepoints', ...
      ... scanner_model_field-strength , ...
      'mr_anat_t1w','mr_func','mr_dwi','mr_fmap',''}; 
      pidpa = 2;
  end
  Treport(pidpa,:) = {
    sub, Vjson.PatientSex, round(Vjson.PatientAge), ...
    Vjson.center, Vjson.StudyID, Vjson.StudyID, ...     
    '',BIDSsubdir,};
  Treport(2:end,:) = sortrows(Treport(2:end,:)); 
  
  cat_io_csv(Preport,Treport,'','',struct('delimiter','\t')); 

  % subreports
end
% =========================================================================
function writePrivateTSV(Vjson,sub,Poutdir,BIDSsubdir,subdir)
% private and participant data
% The private.tsv should contain fields that are removed in the BIDS
% processing such as the real Patient name and his birth data etc. 
% It might be saved in another directory to avoid unwanted uploading?
  Pprivate   = fullfile(Poutdir,subdir,'private.tsv');
  if exist(Pprivate,'file')
    Tprivate = cat_io_csv(Pprivate, '','', struct('delimiter','\t')); 
    pidnum   = find(cellfun(@isnumeric,Tprivate(:,2))); 
    Tprivate(pidnum,2) = cellfun(@num2str,Tprivate(pidnum,2),'UniformOutput',false); 
    pidpa    = find( matches( Tprivate(2:end,1) , sub ) ) + 1;
    if isempty(pidpa) || pidpa<=1, pidpa = size(Tprivate,1) + 1; end
  else
    Tprivate = {'participant_id','PatientID','PatientName', ...
                'PatientSex','PatientAge','PatientWeight', ...
                'PatientBirthDate','AcquisitionDateTime'}; 
    pidpa = 2;
  end
  % private table with possible subject name 
  Tprivate(pidpa,:) = {sub, Vjson.PatientID, Vjson.PatientName, ...
                       Vjson.PatientSex,Vjson.PatientAge,Vjson.PatientWeight', ...
                       Vjson.PatientBirthDate,Vjson.AcquisitionDateTime};
  Tprivate(2:end,:) = sortrows(Tprivate(2:end,:)); 
  cat_io_csv(Pprivate,Tprivate,'','',struct('delimiter','\t')); 

  % subreports
end
% =========================================================================
function writeParticipantTSV(Vjson,sub,Poutdir,BIDSsubdir)
% Create participant file (eg. OpenNeuro) 
  Pparticipants = fullfile(Poutdir,BIDSsubdir,'participants.tsv');
  if exist(Pparticipants,'file')
    Tparticipants = cat_io_csv(Pparticipants, '','', struct('delimiter','\t')); 
    pidpa = find(matches( Tparticipants(2:end,1) , sub )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Tparticipants,1) + 1; end
  else
    Tparticipants = {'participant_id','sex','age'}; 
    pidpa = 2;
  end
  Tparticipants(pidpa,:) = {sub, Vjson.PatientSex, round(Vjson.PatientAge) };
  Tparticipants(2:end,:) = sortrows(Tparticipants(2:end,:)); 

  cat_io_csv(Pparticipants,Tparticipants,'','',struct('delimiter','\t')); 
end
% =========================================================================
function Po = prepNii(Pi,opts,del)
  if opts.gzipi
    waschar = 0; 
    if ischar(Pi)
      waschar = 1; 
      Pi = cellstr(Pi); 
    end
    Po = Pi; 
    for fi = 1:numel(Pi)
      if ~strcmp( spm_file(Pi{fi},'ext'),'nii') 
        Po{fi} = spm_file(Pi{fi},'ext',''); 
        if ~exist(Po{fi},'file') && exist(Pi{fi},'file')
          try
            gunzip(Pi{fi});
          catch
            Po = ''; 
          end
          if exist('del','var') && del
            delete(Pi{fi}); 
          end
        end
      end
    end
    if waschar
      Po = char(Po); 
    end
  else
    Po = Pi; 
  end 
end
% =========================================================================
function Po = prepNiigz(Pi,opts)
  if opts.gzipi
    waschar = 0; 
    if ischar(Pi)
      waschar = 1; 
      Pi = cellstr(Pi); 
    end
    Po = Pi; 
    for fi = 1:numel(Pi)
      if ~cat_io_contains( spm_file(Pi{fi},'ext') ,'gz') 
        Po{fi} = spm_file(Pi{fi},'ext','nii.gz');
        if ~exist(Pi{fi},'file')
          gzip(Pi{fi});
          delete(Pi{fi}); 
        end
      else
        if ~exist(Pi{fi},'file') && exist(spm_file(Pi{fi},'ext',''),'file')
          gzip(spm_file(Pi{fi},'ext',''));
          delete(spm_file(Pi{fi},'ext',''));
        end
      end
    end
    if waschar
      Po = char(Po); 
    end
  else
    Po = Pi; 
  end 
end
% =========================================================================
function Pout = anomize(Pin,opts,datatype) 
  
  if (~strcmp( datatype, 'anat') && opts.anonymize < 2) || opts.anonymize == 0
    Pout = Pin; 
    return
  end

  % prepare gz-input and setup prefix for output
  Pout = spm_file(Pin,'prefix','anon_');
  
  % be lazy if the file already exist
  if exist(Pout,'file') && ~opts.rerun; return; end

  % gunzip raw input file 
  Pin = prepNii(Pin,opts,0);

  % defacing
  try
    V = spm_vol(Pin); 
    if numel(V)>1
      % create average to run defacing on this one
      Y  = spm_read_vols(V);
      Ym = cat_stat_nanmean(Y,4);
      
      %%%%%%%%%% this might bias the realignment!
      Vm = V(1); Vm.fname = spm_file(Vm.fname,'prefix','anon_');
      spm_write_vol(Vm,Ym);
      Pmsk = spm_deface( struct( 'images' , Vm.fname )); 
      Ymsk = spm_read_vols(spm_vol(Pmsk)) > 0; 
      delete(Vm.fname); delete(Pmsk);

      % apply masking
      Va   = V; 
      for vi = 1:numel(Va)
        Va(vi).fname = spm_file(Va(vi).fname,'prefix','anon_');
        Y(:,:,:,vi) = Y(:,:,:,vi) .* Ymsk;
        spm_write_vol(Va(vi),Y(:,:,:,vi));
      end

    else
      spm_deface( struct( 'images' , Pin ));
    end
  catch
    copyfile(Pin,Pout);
    %%% mark as failed?
  end

  % re-zip raw input file
  Pout = prepNiigz(Pout,opts);
  if opts.gzipi, gzip(spm_file(strrep(Pin,'.nii.gz','.nii'))); end

end
% =========================================================================
function matlabbatch = SPMsegment(Pfiles,opts)

  if exist( spm_file(Pfiles, 'prefix', 'l0'), 'file')
    return
  end

  Pfiles = prepNii(Pfiles,opts,0);

  % TPM setting 
  if 0
    PTPM = fullfile(spm('dir'),'TPM','TPM.nii'); ngaus = [1 1 2 3 4 2];
  else
    Ptpm = fullfile(spm('dir'),'TPM','mni0R1p5_TPM7blr.nii'); ngaus = [1 1 1 2 1 1 3];
  end
  Vtpm = spm_vol(Ptpm);

  % SPM segmentation 
  mi = 1; 
  matlabbatch{mi}.spm.spatial.preproc.channel.vols     = {Pfiles}; 
  matlabbatch{mi}.spm.spatial.preproc.channel.biasreg  = 0.001;
  % in general a bit more is better 
  matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm = 45;  % default = 60 
  matlabbatch{mi}.spm.spatial.preproc.channel.write    = [0 1];
  for ci = 1:numel(Vtpm)
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).tpm    = {sprintf('%s,%d',Ptpm,ci)};
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).ngaus  = ngaus(ci);
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).native = [ci<numel(Vtpm) 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).warped = [ci<numel(Vtpm) 0];
  end
  % MRF remove fine anatomical details and it is better to live with random noise/artifacts 
  matlabbatch{mi}.spm.spatial.preproc.warp.mrf     = 0.1; % default = 1 
  matlabbatch{mi}.spm.spatial.preproc.warp.cleanup = 1;
  matlabbatch{mi}.spm.spatial.preproc.warp.reg     = [0 0.0001 0.05 0.005 0.02];
  % we are now in MNI space and this performs better
  matlabbatch{mi}.spm.spatial.preproc.warp.affreg  = 'subj'; % default = 'mni' 
  matlabbatch{mi}.spm.spatial.preproc.warp.fwhm    = 0; 
  matlabbatch{mi}.spm.spatial.preproc.warp.samp    = 3;      % default = 3 
  matlabbatch{mi}.spm.spatial.preproc.warp.write   = [0 1];  % backward forward


  % create label map for quick review
  for wi = 0:1 % subjects/template space
    for li = 0:1 
      mi = mi + 1; 
      for ci = 1:numel(Vtpm)-1
        if li==1, label='l'; else, label='c'; end
        if wi
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.images{ci}(1) = ...
            cfg_dep(sprintf('Segment: wc%d Images',ci), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{ci}, '.','wc', '()',{':'}));
        else
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.images{ci}(1) = ...
            cfg_dep(sprintf('Segment: c%d Images',ci), ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
        end
      end
      if wi
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.prefix = ['\f\f\fw' label '0'];
      else
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.prefix = ['\f\f' label '0'];
      end
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.suffix          = '';
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.outdir          = {''};
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.BIDSdir         = '';
      if li == 0
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = 'i1*2 + i2*3 + i3*1';
      else 
        matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = 'round(i1)*1'; 
        for cii = 2:numel(Vtpm)-1
          matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression = [ ...
            matlabbatch{mi}.spm.tools.cat.tools.mimcalc.expression , ...
            sprintf('+round(i%d)*%d',cii,cii)]; 
        end
      end
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.var             = struct('name', {}, 'value', {});
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.dmtx    = 0;
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.mask    = 0; % masking is not working here
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.interp  = 1 - li;
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.dtype   = 2; % 2
      matlabbatch{mi}.spm.tools.cat.tools.mimcalc.options.coreg   = 0;
    end
  end


  % cleanup 
  %  - remove classes ([w]c1-c#) that were useful for the label map but are
  %    not further required
  mi = mi + 1; 
  for fi = 1:numel(Pfiles)
    for ci = 1:numel(Vtpm)-1
      if fi==1 && ci == 1
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(1) = ...
          cfg_dep(sprintf('Segment: c%d Images',ci), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
      else
        matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(end+1) = ...
          cfg_dep(sprintf('Segment: c%d Images',ci), ...
          substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
          substruct('.','tiss', '()',{ci}, '.','c', '()',{':'}));
      end                
    end
    for ci = 2:numel(Vtpm)-1 % keep GM 
      matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(end+1) = ...
        cfg_dep( sprintf('Segment: wc%d Images',ci), ...
        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
        substruct('.','tiss', '()',{ci}, '.','wc', '()',{':'}));
    end
    matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
  end

  % run batch
  if 1 
    evalc('spm_jobman(''run'',matlabbatch);');  
  else % debugging
    spm_jobman('run',matlabbatch);
  end
end
% =========================================================================
function QR = qualityRating(QM,Pprotocols)
  Pqcprotocols = spm_file(Pprotocols,'prefix','qc','ext','.json');
  if exist(Pqcprotocols,'file'), QCP = cat_io_json(Pqcprotocols); else, QCP = struct(); end

  % prepare output
  FN = fieldnames(QM); 
  for fni = 1:numel(FN)
    QR.(FN{fni}) = nan;
  end

  % quality rating
  FN = fieldnames(QCP); 
  for fni = 1:numel(FN)
    if isfield(QM,FN{fni}) && all(~(isnan(QCP.(FN{fni})))) 
      if numel(QCP.(FN{fni}))==2 && abs(diff(QCP.(FN{fni}))) > 0.001
        QR.(FN{fni}) = max(0,min(1, (QM.(FN{fni}) - QCP.(FN{fni})(1) ) / ( QCP.(FN{fni})(2)*2 - QCP.(FN{fni})(1) ) )) * 5 + .5; 
      else
        QR.(FN{fni}) = 10.5 - 10*(abs(QM.(FN{fni}) - QCP.(FN{fni})(1))<.001);
      end
    else
      QR.(FN{fni}) = nan;
    end
  end
  
  % averaging
  if isempty(FN)
    QR.SQR = nan;
  else
    fc = 2;
    QR.SQR = min(10.5,max(0.5, cat_stat_nanmean( cell2mat(struct2cell(QR)).^fc ).^fc));
  end
end
% =========================================================================
function [Pr,QM] = runQC(P, type, opts, Pprotocols)

  FNQC = {'BSM' 'WSM' 'ISR' 'NSR' 'RES'}; 

  opts.MarkColor  = cat_io_colormaps('marks+',40); 
  col2mark = @(val) opts.MarkColor(min(size(opts.MarkColor,1)-3,max(1,floor( val/9.5 * ...
    size(opts.MarkColor,1)))),:); 

  if opts.gzipi
    Pqc = spm_file(spm_file(P,'ext',''),'prefix','catDCM2NIIqc_','ext','mat');
  else
    Pqc = spm_file(P,'prefix','catDCM2NIIqc_','ext','mat');
  end

  if exist(Pqc,'file')
    [pp,ff,ee] = spm_fileparts(P); 
    Pr = cat_vol_findfiles( pp, ['r' ff ee]);
    load(Pqc,'QM'); QM.SQR = nan; 

  else
    % gunzip 
    try
      P = prepNii(P,opts,0);
    catch
      Pr = P; 
      QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);  
      cat_io_cprintf([0.5 0 0],'QC-failed\n');
      return
    end
    V = spm_vol(P); 
    vx_vol = sqrt(sum(V(1).mat(1:3,1:3).^2));
  
    % reslice in 4D data
    if numel(V) > 1
      Pr = spm_file(P,'prefix','r'); 
    else
      Pr = P; 
    end
   
    Y = single(spm_read_vols(V));
    sig75 = nan(1,size(Y,4)); 
    for di = 1:size(Y,4)
      Ydi = Y(:,:,:,di); 
      sig75(di) = prctile(Ydi(:),75); 
      clear Ydi; 
    end 
    if strcmp(type,'dwi') 
      isepi = sig75 > mean(sig75); 
    end
  
    if ~exist(Pr,'file') 
      %% realignment batch
      clear matlabbatch; 
      V = spm_vol(P);
      switch type
        case 'dwi'
          epiids = find(~isepi); 
        case 'func' 
          epiids = 1:numel(V);
        otherwise
          epiids = 1:numel(V);
      end
      for vi = 1:numel(epiids)
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{vi,1} = ... 
          sprintf('%s,%d',V(vi).fname,epiids(vi)); 
      end
      if strcmp(type,'dwi') 
        epiids = find(isepi); 
        for vi = 1:numel(epiids)
          matlabbatch{1}.spm.spatial.realign.estwrite.data{2}{vi,1} = ... 
            sprintf('%s,%d',V(vi).fname,epiids(vi)); 
        end
      end
      if opts.hrrealign
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality  = 0.95;    
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep      = 1.5;     
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp   = 4;      
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm     = 1;
      else
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality  = 0.8;    
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep      = 4;   
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp   = 1;    
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm     = 2;
      end
      matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm      = 1;
      matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap     = [0 0 0];
      matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight   = '';
      matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which    = [2 0];  
      matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap     = [0 0 0];
      matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask     = 1;
      matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix   = 'r';
      evalc('spm_jobman(''run'',matlabbatch);');  
      movefile( spm_file(P,'ext','.mat') , spm_file(P,'ext','.mat','prefix','rp_'));  
    end
  
    %% 0-none, 1-run for QC, 2-keep
    opts.sliceMotionCor = 1; 
    opts.biasCor        = 1; 
    opts.denoise        = 1; 
    for run = 1%:1 + strcmp(type,'dwi')
      Vr = spm_vol(Pr);
      if strcmp(type,'dwi')
        % sub-set-wise correction
        if run == 1
          Yr = single(spm_read_vols(Vr(~isepi)));
        else
          Yr = single(spm_read_vols(Vr(isepi)));
        end
      else
        Yr = single(spm_read_vols(Vr));
      end
  
      Ym = real(single(cat_stat_nanmean(abs(Yr),4)));     % mean image
      s0 = prctile(Ym(:),75);   % signal intensity 
      Yb = Ym > s0 & (cat_vol_grad(cat_vol_median3(Ym))./Ym < .3);
      s1 = prctile(Ym(Yb(:)),90); 
      Yb = cat_vol_morph(Yb,'lo',2); 
      Yw = cat_vol_approx( abs(Yr(:,:,:,1)) .*  Yb); 
      Yw = cat_vol_smooth3X(Yw,8./vx_vol); 
  
      %% average images
      if strcmp(type,'dwi')
        if run==1, epiids = find(~isepi); else, epiids = find(isepi); end
        Yd = zeros(size(Yr),'single'); 
        Yr = single(spm_read_vols(Vr(epiids)));
        WSM = nan(1,size(Yr,4));
        for vi = 1:size(Yr,4)
          Ya  = Yr(:,:,:,vi) ./ Yw;
    
          if opts.sliceMotionCor 
            %%
            Yw1 = Yw;
            for zi = 1:size(Yr,3)
              if zi == 1
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi+1:zi+1),3); 
              elseif zi == size(Yr,3)
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi-1:2:zi-1),3); 
              else
                Ytmp = Ya(:,:,zi) - mean(Ya(:,:,zi-1:2:zi+1),3); 
              end
              Ytmp = cat_vol_smooth3X( repmat(Ytmp,1,1,3), mean(4./vx_vol(1:2)) );
              Yw1(:,:,zi) = Ytmp(:,:,1);
            end
            Yas = Ya - Yw1 .* abs(Yw1).^.25;
            Yas(Ya==0) = 0; % defacing
            WSM(vi) = cat_stat_nanmean( (Yas(:) - Ya(:)).^2 ).^.5; 
            if opts.sliceMotionCor > 2, Ya = Yas; end
          end
    
          Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end
          
          Vrr = Vr; Vrr(epiids(vi)).fname = spm_file(Vrr(epiids(vi)).fname,'prefix','c'); 
          if opts.biasCor
            spm_write_vol(Vrr(epiids(vi)),Yas .* Yw); 
          else
            spm_write_vol(Vrr(epiids(vi)),Yas .* mean(Yw(:))); 
          end
  
          if run==1
            Yd(:,:,:,vi) = sqrt( (Ya-Yas).^2 * 2 ); % Rician noise
          end
        end
    
        if run==1
          Yn   = mean(Yd,4);
          Yns  = cat_vol_approx(cat_vol_median3(Yn)); 
        end
      end
    end
  
  
    % do measurements
    %Ym  = Ym ./ Yw; % bias corrected
    Ys  = cat_stat_nanstd(Yr,4) ./ Yw; 
    Yss = cat_vol_approx(cat_vol_median3(Ys)); 
    
    % get motion parameters
    Pm  = spm_file(P,'prefix','rp_','ext','.txt');
    if exist(Pm,'file'), rp = load(Pm); else, rp = NaN; end
  
    % final measures
    QM.BSM  = cat_stat_nanmean(cat_stat_nanstd(rp,1).^2).^.5; % average motion (between scan movement)
    QM.ISR  = cat_stat_nanstd(Yw(Yb(:))) ./ s1;   % homogeneity to signal rating
    if strcmp(type,'dwi') && exist('Yns','var') && mean(Yns(:))~=0 
      QM.NSR = cat_stat_nanmean(Yns(Yb(:)));      % noise to signal rating based on the denoising
      QM.WSM = cat_stat_nanmean(WSM.^2).^.5;      % within slice motion 
    elseif isscalar(V)
      [Ya,Ybr]  = cat_vol_resize({Yr./Yw,single(Yb)},'reduceV',1,2,32,'meanm');
      
      % measure variance in background an WM
      Yg     = cat_vol_grad(Ya)./Ya; 
      Ybg    = cat_vol_morph(cat_vol_morph(cat_vol_morph(Yg./Ya>prctile(Yg(:),75),'lc',1),'lo',2),'de',2);
      Ytis   = cat_vol_morph(cat_vol_morph(cat_vol_morph(Yg<prctile(Yg(Ybr(:)>.5),5) & Ya>.8 & Ya<1.2 & Ybr>.5,'lc',1),'lo',1),'de',1);
      [Ygr,Ybgr,Ytisr] = cat_vol_resize({Ya,Ybg,Ytis},'reduceV',1,1,32,'medianm');
      NSRbg  = cat_vol_localstat(Ygr,Ybgr>.9,2,4);  NSRbg  = cat_stat_nanmedian(NSRbg(Ybgr(:)>.9)); 
      NSRtis = cat_vol_localstat(Ygr,Ytisr>.9,2,4); NSRtis = cat_stat_nanmedian(NSRtis(Ytisr(:)>.9)); 
      QM.NSR = cat_stat_nanmean([NSRbg,NSRtis]);       % noise to signal rating based on the approximated 
  
      % estimate denoising difference
      Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end
  
      % average
      QM.NSR = max( QM.NSR , cat_stat_nanmean((Ya(:) - Yas(:)).^2).^.5 );       % noise to signal rating based on the approximated 
      QM.WSM = NaN; 
    else
      QM.NSR = cat_stat_nanmean(Yss(Yb(:)));       % noise to signal rating based on the approximated 
      QM.WSM = NaN; 
    end
    QM.vx_vol = vx_vol; 
    QM.RES    = cat_stat_nanmean(QM.vx_vol.^2).^.5;   % RESolution rating
    QM.SQR    = nan; 
  end

  if opts.verb == 1
    printvals = isempty(Pprotocols); 
    QR = qualityRating(QM,Pprotocols);
    for fni = 1:numel(FNQC)
      if (~printvals && isnan(QR.(FNQC{fni}))) || ...
         ( printvals && isnan(QM.(FNQC{fni})))
        fprintf('    -')
      else
        if printvals % original values
          fprintf('%5.1f',QM.(FNQC{fni}));
        else % ratings 
          cat_io_cprintf(col2mark(QR.(FNQC{fni})),'%5.1f',QR.(FNQC{fni}));
        end
      end
    end
    if isnan(QR.SQR)
      fprintf('    -')
    else
      cat_io_cprintf(col2mark(QR.(FNQC{fni})),'%5.1f',QR.SQR);
    end
    fprintf(' \n');
  else
    cat_io_cprintf('blue','BSM/WSM/ISR/NSR/RES: %4.2f %4.2f %4.2f %4.2f %4.2fmm\n',QM.BSM,QM.WSM,QM.ISR,QM.NSR,QM.RES); 
  end
  save(Pqc,'QM');

  P  = prepNiigz(P ,opts);
  Pr = prepNiigz(Pr,opts);

end
% =========================================================================
function [Pm,Pc0,Pwc1] = segmentanat(Pin,datatype,opts)
  
  %% segment on anonymized data!
  opts.segment = 0; 
  opts.denoise = 0;

  if (strcmp( datatype, 'anat') || opts.segment > 1) 
 
    if opts.denoise
      Pm   = spm_file(Pin,'prefix','msanlm'); 
      Pc0  = spm_file(Pin,'prefix','c0sanlm'); 
      Pwc1 = spm_file(Pin,'prefix','wc1sanlm'); 
    else
      Pm   = spm_file(Pin,'prefix','ms'); 
      Pc0  = spm_file(Pin,'prefix','c0');
      Pwc1 = spm_file(Pin,'prefix','wc1'); 
    end
    if ~exist( Pc0, 'file') 
  
      if opts.denoise && ~exist( spm_file(Pin,'prefix','sanlm'), 'file')
        cat_vol_sanlm(struct('data', Pin,'opts.verb',0));
      end
  
      if opts.denoise
        SPMsegment( spm_file(Pin,'prefix','sanlm'),opts );
      else
        SPMsegment( Pin ,opts);
      end
  
      %% QC 
      qcversion = 'cat_vol_qa201901x';
      %qcversion = 'cat_vol_qa202412';
      if opts.denoise, prefix = 'c0sanlm_'; else, prefix = 'c0'; end
      %%
      Pin2 = prepNii(spm_file({Pin},'prefix',prefix),opts,0);
      cat_vol_qa('p0',Pin2,Pin2,Pin2,'','',...
        struct('prefix',[qcversion '_'],'version',qcversion,'rerun',0,'verb',0) );
  
      %% gzipi
      if opts.gzipi
        prefixes = {'c0','wc0','wc1','l0','wl0','m','y_'};
        for pri=1:numel(prefixes)
          file = spm_file(strrep(Pin,'.nii.gz','.nii'),'prefix',prefixes{pri}) ; 
          if exist(file,'file'), gzip( file ); end
        end
      end
    end
  end
end