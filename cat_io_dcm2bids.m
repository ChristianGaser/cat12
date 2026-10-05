function cat_io_dcm2bids(job)
%cat_io_dcm2bids. Batch definition to convert DICOM to BIDS. 
%
% This batch uses DCM2NIIX to convert DICOM into NIFTI images with JSON  
% sidecars. It stores the converted data and reorganizes the output in BIDS. 
% It allows filtering for specific MRI protocols within a directory that 
% outline relevant MR parameters within a JSON file.
%
% Please check out the CAT sub-directory DCM2BIDS/DZPG3T for an example 
% that includes the definition for structural, functional, and diffusion 
% scans used by the DZPG (German Center for Mental Health, https://www.dzpg.org/).
%
% The batch also applies the SPM anonymizing routine and runs a basic image
% quality control. 
%

% ADD control file that is removed at the end. 
% if it is still there then something went wrong or the process was interupted
% this is also required in case of 4d data!

%  Minor points:    
%  . Try to replace incorrectly introduced German letters in Names?
%    >> maybe better via the csv/json read/write functions
%  . Add simpler test-protocols rather then the DZPG things in CAT as 
%    example?
%  . no acquisition duration in standard protocols !
%    but nice to have otherwise could be done by dime diff.
%
%  * BUGS: 
%     - use minimum age for participant.tsv - (min)age
%     - sites with different entries > cell rather then struct !
%  * feature:
%     - add separate protocols tables by BIDS type?
% 
%  * Separate the QC processing into a separate function/batch
%  * Connect to OHBM DGKN and Esteban 2026


% DEVELOPMENT DOCUMENTATION 
% =========================================================================
% DESIGN QUESTIONS: 
% * Use own subjectIDs ? 
%   No  - to avoid misalignment
%   Yes - to avoid human error - as number using private.tsv >> function 
%       - we will need this at least for the strong anonymizing setting
%       - a subject.tsv would also allow to integrate essential   
%         phenotypical data for groups or time-points 
%
% * Flexibility
%   - tolerance parameter? but what about ordinary variables
%   - fallback option?
%   >> maybe later, first focus on the hard setting
%
% * Maybe a flag for "only complete protocols"? 
% * Maybe a flag to for hard/soft json-para files (use all fields)
%   >> tolerance parameter
%
% * Reports 
%   
%      protocols        anat1, ..., fMRI1...
%   >> study-defined    
%   >> subject-defined
%   >> session-defined
%
% * Overview directories: 
%   * Protocols with subdires
%   * Studies (one-file with subjects and protocol-names)
%      JE2: subjects, images, StudyID
%   *** create a quality file based on the first scans (n>5)?  
%
% * where to save QC and other processed data?
%   - best would be to create an additional derivatives copy directory
%
% TODO: 
% =========================================================================
%
% * Add/extend read-me in the catDCM2NIIdb and catDM2NIIx directories 
%   as well as each result directory that explains what this is!
%
% * Recursive calls / already converted/sorted data (i.e., nii/json import)
%   >> extract BIDS structure & files
%   >> direct import 
%
% * Report Files (sub/scan) ****
%
% * if DZPG protocol ... auto mail creation?
%
% * Anonymizing 
%   * sub: define own sub-ids on level 2 or 3? ******
%   * ses: avoid scan-dates on level 1, use studyID on level 2
%   * json-checks ...
%   * combine the options in one powerful parameter
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
%   - MPM processing (not required so far)
%
%
% * Basic preprocessing, QC and anonymizing of (un)organized data?
% =========================================================================
%   * Preprocessing is a very complex issue that need general consensus! 
%     Yes, but we can make a basic suggestion that can be extended later. 
%     However, this might be irrelevant if things are done by the DZNE!
%     However, it should be an external batch that could be flagged here.
%
%   * Basic preprocessing for data overview not analysis!
%       - T1w/T2w: SPM/CAT
%       - MPMs:    hMRI
%       - dMRI:    diffusion TB
%       - fMRI:    SPM
%
%   * Basic processing to detect neurological outliers?
%     This requires a normative model. 
%     It presents a logical step in case of preprocessing!
%
%   * Basic optimization and inter-modality optimization/harmonization. 
%     This can be seen as additional part to further improve preprocessing, 
%     i.e., after preprocessing is established.
%       1) bias
%       2) denoising
%       3) contrast ("global mean")
%      (4) reorient & BB (rigid registration to MNI-space)
%
%   * How to add the data?
%     >> derivatives/TOOL/...
%
%
% * write report tables
% =========================================================================
%   - Overview for DZPG/Imaging, i.e., on line per study:
%     (Name, ID, number of (in)complete subjects, long/cross-design, modalities, avg. quality score)
%   - Overview for PI, i.e., one file per study/protocol-set with on line per subject:
%     (SubID, ...)
%   - Detailed overview, i.e., one file for all protocol-set
%
%    
%
% * protocol/QC check: 
% =========================================================================
%   - number of images/slices
%   - dcm convert error (missing dcm-files sliced etc.)
%   * add para-fields final DZPG protocol: 
%     - slice timing
%     - # slices
%     - head position ?
%   * QC-test: head-orientation ? **********
%     check AC DC offset [rmse(T) rmse(R)] 
%
%
% * Features
% =========================================================================
%   - parallelization 
%
%
% * ISSUES & BUGS:
% =========================================================================
%  * BUG: handling of multiple runs ... 
%         you can test this by using the localizer/scouts
%         - same scan but different series-number       
%
%
% * TOTEST:
% =========================================================================
%   - get BIDS test cases
%   - QC Tests (MR-ART?, DZPG samples)

 
  def.data                  = {};    % input DCM directories (add JSON/NII input later)
  def.outdir                = {pwd}; % main output directory
  def.subdir                = 'study'; % default
  def.BIDSdir               = 'BIDS';
  
  def.dicts.Pprotocoldirs   = {};    % input DCM protocol directories
  def.dicts.Pcenterdict     = {};    % dictionary for center names (otherwise scanner ID)
  def.dicts.Pstudydict      = {};    % not implemented yet
  def.dicts.Psubjdict       = {};    % not implemented yet

  def.opts.ProtocolFileName = 1;     % replace protocol name by the filenames
                                     % of the evaluation protocols under protocol directories 
  def.opts.gzipi            = 1;     % internal use of nii.gz (save disk space but maybe a bit slower) 
  def.opts.gzipe            = 1;     % external use of nii.gz (save disk space but nonoptimal for SPM processing)
  def.opts.tolerance        = 5;     % tolerance in percent for MR parameters (does not help for ordinal variables)
  def.opts.Pdcm2nii         = getDCM2NIIX; % detect DCM2NIIX installation directory
  def.opts.verb             = 1;     % 0-no, 1-yes
  def.opts.output           = 2;     % 0-no, 1-only json, 2-json+nii (input protocols), 3-full
  %def.opts.studies          = '';    % study filter >> file
  def.opts.subIDform        = 2;     % (1) use only PatientID, i.e. sub-PID 
                                     % (2) add the ScannerID, i.e. sub-SITE-PID 
                                     % (3) add also StudyName, i.e. sub-STUDY-SITE-PID 
                                     % (4) i.e. sub-SITE-STUDY-PID 
  def.opts.anonymize        = 1;     % imaging data: 0-no, 1-light-anat-only,  2-strong-all, 3-skull-strip-light?
                                     % meta data:    0-no, 1-basic-reduced,    2-strong-onlycodeds
  def.opts.deletetmp        = 0;     % remove temporary directory (better not)
  
  def.opts.hrrealign        = 0;     % highres-realignment for further use? difficult to keep this persistent 
  def.opts.rerun            = 0;     % rerun - overwrite existing
  def.opts.tableformat      = 'csv'; % csv/tsv 
  def.opts.longDBnames      = 1;     % 0 - sub-ScannerID-SubjectID/ses-seriesID/snr-seriesNR
                                     % 1 - sub-ScannerID-SubjectID/ses-seriesID-DATE/snr-seriesNR-protocolName
  def.opts.checkgzipi       = 1;     % look for unzipped data and zip it
  def.opts.runqc            = 1;     % run QC (all data)
  def.opts.preprocessing    = 1;     % run segmentation (anat) 
  def.opts.denoise          = 0;     % do denoising
  def.opts.ignoreScouts     = 1;     % remove localizer/scouts ASAP
  def.opts.protocolsubdirs  = 0;     % use additional sub-directories to separate between protocols
  def.opts.BIDSsep          = '-';   % Use extra BIDS-incompatible separator 
  
  if ~exist('job','var'), job = struct(); end
  job = cat_io_checkinopt(job,def); 
  if isempty(job.data), return; end
  
  

  %% ======================================================================
  Pdcmdir      = job.data;
  Pcenterdict  = job.dicts.Pcenterdict;
  Pprodictdirs = job.dicts.Pprotocoldirs;
  %Pstudydict   = job.dicts.Pstudydict;
  %Psubjdict    = job.dicts.Psubjdict;
  Poutdir      = fullfile(job.outdir{1},job.subdir); 
  tol          = job.opts.tolerance;
  if ~exist(Poutdir,'dir'), mkdir(Poutdir); end

  % main DB directory
  Pdbdirnam = 'catDCM2BIDSdb'; 

  % check the gzip-status of the database
  % e.g. in case of user interruptions in previous imports
  checkGzipi(fullfile(Poutdir,Pdbdirnam), job.opts.gzipi, job.opts.checkgzipi); 
  %%%%%% check file status to repair issues with broken files?

  % check BIDS output settings
  if job.opts.output > 0
    if job.opts.protocolsubdirs
      job.BIDSsubdir = strjoin( [ {'BIDS'}; spm_file(Pprodictdirs,'basename') ] , '-');
    else
      job.BIDSsubdir = 'BIDS';
    end
    Popts = fullfile(Poutdir, job.BIDSdir, job.BIDSsubdir,[Pdbdirnam '.mat']);
    [same, job] = checkPreviousSetting(Popts,job);
    if ~same, return; end 
    clear Popts; 
   
    % readme file for each protocol dir
    for pdi = 1:numel(Pprodictdirs)
      txt = {'catDCIM2BIDS BIDS export sub-directory with protocol-(mis)matching specific subdirectories.'};
      cat_io_csv(fullfile(Poutdir, job.BIDSdir, job.BIDSsubdir, 'readme.txt'),txt);
    end

    % readme file for the main BIDS dir
    txt = {'catDCIM2BIDS main BIDS export directory protocol-directory dependent subdirectories. '};
    cat_io_csv(fullfile(Poutdir, job.BIDSdir, 'readme.txt'),txt);
  end
  

  % read site-names and protocol definitions
  sites     = getSites(Pcenterdict);
  protocols = getProtocols(Pprodictdirs);  
  %studies   = getStudies(Pstudydict);  %%%%%%%%%%%%
  %subj      = getSubjects(Psubjdict);  %%%%%%%%%%%%
 

  %%%%% special case of the internal database directories: DCM2NIIX case


  % get all sub-directories
  Pdcmdirs = {}; Pdcmdirs0 = {};
  for di = 1:numel(Pdcmdir)
    Pdcmdirsdi = cat_vol_findfiles(Pdcmdir{di},'*',struct('dirs',1));
    %Pdcmdirsdi(cellfun(@(x) numel(x)>1,strfind(Pdcmdirsdi,Pdbdirnam))) = []; 
    Pdcmdirsdi(cat_io_contains(Pdcmdirsdi,Pdbdirnam)) = []; 
    if job.opts.ignoreScouts
      if di==1, cat_io_cprintf('blue','\n  Skip all localizer and scouts!\n\n'); end
      Pdcmdirsdi(cat_io_contains(lower(Pdcmdirsdi),{'localizer','scout'})) = []; 
    end
    Pdcmdirs   = [Pdcmdirs;  Pdcmdirsdi]; %#ok<AGROW>
    Pdcmdirs0  = [Pdcmdirs0; repmat(Pdcmdir(di),size(Pdcmdirsdi))]; %#ok<AGROW>
  end

  % DB subdirs
  Pmdbdir = fullfile(Poutdir,Pdbdirnam); 
  Ptbldir = fullfile(Poutdir,[Pdbdirnam '-tables']); 
  if ~exist(Pmdbdir,'dir'), mkdir(Pmdbdir); end
  if ~exist(fullfile(Ptbldir,'protocols'),'dir'), mkdir(fullfile(Ptbldir,'protocols')); end
  % readme file
  txt = {
    'catDCIM2BIDS batch main "database" directory. Here we store the imported DICOM images as NIFTI. '; 
    'The images are organized similar to BIDS with subject/session/scan: ';
    '  sub-PatientID/ses-StudyID-StudyDate/ser-SeriesNumber-ProtocolName';
    sprintf(['The sub-directory "+%s" store the original import path, where ' ...
      'only the JSON files stay to detect DICOM reimports from the same source. '],Pdbdirnam);  
    'Be careful with edits!'
  };
  cat_io_csv(fullfile(Pmdbdir,'readme.txt'),txt);
  txt = {
    'Overview directory of files in the database directory. '};
  cat_io_csv(fullfile(Ptbldir,'readme.txt'),txt);
  


  %% basic initialization that might have to be extended
  Vjson = cell(1,numel(Pdcmdirs)); sci = 0;
  sub = Vjson; ses = Vjson; datatype = Vjson; pro = Vjson; acq = Vjson; 
  sn = Vjson; task = Vjson; suffix = Vjson; site = Vjson; scankey = Vjson;
  Panon = Vjson; Pp0 = Vjson; Pwc1 = Vjson; Pdbdir = Vjson;
  Pdbdirpath = Vjson; BIDSpathd = Vjson;
  BIDSpath = Vjson; BIDSdir = Vjson; BIDSfile = Vjson;
  PID = ''; sni = 0; stime = datetime('now'); 
  QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]); 
  for fdiri = 1:numel(Pdcmdirs)
    
    % Convert DCM: 
    % =====================================================================
    Pdirfi = dcm2niix( Pdcmdirs{fdiri} , Pmdbdir, job.opts.Pdcm2nii, job.opts.gzipi, job.opts.rerun); 
%%%%% else gzip or gunzip depending all files
%%%%% special case of the internal database directories: DCM2NIIX case

  
    % get json files
    % =====================================================================
    if job.opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end
    Pjson  = cat_vol_findfiles(Pdirfi,'*.json',struct('depth',1));
    Pnii   = spm_file(Pjson,'ext',niiext); 
    for fscni = 1:numel(Pnii) % for every scan
      sci = sci + 1;

      % get nii2dcm information 
      Vjson{sci} = cat_io_json(Pjson{fscni});
      if isfield( Vjson{sci} , 'ImageTypeText')
        Vjson{sci}.ImageType = unique([Vjson{sci}.ImageType; Vjson{sci}.ImageTypeText ]); 
      end
      Vjson{sci} = assurePatientDCMfields(Vjson{sci});
      fnameparts = strsplit(spm_file(Pjson{fscni},'basename'),'='); 
      Vjson{sci}.ProtocolName = fnameparts{2};


      % create table header
      % ===================================================================
      if ~strcmp(PID, Vjson{sci}.PatientID)
        sni = sni + 1;

        if sci>1
          fprintf('%s\n',repmat('-',1,154)); 
          fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
          stime = datetime('now'); 
        end

        cat_io_cprintf([0 0.2 .8],'\n%-45s',sprintf('%4d) %s', sni, ...
          spm_str_manip(sprintf('sub-%s-%s',Vjson{sci}.DeviceSerialNumber, ...
            strrep(Vjson{sci}.PatientID,'_','')),'a38')));
        if job.opts.output
          if isempty(Pprodictdirs) || isempty(Pprodictdirs{1})
            fprintf('%15s : %-50s%10s %5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              ' ',' ','BSM','WSM','ISR','NSR','RES','SQR'); 
          else
            fprintf('%15s : %-50s%10s %5s%5s%5s%5s%5s%5s\n','dicom-protocol', ...
              'matching-protocol','mismatch','BSM','WSM','ISR','NSR','RES','SQR'); 
          end
        else
          fprintf('%15s : %-50s\n','dicom-protocol','matching-protocol'); 
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
        continue; 
      end
      % ignore 2D data
      if exist(Pnii{fscni},'file')
        Vsz = dir(Pnii{fscni});
        if isempty(Vsz) || ~isfield(Vsz,'bytes')
          cat_io_cprintf([.5 .5 .5],'%60s : ignore 2D data\n',fname1); 
          continue; 
        end
        if Vsz.bytes/1024 < 1000 
          evalc('V = spm_vol(Pnii{fscni});'); 
          if numel(V(1).dim)>2 && any(V(1).dim < 5) 
            cat_io_cprintf([.5 .5 .5],'%60s : ignore 2D data\n',fname1); 
            continue; 
          end
        end
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
      [Pdbdirpath{sci},Pdbdir{sci}] = importScan( ...
        Vjson{sci}, Pjson{fscni}, Pnii{fscni}, Pdcmdirs{fdiri}, ... 
        Pmdbdir, Ptbldir, fnameparts, niiext, job.opts );

      Pjson{fscni} = spm_file(Pjson{fscni},'path', Pdbdirpath{sci});
      Pnii{fscni}  = spm_file(Pnii{fscni}, 'path', Pdbdirpath{sci});

      % test protocols
      [match,pname,mismatchstr,Pnfailedid,Pfname] = ...
        testProtcols(Pjson{fscni},Vjson{sci},Ptbldir,protocols,tol,job,sites);
    
      % site definition 
      site{sci} = setupSites(Vjson{sci},sites);
  
      % study definition %%% need later refinement
      studynam = job.subdir;


      % main BIDS fields (subject, session, datatype, weighting, ...)
      % ===================================================================
      % subject
      bs = def.opts.BIDSsep; 
      switch job.opts.subIDform
        case 1 % sub-PID
          sub{sci} = sprintf('sub-%s', strrep(Vjson{sci}.PatientID,'_','')); %#ok<*SAGROW>
        case 2 % sub-SITE-PID
          sub{sci} = sprintf('sub-%s%s%s', site{sci}, bs, strrep(Vjson{sci}.PatientID,'_','')); %#ok<*SAGROW>
        case 3 % sub-STUDY-PID
          sub{sci} = sprintf('sub-%s%s%s', studynam, bs, strrep(Vjson{sci}.PatientID,'_','')); %#ok<*SAGROW>
        case 4 % sub-SITE-STUDY-PID
          sub{sci} = sprintf('sub-%s%s%s%s%s', site{sci}, bs, studynam, bs, strrep(Vjson{sci}.PatientID,'_','')); %#ok<*SAGROW>
        case 5 % sub-STUDY-SITE-PID
          sub{sci} = sprintf('sub-%s%s%s%s%s', studynam, bs, site{sci}, bs, strrep(Vjson{sci}.PatientID,'_','')); %#ok<*SAGROW>
      end

      % session
      if job.opts.anonymize > 1 % ses-StudyID
        ses{sci}  = sprintf('ses-%s',Vjson{sci}.StudyID);
      else % ses-date
        ses{sci}  = sprintf('ses-%s',fnameparts{3}(1:min(8,numel(fnameparts{3}))));
      end
      if numel(fnameparts{3}) > 8  &&  job.opts.anonymize==0
        ses{sci}  = sprintf('%s%s%s',ses{sci}, bs, fnameparts{3}(9:end)); 
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
          pro{sci}    = strrep(fnameparts{2},'_','-'); 
        else
          BIDSsubdir  = fullfile(job.BIDSdir, job.BIDSsubdir, 'mismatch');
          pro{sci}    = strrep(fnameparts{2},'_','-'); 
        end
      end

      % evaluate protocols
      datatype{sci} = setupDatatype(pro{sci}); 
      series = pro{sci}; 
      FN = {'ProtocolName', 'SeriesDescription', 'SequenceName'};
      for fni = 1:numel(FN)
        if isfield(Vjson{sci},FN{fni}), series = [series ' ' Vjson{sci}.(FN{fni})]; end %#ok<AGROW>
      end
      task{sci}     = setupTask(datatype{sci}, series );
      % %%%%%%%%%%%%%%%%%%%%%%%%%% refine run definition 
      % Instead of the acq-series number it would be nice to have a run variable.
      % I would like to have it all but only necessary cases, e.g. counting
      % from 1 to 2 if there are more equal scans. However, this is typically
      % only clear with the second but not the first scan ...
      sn{sci}       = sprintf('%03.0f',Vjson{sci}.SeriesNumber); % fnameparts{4})); 
      % Count scans of the same subject, session and series (e.g., fieldmap
      % magnitude/phase or multi-echo images). The keys are stored because
      % cleanupVjson later removes the Patient fields from Vjson.
      scankey{sci} = scanKey(Vjson{sci});
      run  = sum( strcmp( scankey(1:sci) , scankey{sci} ) );
      % look ahead to the next file to detect a series with further images
      runs = run;
      if numel(Pjson) > fscni
        Pjsonnext = spm_file(Pjson{fscni+1},'path',Pdirfi); % not yet imported
        if exist(Pjsonnext,'file')
          runs = runs + strcmp( scanKey(cat_io_json(Pjsonnext)) , scankey{sci} );
        end
      end
      if runs > 1
        acq{sci}    = sprintf('acq-%s%s%d%s%s', sn{sci}, bs, run, bs, pro{sci}); 
      else
        acq{sci}    = sprintf('acq-%s%s%s', sn{sci}, bs, pro{sci}); 
      end
      % get suffix 
      [suffix{sci},acq{sci}] = setupSuffix(datatype{sci}, acq{sci}, Vjson{sci}.SeriesDescription);
      

      % define BIDS naming
      % ===================================================================
      BIDSdir{sci}   = fullfile(sub{sci}, ses{sci}, datatype{sci}); 
      BIDSpath{sci}  = fullfile(Poutdir, BIDSsubdir, BIDSdir{sci}); 
      BIDSpathd{sci} = fullfile(Poutdir, BIDSsubdir, 'derivatives', BIDSdir{sci}); 
      BIDSfile{sci}  = sprintf('%s_%s_%s%s_%s.json', ...
        sub{sci}, ses{sci}, acq{sci}, task{sci}, suffix{sci});
      if ~match && ~isempty(protocols) && ~isempty(protocols{1})
        for pfi = 1:numel(Pnfailedid)
          pdir      = fullfile( BIDSpath{sci} , 'matchfiles' );
          if ~exist(pdir,'dir'), mkdir(pdir); end
          acq2      = sprintf('acq-%s%s%s', sn{sci}, bs, pro{sci}); 
          BIDSfile2 = sprintf('%s_%s_%s%s_%s.csv', ...
            sub{sci}, ses{sci}, acq2, task{sci}, suffix{sci});

          try
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','current','required'}; ...
              mismatchstr{ Pnfailedid(pfi) }] ); 
          catch
            cat_io_csv( fullfile( pdir, BIDSfile2) , [{'field','current','required'}; ...
              {mismatchstr(Pnfailedid(pfi),1),'mlt. values','mlt. values'}] ); 
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
        Panon{fscni} = anomize(Pnii{fscni}, job.opts, datatype{sci}); 
        

        % Basic QC with 
        %  - BSM (Between Scan Motion) in 4D data (average correction of the realignment)
        %  - WSM (Within Scan Motion) in 4D data (image variance between scans)
        %  - ISR (Inhomogeneity to Signal Ratio)
        %  - NSR (Noise to Signal Ratio)
        %  - RES (RMS of voxel resolution)
        if job.opts.runqc
          [Prnii{fscni},QM(fscni),Pqc{fscni}] = runQC(spm_file(Panon{fscni}), datatype{sci}, job.opts, Pfname); %#ok<AGROW>
        end
        

        % Basic preprocessing 
        % segment anatomical scan ... what to do otherwise? how to save/fill data?
        if strcmp(datatype{sci},'anat')
          [Pm{fscni}, Pp0{fscni}, Pwc1{fscni}] = segmentanat(Panon{fscni}, datatype{sci}, job.opts); %#ok<AGROW>
        end
      else
        fprintf('\n');
      end
      if job.opts.output == 0, continue; end

      
      
      

      %% create result dir and copy files
      % ===================================================================
      if job.opts.output > 0
        if ~exist(BIDSpath{sci},'dir'), mkdir(BIDSpath{sci}); end
        if job.opts.output > 1
          der = [0 0 0 1 1]; 
          pre = {'','','','catDCM2BIDSqc_','catDCM2BIDSsegus_'};
          ext = {niiext,'.bval','.bvec','.json','.json'}; 
          if job.opts.preprocessing > 1
            % copy also preprocessed files
            der = [der 1 1 1 1]; %#ok<AGROW>
            pre = [pre 'c0' 'm' 'mwc1' 'r']; %#ok<AGROW>
            ext = [ext niiext niiext niiext niiext]; %#ok<AGROW>
          end
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
              % in case of niftis we try to avoid to copy to save time  
                if job.opts.gzipe
                  if ~exist(spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei}),'file')
                    if ~job.opts.gzipi
                      gzip( file )
                      movefile( [file '.gz'], spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei})); 
                    else
                      copyfile( file, spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei})); 
                    end
                  end
                else
                  if ~exist(spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii','prefix',pre{ei}),'file')
                    if job.opts.gzipi
                      gunzip( file )
                      movefile( [file(1:end-3)], spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii','prefix',pre{ei})); 
                    else
                      copyfile( file, spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii','prefix',pre{ei})); 
                    end
                  end
                end
              else
                copyfile( file, spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext',ext{ei},'prefix',pre{ei})); 
              end
            end
          end
        end

        % anonymize and save json file
        Vjsonp = cleanupVjson(Vjson{sci},'protocol'); 
        cat_io_json( spm_file( fullfile(BIDSpath{sci},BIDSfile{sci}),'ext','.json'), Vjsonp); 


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
      writeParticipantTSV(Vjson{sci}, sub{sci}, Poutdir, BIDSsubdir, job.opts.anonymize );
      writeSubjectTSV(Vjson{sci}, sub{sci}, ses{sci}, Poutdir, BIDSsubdir, job.opts.anonymize); % session-specific data
      %writeScanReportTSV(Vjson{sci},sub{sci},site{sci},Poutdir,BIDSsubdir,'report'); %%%%%%%%%%
      BIDSsubdirprivate = strrep(BIDSsubdir, fullfile(job.BIDSdir,job.BIDSsubdir),fullfile([job.BIDSdir '-private'],job.BIDSsubdir)); 
      writePrivateTSV(Vjson{sci}, sub{sci}, fullfile(Poutdir,BIDSsubdirprivate)); 

    end
  end 
  if sci>1
    fprintf('%s\n',repmat('-',1,154)); 
    fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
  end
  fprintf('\nDCM2BIDS - import done.\n')

  if job.opts.output > 0  &&  isfield(job,'BIDSsubdir')
    % create final report form result dir
    %spm_file(Pprodictdirs,'basename')
    % study||protocol||subjectID|sex|minage|maxage||#anat/sub|#dwi/sub|#func/sub||aQR|dQR|fQR|SQR||Vgm|Vwm|Vcsf|TIV|mnFA|mn...
    writeSubReportTSV(Poutdir,BIDSsubdir,'report'); 
    % study||protocol||#subjects|%male|minage|maxage||#anat/sub|#dwi/sub|#func/sub||aQR|dQR|fQR|SQR||Vgm|Vwm|Vcsf|TIV|mnFA|mn...

    %% reports
    % ===================================================================
    % basic field: 
    %  - sub, age-range, age(mn+sd), #MRI (session/subject)
    % quality-ratings (5): 
    %  - anat/func/dwi/fmap/total
    % subject-measures: 
    %  - anat:     TIV, rGMV, rWMV, rCSFV, GMT,
    %  - func-rs:  network conectivity
    %  - dwi:      AD,FD, ...
    %
    Psubjecttable = fullfile(Poutdir,[job.BIDSdir '-report'],sprintf('report_%s.%s', job.BIDSsubdir, job.opts.tableformat));
    if exist(Psubjecttable,'file'), delete(Psubjecttable); end
    if job.opts.protocolsubdirs
      Pdirs = cat_vol_findfiles( fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir),'*',struct('depth',1,'dirs',1)); 
    else
      Pdirs = {fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir)}; 
    end
    for pdi = 1:numel(Pdirs)
      Protocol = spm_file(Pdirs{pdi},'basename'); 
      if job.opts.protocolsubdirs
        Pparticipants = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,Protocol,'participants.tsv');
      else
        Pparticipants = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,'participants.tsv');
      end
      Thdr = {'project', 'protocol','#subjects', '#sessions/subject', ...
        '#anat/subjects', '#dwi/subjects', '#func/subjects', ...
        '%males', 'mean(age)', 'std(age)', 'min(age)', 'max(age)' }; 
      if ~exist(Pparticipants,'file')
        Tparticipants = Thdr;
      else
        Tparticipants = cat_io_csv(Pparticipants,'','',struct('delimiter','\t'));
      end
      
      % look for available data
      if ~exist(Pdirs{pdi},'dir'), continue; end
      subjects = cat_vol_findfiles( Pdirs{pdi}, 'sub-*', struct('depth',1,'dirs',1));
      sessions = cat_vol_findfiles( Pdirs{pdi}, 'ses-*', struct('depth',2,'dirs',1));
      anat     = cat_vol_findfiles( Pdirs{pdi}, 'anat',  struct('depth',3,'dirs',1));
      func     = cat_vol_findfiles( Pdirs{pdi}, 'func',  struct('depth',3,'dirs',1));
      dwi      = cat_vol_findfiles( Pdirs{pdi}, 'dwi',   struct('depth',3,'dirs',1));

      % add row
      if size(Tparticipants,1)>1
        Tnewrow = {job.subdir, Protocol, numel(subjects), numel(sessions)/numel(subjects), ...
          numel(anat)/numel(subjects), numel(dwi)/numel(subjects), numel(func)/numel(subjects), ...
          mean(cellfun(@(x) strcmp(x,'M'), Tparticipants(2:end,2))), mean(cell2mat(Tparticipants(2:end,3))), ...
          std(cell2mat(Tparticipants(2:end,3))), min(cell2mat(Tparticipants(2:end,3))), max(cell2mat(Tparticipants(2:end,3))), ...
          }; 
      end
      
      %
      updateTable(Psubjecttable,Thdr,Tnewrow,pdi,0);
    % create/extend BIDS csv files ???
    end
  
  %%
    % handling GZIP in output directory
    gzipunzipOutputdata(Poutdir,job.opts)
  end

end
% =========================================================================
function key = scanKey(V)
%scanKey. Subject/session/series identifier to count images of one series.
  if ~isfield(V,'SeriesNumber') || isempty(V.SeriesNumber), key = ''; return; end
  key = sprintf('%d',V.SeriesNumber);
  FN  = {'PatientID','StudyInstanceUID'};
  for fni = 1:numel(FN)
    if isfield(V,FN{fni}), key = [key '|' char(string(V.(FN{fni})))]; end %#ok<AGROW>
  end
end
% =========================================================================
function [Pdbdirpath,Pdbdir] = importScan(Vjson, Pjson, Pnii, Pdcmdirs, Pmdbdir, Ptbldir, fnameparts, niiext, opts )
%

%datetime( Vjson.AcquisitionDateTime , 'Format','yyyyMMdd-hhmmss'))

  % DB lists 
  Pstudytable    = spm_file(fullfile(Ptbldir,'studies'),  'ext', opts.tableformat); 
  Psubjecttable  = spm_file(fullfile(Ptbldir,'subjects'), 'ext', opts.tableformat); 
  Pprotocoltable = spm_file(fullfile(Ptbldir,'protocols'),'ext', opts.tableformat); 

  if opts.longDBnames
    Pdbdir = fullfile( ...
      sprintf('sub-%s-%s',   Vjson.DeviceSerialNumber,Vjson.PatientID), ...
      sprintf('ses-%s-%s',   Vjson.StudyID, datetime( Vjson.ScanDate , 'Format','yyyyMMdd')), ...
      sprintf('snr-%04d-%s', Vjson.SeriesNumber, Vjson.ProtocolName)); 
  else
    Pdbdir = fullfile( ...
      sprintf('sub-%s-%s',   Vjson.DeviceSerialNumber,Vjson.PatientID), ...
      sprintf('ses-%s',      Vjson.StudyID), ...
      sprintf('snr-%04d',    Vjson.SeriesNumber)); 
  end

  % update 
  Pdbdirpath = fullfile(Pmdbdir, Pdbdir); 
  if ~exist( spm_file(Pnii,'path',Pdbdirpath), 'file') && ~exist(Pnii,'file')
    dcm2niix( Pdcmdirs , Pmdbdir, opts.Pdcm2nii, opts.gzipi, 1); 
    [Pdbdirpath,Pdbdir] = importScan( Vjson, Pjson, Pnii, Pdcmdirs, Pmdbdir, Ptbldir, fnameparts, niiext, opts );
    return
  end
  Pconvj = spm_file(Pjson,'path', Pdbdirpath); 
  if exist(Pconvj,'file'), Pjson = Pconvj; end
  Pconvn = spm_file(Pjson,'ext',niiext,'path', Pdbdirpath); 
  if exist(Pconvn,'file'), Pjson = Pconvn; end
  if isfield(Vjson,'ProcedureStepDescription')
    ProcedureStepDescription = getFileString(Vjson.ProcedureStepDescription);
  else
    ProcedureStepDescription = getFileString(Vjson.SeriesDescription);
  end
  % SCAN-LIST:
  %  - to avoid double entries ... but what to do ...
  %  - a flag (keep old / take new)
%%%%%%%  * this grows too strong >> internal variable + mat     
  if 0
    Thdr    = {'SeriesInstanceUID','DBpath'}; %,'IMPORTpath'}; 
    Tnewrow = {Vjson.SeriesInstanceUID, Pdbdir}; %, Pjson}; 
    updateTable(Pscantable,Thdr,Tnewrow,1,opts.rerun);
  end

  % SESSION-LIST:
  %  - to avoid double entries ... but what to do ...
  %  - a flag (keep old / take new)
%%%%%%%  * this one is usesless      
  if 0
    Thdr    = {'StudyInstanceUID','DBpath','ProcedureStepDescription'}; %,'IMPORTpath'}; 
    Tnewrow = {Vjson.StudyInstanceUID, spm_fileparts(Pdbdir), Vjson.ProcedureStepDescription}; %, Pjson}; 
    updateTable(Psessiontable,Thdr,Tnewrow,1,opts.rerun);
  end

  % SUBJECT-LIST ("constant" data):
%%%%%%%  * this one is nice but a session counter would be nice 
%%%%%%%  * maybe the scan period 
%%%%%%%  * import path (not robust :/)
  Thdr    = {'PatientID','PatientName', 'PatientSex','PatientBirthDate'}; 
  Tnewrow = {Vjson.PatientID, Vjson.PatientName, ...
             Vjson.PatientSex, Vjson.PatientBirthDate}; 
  updateTable(Psubjecttable,Thdr,Tnewrow,1,opts.rerun);
  
  % STUDY-LIST (defined by ProcedureStepDescription with additional list of subjects):
  % - ProcedureStepDescription with subjects
  Pstudysubtable = spm_file(fullfile(Ptbldir,'study_subject_lists',ProcedureStepDescription),'ext',opts.tableformat); 
  Thdr    = {'Subjects'}; 
  Tnewrow = {Vjson.PatientID}; 
  nsub    = updateTable(Pstudysubtable,Thdr,Tnewrow,1,opts.rerun);
  % - study-list with number of subjects
  Thdr    = {'ProcedureStepDescription','nSubjects'}; 
  Tnewrow = {ProcedureStepDescription,nsub}; 
  updateTable(Pstudytable,Thdr,Tnewrow,1,opts.rerun);
% with counters for anat, func, fmap, dwi ?      
  
  % PROTOCOLS-LIST:
  datatype0         = setupDatatype(Vjson.ProtocolName); 
  ProtocolName      = getFileString(strrep(fnameparts{2},'_','-'));
  %Pprotocolsubtable = spm_file(fullfile(Ptbldir,ProcedureStepDescription,datatype0,ProtocolName),'ext',opts.tableformat); 
  Pprotocolsubtable = spm_file(fullfile(Ptbldir,'protocols',datatype0,ProtocolName),'ext',opts.tableformat); 

  % save shorted protocol as json to be used as filter
  Vjson0 = cleanupVjson(Vjson,'protocol'); 
  Pjson0 = spm_file( Pprotocolsubtable,'ext','.json'); 
  ProtocolName0 = ProtocolName; 
  if ~exist(Pjson0,'file')
    cat_io_json( Pjson0, Vjson0); 
  else
    %Pjsons = cat_vol_findfiles( fullfile(Ptbldir,ProcedureStepDescription,datatype0)  , '*.json');
    Pjsons = cat_vol_findfiles( fullfile(Ptbldir,'protocols',datatype0)  , '*.json');
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
  Tnewrow = {Vjson.PatientID}; 
  nsub    = updateTable(Pprotocolsubtable,Thdr,Tnewrow,1,opts.rerun);

  % - study-list with number of subjects
  Thdr    = {'ProtocolName','nSubjects','ProcedureStepDescription', ...
             'MRAcquisitionType','ScanningSequence', 'SequenceVariant', ...
             'RepetitionTime', 'EchoTime', ...
             'FlipAngle', 'SliceThickness', 'ImageType'}; 
  Vjson0 = Vjson; 
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
  updateTable(Pprotocoltable,Thdr,Tnewrow,1,opts.rerun);


  % import files
  Pdbdirpath = fullfile(Pmdbdir, Pdbdir); 
  if ~exist(Pdbdirpath,'dir') || ~exist( spm_file(Pnii,'path',Pdbdirpath), 'file')
    if ~exist(Pdbdirpath,'dir'), mkdir(Pdbdirpath); end
    ext = {niiext,'.bval','.bvec'};
    for ei = 1:numel(ext)
      if exist(spm_file(Pjson,'ext',ext{ei}),'file') 
        if ~exist( spm_file(Pjson,'ext',ext{ei},'path', Pdbdirpath),'file')
          movefile( spm_file(Pjson,'ext',ext{ei}) , Pdbdirpath );
        else
          delete( spm_file(Pjson,'ext',ext{ei}) );
        end
      end
    end
    if ~exist( spm_file(Pjson,'path',Pdbdirpath), 'file')
      copyfile( Pjson , Pdbdirpath );
    end
  end

end
% =========================================================================
function checkGzipi(Pdbdir, gzipi, run)
% gzipunzipOutputdata. Assure the internal nifti gz-status defined by gzipi
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
      T = T & isequal(S1.(FN1{fni}),S2.(FN1{fni})); 
    elseif ischar( S1.(FN1{fni}) ) && ischar( S2.(FN1{fni}) ) 
      T = T & (strcmp(S1.(FN1{fni}),S2.(FN1{fni})));
    elseif isstruct( S1.(FN1{fni}) ) && isstruct( S2.(FN1{fni}) ) 
      T = T & structEqual(S1.(FN1{fni}),S2.(FN1{fni}));
    elseif iscellstr( S1.(FN1{fni}) ) && iscellstr( S2.(FN1{fni}) )   %#ok<ISCLSTR>
      T = T & strcmp(char(join(S1.(FN1{fni}))),char(join(S2.(FN1{fni})))); 
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
    if ~exist(fileparts(Ptable),'dir'), mkdir(fileparts(Ptable)); end
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
function [same,job] = checkPreviousSetting(Popts,job)
%checkPreviousSetting. Test if settings are identical to previous runs.
  if exist(Popts,'file')
    load(Popts,'opts','dicts');
    RMFN   = setxor(fieldnames(job.opts),{'ProtocolFileName','anonymize','gzipe','tolerance','subIDform'});
    jopts0 = rmfield(job.opts, intersect(fieldnames(job.opts),RMFN));  
    opts0  = rmfield(opts    , intersect(fieldnames(opts),RMFN));  
    same   = structEqual(jopts0,opts0) && structEqual(job.dicts,dicts);
    if ~same
      cat_io_cprintf('err', ...
        ['  Error the underlying structure of the BIDS directory does not fit to the current settings. \n', ...
         '  Use previous parameters or or change the output directory! \n\n']);
      if ~structEqual(jopts0,opts0)
        cat_io_cprintf('blue','    Old opts:\n'); 
        disp(opts); 
      end
      if ~structEqual(job.dicts,dicts)
        cat_io_cprintf('blue','    Old dicts:\n'); 
        disp(dicts); 
      end
      
      % Request user interaction 
      p = spm_input('Different BIDS setup detected. Select how to go on.',1,'m', ...
        {'Stop processing to update settings','Use old BIDS settings','Replace old BIDS directory'},[0 1 2],1); 
      if p == 1
        same      = 1;
        job.opts  = opts; 
        job.dicts = dicts; 
      elseif p == 2
        spm_figure('Clear',spm_figure('FindWin','Interactive'));
        p = spm_input('Realy replace old directory?',1,'Yes|No',[1 0],2); 
        if p
          same = 1; 
          rmdir(spm_file(Popts,'path'),'Recursive',true)
          mkdir(spm_file(Popts,'path')); 
          opts = job.opts; dicts = job.dicts; 
          save(Popts,'opts','dicts');
        end
      end
    end  
  else
    if ~exist(spm_file(Popts,'path'),'dir'), mkdir(spm_file(Popts,'path')); end
    opts = job.opts; dicts = job.dicts; 
    save(Popts,'opts','dicts');
    same = 1; 
  end
end
% =========================================================================
function gzipunzipOutputdata(Poutdir,opts)
%gzipunzipOutputdata. Assure the internal .nii.gz status defined by gzipi
  BIDSsubdirs = cat_vol_findfiles( Poutdir , 'BIDS*', struct('dirs',1)); 
  if opts.gzipi == 1  &&  opts.gzipe == 0  &&  ~isempty(BIDSsubdirs)
    % gunzip all files in the result directory
    fprintf('DCM2BIDS - gunzip BIDS output data')
    for bi = 1:numel(BIDSsubdirs)
      Pniigz = cat_vol_findfiles( BIDSsubdirs{bi} , '*.nii.gz' ); 
      for fi = 1:numel(Pniigz), gunzip(Pniigz{fi}); delete(Pniigz{fi}); end
    end
    fprintf('done.\n')
  elseif opts.gzipi == 0  &&  opts.gzipe == 1  &&  ~isempty(BIDSsubdirs)
    % gzip all files in the result directory
    fprintf('DCM2BIDS - gzip BIDS output data')
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
    %%%%%%% maybe use input path or include it in CAT?
    error('cat_io_dcm2bids:noDcm2niix', ...
      'Cannot find dcm2niix. Please install from: \n  %s\n\n', ... 
      spm_file('https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage' , ...
        'link','https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage')); 
  end
end
% =========================================================================
function tmpdir = dcm2niix( Pdcmdirs , Poutdir, Pdcm2niix, gzipi, rerun)
%dcm2niix. Import and convert DICOM data 
  tmpdir = fullfile(Poutdir,'catDCM2BIDSimportpath',Pdcmdirs); 
  if ~exist(tmpdir,'dir') || rerun
    if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
  
    if gzipi, gz = 'y'; else, gz = 'n'; end

    % convert DCM in directory
    % - use = as more unique separator 
    cmd = sprintf('%s -f "%%f=%%p=%%t=%%s" -p n -z %s -ba n -o "%s" "%s"', ...
      Pdcm2niix, gz, tmpdir,  Pdcmdirs); 
    [status,cmdout] = system(cmd); %#ok<ASGLU>
  end

  if gzipi
    P = cat_vol_findfiles( tmpdir , '*.nii' ,struct('depth',1)); 
    for fi=1:numel(P), gzip(P{fi}); delete(P{fi}); end
  else
    P = cat_vol_findfiles( tmpdir , '*.nii.gz' ,struct('depth',1)); 
    for fi=1:numel(P), gunzip(P{fi}); delete(P{fi}); end
  end
end
% =========================================================================
function [match,pname,mismatchstr,Pnfailedid,Pfname] = testProtcols(Pjson,Vjson,Ptbldir,protocols,tol,job,sites)
% test protocols

  if isempty(protocols) || isempty(protocols{1,1})
  % no given protocols
    pmatchs         = inf; 
    pmatch          = 0; 
    mismatchstr{1}  = cell(0,3); 
    mismatchcnt(1)  = 0; 
    
  else
    pmatch      = ones(1,size(protocols,1)); 
    pmatchs     = zeros(1,size(protocols,1)); 
    mismatchstr = cell(1,size(protocols,1)); 
    mismatchcnt = zeros(1,size(protocols,1)); 
    for pri = 1:size(protocols,1)
      mismatchstr{pri} = cell(0,3); 
      FNpi = fieldnames(protocols{pri,3});
     % FNpi(cat_io_contains(FNpi,{'ConsistencyInfo','PulseSequenceDetails', ...
     %   'ImageComments','SequenceName','ProtocolName','SeriesDescription'})) = []; 
      
      pmatchs(pri) = numel(FNpi); 
      for fni = 1:numel(FNpi)
        if isfield( Vjson, FNpi{fni} ) 
          if strcmp(FNpi{fni}(1:2),'x_'), continue; end % comment
          
          if ( islogical( Vjson.(FNpi{fni})) || isnumeric( Vjson.(FNpi{fni})) ) && ...
             ( islogical( protocols{pri,3}.(FNpi{fni})) || isnumeric( protocols{pri,3}.(FNpi{fni})) )
            if numel( Vjson.(FNpi{fni}) ) == numel( protocols{pri,3}.(FNpi{fni}))
              matstr =  all( Vjson.(FNpi{fni}) >= protocols{pri,3}.(FNpi{fni})*(1-tol/100) & ...
                             Vjson.(FNpi{fni}) <= protocols{pri,3}.(FNpi{fni})/(1-tol/100)); 
              %matstr = max(0,min(1, max(0,abs( Vjson.(FNpi{fni}) - protocols{pri,3}.(FNpi{fni}) ) - (protocols{pri,3}.(FNpi{fni})*tol/100) ))); 
              mismatchcnt(pri) = mismatchcnt(pri) + (1-matstr) * .5; 
              pmatchn = matstr;
            else
              mismatchcnt(pri) = mismatchcnt(pri) + 0.5;
              pmatchn = 0; 
            end
          elseif ischar( Vjson.(FNpi{fni}) ) && ischar( protocols{pri,3}.(FNpi{fni}) )
            matstr = strcmp( Vjson.(FNpi{fni}) , protocols{pri,3}.(FNpi{fni}) );
            mismatchcnt(pri) = mismatchcnt(pri) + (1-matstr);
            pmatchn = matstr; 
          elseif iscell( Vjson.(FNpi{fni}) ) && iscell( protocols{pri,3}.(FNpi{fni}) )
            %pmatchn = strcmp( char(sort(protocols{pri,3}.(FNpi{fni}(:)))) , char(sort(Vjson.(FNpi{fni}(:)))) ); 
            matstr  = 0; 
            for pmi = 1:numel( Vjson.(FNpi{fni}) )
               matstr = matstr + ...
                (1 - max([0;cat_io_contains(protocols{pri,3}.(FNpi{fni}),Vjson.(FNpi{fni})(pmi))])); 
            end
            for pmi = 1:numel( protocols{pri,3}.(FNpi{fni}) )
              matstr = matstr + ...
                (1 - max([0;cat_io_contains(Vjson.(FNpi{fni}),protocols{pri,3}.(FNpi{fni})(pmi))])); 
            end
            mismatchcnt(pri) = mismatchcnt(pri) + max(0,matstr); % - tol/2);
            pmatchn = mismatchcnt(pri) < 1;
          else
            mismatchcnt(pri) = mismatchcnt(pri) + 1;
            pmatchn = 0; 
            % need refinement ! ... image type             
          end
          pmatch(pri) = pmatch(pri) & pmatchn; 
  
          if ~pmatchn
            if iscell( Vjson.(FNpi{fni}) ) 
              tstr1 = char(join(sort(Vjson.(FNpi{fni}))));
              tstr2 = char(join(sort(protocols{pri,3}.(FNpi{fni})))); 
            else
              if isscalar(Vjson.(FNpi{fni}))
                tstr1 = Vjson.(FNpi{fni});
              elseif ischar(Vjson.(FNpi{fni})) || isstring(Vjson.(FNpi{fni}))
                tstr1 = join(Vjson.(FNpi{fni}),1); 
              else
                tstr1 = sprintf('%dx%d %s', size(Vjson.(FNpi{fni}),1), ...
                  size(Vjson.(FNpi{fni}),2), class(Vjson.(FNpi{fni}))); 
              end
              if isscalar(protocols{pri,3}.(FNpi{fni}))
                tstr2 = protocols{pri,3}.(FNpi{fni});
              elseif ischar(protocols{pri,3}.(FNpi{fni})) || isstring(protocols{pri,3}.(FNpi{fni}))
                tstr2 = join(protocols{pri,3}.(FNpi{fni}),1); 
              else
                tstr2 = sprintf('%dx%d %s', size(protocols{pri,3}.(FNpi{fni}),1), ...
                  size(protocols{pri,3}.(FNpi{fni}),2), class(protocols{pri,3}.(FNpi{fni}))); 
              end
            end
            mismatchstr{pri} = [ mismatchstr{pri}; {FNpi{fni} tstr1 tstr2 }]; 
          end

          pmatch(pri) = pmatch(pri) & pmatchn; 
        end
      end
    end
  end
%% always print session 

  % display progress
  % =======================================================================
  fname1     = sprintf('%3d) %s', Vjson.SeriesNumber, spm_str_manip(strrep(sprintf('%s_%s_%s', ...
                strrep(Vjson.PatientID, strrep(char(Vjson.ScanDate),'-',''),''), ...
                strrep(char(Vjson.ScanDate),'-',''), Vjson.ProtocolName),'__','_'),'l55'));
  if isempty(protocols) || isempty(protocols{1,1})
  % no given protocols: no match but use the DICOM protocol name
    match      = 0;
    Pnfailedid = [];
    pname      = Vjson.ProtocolName;
    Pfname     = '';
    if job.opts.verb
      fprintf('%60s : ',fname1);
      datatype = setupDatatype(pname);
      cat_io_cprintf([0 0 0.5],sprintf('%-50s%10s ', [datatype filesep pname], ''));
    end
    return
  elseif all( mismatchcnt > 0.05 )
  % unknown protocol
    Pnfailed   = mismatchcnt; %cellfun(@(x) size(x,1),mismatchstr);
    Pnfailedid = find(Pnfailed == min(Pnfailed));% & Pnfailed < 6);
    pname      = Vjson.ProtocolName;
    Pfname     = '';
    if job.opts.verb
      fprintf('%60s : ',fname1);
      pname0 = spm_file(protocols{Pnfailedid(1),2},'basename');
      if min(Pnfailed) < 1
        % very close protocol
        Pfname     = protocols{Pnfailedid(1),2}; 
        cat_io_cprintf([.5 .5 0],sprintf('%-50s%10s ',pname0,'~'));
      elseif min(Pnfailed) < 6 
      % close protocol
        Pfname     = protocols{Pnfailedid(1),2}; 
        cat_io_cprintf([1 .5 0],sprintf('%-50s%10s ', ...
          pname0, sprintf('%2.0f/%2.0f',min(Pnfailed), numel(Pnfailed))));
      else
        datatype = setupDatatype(pname);
        cat_io_cprintf([.7 0 0],sprintf('%-50s%10s ', ...
          sprintf('Unknown %s protocol',datatype), ...
          sprintf('%2.0f/%2.0f',min(Pnfailed), numel(Pnfailed))));
      end
    end
  else
    Pnfailedid = {}; 
    %pname  = spm_file(protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2},'basename'); 
    %Pfname = protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2}; 
    pname  = spm_file(protocols{find(mismatchcnt == min(mismatchcnt),1,'first'),2},'basename'); 
    Pfname = protocols{find(mismatchcnt == min(mismatchcnt),1,'first'),2}; 
    if job.opts.verb 
      pname0 = spm_str_manip( pname, 'l50');
      fprintf('%60s : ',fname1); 
      cat_io_cprintf([0 .5 0],sprintf('%-50s%10s ',pname0,''));
    end
  end
  % match = 0 for no fitting protocol, otherwise the index of the fitting protocol
  [mincnt,matchid] = min(mismatchcnt);
  match = (mincnt <= 0.05) * matchid;

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
  elseif max(0,1 - min(mismatchcnt)) > .95
    proSubdir = 'conform';
  elseif max(0,1 - min(mismatchcnt)) > .5
    proSubdir = 'accepted';
  elseif min(Pnfailed) < 6
    proSubdir = 'semiconform';
  else    
    proSubdir = 'nonconform'; 
  end

  % evaluate protocols
  proname = strrep(strrep(Vjson.ProtocolName,'_','-'),' ','-'); 
  datatype = setupDatatype(proname); 

  % study


  site     = setupSites(Vjson,sites);
  Pprodir  = fullfile(Ptbldir,'protocols-details',datatype,proSubdir,proname);
  if ~exist(Pprodir,'dir'), mkdir(Pprodir); end

  Vjsonc   = cleanupVjson(Vjson); 
  if ~match
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
        cat_io_json(Pprofile,Vjson3); clear Vjson3;

        Pprofiled = fullfile(Pprosubdir,sprintf('%s-%s-diff.csv', ...
          spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        cat_io_csv(Pprofiled,mismatchstr{Pnfailedid(fi)});

        Pqc = spm_file(protocols{Pnfailedid(fi),2},'prefix','qc');
        if exist(Pqc,'file')
          Pprofileqc = fullfile(Pprosubdir,sprintf('qc%s-%s.json', ...
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
  FNM = {}; for fni=1:numel(FN), FNM = [FNM; FN{fni}]; end; FNM = unique(FNM); %#ok<AGROW>
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
  if isempty(Pprodictdirs) || isempty(Pprodictdirs{1}), return; end
  for di = 1:numel(Pprodictdirs)
    if isempty(Pprodictdirs{di}), continue; end  
    if ~exist(Pprodictdirs{di},'dir')
      error('cat_io_dcm2bids:protocoldir','Protocol directory %d  "%s" does not exist.\n',di,Pprodictdirs{di}); 
    end
    Pprodictsubdirs = cat_vol_findfiles(Pprodictdirs{di},'*.json');
    Pprodictsubdirs(cat_io_contains(Pprodictsubdirs,[filesep 'qc'])) = []; 
    if isempty(Pprodictsubdirs) % 'cat_io_dcm2bids:noProtocols',
      cat_io_cprintf('warn','No protocols found in:\n  %s\n',Pprodictdirs{di});
    end
    for fi = 1:numel(Pprodictsubdirs)
      pdi = pdi + 1; 
      protocols{pdi,1} = spm_file(Pprodictdirs{di},'basename'); 
      protocols{pdi,2} = Pprodictsubdirs{fi}; 
      protocols{pdi,3} = cat_io_json(Pprodictsubdirs{fi});
    end
  end
  protocols = protocols(1:max(1,pdi),:); % remove empty rows of directories without protocols
  if isempty(protocols{1}) % 'cat_io_dcm2bids:noProtocols',
    cat_io_cprintf('warn','No protocols found in any protocol directory:\n  %s\n', ...
      strjoin(Pprodictdirs(:)',sprintf('\n  ')));
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
        'SeriesInstanceUID'; 'StudyInstanceUID'; 'StudyID'; 
        'InstitutionName'; 'InstitutionalDepartmentName'; 'InstitutionAddress'; 
        'StationName'; 
        ... 'ProcedureStepDescription'; 
        'BodyPartExamined';
        'AcquisitionTime'; 
        ... 'AcquisitionDateTime'
        'ShimSetting'; 'TxRefAmp'; 
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

  
  %% text fields
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

  % (re)estimate age
  if isfield(V,'PatientBirthDate') && isfield(V,'ScanDate') 
    ScanDate     = datetime(V.ScanDate,'Format','uuuu-MM-dd');
    BirthDate    = datetime(V.PatientBirthDate,'Format','uuuu-MM-dd');
    V.PatientAge = char(duration( ScanDate - BirthDate, 'Format','y'));
    V.PatientAge = str2double(V.PatientAge(1:end-4)); 
  elseif ~isfield(V,'PatientAge') 
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
% first level
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
    cat_io_cprintf('red', sprintf('\n  Unknown BIDS datatype (anat/func/...) for protocol "%s"\n', lower(pro)) ) 
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
function writePrivateTSV(Vjson,sub,Poutdir)
% private and participant data
% The private.tsv should contain fields that are removed in the BIDS
% processing such as the real Patient name and his birth data etc. 
% It might be saved in another directory to avoid unwanted uploading?
  Pprivate   = fullfile(Poutdir,sprintf('private_participants.tsv'));
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
function writeParticipantTSV(Vjson,sub,Poutdir,BIDSsubdir,anon)
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
  Tparticipants(pidpa,:) = {sub, Vjson.PatientSex, round(Vjson.PatientAge,2-anon) };
  %%%%%% here we could add subject (but not session) depending data
  Tparticipants(2:end,:) = sortrows(Tparticipants(2:end,:)); 

  cat_io_csv(Pparticipants,Tparticipants,'','',struct('delimiter','\t')); 
end
% =========================================================================
function writeSubjectTSV(Vjson,sub,ses,Poutdir,BIDSsubdir,anon)
% Create subject-wise file with session-specific data 
  Psubject = fullfile(Poutdir,BIDSsubdir,sub,sprintf('%s.tsv',sub));
  if exist(Psubject,'file')
    Tsubject = cat_io_csv(Psubject, '','', struct('delimiter','\t')); 
    pidpa = find(matches( Tsubject(2:end,1) , ses )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Tsubject,1) + 1; end
  else
    Tsubject = {'ses_id','age','weight','anat','dwi','func','SQR'}; 
    pidpa = 2;
  end
  %
  if exist( fullfile(Poutdir,BIDSsubdir,'derivatives','catDCM2BIDS',sub),'dir')
    Pqc = cat_vol_findfiles( fullfile(Poutdir,BIDSsubdir,'derivatives','catDCM2BIDS',sub), ...
      'catDCM2BIDSqc_*.json',struct('depth',3));
  else
    Pqc = cat_vol_findfiles( fullfile(Poutdir,BIDSsubdir,sub) , 'catDCM2BIDSqc_*.json',struct('depth',3));
  end
  qc = cat_io_json(Pqc); Pqc{end+1} = 'anat/nan';

  try
    SQR = nan(size(qc)); for si=1:numel(qc), SQR(si) = qc(si).qualityrating.SQR; end; SQR(end+1) = nan; 
    anat = cat_io_contains( spm_file(spm_file(Pqc,'path'),'basename'),'anat');
    dwi  = cat_io_contains( spm_file(spm_file(Pqc,'path'),'basename'),'dwi');
    func = cat_io_contains( spm_file(spm_file(Pqc,'path'),'basename'),'func');
    %%%% here one could add other subject-session specific data by a table 
    Tsubject(pidpa,:) = {ses, round(Vjson.PatientAge,2-anon), round(Vjson.PatientWeight), ...
      cat_stat_nanmean(SQR(anat)), cat_stat_nanmean(SQR(dwi)), cat_stat_nanmean(SQR(func)), ...
      cat_stat_nanmean(SQR) };
  catch
    Tsubject(pidpa,:) = {ses, round(Vjson.PatientAge,2-anon), round(Vjson.PatientWeight), ...
      nan, nan, nan, nan };
  end
  Tsubject(2:end,:) = sortrows(Tsubject(2:end,:)); 
  
  % write
  cat_io_csv(Psubject,Tsubject,'','',struct('delimiter','\t')); 
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
      if strcmp( spm_file(Pi{fi},'ext'),'gz') 
        Po{fi} = spm_file(Pi{fi},'ext',''); 
      end
      if ~exist(Po{fi},'file') && exist(Pi{fi},'file')
        try
          gunzip(Pi{fi});
        catch
          % corrupted files
          if exist(Pi{fi},'file'), delete(Pi{fi}); end
          if exist(spm_file(Po{fi},'ext','json'),'file'), delete(spm_file(Po{fi},'ext','json')); end
        end
        if exist('del','var') && del
          delete(Pi{fi}); 
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
        if ~exist(Po{fi},'file')
          gzip(Pi{fi});
          delete(Pi{fi}); 
        end
      else
        if ~exist(Po{fi},'file') && exist(spm_file(Pi{fi},'ext',''),'file')
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
%SPMsegment. Run SPM segmentation for anatomical data

  if exist( spm_file(Pfiles, 'prefix', 'l0'), 'file')
    return
  end

  Pfiles = cellstr(prepNii(Pfiles,opts,0));

  % TPM setting 
  if exist(fullfile(spm('dir'),'TPM','mni0R1p5_TPM7blr.nii'),'file')
    Ptpm = fullfile(spm('dir'),'TPM','mni0R1p5_TPM7blr.nii'); ngaus = [1 1 1 2 1 1 3];
  else
    Ptpm = fullfile(spm('dir'),'TPM','TPM.nii'); ngaus = [1 1 2 3 4 2];
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
    matlabbatch{mi}.spm.spatial.preproc.tissue(ci).warped = [ci<numel(Vtpm) ci<2]; % unmod mod
  end
  % MRF remove fine anatomical details and it is better to live with random noise/artifacts 
  matlabbatch{mi}.spm.spatial.preproc.warp.mrf     = 0.1; % default = 1 
  matlabbatch{mi}.spm.spatial.preproc.warp.cleanup = 1;
  matlabbatch{mi}.spm.spatial.preproc.warp.reg     = [0 0.0001 0.05 0.005 0.02];
  % we are now in MNI space and this performs better
  matlabbatch{mi}.spm.spatial.preproc.warp.affreg  = 'subj'; % default = 'mni' 
  matlabbatch{mi}.spm.spatial.preproc.warp.fwhm    = 0; 
  matlabbatch{mi}.spm.spatial.preproc.warp.samp    = 6;      % default = 3 
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
  %spm_jobman('run',matlabbatch); % just for debugging
  evalc('spm_jobman(''run'',matlabbatch);');  
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
  if isfield(QCP,'RES'), QCP.vx_vol = QCP.RES; FN = [FN; intersect(fieldnames(QM),'vx_vol')]; end
  for fni = 1:numel(FN)
    if isfield(QM,FN{fni}) && all(~(isnan(QCP.(FN{fni})))) && all(~isnan(QM.(FN{fni})))
      if numel(QCP.(FN{fni}))==2 && abs(diff(QCP.(FN{fni}))) > 0.001
        for i=1:numel(QM.(FN{fni}))
          QR.(FN{fni})(i) = max(0.5,min(10.5, (QM.(FN{fni})(i) - QCP.(FN{fni})(1) ) / ...
            ( QCP.(FN{fni})(2)*5/6 - QCP.(FN{fni})(1) ) * 10 + 1)); 
        end
        QR.(FN{fni}) = mean( QR.(FN{fni}).^2 )^.5;  
      else
        QR.(FN{fni}) = 10.5 - 9.5*(abs(QM.(FN{fni}) - QCP.(FN{fni})(1))<.001);
      end
    else
      QR.(FN{fni}) = nan;
    end
  end
  
  % averaging
  if isfield(QR,'vx_vol'), QR0 = rmfield(QR,'vx_vol'); end
  if isempty(FN)
    QR.SQR = nan;
  else
    fc = 2;
    QR.SQR = min(10.5,max(0.5, cat_stat_nanmean( cell2mat(struct2cell(QR0)).^fc ).^(1/fc)));
  end
end
% =========================================================================
function [Pr,QM,Pqc] = runQC(P, type, opts, Pprotocols)

  FNQC = {'BSM' 'WSM' 'ISR' 'NSR' 'RES'}; 

  opts.MarkColor  = cat_io_colormaps('marks+',40); 
  col2mark = @(val) opts.MarkColor(min(size(opts.MarkColor,1)-3,max(1,floor( val/9.5 * ...
    size(opts.MarkColor,1)))),:); 

  if opts.gzipi
    Pqc = spm_file(spm_file(strrep(P,'.nii.gz','.nii'),'ext',''),'prefix','catDCM2BIDSqc_','ext','mat');
  else
    Pqc = spm_file(strrep(P,'.nii.gz','.nii'),'prefix','catDCM2BIDSqc_','ext','mat');
  end
  Pqcj = spm_file(Pqc,'ext','json');

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
   
    if V(1).dim(3) < 5
      % spectroscopy overview image
      QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]);  
      cat_io_cprintf('blue','spectroscopy preview?\n'); 
      return
    end

    try
      Y = single(spm_read_vols(V));
    catch 
      %if cat_io_contains(spm_file(P,'basename'),'anan_')
      %if opts.gzipi && exist([P '.gz'],'file'), delete([P '.gz']); else; delete(P); end
      rmdir(fileparts(P),'Recursive',true);
      error('cat_io_dcm2bids:runQC','Remove "%s" to reimport. Please rerun!\n',P); 
      %end
    end
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
  

    %% main QC estimation 
    %  - special options: 0-none, 1-run for QC, 2-keep
    %  - this needs further work to run quick or use data permanently 
    %  
    opts.sliceMotionCor = 1; 
    opts.biasCor        = 1; 
    opts.denoise        = 1; 
    opts.tlim           = min(size(Y,4),8); 
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
 
      % basic segmentation of object/background
      % Ym  .. mean image
      % Ybb .. boundary map
      % Yb  .. brain/head/object mask - non-noise area
      % Yw  .. bias map
      % Yg  .. gradient/edge/noise map
      % g0  .. threshold in Yg
      % s1  .. signal intensity in Ym (for GM-WM intensity)
      Ym  = real(single(cat_stat_nanmean(Yr,4)));
      bbs = min(4,size(Ym)/4); 
      Ybb = true(size(Ym)); Ybb(bbs(1):end-bbs(1),bbs(2):end-bbs(2),bbs(3):end-bbs(3)) = false; 
      if  nnz(Ym(:)<0) ./ numel(Ym)  < .3 %~strcmp(type,'fmap')  &&  (
        %% typical image with mostly positive values and high intensity object
        Yg  = cat_vol_grad(Ym) ./ Ym; % this function has issues with negative non-noise structures
        g0  = prctile(Yg(:),10) * 2;
        Yb  = ~cat_vol_morph(Yg > g0*2,'ldc',2); 
        g0  = prctile(Yg(~Yb(:)),10) * 2;
        s0  = prctile(Ym(Yg(:) < g0 & Ym(:) > prctile(Ym(:),80) & Ym(:) < prctile(Ym(:),95)),80); 
        Yo  = Yg < g0  &  Ym > s0*.4  &  Ym < s0*1.5;
        Yw  = cat_vol_smooth3X(cat_vol_approx( abs(Yr(:,:,:,1)) .* Yo),4); 
        Yb  = cat_vol_morph(Yg .* (Ym./Yw) > g0,'ldc',2); 
        Ynr = cat_vol_localstat(Ym./Yw,Yb,1,4); nr = cat_stat_nanmean(Ynr(Yb(:)))*2; 
        Ywm = cat_vol_morph(Yg<g0*2 & (Ym./Yw)>.8-nr & (Ym./Yw)<1.2+nr,'ldo',0);
        s1  = prctile(Ym(Ywm(:)),90); 
      else
        %% typical fieldmap with positive and negative values
        %  - here the bias is the information (so no correction) surrounded 
        %    by heavy noise that defines the signal intensity 
        Yg  = cat_vol_localstat(Ym,true(size(Ym)),2,4) ./ ...
              cat_vol_localstat(cat_vol_smooth3X(abs(Ym),1),true(size(Ym)),1,4);
        g0  = prctile(Yg(:),10) * 2;
        Yb  = cat_vol_morph( cat_vol_morph(Yg < g0 & ~Ybb,'ldo',1), 'ldc', 4); 
        s1  = prctile(Ym(~Yb(:)),90);
        Yw  = ones(size(Yg)) * s1;
      end  
      clear Yg g0 Ybb 
      


      %% 4D data data evaluation
      if ~strcmp(type,'anat') %strcmp(type,'dwi') &&  strcmp(type,'func')
        if strcmp(type,'dwi') 
          % in case of dwi, we need to split between the EPI images and
          % direction weighted scans
          if run==1, epiids = find(~isepi); else, epiids = find(isepi); end
        else
          epiids = 1:opts.tlim;
        end

        Yd  = zeros(size(Yr),'single'); 
        Yr  = single(spm_read_vols(Vr(epiids)));
        WSM = nan(1,size(Yr,4));
        for vi = 1:min(size(Yr,4),opts.tlim)
        % for each time-point / direction apply the general bias correction
          Ya  = Yr(:,:,:,vi) ./ Yw;
    

          % If the data was realigned, we can correct for slice-wise motion
          % artifacts and interpret this as within-slice motion (WSM).
          if opts.sliceMotionCor 
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
            Yas(Ya==0) = 0; % apply defacing
            WSM(vi) = cat_stat_nanmean( (Yas(:) - Ya(:)).^2 ).^.5; 
            if opts.sliceMotionCor > 2, Ya = Yas; end
          end
    

          % Denoising of a single slice to quantify the amount of noise 
          % in the difference image
          Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end
          Vrr = Vr; Vrr(epiids(vi)).fname = spm_file(Vrr(epiids(vi)).fname,'prefix','c'); 
          if opts.biasCor
            spm_write_vol(Vrr(epiids(vi)),Yas .* Yw); 
          else
            spm_write_vol(Vrr(epiids(vi)),Yas .* mean(Yw(:))); 
          end
          %if run==1 % why only run 1?
          Yd(:,:,:,vi) = sqrt( (Ya-Yas).^2 * 2 ); % Rician noise
          %end
        end
    
        if run==1 % why only run 1?
          Yn   = mean(Yd(:,:,:,1:min(size(Yr,4),opts.tlim)), 4);
          Yns  = cat_vol_approx(cat_vol_median3(Yn)); 
        end
      end
    end
    

    % do measurements
    %Ym  = Ym ./ Yw; % bias corrected
    Ys  = cat_stat_nanstd(Yr,4) ./ Yw; 
    try
      Yss = cat_vol_approx(cat_vol_median3(Ys)); 
    catch
      Yss = cat_vol_approx(smooth3(Ys)); 
    end

    % get motion parameters
    Pm  = spm_file(P,'prefix','rp_','ext','.txt');
    if exist(Pm,'file'), rp = load(Pm); else, rp = NaN; end
  
    % final measures
    QM.BSM  = cat_stat_nanmean(cat_stat_nanstd(rp,1).^2).^.5; % average motion (between scan movement)
    QM.ISR  = cat_stat_nanstd(Yw(Yb(:))) ./ s1;   % homogeneity to signal rating
    if exist('Yns','var') && mean(Yns(:))~=0 % (strcmp(type,'dwi') || strcmp(type,'func')) && 
      QM.NSR = cat_stat_nanmean(Yns(Yb(:)));      % noise to signal rating based on the denoising
      QM.WSM = cat_stat_nanmean(WSM(1:min(numel(WSM),opts.tlim)).^2).^.5;      % within slice motion 
    

    elseif   strcmp(type,'anat')  &&  ( nnz(Ym(:)<0) ./ numel(Ym)  < .3 ) 
      [Ya,Ybr]  = cat_vol_resize({mean(Yr,4)./Yw,single(Yb)},'reduceV',vx_vol,2,32,'meanm');
      
      % measure variance in background and foreground
      Yg     = abs(cat_vol_grad(Ya+min(Ya(:)))) ./ abs(Ya+min(Ya(:))); 
      Ytis   = Yg<prctile(Yg(Ybr(:)>.5),50) & Ya>.5 & Ya<1.5 & Ybr>.5; 
      Ytis   = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ytis,'l',1),'lc',1),'de',1);
      Ytis   = Yg<prctile(Yg(Ytis(:)>.5),50) & Ya>.5 & Ya<1.5 & Ybr>.5; 
      Ytis   = cat_vol_morph(Ytis,'ldo',1);
      [Ygr,Ybgr,Ytisr] = cat_vol_resize({Ya,~Ybr,Ytis},'reduceV',1,2,32,'meanm');
      NSRbg  = cat_vol_localstat(Ygr,Ybgr>.9,2,4);  NSRbg  = cat_stat_nanmedian(NSRbg(Ybgr(:)>.9)); 
      NSRtis = cat_vol_localstat(Ygr,Ytisr>.9,2,4); NSRtis = cat_stat_nanmedian(NSRtis(Ytisr(:)>.9)); 
      QM.NSR = min([NSRbg,NSRtis]);       % noise to signal rating based on the approximated 
  
      %% estimate denoising difference
      Yas = Ya + 0; if opts.denoise, cat_sanlm(Yas,1,3); end
  
      % average
      QM.NSR = QM.NSR; % , cat_stat_nanmean((Ya(:) - Yas(:)).^2).^.5 );       % noise to signal rating based on the approximated 
      QM.WSM = cat_stat_nanmean((Ya(:) - Yas(:)).^2).^.5; 
    else
      QM.NSR = cat_stat_nanmean(Yss(Yb(:)));       % noise to signal rating based on the approximated bias field
      if exist('WSM','var') 
        QM.WSM = cat_stat_nanmean(WSM(1:opts.tlim).^2).^.5;
      else
        QM.WSM = NaN; 
      end
    end
    QM.vx_vol = vx_vol; 
    QM.RES    = cat_stat_nanmean(QM.vx_vol.^2).^.5;   % RESolution rating
    QM.SQR    = nan; 
  end

  QR = qualityRating(QM,Pprotocols);
  printvals = isempty(Pprotocols); 
  for fni = 1:numel(FNQC)
    if (~printvals && isnan(QR.(FNQC{fni}))) || ...
       ( printvals && isnan(QM.(FNQC{fni})))
      fprintf('    -')
    else
      if printvals % original values
        fprintf('%5.2f',QM.(FNQC{fni}));
      else % ratings 
        cat_io_cprintf(col2mark(QR.(FNQC{fni})),'%5.1f',QR.(FNQC{fni}));
      end
    end
  end
  if isnan(QR.SQR)
    fprintf('    -')
  else
    cat_io_cprintf(col2mark(QR.SQR),'%5.1f',QR.SQR);
  end
  fprintf(' \n');

  save(Pqc,'QM');
  cat_io_json(Pqcj,struct('qualitymeasures',QM,'qualityrating',QR)); 

  %P = prepNiigz(P ,opts);
  Pr = prepNiigz(Pr,opts);

end
% =========================================================================
function [Pm,Pc0,Pwc1,Pseg] = segmentanat(Pin,datatype,opts)
  
  %% segment on anonymized data!
  Pm = ''; Pc0 = ''; Pwc1 = ''; Pseg = ''; 

  if ~opts.preprocessing, return; end

  if (strcmp( datatype, 'anat') || opts.preprocessing > 1) 
    if opts.denoise, prefn = 'sanlm_'; else, prefn = ''; end

    Pm   = spm_file(Pin,'prefix',['m' prefn]); 
    Pc0  = spm_file(Pin,'prefix',['c0' prefn]); 
    Pwc1 = spm_file(Pin,'prefix',['wc1' prefn]); 
    Pmat = spm_file(strrep(Pin,'.nii.gz','.nii'),'prefix',prefn,'suffix','_seg8','ext','mat'); 
    Pseg = spm_file(strrep(strrep(Pin,'.nii.gz','.nii'),[filesep 'anon_'], filesep), ...
      'prefix','catDCM2BIDSsegus_','ext','json');

    if ~exist( Pc0, 'file') || opts.rerun
      % denoising
      if opts.denoise && ~exist( spm_file(Pin,'prefix',prefn), 'file')
        cat_vol_sanlm(struct('data', {{Pin}},'verb',0,'prefix',prefn));
      end
  

      % run SPM segmentation  
      if opts.denoise, Pin2 = spm_file(Pin,'prefix',prefn); else, Pin2 = Pin; end
      SPMsegment( Pin2 ,opts);


      % CAT QC 
      if opts.runqc > 1 
        qcversion = 'cat_vol_qa201901x';
        if opts.denoise, prefix = ['c0' prefn]; else, prefix = 'c0'; end
        Pin2 = prepNii(spm_file({Pin},'prefix',prefix),opts,0);
        cat_vol_qa('p0',Pin2,Pin2,Pin2,'','',...
          struct('prefix',[qcversion '_'],'version',qcversion,'rerun',opts.rerun,'verb',0) );
      end
    
   
      % evaluate segmentation  
      seg8 = load(Pmat);
      if strcmp(spm_file(Pc0,'ext'),'gz') && exist(Pc0,'file')
        try
          evalc('V = spm_vol( Pc0 );'); 
        catch
          SPMsegment( Pin2 ,opts);
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        end
      elseif exist(spm_file(Pc0,'ext',''),'file')
        try 
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        catch
          SPMsegment( Pin2 ,opts);
          evalc('V = spm_vol( spm_file(Pc0,''ext'','''' ));'); 
        end
      end
      Y      = spm_read_vols(V); 
      vx_vol = sqrt(sum(V(1).mat(1:3,1:3).^2));

      % tissue volumes
      spmus.TIV   = nnz(Y(:)>0.5) * prod(vx_vol) / 1000;
      spmus.aGMV  = nnz(round(Y(:))==2) * prod(vx_vol) / 1000; 
      spmus.aWMV  = nnz(round(Y(:))==3) * prod(vx_vol) / 1000; 
      spmus.aCSFV = nnz(round(Y(:))==1) * prod(vx_vol) / 1000; 
      spmus.rGMV  = spmus.aGMV  ./ spmus.TIV; 
      spmus.rWMV  = spmus.aWMV  ./ spmus.TIV; 
      spmus.rCSFV = spmus.aCSFV ./ spmus.TIV; 
      
      % tissue intensities
      spmus.iGM   = seg8.mn(seg8.lkp==1) * seg8.mg(seg8.lkp==1);
      spmus.iWM   = max(seg8.mn(seg8.lkp==2));
      spmus.iCSF  = min(seg8.mn(seg8.lkp==3));

      % QC like parameters
      %   ll  = log-likelihood
      %   NCR = noise-to-contrast-ratio as minimum brain tissue variance 
      %         divided by the average tissue contrast
      spmus.qc.TPMll  = seg8.ll;
      spmus.qc.NCR    = min( shiftdim(seg8.vr(seg8.lkp(:)<4).^.5) ) ./ ...
                        mean( [ abs(spmus.iGM-spmus.iWM) abs(spmus.iGM-spmus.iCSF) ...
                                abs(spmus.iWM-spmus.iCSF)] * 2 * 3); 
      cat_io_json(Pseg,spmus);
      
    
      % gzipi
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