function out = cat_io_dcm2bids(job)
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
% The import first reads only the DICOM headers (dcm2niix -b o) to store 
% the JSON sidecars of new scans with their DICOM directory in the internal
% database (catDCM2BIDSdb), whereas the images are converted only if they 
% are required (not for the JSON-only outputs). Besides DICOM directories, 
% directories of the database of the output (catDCM2BIDSdb, its sub-*, 
% ses-*, or snr-* directories) can be selected to export already imported
% data again (see cat_io_dcm2bids_import and cat_io_dcm2bids_importDCM).
%
% The functions are organized in the following files, where the functions
% of other files are called with their name as action, e.g.
% cat_io_dcm2bids_db('getDBdir',V):
%   cat_io_dcm2bids_defaults   .. default settings (job structure)
%   cat_io_dcm2bids_import     .. input sources and import stage
%   cat_io_dcm2bids_importDCM  .. DICOM import and image conversion (dcm2niix)
%   cat_io_dcm2bids_importBIDS .. import of BIDS datasets
%   cat_io_dcm2bids_db         .. internal database (layout, scan identity,
%                                 import information)
%   cat_io_dcm2bids_protocols  .. sites, protocol test, and overview tables
%   cat_io_dcm2bids_bids       .. BIDS naming and output files
%   cat_io_dcm2bids_pp         .. anonymization, segmentation, and registration
%   cat_io_dcm2bids_overlap    .. brain coverage of the scans in MNI space
%   cat_io_dcm2bids_qc         .. quality control
%   cat_io_dcm2bids_report     .. report tables
%   cat_io_dcm2bids_render     .. rendering of the BIDS images
%   cat_io_dcm2bids_helper     .. shared helper functions (text/JSON files,
%                                 paths, NIfTI)
%   cat_io_structEqual         .. comparison of structures (general CAT function)
%
%


  % output structure with the BIDS datatypes and suffixes (also used by the 
  % batch dependencies, i.e. out = cat_io_dcm2bids('outputs'))
  out = cat_io_dcm2bids_bids('BIDSoutputs'); 
  if exist('job','var') && ischar(job) && strcmp(job,'outputs'), return; end
 
  % default settings (see cat_io_dcm2bids_defaults)
  def = cat_io_dcm2bids_defaults; 

  if ~exist('job','var'), job = struct(); end
  job = cat_io_checkinopt(job,def);
  if isempty(job.data), return; end

  % GUI render choice: "norender" or "render" with the (expert) subfields
  if isfield(job.opts.render,'norender')
    job.opts.render     = rmfield(job.opts.render,'norender');
    job.opts.render.run = 0;
  elseif isfield(job.opts.render,'render')
    % direct copy, as cat_io_updateStruct ignores empty fields (e.g. no 
    % overview slices) and does not replace a matrix by a structure array
    rjob = job.opts.render.render; 
    job.opts.render = rmfield(job.opts.render,'render'); 
    fn = fieldnames(rjob); 
    for fni = 1:numel(fn), job.opts.render.(fn{fni}) = rjob.(fn{fni}); end
    job.opts.render.run = 1;
  end
  
  

  %% ======================================================================
  Pdcmdir      = job.data;
  Pcenterdict  = job.dicts.Pcenterdict;
  Pprodictdirs = job.dicts.Pprotocoldirs;
  %Pstudydict   = job.dicts.Pstudydict;
  %Psubjdict    = job.dicts.Psubjdict;
  Poutdir      = fullfile(job.outdir{1},job.subdir); 
  tol          = job.opts.tolerance;

  % main DB directory
  Pdbdirnam = 'catDCM2BIDSdb'; 
  Pmdbdir   = fullfile(Poutdir,Pdbdirnam); 

  % input: DICOM directories (import) and/or directories of the database of
  % this output (export of imported data), see cat_io_dcm2bids_import (before
  % any output is written, as directories of other databases give an error)
  Psrc = cat_io_dcm2bids_import('getInputSources',Pdcmdir, Pmdbdir, job.opts.ignoreScouts); 
  if ~exist(Poutdir,'dir'), mkdir(Poutdir); end

  % check the gzip-status of the database
  % e.g. in case of user interruptions in previous imports
  cat_io_dcm2bids_db('checkGzipi',fullfile(Poutdir,Pdbdirnam), job.opts.gzipi, job.opts.checkgzipi); 
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
      txt = {'catDCM2BIDS BIDS export sub-directory with protocol-(mis)matching specific subdirectories.'};
      cat_io_csv(fullfile(Poutdir, job.BIDSdir, job.BIDSsubdir, 'readme.txt'),txt);
    end

    % readme file for the main BIDS dir
    txt = {'catDCM2BIDS main BIDS export directory protocol-directory dependent subdirectories. '};
    cat_io_csv(fullfile(Poutdir, job.BIDSdir, 'readme.txt'),txt);
  end
  

  % read site-names and protocol definitions
  sites     = cat_io_dcm2bids_protocols('getSites',Pcenterdict);
  protocols = cat_io_dcm2bids_protocols('getProtocols',Pprodictdirs);  
  %studies   = getStudies(Pstudydict);  %%%%%%%%%%%%
  %subj      = getSubjects(Psubjdict);  %%%%%%%%%%%%
 

  % DB subdirs
  Ptbldir = fullfile(Poutdir,[Pdbdirnam '-tables']); 
  if ~exist(Pmdbdir,'dir'), mkdir(Pmdbdir); end
  if ~exist(fullfile(Ptbldir,'protocols'),'dir'), mkdir(fullfile(Ptbldir,'protocols')); end
  % readme file
  txt = {
    'catDCM2BIDS batch main "database" directory. Here we store the imported DICOM data as JSON/NIFTI. '; 
    'The scans are organized similar to BIDS with subject/session/scan: ';
    '  sub-DeviceSerialNumber-PatientID/ses-StudyID-ScanDate/snr-SeriesNumber-ProtocolName';
    ['Each scan has the JSON sidecar of dcm2niix "folder=protocol=time=series.json", the ' ...
      'import information "catDCM2BIDSimport_*.json" with the original DICOM directory, and ' ...
      'the image (converted only if required, i.e., not for JSON-only outputs). ']; 
    ''
    ['Directories of this database (the database itself, sub-*, ses-*, snr-*) can be ' ...
      'selected as input to export the imported data again (e.g. with other protocols). ']; 
    ''
    'Be careful with edits!'
  };
  cat_io_csv(fullfile(Pmdbdir,'readme.txt'),txt);
  txt = {
    'Overview directory of files in the database directory.'
    'Can be deleted in case of problems to be recreated in the next run. '};
  cat_io_csv(fullfile(Ptbldir,'readme.txt'),txt);
  


  % import the DICOM headers of new scans (fast, dcm2niix header only) and
  % the scans of BIDS datasets into the database, get the scan directories
  % of all selected scans (DICOM, BIDS, and database input), and convert
  % (or copy) their images if required (not for the JSON-only outputs 0
  % and 1), see cat_io_dcm2bids_import
  [Psnr, convmsg] = cat_io_dcm2bids_import('importSources', Psrc, Pmdbdir, job.opts);


  %% basic initialization that might have to be extended
  Vjson = cell(1,numel(Psnr)); sci = 0;
  sub = Vjson; ses = Vjson; datatype = Vjson; pro = Vjson; acq = Vjson; 
  sn = Vjson; task = Vjson; suffix = Vjson; site = Vjson; scankey = Vjson; scanid = Vjson;
  Panon = Vjson; Pp0 = Vjson; Pwc1 = Vjson; Pdbdir = Vjson;
  Pdbdirpath = Vjson; BIDSpathd = Vjson; Pdbnii = Vjson;
  BIDSpath = Vjson; BIDSdir = Vjson; BIDSfile = Vjson; BIDSsub = Vjson; Pprot = Vjson;
  PID = ''; sni = 0; stime = datetime('now'); bidsroots = {}; convver = {}; 
  srcmap   = containers.Map('KeyType','char','ValueType','any'); % new names of BIDS inputs
  intended = {}; % exported JSONs with IntendedFor references
  partdata = containers.Map('KeyType','char','ValueType','any'); % participant data
  if isfield(job,'BIDSsubdir') % BIDS output directory (protocol specific in the loop)
    BIDSsubdir = fullfile(job.BIDSdir, job.BIDSsubdir); 
  end
  QM = struct('NSR',[],'ISR',[],'RES',[],'BSM',[],'WSM',[],'vx_vol',[],'SQR',[]); 
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
      fnameparts = strsplit(spm_file(Pjson{fscni},'basename'),'='); 
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
      % derived images (also ignored by dcm2niix -i y), e.g. for imports 
      % without this option 
      if job.opts.ignoreScouts && isfield(Vjson{sci},'ImageType') && ...
          any(strcmpi(cellstr(Vjson{sci}.ImageType),'DERIVED'))
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
      cat_io_dcm2bids_protocols('updateDBtables',Vjson{sci}, Pnii{fscni}, Pdbdir{sci}, Ptbldir, fnameparts, sites, job.opts); 

      % test protocols
      [match,pname,mismatchstr,Pnfailedid,Pfname] = ...
        cat_io_dcm2bids_protocols('testProtocols',Pjson{fscni},Vjson{sci},Ptbldir,protocols,tol,job,sites);
      Pprot{sci} = Pfname; % (closest) protocol file, e.g. for the coverage mask
      if job.opts.output == 0, fprintf('\n'); continue; end % only the overview tables
    
      % site definition 
      site{sci} = cat_io_dcm2bids_protocols('setupSites',Vjson{sci},sites);
  
      % study definition %%% need later refinement
      studynam = cat_io_dcm2bids_bids('bidsLabel',job.subdir);


      % main BIDS fields (subject, session, datatype, weighting, ...)
      % ===================================================================
      % subject
      bs   = job.opts.BIDSsep; 
      NPID = cat_io_dcm2bids_bids('bidsLabel',Vjson{sci}.PatientID); 
      switch job.opts.subIDform
        case 1 % sub-PID
          sub{sci} = sprintf('sub-%s', strrep(NPID,'_','')); %#ok<*SAGROW>
        case 2 % sub-SITE-PID
          sub{sci} = sprintf('sub-%s%s%s', cat_io_dcm2bids_bids('bidsLabel',site{sci}), bs, strrep(NPID,'_','')); %#ok<*SAGROW>
        case 3 % sub-STUDY-PID
          sub{sci} = sprintf('sub-%s%s%s', studynam, bs, strrep(NPID,'_','')); %#ok<*SAGROW>
        case 4 % sub-SITE-STUDY-PID
          sub{sci} = sprintf('sub-%s%s%s%s%s', site{sci}, bs, studynam, bs, strrep(NPID,'_','')); %#ok<*SAGROW>
        case 5 % sub-STUDY-SITE-PID
          sub{sci} = sprintf('sub-%s%s%s%s%s', studynam, bs, site{sci}, bs, strrep(NPID,'_','')); %#ok<*SAGROW>
      end
      
      % session
      if job.opts.anonymize > 1 % ses-StudyID
        ses{sci}  = sprintf('ses-%s',cat_io_dcm2bids_bids('bidsLabel',Vjson{sci}.StudyID));
      elseif isempty(regexp(fnameparts{3},'^[0-9]+$','once')) 
        % session label (e.g. BIDS input without scan date)
        ses{sci}  = sprintf('ses-%s',cat_io_dcm2bids_bids('bidsLabel',fnameparts{3}));
      else % ses-date
        ses{sci}  = sprintf('ses-%s',fnameparts{3}(1:min(8,numel(fnameparts{3}))));
        if numel(fnameparts{3}) > 8  &&  job.opts.anonymize==0
          ses{sci}  = sprintf('%s%s%s',ses{sci}, bs, fnameparts{3}(9:end)); 
        end
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
          pro{sci}    = strrep(fnameparts{2},'_',job.opts.BIDSsep); 
        else
          BIDSsubdir  = fullfile(job.BIDSdir, job.BIDSsubdir, 'mismatch');
          pro{sci}    = strrep(fnameparts{2},'_',job.opts.BIDSsep); 
        end
        pro{sci} = cat_io_dcm2bids_bids('bidsLabel',pro{sci});
      end
      BIDSsub{sci} = BIDSsubdir; 

      % evaluate protocols
      datatype{sci} = cat_io_dcm2bids_bids('setupDatatype',pro{sci}); 
      series = pro{sci}; 
      FN = {'ProtocolName', 'SeriesDescription', 'SequenceName'};
      for fni = 1:numel(FN)
        if isfield(Vjson{sci},FN{fni}), series = [series ' ' Vjson{sci}.(FN{fni})]; end %#ok<AGROW>
      end
      task{sci}     = cat_io_dcm2bids_bids('setupTask',datatype{sci}, series );
      if isbids % BIDS input: datatype and task of the source dataset
        datatype{sci} = Isrc.Datatype; 
        if isfield(Isrc.Entities,'task'), task{sci} = ['_task-' Isrc.Entities.task]; else, task{sci} = ''; end
      end
      % %%%%%%%%%%%%%%%%%%%%%%%%%% refine run definition 
      % Instead of the acq-series number it would be nice to have a run variable.
      % I would like to have it all but only necessary cases, e.g. counting
      % from 1 to 2 if there are more equal scans. However, this is typically
      % only clear with the second but not the first scan ...
      sn{sci}       = sprintf('%03.0f',Vjson{sci}.SeriesNumber); % fnameparts{4})); 
      % Count scans of the same subject, session and series (e.g., fieldmap
      % magnitude/phase or multi-echo images). The keys are stored because
      % cleanupVjson later removes the Patient fields from Vjson.
      scankey{sci} = cat_io_dcm2bids_db('scanKey',Vjson{sci});
      run  = sum( strcmp( scankey(1:sci) , scankey{sci} ) );
      % look ahead to the next file to detect a series with further images
      runs = run;
      if numel(Pjson) > fscni
        Pjsonnext = spm_file(Pjson{fscni+1},'path',Pdirfi); % not yet imported
        if exist(Pjsonnext,'file')
          % count only further images of the series but no repeated scans 
          % (e.g. the same series converted from another directory)
          Vjsonnext = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pjsonnext)); 
          fnamenext = strsplit(spm_file(Pjsonnext,'basename'),'='); 
          runs = runs + ( strcmp( cat_io_dcm2bids_db('scanKey',Vjsonnext) , scankey{sci} ) && ...
            ~any(strcmp( scanid(1:sci) , cat_io_dcm2bids_db('getScanID',Vjsonnext,fnamenext) )) );
        end
      end
      if runs > 1
        acq{sci}    = sprintf('acq-%s%s%d%s%s', sn{sci}, bs, run, bs, pro{sci}); 
      else
        acq{sci}    = sprintf('acq-%s%s%s', sn{sci}, bs, pro{sci}); 
      end
      acq{sci} = cat_io_dcm2bids_bids('bidsLabel',acq{sci},1);
      % get suffix 
      [suffix{sci},acq{sci}] = cat_io_dcm2bids_bids('setupSuffix',datatype{sci}, acq{sci}, Vjson{sci}.SeriesDescription);
      if isbids, suffix{sci} = Isrc.Suffix; end % BIDS input: suffix of the source dataset
      

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

      %% (anonymized) input file for internal processing (not zipped!) 
      if job.opts.output > 1
        Panon{fscni} = cat_io_dcm2bids_pp('anonymize',Pnii{fscni}, job.opts, datatype{sci}); 
        

        % Basic QC with 
        %  - BSM (Between Scan Motion) in 4D data (average correction of the realignment)
        %  - WSM (Within Scan Motion) in 4D data (image variance between scans)
        %  - ISR (Inhomogeneity to Signal Ratio)
        %  - NSR (Noise to Signal Ratio)
        %  - RES (RMS of voxel resolution)
        if job.opts.runqc
          [Prnii{fscni},QM(fscni),Pqc{fscni}] = cat_io_dcm2bids_qc(spm_file(Panon{fscni}), datatype{sci}, job.opts, Pfname); %#ok<AGROW>
        else
          fprintf('\n'); % end of the scan line (otherwise ended by the QC values)
        end
        

        % Basic preprocessing 
        % segment anatomical scan ... what to do otherwise? how to save/fill data?
        if strcmp(datatype{sci},'anat')
          [Pm{fscni}, Pp0{fscni}, Pwc1{fscni}] = cat_io_dcm2bids_pp('segmentanat',Panon{fscni}, datatype{sci}, job.opts); %#ok<AGROW>
        end
      else
        fprintf('\n');
      end

      
      
      

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
                  if 1%~exist(spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei}),'file')
                    if ~job.opts.gzipi
                      gzip( file )
                      movefile( [file '.gz'], spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei})); 
                    else
                      copyfile( file, spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii.gz','prefix',pre{ei})); 
                    end
                  end
                else
                  if 1%~exist(spm_file( fullfile(BIDSpath0,BIDSfile{sci}),'ext','.nii','prefix',pre{ei}),'file')
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
        Vjsonp = cat_io_dcm2bids_protocols('cleanupVjson',Vjson{sci},'bids'); 
        if ~isempty(task{sci}) && ~isfield(Vjsonp,'TaskName') % required by BIDS
          Vjsonp.TaskName = regexprep(task{sci},'^_task-',''); 
        end
        % image dimensions and voxel size (if the image is available)
        if isfield(Vjson{sci},'Dimensions'), Vjsonp.Dimensions = Vjson{sci}.Dimensions; end
        if isfield(Vjson{sci},'VoxelSize'),  Vjsonp.VoxelSize  = Vjson{sci}.VoxelSize;  end
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
      %cat_io_dcm2bids_report('writeScanReportTSV',Vjson{sci},sub{sci},site{sci},Poutdir,BIDSsubdir,'report'); %%%%%%%%%%

    end

    % end of a session (the scan directories are sorted, i.e., the scans of a
    % session are processed one after another): affine registration of the 
    % session to MNI space with the orientation rating of the session and 
    % the coverage rating of its scans (see writeSessionAffines in 
    % cat_io_dcm2bids_pp) and the session data of the subject (sessions.tsv) 
    if job.opts.output > 0 && (fdiri == numel(Psnr) || ~strcmp(fileparts(Psnr{fdiri+1}),fileparts(Pdirfi)))
      sesids = find(strcmp(cellfun(@fileparts,Pdbdirpath(1:sci),'UniformOutput',false),fileparts(Pdirfi)) & ...
        ~cellfun('isempty',BIDSsub(1:sci))); 
      SOR = nan; 
      if job.opts.output > 1 && ~isempty(sesids)
        S = cat_io_dcm2bids_pp('writeSessionAffines',Pdbdirpath(sesids), Pdbnii(sesids), datatype(sesids), ...
          suffix(sesids), BIDSpathd(sesids), sub(sesids), ses(sesids), job.opts, BIDSfile(sesids), Pprot(sesids)); 
        if ~isempty(S), SOR = S(1).SOR; end
      end
      for bsd = unique(BIDSsub(sesids)) % each BIDS (protocol) directory of the session
        id = sesids(find(strcmp(BIDSsub(sesids),bsd{1}),1,'last')); 
        cat_io_dcm2bids_bids('writeSubjectTSV',Vjson{id}, sub{id}, ses{id}, Poutdir, bsd{1}, job.opts.anonymize, SOR); 
      end
    end
  end 
  if sci>1
    fprintf('%s\n',repmat('-',1,154)); 
    fprintf('%154s\n',['duration: ' char(duration(datetime('now') - stime))]); 
  end
  fprintf('\nDCM2BIDS - import done.\n')

  % IntendedFor references of BIDS inputs with the new file names
  cat_io_dcm2bids_bids('remapIntendedFor',intended, srcmap); 

  % participants.tsv and private table (once for all participants)
  if partdata.Count > 0
    P = partdata.values; P = [P{:}]; 
    for bsd = unique({P.BIDSsubdir})
      Pb = P(strcmp({P.BIDSsubdir},bsd{1})); 
      cat_io_dcm2bids_bids('writeParticipantTSV',Pb, fullfile(Poutdir,bsd{1},'participants.tsv'), job.opts.anonymize); 
      BIDSsubdirprivate = strrep(bsd{1}, fullfile(job.BIDSdir,job.BIDSsubdir),fullfile([job.BIDSdir '-private'],job.BIDSsubdir)); 
      cat_io_dcm2bids_bids('writePrivateTSV',Pb, fullfile(Poutdir,BIDSsubdirprivate,'private_participants.tsv')); 
    end
  end

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
  
  if job.opts.output > 0  &&  isfield(job,'BIDSsubdir')

    % BIDS directories (protocol-specific subdirectories)
    if job.opts.protocolsubdirs
      Pdirs = cat_vol_findfiles( fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir),'*',struct('depth',1,'dirs',1));
    else
      Pdirs = {fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir)};
    end

    % create final report from result dir (see cat_io_dcm2bids_report)
    %spm_file(Pprodictdirs,'basename')
    cat_io_dcm2bids_report('writeSubReportTSV',Poutdir,BIDSsubdir,'report');
    cat_io_dcm2bids_report('writeReportTable',Poutdir,Pdirs,job);

    % dataset description, README, and participants description of each 
    % BIDS directory
    sources = cat_io_dcm2bids_bids('getBIDSsources',unique(bidsroots)); 
    for pdi = 1:numel(Pdirs)
      if exist(Pdirs{pdi},'dir'), cat_io_dcm2bids_bids('writeDatasetFiles',Pdirs{pdi}, job, sources, unique(convver)); end
    end
  
  %%
    % handling GZIP in output directory
    cat_io_dcm2bids_bids('gzipunzipOutputdata',Poutdir,job.opts)

    % render one slice per scan for each datatype/subtype
    if job.opts.render.run && job.opts.output > 1
      fprintf('\nDCM2BIDS - Render BIDS NIFTIs.\n')

      cat_io_dcm2bids_render( fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir), ...
        fullfile(Poutdir,[job.BIDSdir '-report'],'render',job.BIDSsubdir), job.opts.render);
      
      fprintf('\nDCM2BIDS - Render BIDS NIFTIs done.\n')
    end

    % converted raw BIDS images for the batch dependencies
    out = cat_io_dcm2bids_bids('getBIDSoutputs', fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir), out ); 
  end

end
% =========================================================================
function [same,job] = checkPreviousSetting(Popts,job)
%checkPreviousSetting. Test if settings are identical to previous runs.
  if exist(Popts,'file')
    load(Popts,'opts','dicts');
    % only the fields that define the BIDS structure (independent of other, 
    % e.g. removed or new, fields of the current or previous options)
    KEEP   = {'ProtocolFileName','anonymize','gzipe','tolerance','subIDform'}; 
    jopts0 = rmfield(job.opts, setdiff(fieldnames(job.opts),KEEP));  
    opts0  = rmfield(opts    , setdiff(fieldnames(opts),KEEP));  
    same   = cat_io_structEqual(jopts0,opts0) && cat_io_structEqual(job.dicts,dicts);
    if ~same
      cat_io_cprintf('err', ...
        ['  Error the underlying structure of the BIDS directory does not fit to the current settings. \n', ...
         '  Use previous parameters or change the output directory! \n\n']);
      if ~cat_io_structEqual(jopts0,opts0)
        cat_io_cprintf('blue','    Old opts:\n'); 
        disp(opts); 
      end
      if ~cat_io_structEqual(job.dicts,dicts)
        cat_io_cprintf('blue','    Old dicts:\n'); 
        disp(dicts); 
      end
      
      % Request user interaction 
      p = spm_input('Different BIDS setup detected. Select how to go on.',1,'m', ...
        {'Stop processing to update settings','Use old BIDS settings','Replace old BIDS directory'},[0 1 2],1); 
      if p == 1
        same      = 1;
        for fni = 1:numel(KEEP) % only the BIDS-defining fields
          if isfield(opts,KEEP{fni}), job.opts.(KEEP{fni}) = opts.(KEEP{fni}); end
        end
        job.dicts = dicts; 
      elseif p == 2
        spm_figure('Clear',spm_figure('FindWin','Interactive'));
        p = spm_input('Really replace old directory?',1,'Yes|No',[1 0],2); 
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
