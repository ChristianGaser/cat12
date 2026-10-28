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
% are required (not for the JSON-only outputs, directly after the header of
% each new DICOM directory, i.e. its files are read once). Besides DICOM directories, 
% directories of the database of the output (catDCM2BIDSdb, its sub-*, 
% ses-*, or snr-* directories) can be selected to export already imported
% data again (see cat_io_dcm2bids_import and cat_io_dcm2bids_importDCM).
%
% The anonymized DICOM export (opts.dicomexport) writes the DICOM files of
% the imported scans with cleaned headers and defaced images as zip files
% of the sessions (catDCM2BIDSdb-dicom), e.g. to transfer the data to a
% central site that runs the batch again (see cat_io_dcm2bids_export).
% THE EXPORT IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE
% HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED.
%
% The BIDS and DICOM exports are extended by each run (opts.exportmode 1),
% or contain only the scans of this run (2, the directories <BIDS>-selected
% and catDCM2BIDSdb-dicom-selected are written again by each run), whereas
% the database is always extended. A BIDS directory remembers its settings,
% protocol definitions, and naming version, i.e. the user is asked to
% replace it if they changed, as renamed scans would keep their old files
% (see checkPreviousSetting).
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
%   cat_io_dcm2bids_phenotypic .. phenotypic data of the sessions (user tables)
%   cat_io_dcm2bids_parallel   .. subject packages and their (parallel) processing
%   cat_io_dcm2bids_process    .. processing of the scans of a package (protocol
%                                 test, BIDS export, QC, session preprocessing)
%   cat_io_dcm2bids_pp_anat    .. anonymization, segmentation, registration, 
%                                 and shared preprocessing functions (e.g. 
%                                 distortion correction, preprocessing level)
%   cat_io_dcm2bids_pp_dwi     .. DTI measures of DWI scans
%   cat_io_dcm2bids_pp_rsfmri  .. quality and connectivity measures of rs-fMRI
%   cat_io_dcm2bids_overlap    .. brain coverage of the scans in MNI space
%   cat_io_dcm2bids_qc         .. quality control
%   cat_io_dcm2bids_tables     .. tables of the database
%   cat_io_dcm2bids_export     .. anonymized DICOM export (cleaned headers,
%                                 defaced images)
%   cat_io_dcm2bids_anonymize  .. numbers and date shifts of the anonymization 
%                                 level 2
%   cat_io_dcm2bids_report     .. report tables
%   cat_io_dcm2bids_reportfig  .. report figure (PDF)
%   cat_io_dcm2bids_stats      .. statistics of the report figure
%   cat_io_dcm2bids_render     .. rendering of the BIDS images
%   cat_io_dcm2bids_readme     .. readme files of the output directories
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
  if isfield(job,'opts') && isfield(job.opts,'render') && ~isstruct(job.opts.render)
    job.opts = rmfield(job.opts,'render'); 
    job.opts.render = def.opts.render; 
  end
  job = cat_io_checkinopt(job,def);
  if isempty(job.data), return; end

  % anonymization level (0 - no, 1 - pseudonymized, 2 - anonymized) with 
  % its BIDS directory (e.g. BIDS-pseudonymized, also of the private and
  % report directories), see cat_io_dcm2bids_anonymize
  job.opts.anonymize = min(2,max(0,round(job.opts.anonymize))); 
  levelname   = {'','-pseudonymized','-anonymized'}; 
  job.BIDSdir = [job.BIDSdir levelname{job.opts.anonymize + 1}];
  % disclaimer of the anonymization level: always without anonymization
  % (personal data), otherwise only with phenotypic data (the user has to
  % check these data, see the readme files for the disclaimer of all levels)
  disclaimer = job.opts.anonymize == 0 || ...
    any(~cellfun('isempty',strtrim(cellstr(job.dataset.phenotypic))));
  if disclaimer, cat_io_cprintf('err','\n%s\n',cat_io_dcm2bids_readme('disclaimer',job.opts.anonymize)); end

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

  % main DB directory
  Pdbdirnam = 'catDCM2BIDSdb'; 
  Pmdbdir   = fullfile(Poutdir,Pdbdirnam); 

  % output within a DataLad dataset: files saved by DataLad are locked 
  % (read-only), i.e. the database cannot be updated by later runs
  Pds = cat_io_dcm2bids_helper('absPath',Poutdir); 
  while ~exist(fullfile(Pds,'.datalad'),'dir') && ~strcmp(fileparts(Pds),Pds), Pds = fileparts(Pds); end
  if exist(fullfile(Pds,'.datalad'),'dir')
    cat_io_cprintf('warn',['  The output directory is part of the DataLad dataset "%s". Files saved by DataLad are ' ...
      'locked (read-only), i.e. use "datalad unlock" (or "datalad run") before a run or keep the output outside ' ...
      'of the dataset.\n'],Pds); 
  end

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
    % export of only the selected scans: own BIDS directory (with its private
    % and report directories) that is written again by each run
    if job.opts.exportmode == 2
      job.BIDSsubdir = [job.BIDSsubdir '-selected'];
      Psel = [fullfile(Poutdir, {job.BIDSdir, [job.BIDSdir '-private']}, job.BIDSsubdir), ...
        fullfile(Poutdir, [job.BIDSdir '-report'], {'tables','regions','render','protocols'}, job.BIDSsubdir)];
      for si = 1:numel(Psel), if exist(Psel{si},'dir'), rmdir(Psel{si},'s'); end, end
      D = dir(fullfile(Poutdir, [job.BIDSdir '-report'], ['report_' job.BIDSsubdir '*']));
      for si = 1:numel(D), delete(fullfile(D(si).folder,D(si).name)); end
      clear Psel D si;
    end
    Popts = fullfile(Poutdir, job.BIDSdir, job.BIDSsubdir,[Pdbdirnam '.mat']);
    [same, job] = checkPreviousSetting(Popts,job);
    if ~same, return; end 
    clear Popts; 
  end
  

  % read site-names and protocol definitions
  sites     = cat_io_dcm2bids_protocols('getSites',Pcenterdict);
  protocols = cat_io_dcm2bids_protocols('getProtocols',Pprodictdirs);  
  %studies   = getStudies(Pstudydict);  %%%%%%%%%%%%
  %subj      = getSubjects(Psubjdict);  %%%%%%%%%%%%
 

  % DB subdirs
  Ptbldir = fullfile(Poutdir,[Pdbdirnam '-tables']); 
  if ~exist(Pmdbdir,'dir'), mkdir(Pmdbdir); end
  cat_io_dcm2bids_tables('initTables',Ptbldir); % remove tables of older versions
  if ~exist(fullfile(Ptbldir,'protocols'),'dir'), mkdir(fullfile(Ptbldir,'protocols')); end
  % readme files of the output directory and its subdirectories (written
  % again at the end, see cat_io_dcm2bids_readme)
  cat_io_dcm2bids_readme('writeReadmes',Poutdir,job);
  


  % import the DICOM headers of new scans (fast, dcm2niix header only) and
  % the scans of BIDS datasets into the database and get the scan 
  % directories of all selected scans (DICOM, BIDS, and database input), 
  % see cat_io_dcm2bids_import (the images of new scans are converted by the
  % import for outputs with images, missing images for each subject package,
  % see below)
  Psnr = cat_io_dcm2bids_import('importSources', Psrc, Pmdbdir, job.opts);

  % site abbreviations of the scanners: the site dictionary and for all
  % other scanners the registry of the database (i.e. no serial numbers in
  % the BIDS labels, see siteRegistry in cat_io_dcm2bids_protocols)
  sites = cat_io_dcm2bids_protocols('siteRegistry',Pmdbdir,job.opts,sites); 

  % numbers and date offsets of the subjects for the anonymization level 2
  % (see cat_io_dcm2bids_anonymize)
  if job.opts.anonymize > 1, anon = cat_io_dcm2bids_anonymize('registry',Pmdbdir); else, anon = []; end

  % level 2: scanners of this run without siteID (site label only with the
  % manufacturer and field strength, see setupSites in cat_io_dcm2bids_protocols)
  if job.opts.anonymize > 1
    ser = unique(cellfun(@(p) char(regexp(spm_file(fileparts(fileparts(p)),'basename'),'(?<=^sub-).+?(?=-)','match','once')), ...
      Psnr,'UniformOutput',false)); 
    nosid = ser(cellfun(@(x) isempty(cat_io_dcm2bids_protocols('siteID',x,sites)) && ... % (not for the site label of an anonymized export)
      isempty(regexp(x,'^[A-Za-z0-9]{0,16}[A-Z]\d{2}T$','once')),ser));
    if ~isempty(nosid)
      cat_io_cprintf('warn',['  Anonymization level 2: %d scanner(s) without siteID, i.e. their site labels contain only the ' ...
        'manufacturer and field strength (e.g. SI30T) and the subject labels are only unique within this database. ' ...
        'Define the "siteID" of the scanners in the site dictionary (or in catDCM2BIDSdb/catDCM2BIDSsites.json) to ' ...
        'pool the anonymized data of several sites.\n'],numel(nosid)); 
    end
  end

  % phenotypic data of the sessions (user tables matched to the sessions of
  % the database, exported to the sessions files, see cat_io_dcm2bids_phenotypic)
  cat_io_dcm2bids_phenotypic('importTables',job.dataset.phenotypic,job.dataset.phenotypickeys,Pmdbdir,Ptbldir);

  % age limit of the sessions (opts.agelimit): sessions outside the limit
  % or without age (e.g. age 0 of an undefined age) are not processed and
  % are listed in the table excluded_sessions of the database (see ageLimit)
  Psnr = ageLimit(Psnr,Pmdbdir,Ptbldir,job.opts);

  % processing of the scans in subject packages (all scans of a subject,
  % see cat_io_dcm2bids_parallel): conversion of the images if required 
  % (not for the JSON-only outputs 0 and 1), protocol test, BIDS export, 
  % QC, and preprocessing of the sessions (see cat_io_dcm2bids_process), 
  % by job.nproc parallel processes (0 - in this process)
  % work directory of the unzipped images on a faster disk (e.g. the 
  % internal SSD for a database on an external hard disk, see workSetup in
  % cat_io_dcm2bids_helper)
  if job.opts.output > 1
    job.opts.workroot = cat_io_dcm2bids_helper('workSetup',Pmdbdir,job.opts,job.nproc); 
    % rough estimate of the disk space of this run (warning for disks with
    % less free space, e.g. the system disk, see diskCheck)
    cat_io_dcm2bids_helper('diskCheck',Psnr,Pmdbdir,fullfile(Poutdir,job.BIDSdir),job.opts,job.nproc); 
  end
  C = struct('job',job,'Poutdir',Poutdir,'Pmdbdir',Pmdbdir,'Ptbldir',Ptbldir, ...
    'Pprodictdirs',{Pprodictdirs},'sites',{sites},'protocols',{protocols},'anon',anon); 
  R = cat_io_dcm2bids_parallel('run',C,Psnr,job.nproc); clear C
  fprintf('\nDCM2BIDS - import done.\n')

  % subject table and protocol variants of the scans of the database in the
  % order of the scans (after all packages, as parallel processes would 
  % write these files at the same time, see scanRecord in 
  % cat_io_dcm2bids_process)
  for vi = 1:numel(R.dbV)
    cat_io_dcm2bids_protocols('updateDBtables',R.dbV{vi}, Ptbldir, sites, job.opts); 
  end
  for ri = 1:numel(R.rec)
    if ~isempty(R.rec{ri}) && ~isempty(R.varg{ri})
      v = R.varg{ri}; 
      R.rec{ri}.Variant = cat_io_dcm2bids_protocols('protocolVariant',Ptbldir,v.PT,v.ProtocolName, ...
        v.datatype,v.type,job.opts.tolerance,job.opts); 
    end
  end

  % tables of the database (scans, sites, protocols, see cat_io_dcm2bids_tables)
  rec = [R.rec{~cellfun('isempty',R.rec)}]; 
  cat_io_dcm2bids_tables('writeTables',Ptbldir,rec,sites,job.opts);

  % anonymized DICOM export of the scans of this run (catDCM2BIDSdb-dicom,
  % see cat_io_dcm2bids_export)
  if job.opts.dicomexport > 0
    try
      cat_io_dcm2bids_export('exportDICOM',R.Pdbjson,R.datatype,R.pstatus,R.site, ...
        R.sub,R.ses,R.BIDSfile,Pmdbdir,job.opts);
    catch e
      cat_io_cprintf('err','  DICOM export failed: %s (%s, line %d)\n',e.message,e.stack(1).name,e.stack(1).line);
    end
  end

  % IntendedFor references of BIDS inputs with the new file names
  cat_io_dcm2bids_bids('remapIntendedFor',R.intended, R.srcmap); 

  % participants.tsv and private table (once for all participants)
  if R.partdata.Count > 0
    P = R.partdata.values; P = [P{:}]; 
    for bsd = unique({P.BIDSsubdir})
      Pb = P(strcmp({P.BIDSsubdir},bsd{1})); 
      cat_io_dcm2bids_bids('writeParticipantTSV',Pb, fullfile(Poutdir,bsd{1},'participants.tsv'), job.opts.anonymize); 
      BIDSsubdirprivate = strrep(bsd{1}, fullfile(job.BIDSdir,job.BIDSsubdir),fullfile([job.BIDSdir '-private'],job.BIDSsubdir)); 
      cat_io_dcm2bids_bids('writePrivateTSV',Pb, fullfile(Poutdir,BIDSsubdirprivate,'private_participants.tsv')); 
    end
  end

  
  if job.opts.output > 0  &&  isfield(job,'BIDSsubdir')

    % BIDS directories (protocol-specific subdirectories)
    if job.opts.protocolsubdirs
      Pdirs = cat_vol_findfiles( fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir),'*',struct('depth',1,'dirs',1));
    else
      Pdirs = {fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir)};
    end

    % report tables of the BIDS directory (scans, sessions, project, see 
    % cat_io_dcm2bids_report) and the generalized MRI protocols of its scans
    cat_io_dcm2bids_report('writeReportTables',Poutdir,job);
    cat_io_dcm2bids_bids('writeProtocols',fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir), ...
      fullfile(Poutdir,[job.BIDSdir '-report'],'protocols',job.BIDSsubdir),Ptbldir,job.opts.anonymize);
    if job.opts.reportfig.run % diagrams of the report tables (see cat_io_dcm2bids_reportfig)
      try
        cat_io_dcm2bids_reportfig('writeReportFigure',Poutdir,job);
      catch e
        cat_io_cprintf('err','  Report figure failed: %s (%s, line %d)\n',e.message,e.stack(1).name,e.stack(1).line);
      end
    end

    % dataset description, README, and participants description of each 
    % BIDS directory
    sources = cat_io_dcm2bids_bids('getBIDSsources',unique(R.bidsroots)); 
    for pdi = 1:numel(Pdirs)
      if exist(Pdirs{pdi},'dir'), cat_io_dcm2bids_bids('writeDatasetFiles',Pdirs{pdi}, job, sources, unique(R.convver)); end
    end
  
  %%
    % handling GZIP in output directory
    cat_io_dcm2bids_bids('gzipunzipOutputdata',Poutdir,job.opts)

    % render one slice per scan for each datatype/subtype
    if job.opts.render.run && job.opts.output > 1
      fprintf('\nDCM2BIDS - Render BIDS NIFTIs.\n')

      % (the subject and session reports of the parallel processes are 
      % already rendered, see cat_io_dcm2bids_parallel)
      ropts = job.opts.render; if R.rendered, ropts.subjects = 0; end
      ropts.dbdir = Pmdbdir; % (derivatives of the database, see opts.preprocessing = 1)
      cat_io_dcm2bids_render( fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir), ...
        fullfile(Poutdir,[job.BIDSdir '-report'],'render',job.BIDSsubdir), ropts);
      
      fprintf('\nDCM2BIDS - Render BIDS NIFTIs done.\n')
    end

    % converted raw BIDS images for the batch dependencies
    out = cat_io_dcm2bids_bids('getBIDSoutputs', fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir), out ); 
  end

  % readme files of the new directories of this run
  cat_io_dcm2bids_readme('writeReadmes',Poutdir,job);
  if disclaimer, cat_io_cprintf('err','\n%s\n',cat_io_dcm2bids_readme('disclaimer',job.opts.anonymize)); end
  if R.failed > 0
    cat_io_cprintf('err','DCM2BIDS - %d subject package(s) failed (see above).\n',R.failed); 
  end
end
% =========================================================================
function [same,job] = checkPreviousSetting(Popts,job)
%checkPreviousSetting. Test if settings are identical to previous runs.
  if exist(Popts,'file')
    load(Popts,'opts','dicts');
    % only the fields that define the BIDS structure (independent of other, 
    % e.g. removed or new, fields of the current or previous options)
    KEEP   = {'ProtocolFileName','anonymize','gzipe','tolerance','subIDform','sesunit'}; 
    KEEP   = intersect(KEEP,intersect(fieldnames(job.opts),fieldnames(opts))); % fields of both settings
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
          rmdir(spm_file(Popts,'path'), 's');
          mkdir(spm_file(Popts,'path'));
          saveSetting(Popts,job);
        end
      end
    end

    % changed protocol definitions or BIDS naming (see namingVersion in
    % cat_io_dcm2bids_bids): scans can get other names, whereas the files of
    % their former names remain in the BIDS directory
    if same
      S = load(Popts); prot = protocolSetup(job.dicts.Pprotocoldirs);
      nver = cat_io_dcm2bids_bids('namingVersion');
      if ~isfield(S,'protocols'), S.protocols = cell(0,2); end
      if ~isfield(S,'namingver'), S.namingver = 1; end % (BIDS directories without version)
      chg = setxor(S.protocols(:,1),prot(:,1)); % added or removed protocols
      for pri = reshape(find(ismember(prot(:,1),S.protocols(:,1))),1,[])
        if ~isequaln(prot{pri,2},S.protocols{strcmp(S.protocols(:,1),prot{pri,1}),2}), chg{end+1} = prot{pri,1}; end %#ok<AGROW>
      end
      if ~isempty(chg) || S.namingver ~= nver
        cat_io_cprintf('err',['  The protocol definitions or the BIDS naming changed since the last export to this ' ...
          'BIDS directory, i.e. scans can get other names, whereas the files of their former names remain.\n']);
        if S.namingver ~= nver, cat_io_cprintf('blue','    BIDS naming version %d (old %d)\n',nver,S.namingver); end
        if ~isempty(chg), cat_io_cprintf('blue','    Changed protocols:\n      %s\n',strjoin(sort(chg),sprintf('\n      '))); end
        p = spm_input('Changed protocols or BIDS naming. Select how to go on.',1,'m', ...
          {'Stop processing','Keep BIDS directory (old names remain)','Replace old BIDS directory'},[0 1 2],1);
        if p == 1
          saveSetting(Popts,job);
        elseif p == 2
          spm_figure('Clear',spm_figure('FindWin','Interactive'));
          if spm_input('Really replace old directory?',1,'Yes|No',[1 0],2)
            rmdir(spm_file(Popts,'path'), 's');
            mkdir(spm_file(Popts,'path'));
            saveSetting(Popts,job);
          else
            same = 0;
          end
        else
          same = 0;
        end
      end
    end
  else
    if ~exist(spm_file(Popts,'path'),'dir'), mkdir(spm_file(Popts,'path')); end
    saveSetting(Popts,job);
    same = 1;
  end
end
% =========================================================================
function saveSetting(Popts,job)
%saveSetting. Settings of the BIDS directory (see checkPreviousSetting):
%  options, dictionaries, protocol definitions, and the BIDS naming version.
  opts = job.opts; dicts = job.dicts;
  protocols = protocolSetup(job.dicts.Pprotocoldirs);
  namingver = cat_io_dcm2bids_bids('namingVersion');
  save(Popts,'opts','dicts','protocols','namingver');
end
% =========================================================================
function prot = protocolSetup(Pprodictdirs)
%protocolSetup. Normalized protocol definitions (without QC definitions,
%  see getProtocols in cat_io_dcm2bids_protocols) with their name relative
%  to the protocol directory (e.g. DZPG3T/anat/DZPG-anat-T2w-SPC.json),
%  i.e. independent of the location and format of the files.
  prot = cell(0,2);
  if isempty(Pprodictdirs) || all(cellfun('isempty',cellstr(Pprodictdirs))), return; end
  protocols = cat_io_dcm2bids_protocols('getProtocols',cellstr(Pprodictdirs));
  for pri = 1:size(protocols,1)
    if isempty(protocols{pri,2}), continue; end
    rel = regexp(protocols{pri,2},['(?<=[\\/])' regexptranslate('escape',protocols{pri,1}) '[\\/].*$'],'match','once');
    if isempty(rel), rel = spm_file(protocols{pri,2},'filename'); end
    prot(end+1,:) = {strrep(rel,'\','/'), protocols{pri,3}}; %#ok<AGROW>
  end
end
% =========================================================================
function Psnr = ageLimit(Psnr,Pmdbdir,Ptbldir,opts)
%ageLimit. Scan directories Psnr of the sessions within the age limit
%  opts.agelimit ([min max] in years, the age of the first sidecar of the
%  session with age, see assurePatientDCMfields in cat_io_dcm2bids_db).
%  Sessions outside the limit or without age (e.g. age 0, the undefined age
%  of many DICOMs, for a lower limit > 0) are removed with a warning and are
%  listed in the table excluded_sessions.<fmt> of the database tables
%  (cumulative, i.e. also the excluded sessions of former runs that are not
%  part of this run, see excluded_sessions.mat). An empty limit or [-Inf Inf]
%  keeps all sessions (also sessions without age).
  lim = double(opts.agelimit); 
  if numel(lim) ~= 2 || all(isinf(lim)), return; end
  Pses = cellfun(@fileparts,Psnr,'UniformOutput',false); 
  [uses,~,si] = unique(Pses,'stable'); 
  keep = true(numel(uses),1); 
  hdr  = {'Session','PatientID','PatientSex','ScanDate','PatientAge','Reason','Scans'}; T = cell(0,numel(hdr)); 
  for ui = 1:numel(uses)
    V = sessionAge(Psnr(si == ui)); 
    keep(ui) = V.PatientAge >= lim(1) && V.PatientAge <= lim(2); % (false for NaN)
    if ~keep(ui)
      if isnan(V.PatientAge),     reason = 'no age'; 
      elseif V.PatientAge < lim(1), reason = sprintf('age < %g',lim(1)); 
      else,                       reason = sprintf('age > %g',lim(2)); 
      end
      T(end+1,:) = {uses{ui}(numel(Pmdbdir)+2:end), V.PatientID, V.PatientSex, char(V.ScanDate), ...
        V.PatientAge, reason, nnz(si == ui)}; %#ok<AGROW>
    end
  end
  Psnr = Psnr(keep(si));
  if ~isempty(T)
    cat_io_cprintf('warn',['\n  Age limit [%g %g] years: %d of %d sessions are not processed (%d without age, ' ...
      'see opts.agelimit and %s):\n'],lim,size(T,1),numel(uses),nnz(strcmp(T(:,6),'no age')), ...
      fullfile(Ptbldir,['excluded_sessions.' opts.tableformat]));
    for ti = 1:min(size(T,1),50)
      if isnan(T{ti,5}), age = ''; else, age = sprintf(' (%g)',T{ti,5}); end
      cat_io_cprintf('warn','    %-70s %s%s\n',T{ti,1},T{ti,6},age);
    end
    if size(T,1) > 50, cat_io_cprintf('warn','    ... and %d further sessions\n',size(T,1) - 50); end
  end

  % cumulative table: the excluded sessions of former runs without the 
  % sessions of this run (that are tested again, e.g. for another limit)
  Pmat = fullfile(Ptbldir,'excluded_sessions.mat'); 
  if exist(Pmat,'file')
    S = load(Pmat,'T'); 
    if ~isempty(S.T), S.T = S.T(~ismember(S.T(:,1),cellfun(@(p) p(numel(Pmdbdir)+2:end),uses,'UniformOutput',false)),:); end
    T = [S.T; T];
  end
  if ~isempty(T) || exist(Pmat,'file')
    [~,so] = sort(T(:,1)); T = T(so,:);
    save(Pmat,'T');
    cat_io_dcm2bids_helper('writeTable',fullfile(Ptbldir,['excluded_sessions.' opts.tableformat]),[hdr; T],opts.tableformat);
  end
end
% =========================================================================
function V = sessionAge(Psnr)
%sessionAge. Sidecar fields of the first scan of the scan directories Psnr
%  (of one session) with age (or of the last readable scan without age). 
  V = struct('PatientID','NA','PatientSex','NA','ScanDate','','PatientAge',nan); 
  for di = 1:numel(Psnr)
    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{di}); 
    for ji = 1:numel(Pjson)
      try
        V = cat_io_dcm2bids_db('assurePatientDCMfields',cat_io_json(Pjson{ji})); 
      catch % (unreadable sidecar)
        continue
      end
      if ~isnan(V.PatientAge), return; end
    end
  end
end
