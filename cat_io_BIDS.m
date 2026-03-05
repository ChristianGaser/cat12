function out = cat_io_BIDS(varargin)
%cat_io_BIDSdata. Function to analyse and create a BIDS file structure.
% Because each file can define a fully different BIDS case a structure 
% for each file is created with similar entries. In case of processing 
% one element can be handled independent. To evaluate the full sample 
% use matlab brackets. For instance to get the surface result dirs of
% all files you can use  [{BIDS(fi).surfdir}]'  , where the chars are 
% first converted to cell, then grouped and then transposed. 
%
% The function is typically called in the preprocessing e.g. cat_run
% and added to the job structure to be used together with the job.subj
% selector. 
%
% The function is designed to be combined with spm_file and therefore 
% kept simple:
%
%  BIDS = cat_io_BIDS( files, job )
%  path = cat_io_BIDS( BIDS(1), BIDSfield )
%  file = cat_io_BIDS( BIDS(1), BIDSfield, spm_file_para_pairs, ...) 
%
%  files          .. input files  or  call of test dataset 
%
%  job            .. CAT preprocessing job variable
%   .output       .. CAT GUI structure (This is prefered to extopts!)
%   .extopts      .. CAT extopts struture partially defined in cat_defaults.m
%    .BIDS_folder .. BIDS result folder specified in cat_defaults
%                    (eg. ../derivatives/CAT_VER)
%    .subfolder   .. Use of CAT result folder specified in cat_defaults.
%                    Subfolders are never used in case of BIDS but else.
%                    0 - none
%                    1 - use catsubfolder in non-BIDS (default)   
%    .logdircase  .. 1 - main common directory of all inputs
%                    2 - current working directory (default)
%                    3 - first input directory 
%    .resdircase  .. 0 - noBIDS, allways maindir/CAT subdirs
%                    1 - BIDSdir only for BIDS input (default)
%                        no relative paths, i.e., ingnore "../"
%                    2 - BIDSdir allways
%                        with relative paths, i.e., use leading "../"
%                    3 - BIDSdir as pure relative dir independed of BIDS  
%    .createdirs  .. 0 - no (default)
%                    1 - yes, dependent on job.output, e.g. ceate surf 
%                        if surfaces are written
%    .mkBIDSdir   .. create (BIDS) result directories (0-no, default; 1-yes)
%    .verboseBIDS .. be verbose and print out (0-no,1-yes-warnings,2-yes-full)
%
%
%  BIDS(:)        .. full BIDS structure for each input file of files
%   .file         .. cleaned up input file
%
%   .SUB          .. BIDS subject directory or filename if noBIDS
%   .SES          .. BIDS session directory or empty
%   .ANA          .. BIDS anat directory or empty
%   .RUN          .. BIDS run filename specifier for rescans
%   .MOD          .. MRI modality (T1w,T2w,PDw,FLAIR)
%
%   .isBIDS       .. definition if the file/setup requires BIDS
%   .is[SUB|SES|ANA|RUN|MOD] .. availability of BIDS entry 
%
%   .postdirs     .. BIDS subdirectories, eg. sub*/ses*/anat
%   .predirs      .. addition relative directories ahead of the main 
%                    BIDS directory (see out.main_npredirs)
%
%   .rawdir       .. path of the raw BIDS input folder 
%                    (e.g, the parent of sub in case of BIDS)
%   .pdir         .. parent directory of input file or the sub directory 
%                    in case of BIDS
%   .resdir       .. path of the new BIDS derivative folder 
%   .tmpdir       .. force subdirectory 
%   .[label|mri|surf|report|err]dir 
%                 .. full path to specific result directories 
%                    (dependent on the option they are the same as resdir)
%
%
% Examples: 
%  1) Call testfunction to get a test-set of BIDS / nonBIDS files:
%       files = cat_io_BIDS( 'testfiles' );
%       BIDS  = cat_io_BIDS( files ); 
%       [{BIDS(fi).surfdir}]'
%
%     Predefined self tests: 
%       BIDS = cat_io_BIDS( 'selftest_long' );
%       BIDS = cat_io_BIDS( 'selftest_nonbids_auto' );
%
%  2) Get full BIDS structure for the first 3 test files
%       BIDS = cat_io_BIDS( files(1:3) )
%
%  3) Create the surface filename for the fifth file: 
%      spm_file( files(5), 'path', cat_io_BIDS( files(5), '', 'surfdir', ... 
%         'prefix', 'lh.central.', 'suffix', '.topofix','ext','')
%      cat_io_BIDS( files(5), 'surfdir', 'prefix', 'lh.central.', 'suffix', '.topofix','ext','')
%
%  4) Run surface filename creation for the different BIDS result settings: 
%      for i = 0:3, cat_io_BIDS( i , struct('extopts',struct('verboseBIDS',1))); end
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  %#ok<*AGROW>

  % == check input ==
  if nargin < 1
    help cat_io_BIDS; 
  
  elseif nargin == 1
    % basic testfunction 
    out = testfunction( varargin{1} ); 
  
  else
    if isstruct( varargin{1} ) && ~isstruct( varargin{2})
      % selection of specific field or files usign spm_file
      %   path = cat_io_BIDS( BIDS(1), BIDSfield )
      %   file = cat_io_BIDS( BIDS(1), BIDSfield, spm_file_para_pairs, ...) 
      out = getBIDS(varargin{:});

    elseif isstruct( varargin{2} ) && ... 
      ( iscell( varargin{1} ) || ischar( varargin{1} ) ) 
      % BIDS creation/setup
      out = setupBIDS(varargin{1:2});
      
      if nargin > 2
      % spm file setup: 
      %   file = cat_io_BIDS( BIDS(1), BIDSfield, spm_file_para_pairs, ...) 
        out = getBIDS(out, varargin{3:end});
      end

    else
      error('cat_io_BIDS: Incorrect use.'); 
      
    end
  end
end
%==========================================================================
function BIDS = setupBIDS(files,job)
%setupBIDS. Creation of BIDS structure based on file and parameter input.

  % guaranty cell strings and remove SPM image number
  files = cellstr(files); 
  files = spm_file(files, 'number', ''); 
  

  % prepare/complete job variable
  job = updateJob(job);
  

  % == main relative (BIDS) result directory ==
  % setup variables for subfolders (e.g, mri/surf) and 
  % the relative BIDSfolder_rel (e.g., ../derivatives/CAT0815)
  for fi = 1:numel(files) 
    % get function stack to check for longitudinal call
    spmstack = dbstack; 

    % BIDS parameters/settings
    BIDS(fi).main_BIDSfolder  = job.extopts.bids_folder;
    BIDS(fi).main_resdircase  = job.extopts.resdircase;
    BIDS(fi).main_subfolders  = job.extopts.subfolders;
    BIDS(fi).main_logdircase  = job.extopts.logdircase;
    BIDS(fi).main_mkBIDSdir   = job.extopts.mkBIDSdir; 
    BIDS(fi).main_verboseBIDS = job.extopts.verboseBIDS; 
    BIDS(fi).main_catfolders  = job.extopts.catfolders;
    BIDS(fi).main_long        = any( cat_io_contains( {spmstack.file} , 'cat_long_multi_run') );  

    % In the default BIDS case, no-relative path is allowed!
    if BIDS(fi).main_resdircase==1 
      if contains( BIDS(fi).main_BIDSfolder , ['..' filesep] ) && BIDS(fi).main_verboseBIDS
        cat_io_cprintf('warn','Warning:  Relative BIDS folder are not supported for "autobids" (resdircase=1)! \n')
      end
      BIDS(fi).main_BIDSfolder = strrep(BIDS(fi).main_BIDSfolder,['..' filesep],''); 
    end

  
    % count the number of wished recursive directories that we will have to add
    % after the BIDS result folder together with the BIDS subject subfolders
    BIDS(fi).main_npredirs = max([0,strfind(BIDS(fi).main_BIDSfolder,['..' filesep])]);

    BIDS(fi).file  = files{fi};

    % Use strsplit to handle the directories. On Mac/Linux the first /PATH 
    % result in an empty entry that is essential to create the full path /PATH. 
    % Otherwise, spm_file interpretate the PATH as relative and add the current
    % working directory!
    sdirs = strsplit(files{fi},filesep); 
    

    % Detect BIDS directories by looking for key words used in to define
    % the subject, session, and the anatomic directory.
    [BIDS(fi).SUB, BIDS(fi).SES, BIDS(fi).ANA, BIDS(fi).RUN, BIDS(fi).MOD, ...
     BIDS(fi).isSUB, BIDS(fi).isSES, BIDS(fi).isANA, BIDS(fi).isRUN, BIDS(fi).isMOD, ...
     BIDS(fi).postdirs, BIDS(fi).npostdirs] = getBIDSpostdirs0(sdirs);
  
    % define if output can or has to be BIDS (special case for longitudinal 
    % pipeline that is call as nonBIDS after preparing the BIDSdirs directly)
    BIDS(fi).isBIDS = BIDS(fi).isSUB & ( BIDS(fi).main_long | ~isempty(BIDS(fi).main_BIDSfolder)) & BIDS(fi).main_resdircase ~= 3; 

    % Avoid multiple derivatives and remove the leading part. In principle,
    % resdircase==0 should ignore such cases but as we had issues with
    % double derivatives in the longitudinal pipeline it is now just one
    if 1 % BIDS(fi).main_resdircase
      sdirs = clearDoubleDerivatives(sdirs, BIDS(fi).main_BIDSfolder, BIDS(fi).npostdirs);
    end

    % get preceding directories
    subhome = max(1,numel(sdirs) - BIDS(fi).npostdirs); 
    BIDS(fi).rawdir = [ repmat( filesep, isempty(sdirs{1}), isempty(sdirs{1})) ...
      fullfile( sdirs{ 1:max( 1, subhome - BIDS(fi).main_npredirs - 1 ) } )]; 
    if BIDS(fi).main_npredirs > 0  &&  numel(sdirs) > 1
      BIDS(fi).predirs = fullfile( sdirs{ max( 2, subhome - BIDS(fi).main_npredirs) : max(2,subhome - 1) } ); 
    else
      BIDS(fi).predirs= '';
    end
    clear subhome

    % no preciding directories in this case
    if BIDS(fi).main_resdircase==1 && ~BIDS(fi).isBIDS
       BIDS(fi).rawdir = fileparts(BIDS(fi).file);
      if isempty(BIDS(fi).rawdir), BIDS(fi).rawdir = pwd; end
      BIDS(fi).predirs = '';
    end

    [~,ff,ee] = fileparts(BIDS(fi).rawdir); BIDS(fi).pdir = [ff ee]; 
    
    % define main resultdir and subdirs
    if BIDS(fi).main_resdircase<1 || ( BIDS(fi).main_resdircase==1 && ~BIDS(fi).isBIDS )
      BIDS(fi).BIDSdir = '';
    else
      BIDS(fi).BIDSdir = fullfile( strrep( BIDS(fi).main_BIDSfolder , ...
       ['..' filesep], ''), BIDS(fi).predirs ); 
    end

    % create path to main result directory resdir and further cat sub directories mri/surf/...
    BIDS(fi).resdir = fullfile( BIDS(fi).rawdir , BIDS(fi).BIDSdir , BIDS(fi).postdirs);  
    for sfi = 1:numel(BIDS(fi).main_catfolders )
      if BIDS(fi).main_subfolders  &&  ~BIDS(fi).isBIDS 
        BIDS(fi).([BIDS(fi).main_catfolders{sfi} 'dir']) = fullfile( BIDS(fi).rawdir , ...
          BIDS(fi).BIDSdir , BIDS(fi).postdirs, BIDS(fi).main_catfolders{sfi} ); 
      else
        BIDS(fi).([BIDS(fi).main_catfolders{sfi} 'dir']) = fullfile( BIDS(fi).rawdir , ...
          BIDS(fi).BIDSdir , BIDS(fi).postdirs ); 
      end
    end      
  end


  % == final directory setup & creation ==
  BIDS = prepare_logpath(BIDS); 
  create_resultdirs(BIDS,job); 


  % == verbose debugging output ==
  if job.extopts.verboseBIDS > 1
    out = BIDS; 
    if numel(out,1)
      fprintf('\nVerbose cat_io_BIDS:\n')
      if iscell(out),       for fi = 1:size(out,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out{fi}); end
      elseif ischar(out),   for fi = 1:size(out,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out(fi,:)); end
      elseif isstruct(out), for fi = 1:numel(out),  cat_io_cprintf('blue','%4d)  %s\n', fi, out(fi).surfdir); end
      end
    else
      
      if iscell(out),     cat_io_cprintf('blue','\n  Verbose cat_io_BIDS: %s\n', out{1}); 
      elseif ischar(out), cat_io_cprintf('blue','\n  Verbose cat_io_BIDS: %s\n', out);
      elseif isstruct(out) 
        cat_io_cprintf('blue','\nVerbose cat_io_BIDS:\n');
        disp( out ); 
      end
    end
    fprintf('\n'); 
  end
end
%==========================================================================
function out = getBIDS(BIDS,subfield,varargin)
%getBIDS. Output (and manipulate) of a given subfield using spm_file.  

  if isfield( BIDS , subfield )
    if ischar( BIDS(1).(subfield) ) || iscell( BIDS(1).(subfield) ) 
      out = ({BIDS(:).(subfield)})';
      if nargin > 3 % call spm_file
        % add file name
        [~,ff,ee] = fileparts( {BIDS(:).file}' ); 
        out = spm_file( fullfile(out,strcat(ff,ee)) , varargin{:} ); 
      end
      if isscalar( out ), out = char(out); end
      
    else 
      out = ([BIDS(:).(subfield)])';
    end
  else
    error('cat_io_BIDS:unkownField','Unknown field "%s". \n', subfield)
  end
end
%==========================================================================
function BIDS = prepare_logpath(BIDS)
%prepare_logdir. Create a common log directory.

  switch BIDS(1).main_logdircase
    case 1
      % one main directory for all cases
      if numel(sfiles) > 1
        [~,S] = spm_str_manip( {BIDS(:).files} ,'C'); 
        logdir = fullfile(S.s, BIDS(1).main_BIDSfolder, 'log');
      else
        logdir = fullfile( {BIDS(1).resdir} , 'log');
      end
      [~,ff,ee] = fileparts(logdir); plogdir = [ff ee '_'];
    case 3
      % write into first subject directory
      logdir = fullfile(BIDS(1).resdir, 'log');
      [~,ff,ee] = fileparts(logdir); plogdir = [ff ee '_'];
    otherwise
      % write into current directory
      logdir = fullfile(pwd, 'log');
      [~,ff,ee] = fileparts(logdir); plogdir = [ff ee '_'];
  end
        
  if cat_io_contains(logdir, spm('dir'))
    % in case of the SPM directory we file this in a dataspecific subdirectory
    logdir = fullfile(spm('dir'),'toolbox','CAT', 'logs', char(datetime('now','Format','yyyyMMdd')) );
    plogdir = '';
    
    % give some feedback if directories are created (only then otherwise the warning comes with every call of this function)
    if BIDS(1).main_mkBIDSdir && BIDS(1).main_verboseBIDS && ~exist(logdir,'dir')
      cat_io_cprintf('warn','cat_io_BIDS:  Creation of log directory in SPM/CAT path redirected to: \n  %s\n', logdir)
    end
  end  

  % == create log and main result directory == 
  if ~exist(logdir,'dir')  &&  BIDS(1).main_mkBIDSdir 
    try
      mkdir(logdir); 
    catch
      main_logpath0 = logdir;
      logdir  = fullfile(BIDS(1).rawBIDSpath, BIDS(1).main_ressubdir, 'log');
      cat_io_cprintf('warning','cat_io_BIDS: Creation of log directory "%s" failed. \nUse path of first subject directory "%s".\n', ...
        main_logpath0, logdir)
    end
  end

  for i=1:numel(BIDS)
    BIDS(i).logdir  = logdir; 
    BIDS(i).plogdir = plogdir; 
  end
end
%==========================================================================
function create_resultdirs(BIDS,job)
%create_resultdirs. Create result directories.

  % for longitudinal we need update
  if BIDS(1).main_mkBIDSdir == 0; return; end

  % check if volume or surface atlas data is written 
  writeVatlases = isfield(job,'output') && isfield(job.output,'ROImenu') && isfield(job.output.ROImenu,'atlases') && ...
                  any(cell2mat(struct2cell(rmfield(job.output.ROImenu.atlases,'ownatlas')))); 
  writeSatlases = isfield(job,'output') && isfield(job.output,'sROImenu') && isfield(job.output.sROImenu,'atlases') && ...
                  any(cell2mat(struct2cell(rmfield(job.output.sROImenu.atlases,'ownatlas')))); 

  % check if any volume map has to be written
  writeVols = 0; 
  if isfield(job,'output')
    VFN = {'GM','WM','CSF','ct','pp','WMH','SL','TMPC','atlas','label','labelnative','bias','las','jocobianwarped','warps','rmat'}; 
    for VFNi = 1:numel(VFN)
      writeVols = writeVols || ( isfield(job.output,VFN{VFNi}) && any(cell2mat(struct2cell(job.output.GM))) );
    end
  end

  % check surface output
  writeSurfs = isfield(job,'output') && isfield(job.output,'surface') && job.output.surface;

  % force writing
  if BIDS(1).main_mkBIDSdir > 1 
     writeVatlases = 1; 
     writeVols     = 1;
     writeSurfs    = 1;
  end

  % create data specific outputs 
  for fi = 1:numel(BIDS)
    % this is not in the catdirlist
    if ~exist(BIDS(fi).resdir,'dir') 
      try
        mkdir(BIDS(fi).resdir); 
      catch
        cat_io_cprintf('err','Error: Cannot create directory "%s".', BIDS(fi).resdir);
        continue
      end
    end

    % mri/surf/label depending on if output is writen
    if ~exist(BIDS(fi).labeldir,'dir') && ( writeVatlases || writeSatlases || BIDS(1).main_mkBIDSdir > 1 )
      mkdir(BIDS(fi).labeldir); 
    end
    if ~exist(BIDS(fi).mridir,'dir') && ( writeVols || BIDS(1).main_mkBIDSdir > 1 )
      mkdir(BIDS(fi).mridir); 
    end
    if ~exist(BIDS(fi).surfdir,'dir') && ( writeSurfs || BIDS(1).main_mkBIDSdir > 1 )
      mkdir(BIDS(fi).surfdir); 
    end

    % create other listed directories besides the special cases above and 
    % the tmp and err dir that should be created only if needed
    subfolders = setdiff( BIDS(fi).main_catfolders , {'mri','surf','label','tmp','err'} ); 
    for si = 1:numel(subfolders)
      folder = BIDS(fi).( [ subfolders{si} 'dir' ] );
      if ~exist(folder,'dir')
        mkdir(folder); 
      end
    end
  end
end
%==========================================================================
function job = updateJob(job)
%updateJob. Check all varialbes. 
% job = updateJob(job)

  if isempty(job) && ~isstruct(job), clear job; job = struct(); end

  % update by GUI output variable
  if isfield(job,'output') && isfield(job.output,'BIDS')
    if isfield(job.output.BIDS,'BIDSno') 
      job.extopts.bids_folder = '';
      job.extopts.resdircase  = 0; 
    elseif isfield(job.output.BIDS,'BIDSyes') 
      job.extopts.bids_folder = job.output.BIDS.BIDSyes.BIDSfolder;
      job.extopts.resdircase  = 1; 
    elseif isfield(job.output.BIDS,'BIDSrel')
      job.extopts.bids_folder = job.output.BIDS.BIDSrel.BIDSfolder;
      job.extopts.resdircase  = 2;
    elseif isfield(job.output.BIDS,'relative')  
      job.extopts.bids_folder = job.output.BIDS.relative.BIDSfolder;
      job.extopts.resdircase  = 3;
    end
  end
  
  % basic defaults
  def.extopts             = cat_get_defaults('extopts');
  def.extopts.bids_folder = cat_get_defaults('extopts.bids_folder');  % already defined
  def.extopts.resdircase  = 1; % 0-noBIDS, 1-onlyBIDS, 2-allways
  def.extopts.logdircase  = 2; 
  def.extopts.mkBIDSdir   = 1; % create the result directories: 0-no,1-yes,needed,2-all
  def.extopts.verboseBIDS = 1; % 0-no, 1-yes-warings, 2-yes-full
  def.extopts.catfolders  = {'label','mri','surf','report','err','tmp'}; % basic catfolder definition
  
  job = cat_io_checkinopt(job,def); 

  % addapt for OS in all cases
  job.extopts.bids_folder = cat_io_strrep( job.extopts.bids_folder , {'/','\'}, filesep ); 

end
%==========================================================================
function sdirs = clearDoubleDerivatives(sdirs, ressubfolder, subdirs)
%clearDoubleDerivatives. Remove double derivatives directories. 
% sdirs = clearDoubleDerivatives(sdirs)

  % keep the first / 
  files = [repmat( filesep, isempty(sdirs{1}), isempty(sdirs{1}))  fullfile(sdirs{:}) ];
  files = strrep( files , strrep(ressubfolder ,['..' filesep],''),'');
  sdirs = strsplit(files,filesep);  

  dev = find( cellfun( @(x) strcmp(x,'derivatives'),sdirs) == 1 );
  if ~isempty(ressubfolder) && ~isempty(dev)
    % ../derivatives/d1/derivatives/d2/.. >> ../derivatives/d2/..
    sdirs(dev(1):end-1-subdirs) = []; 
  end
end
%==========================================================================
function [SUB, SES, ANA, RUN, MOD, isSUB, isSES, isANA, isRUN, isMOD, ...
  BIDSsubpath, BIDSsubpathdepth] = getBIDSpostdirs0(sdirs)
%getpostdirs. Get information about possible subject, session, and anat directories.  
% Detect BIDS directories by looking for key words used in to define
% the subject, session, and the anatomic directory.
%
%   [SUB, SES, ANA, RUN, MOD,  isSUB, isSES, isANA, isRUN, isMOD, .. 
%     BIDSsubpath, BIDSsubpathdepth] = getpostdirs(sdirs)
%

  wlist = {'_t1w','_t2w','_pdw','_flair'};

  % string variables
  SUB = spm_str_manip(sdirs{end},'r');  % use filename in case of non-BIDS (check for uniqueness?)
  SES = ''; 
  ANA = ''; 
  RUN = ''; 
  MOD = ''; 
  
  % boolean variables
  isANA = 0; 
  isSUB = 0; 
  isSES = 0;
  isRUN = cat_io_contains(lower(sdirs{end}),{'_run-'});
  isMOD = cat_io_contains(lower(sdirs{end}),wlist);

  % update boolean variables based on dirs
  if numel(sdirs) > 1
  % get the subpath with sub-*[/ses-*][/anat]
    isANA = cat_io_contains(lower(sdirs{end-1}(1:min(4,numel(sdirs{end-1})))),'anat');
    isSES = any(cat_io_contains(lower(sdirs{end-isANA-1}(1:min(4,numel(sdirs{end-1-isANA})))),{'ses-','ses_','sess-','sess_'}));
    isSUB = any(cat_io_contains(lower(sdirs{end-isANA-1-isSES}(1:min(4,numel(sdirs{end-1-isANA-isSES})))),{'sub-','sub_'})); 
   
    if isANA, ANA = sdirs{end-1}; end
    if isSES, SES = sdirs{end-1-isANA}; end
    if isSUB, SUB = sdirs{end-1-isANA-isSES}; end

    % more tolerant defintion
    if ~isSES && ~isSUB
      subcan = cat_io_contains(lower(sdirs(1:end-1)),{'sub-','sub_'});
      for si = sort(find(subcan),'des')
        isSUB = (numel(sdirs) - si) * any(cat_io_contains(lower(sdirs{si}(1:min(4,numel(sdirs{si})))),{'sub-','sub_'})); 
        if isSUB>0, break; end
      end
    end
  end

  % reruns
  filename = spm_str_manip( sdirs{end} , 'tr'); 
  if isRUN
    fparts = strsplit(filename,'_'); 
    subid  = find( cat_io_contains(lower(sdirs),{'sub-','sub_','ses-','ses_','sess-','sess_'})==1,1,'last');
    runid  = find( cat_io_contains(lower(sdirs{min(numel(sdirs),subid+1):end}),'run-')==1,1,'last');
    RUN    = fparts{runid};
  end

  % weighting
  if isMOD
    fparts = strsplit(filename,'_'); 
    modid  = cat_io_contains(lower(fparts),strrep(wlist,'_',''));
    MOD    = strcat(fparts{modid});
  end

  % directory information 
  BIDSsubpathdepth = isANA + isSES + isSUB; 
  if BIDSsubpathdepth > 0
    BIDSsubpath = fullfile( sdirs{ end - BIDSsubpathdepth : end - 1} );
  else
    BIDSsubpath = ''; 
  end
end 
%==========================================================================
function files = adaptWindows(files)
  if ispc
      files  = cat_io_strrep(files, {'/Users/','/'},{'C:\Users\','\'});
      % replace the starting / by 
      cfiles = char(files); 
      corr   = cfiles(:,1)=='\'; 
      files(corr)  = cellstr( strcat( 'C:' , files(corr) ) );
  end
end
%==========================================================================
function out = testfunction( testcase )
%testfunction. Create a set of test files. 

  if ~exist('testcase','var')
    testcase = 1;
  end

  if ~isnumeric(testcase)
  % get testfiles or run self test
    switch lower(testcase)
      case 'testfiles'
        out = getTestFiles;
      case 'testfiles_long'
        out = getLongitudinalTestFiles;
      case {'test_long','selftest_long'}
        out = runLongitudinalPathSelftest;
      case {'test_nonbids_auto','selftest_nonbids_auto'}
        out = runNonBIDSAutoSelftest();
    end
  else
  % run test for the different BIDS resdircase outputs

    % gest general testcase
    files = getTestFiles;

    % run it
    clear job
    job.extopts.resdircase  = max(0,min(3,testcase));
    job.extopts.verboseBIDS = 1; 
    job.extopts.mkBIDSdir   = 0; 
    switch job.extopts.resdircase
      case 0, job.output.BIDS.BIDSno              = struct();
      case 1, job.output.BIDS.BIDSyes.BIDSfolder  = fullfile('../derivatives','CATdefbids');
      case 2, job.output.BIDS.BIDSrel.BIDSfolder  = fullfile('../derivatives','CATrelbids');
      case 3, job.output.BIDS.relative.BIDSfolder = fullfile('../derivatives','CATreldirs');
    end    

    out = cat_io_BIDS(files,job); 
  end
end
%==========================================================================
function files = getTestFiles

  files = {
    ... BIDS
    '/Users/tomcat/BIDS_A/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii'                 % sub-ses-anat                                                
    '/Users/tomcat/BIDS_B/group-01/sub-01/ses-BL/anat/sub-01_ses-01_T1w.nii'        % additional group dir                                        
    '/Users/tom cat/BIDSr/group-01/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii' % additional rescan field                                    
    '/Users/tomcat/BIDS_C/sub-AA/BL/sub-01_T1w.nii'                                 % incorrect session + insufficient long filename + no anat                                                      
    '/Users/tomcat/BIDS_D/sub-01/ses-00/anat2/T1w.nii'                              % bad filename + multiple anats                                                    
    '/Users/tomcat/BIDS_E/sub-01/T1w.nii'                                           % insufficient filename                                                       
    '/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii'                              % main path                         
    ... no BIDS
    '/Users/tomcat/noBIDS_A/sub-01_T1.nii'                                                                                 
    '/Users/tomcat/noBIDS_B/mdir1/mdir2/sub-01_T1w.nii'                                                                           
    '/Users/tomcat/noBIDS_C/mdir3/mdir4/mdir5/sub-01_T1w.nii'                       % just some more directories                                                                    
    '/sub-01_ses-01_run-01_T1w.nii'                                     
    ... previous derivatives BIDS
    '/Users/tomcat/BIDSD_A/derivatives/CAT26.0.rc1_3074/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                                
    '/Users/tomcat/BIDSD_B/Project_J_derivatives/derivatives/CAT26.0.rc1_3074/derivatives/CAT26.0.rc1_3074/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                
    '/Users/tomcat/BIDSD_C/derivatives/CAT26.0.rc1_3074/group-01/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii'                           
    '/derivatives/CAT26.0.rc1_3074/sub-01/ses-01/anat/lh.thickness.sub-01_ses-01_run-01_T1w'                                     
    ... previous derivatives no BIDS
    '/Users/tomcat/DnoBIDSD_A/Project_G/mdir3/mdir4/derivatives/CAT12.9_2565/mri/p0sub-01_T1w.nii'                                            
    '/derivatives/DnoBIDSD_B/surf/lh.thickness.T1w'                                     
    '/noBIDSD_B/derivatives/CAT26.0.rc1_3074_long/derivatives/CAT26.0.rc1_3074_long/HC/M00/sub1.nii'
    ... longitudinal BIDS
    ...'/Users/tomcat/BIDSTEST/Project_LongBIDS/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii'
    ...'/Users/tomcat/BIDSTEST/Project_LongBIDS/derivatives/CAT12.9_2565/sub-01/ses-01/anat/p0sub-01_ses-01_T1w.nii' ... by longjob settings noBIDS!
    ...'/Users/tomcat/BIDSTEST/Project_LongBIDS/derivatives/CAT12.9_2565/sub-01/ses-01/anat/lh.central.sub-01_ses-01_T1w.gii'
    ... longitidunal 
    ...'/Users/tomcat/BIDSTEST/Project_Long/BL/T1w.nii'
    ...'/Users/tomcat/BIDSTEST/Project_Long/BL/derivatives/CAT12.9_2565/mrip0sub-01_ses-01_T1w.nii' ... by longjob settings noBIDS!
    ...'/Users/tomcat/BIDSTEST/Project_Long/BL/derivatives/CAT12.9_2565/surf/lh.central.sub-01_ses-01_T1w.gii'
  };

  files = adaptWindows(files);
end
%==========================================================================
function files = getLongitudinalTestFiles
  files = {
    '/Users/tomcat/BIDSTEST/Project_Long/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii'
    '/Users/tomcat/BIDSTEST/Project_Long/sub-01/ses-02/anat/sub-01_ses-02_T1w.nii'
    '/Users/tomcat/BIDSTEST/Project_Long/derivatives/CAT12.9_2565/sub-01/ses-01/anat/mri/p0sub-01_ses-01_T1w.nii'
    '/Users/tomcat/BIDSTEST/Project_Long/derivatives/CAT12.9_2565/sub-01/ses-02/anat/mri/p0sub-01_ses-02_T1w.nii'
    '/Users/tomcat/BIDSTEST/Project_Long/derivatives/CAT12.9_2565/sub-01/ses-01/anat/surf/lh.central.sub-01_ses-01_T1w.gii'
  };

  files = adaptWindows(files);
end
%==========================================================================
function out = runLongitudinalPathSelftest
  files = getLongitudinalTestFiles;

  clear job
  job.extopts.mkBIDSdir = 0;
  job.extopts.verbBIDS  = 0;
  job.output.BIDS.BIDSyes.BIDSfolder = fullfile('../derivatives','CATlongtest');

  BIDS = cat_io_BIDS(files,job);
  resdirs = {BIDS(:).resdir}';
  baddup  = [filesep 'derivatives' filesep 'derivatives' filesep];
  expectedroot = [filesep 'derivatives' filesep 'CATlongtest' filesep];

  assert(~any(contains(resdirs,baddup)), 'cat_io_BIDS:selftest_long:DoubleDerivatives', ...
    'Detected doubled derivatives path in longitudinal BIDS mapping.');
  assert(all(contains(resdirs,expectedroot)), 'cat_io_BIDS:selftest_long:WrongRoot', ...
    'Longitudinal BIDS mapping did not keep a single derivatives root.');
  assert(all(contains(resdirs,[filesep 'sub-01' filesep])), 'cat_io_BIDS:selftest_long:MissingSubject', ...
    'Longitudinal BIDS mapping lost the subject hierarchy.');
  assert(all(contains(resdirs,[filesep 'anat'])), 'cat_io_BIDS:selftest_long:MissingAnat', ...
    'Longitudinal BIDS mapping lost the anat hierarchy.');

  out = BIDS;
end
%==========================================================================
function out = runNonBIDSAutoSelftest()
  files = {
    '/Volumes/UltraMax/ADNI50/AD25_00/ADNI_005_S_0221.nii'
    '/Volumes/UltraMax/ADNI50/AD25_12/ADNI_005_S_0221.nii'
    '/Volumes/UltraMax/ADNI50/AD25_24/ADNI_005_S_0221.nii'
  };
  files = adaptWindows(files);

  clear job
  job.extopts.mkBIDSdir = 0;
  job.extopts.verbBIDS  = 0;
  job.output.BIDS.BIDSyes.BIDSfolder = fullfile('../derivatives','CATautobids');

  BIDS = cat_io_BIDS(files,job);
  resdirs = {BIDS(:).resdir}';
  srcdirs = cellfun(@fileparts,files,'UniformOutput',false)';

  assert(~any(contains(resdirs,[filesep 'derivatives' filesep])), 'cat_io_BIDS:selftest_nonbids_auto:UnexpectedDerivatives', ...
    'Auto-BIDS fallback for non-BIDS input still writes to derivatives.');
  samefallback = true;
  for fi = 1:numel(files)
    samefallback = samefallback && strcmp(resdirs{fi},srcdirs{fi});
  end
  assert(samefallback, 'cat_io_BIDS:selftest_nonbids_auto:WrongFallback', ...
    'Auto-BIDS fallback for non-BIDS input does not keep original directory.');


  clear job
  job.extopts.mkBIDSdir = 0;
  job.extopts.verbBIDS  = 0;
  job.output.BIDS.BIDSrel.BIDSfolder = fullfile('../derivatives','CATforced');
  BIDSforced = cat_io_BIDS(files,job);
  resdirsforced = {BIDSforced(:).resdir}';
  assert(all(contains(resdirsforced,[filesep 'derivatives' filesep 'CATforced' filesep])), 'cat_io_BIDS:selftest_nonbids_auto:ForcedModeBroken', ...
    'Forced derivatives mode (BIDSrel) no longer writes derivatives for non-BIDS input.');


  % idempotency for longitudinal-style reruns: if input is already copied
  % into derivatives/CATforced, a second mapping must keep the same resdir
  %rerunFiles = cellfun(@(f,rd) fullfile(rd,[spm_str_manip(f,'r') '.nii']), files, resdirsforced, 'UniformOutput', false);
  rerunFiles = cellfun(@(f,rd) fullfile(rd,[spm_str_manip(f,'tr') '.nii']), files, resdirsforced, 'UniformOutput', false);
  BIDSrerun = cat_io_BIDS(rerunFiles,job);
  resdirsrerun = {BIDSrerun(:).resdir}';

  sameresdir = true;
  for fi = 1:numel(files)
    sameresdir = sameresdir && strcmp(resdirsrerun{fi},resdirsforced{fi});
  end
  assert(sameresdir, 'cat_io_BIDS:selftest_nonbids_auto:ForcedModeNested', ...
    'Repeated non-BIDS forced-derivatives mapping creates nested derivatives paths.');
  assert(~any(contains(resdirsrerun,[filesep 'derivatives' filesep 'CATforced' filesep 'derivatives' filesep 'CATforced' filesep])), ...
    'cat_io_BIDS:selftest_nonbids_auto:ForcedModeDoubleDerivatives', ...
    'Repeated non-BIDS forced-derivatives mapping creates double derivatives folders.');

  out = BIDS;
end
%==========================================================================
