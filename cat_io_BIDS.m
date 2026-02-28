function out = cat_io_BIDS(varargin)
%cat_io_BIDSdata. Function to analyse and create a BIDS file structure.
% Because each file can define a fully different BIDS case a structure 
% for each file is created with similar entries. In case of processing 
% one element can be handled independent. To evaluate the full sample 
% use matlab brackets. For instance to get the surface result dirs of
% all files you can use  [{BIDS(fi).surfdir}]'  , where the chars are first
% converted to cell, then grouped and then transposed. 
%
% The function is designed to be combined with spm_file and therefore 
% kept simple.
%  file = cat_io_BIDS(BIDS(1),action,para, ...) 
%
%  BIDS = cat_io_BIDS(files, job);  
%
%  files          .. input files  or  call of test dataset 
%
%  job            .. CAT preprocessing job variable
%   .output       .. CAT GUI structure (This is prefered to extopts!)
%   .extopts      .. CAT extopts struture partially defined in cat_defaults.m
%    .BIDS_folder .. BIDS result folder specified in cat_defaults
%                    (eg. ../derivatives/CAT_VER)
%    .subfolder   .. use of CAT result folder specified in cat_defaults 
%                    0 - none
%                    1 - use catsubfolder in non-BIDS (default)   
%    .logdircase  .. 1 - main common directory of all inputs
%                    2 - current working directory (default)
%                    3 - first input directory 
%    .resdircase  .. 0 - noBIDS
%                    1 - BIDSdir only for BIDS input (default)
%                        no relative paths, i.e., ingnore "../"
%                    2 - BIDSdir allways
%                        with relative paths, i.e., use leading "../"
%                    3 - BIDSdir as pure relative dir independed of BIDS  
%    .createdirs  .. 0 - no
%                    1 - yes, dependent on job.output, e.g. ceate surf 
%                        if surfaces are written
%    .mkBIDSdir   .. create (BIDS) result directories (0-no, default; 1-yes)
%    .verbBIDS    .. be verbose and print out (0-no,1-yes)
%
%
%  BIDS(:)         .. full BIDS structure for each input file of files
%
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
%  2) Get full BIDS structure for the first 3 test files
%       BIDS = cat_io_BIDS( files(1:3) )
%
%  3) Create the surface filename for the fifth file: 
%      spm_file( file(5), 'path', cat_io_BIDS( files(5), '', 'surfpath', 'prefix', 'lh.central.', 'suffix', '.topofix','ext','')
%
%  4) Run surface filename creation for the different BIDS result settings: 
%      for i = 0:3, cat_io_BIDS( i , struct('extopts',struct('verbBIDS',1))); end
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
    if isstruct( varargin{1} )
      out = getBIDS(varargin{:});

    elseif isstruct( varargin{2} ) && ... 
      ( iscell( varargin{1} ) || ischar( varargin{1} ) ) 
      % BIDS setup and test functions
      out = setupBIDS(varargin{1:2});
      
      if nargin > 2
        out = getBIDS(out, varargin{3:end});
      end

    else
      error('cat_io_BIDS: Incorrect use'); 
      
    end
  end
end

function BIDS = setupBIDS(files,job)
  % guaranty cell strings and remove SPM image number
  files = cellstr(files); 
  files = spm_file(files, 'number', ''); 
  

  % prepare/complete job variable
  job = updateJob(job);
  

  % == main relative (BIDS) result directory ==
  % setup variables for subfolders (e.g, mri/surf) and 
  % the relative BIDSfolder_rel (e.g., ../derivatives/CAT0815)
  for fi = 1:numel(files) 
    % parameters
    BIDS(fi).main_BIDSfolder = job.extopts.bids_folder;
    BIDS(fi).main_resdircase = job.extopts.resdircase;
    BIDS(fi).main_subfolders = job.extopts.subfolders;
    BIDS(fi).main_logdircase = job.extopts.logdircase;
    BIDS(fi).main_mkBIDSdir  = job.extopts.mkBIDSdir; 
    BIDS(fi).main_warning    = job.extopts.BIDSwarning; 
    BIDS(fi).main_catfolders = job.extopts.catfolders; 
    

    % update by GUI output variable
    if isfield(job,'output') && isfield(job.output,'BIDS')
      if isfield(job.output.BIDS,'BIDSno') 
        BIDS(fi).main_BIDSfolder = '';
        BIDS(fi).main_resdircase = 0; 
      elseif isfield(job.output.BIDS,'BIDSyes') 
        BIDS(fi).main_BIDSfolder = job.output.BIDS.BIDSyes.BIDSfolder;
        BIDS(fi).main_resdircase = 1; 
      elseif isfield(job.output.BIDS,'BIDSrel')
        BIDS(fi).main_BIDSfolder = job.output.BIDS.BIDSrel.BIDSfolder;
        BIDS(fi).main_resdircase = 2;
      elseif isfield(job.output.BIDS,'relative')  
        BIDS(fi).main_BIDSfolder = job.output.BIDS.relative.BIDSfolder;
        BIDS(fi).main_resdircase = 3;
      end
    end


    % In the default BIDS case, no-relative path is allowed!
    if BIDS(fi).main_resdircase==0 && BIDS(fi).main_warning
      if contains( BIDS(fi).main_BIDSfolder , ['..' filesep] )
        cat_io_cprintf('warn','Warning:  Relative BIDS folder are not supported for "autobids" (resdircase=1)! \n')
      end
      BIDS(fi).main_BIDSfolder = strrep(BIDS(fi).main_BIDSfolder,['..' filesep],''); 
    end

  
    % count the number of wished recursive directories
    % that we will have to add after the BIDS result folder together with 
    % the BIDS subject subfolders
    BIDS(fi).main_npredirs = max([0,strfind(BIDS(fi).main_BIDSfolder,['..' filesep])]);
  
    BIDS(fi).file  = files{fi};
    sdirs = strsplit(files{fi},filesep); 
    

    % Detect BIDS directories by looking for key words used in to define
    % the subject, session, and the anatomic directory.
    [BIDS(fi).SUB, BIDS(fi).SES, BIDS(fi).ANA, BIDS(fi).RUN, BIDS(fi).MOD, ...
     BIDS(fi).isSUB, BIDS(fi).isSES, BIDS(fi).isANA, BIDS(fi).isRUN, BIDS(fi).isMOD, ...
     BIDS(fi).postdirs, BIDS(fi).npostdirs] = getBIDSpostdirs(sdirs);
  
    % define if output can or has to be BIDS
    BIDS(fi).isBIDS = BIDS(fi).isSUB & ~isempty(BIDS(fi).main_BIDSfolder) & BIDS(fi).main_resdircase ~= 3; 

    % avoid multiple derivatives and remove the leading part
    if BIDS(fi).isBIDS
      sdirs = clearDoubleDerivatives(sdirs, BIDS(fi).main_BIDSfolder, BIDS(fi).npostdirs);
    end

    % detect extra directories
    subhome = max(1,numel(sdirs) - BIDS(fi).npostdirs); 
    BIDS(fi).rawBIDSpath = [ repmat( filesep, isempty(sdirs{1}), isempty(sdirs{1})) ...
      fullfile( sdirs{ 1:max( 1, subhome - BIDS(fi).main_npredirs - 1 ) } )]; 
    if BIDS(fi).main_npredirs > 0  &&  numel(sdirs) > 1
      BIDS(fi).predirs = fullfile( sdirs{ max( 2, subhome - BIDS(fi).main_npredirs) : max(2,subhome - 1) } ); 
    else
      BIDS(fi).predirs= '';
    end
    [~,ff,ee] = fileparts(BIDS(fi).rawBIDSpath); BIDS(fi).pdir = [ff ee]; 

    
    % define main resultdir and subdirs
    BIDS(fi).BIDSdir = fullfile( strrep( BIDS(fi).main_BIDSfolder , ...
      ['..' filesep], ''), BIDS(fi).predirs ); 
    BIDS(fi).resBIDSpath = fullfile( BIDS(fi).rawBIDSpath , ...
      strrep( BIDS(fi).main_BIDSfolder ,['..' filesep], '') , BIDS(fi).predirs ); 
   
    BIDS(fi).resdir = fullfile( BIDS(fi).rawBIDSpath , BIDS(fi).BIDSdir , BIDS(fi).postdirs);  
    for sfi = 1:numel(BIDS(fi).main_catfolders)
      % full path
      if BIDS(fi).main_subfolders && ~BIDS(fi).isBIDS
        BIDS(fi).([BIDS(fi).main_catfolders{sfi} 'dir']) = fullfile( BIDS(fi).rawBIDSpath , ...
          BIDS(fi).BIDSdir , BIDS(fi).postdirs, BIDS(fi).main_catfolders{sfi} ); 
      else
        BIDS(fi).([BIDS(fi).main_catfolders{sfi} 'dir']) = fullfile( BIDS(fi).rawBIDSpath , ...
          BIDS(fi).BIDSdir , BIDS(fi).postdirs ); 
      end
    end      
  end


  % == final directory setup & creation ==
  BIDS = prepare_logpath(BIDS); 
  create_resultdirs(BIDS,job); 


  % == debugging ==
  if job.extopts.verbBIDS 
    out = BIDS; 
    if size(out,1) > 1
      fprintf('\nVerbose cat_io_BIDS:\n')
      if iscell(out),       for fi = 1:size(out,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out{fi}); end
      elseif ischar(out),   for fi = 1:size(out,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out(fi,:)); end
      elseif isstruct(out), for fi = 1:size(out.resBIDSpath,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out.resBIDSpath{fi}); end
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
    if BIDS(1).main_mkBIDSdir && BIDS(1).main_warning
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
      cat_io_cprintf('warning','cat_io_BIDS: Creation of log directory "%s" failed. \nUse path of first subject "%s".\n', ...
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
      mkdir(BIDS(fi).resdir); 
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

  def.extopts             = cat_get_defaults('extopts');
 %def.extopts.bids_folder = cat_get_defaults('extopts.bids_folder');  % already defined
  def.extopts.resdircase  = 1; % 0-noBIDS, 1-onlyBIDS, 2-allways
  def.extopts.logdircase  = 2; 
  def.extopts.mkBIDSdir   = 1; 
  def.extopts.verbBIDS    = 0; 
  def.extopts.BIDSwarning = 0; 
  % basic catfolder definition
  def.extopts.catfolders  = {'label','mri','surf','report','err','tmp'};
  
  job = cat_io_checkinopt(job,def); 
  
  % addapt for OS in all cases
  job.extopts.bids_folder = cat_io_strrep( job.extopts.bids_folder , {'/','\'}, filesep ); 

end
%==========================================================================
function sdirs = clearDoubleDerivatives(sdirs, ressubfolder, subdirs)
%clearDoubleDerivatives. Remove double derivatives directories. 
% sdirs = clearDoubleDerivatives(sdirs)

  dev = find( cellfun( @(x) strcmp(x,'derivatives'),sdirs) == 1 );
  if ~isempty(ressubfolder) && ~isempty(dev)
    % ../derivatives/d1/derivatives/d2/.. >> ../derivatives/d2/..
    sdirs(dev(1):end-1-subdirs) = []; 
  end
end
%==========================================================================
function [SUB, SES, ANA, RUN, MOD, isSUB, isSES, isANA, isRUN, isMOD, ...
  BIDSsubpath, BIDSsubpathdepth] = getBIDSpostdirs(sdirs)
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
function out = testfunction( testcase )
%testfunction. Create a set of test files. 

  if ~exist('testcase','var')
    testcase = 1;
  else
    if isnumeric(testcase)
      testcase = max(0,min(3,testcase));
    end
  end


  files = {
    ... full BIDS: sub-ses-anat
    '/Users/tomcat/BIDSTEST/Project_A_fulllong/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii'                                                    
    '/Users/tomcat/BIDSTEST/Project_A_fulllong/sub-01/ses-02/anat/sub-01_ses-02_T1w.nii'                                                   
    ... no BIDS
    '/Users/tomcat/BIDSTEST/Project_B_noDirs/sub-01_T1.nii'                                                                                 
    '/Users/tomcat/BIDSTEST/Project_B_noDirs/sub-02_T1.nii'     
    ... some BIDS: incorrect session + bad filename + no anat
    '/Users/tomcat/BIDSTEST/Project_C_badses_noAnat/sub-AA/BL/sub-01_T1w.nii'                                                         
    '/Users/tomcat/BIDSTEST/Project_C_badses_noAnat/sub-AA/FU/sub-01_T1w.nii'        
    ... BIDS: bad filename + multiple anats
    '/Users/tomcat/BIDSTEST/Project_D_badfname_multiAnat/sub-01/ses-00/anat/T1w.nii'                                                        
    '/Users/tomcat/BIDSTEST/Project_D_badfname_multiAnat/sub-01/ses-00/anat2/T1w.nii'       
    ... BIDS: 
    '/Users/tomcat/BIDSTEST/Project_E_sameName/sub-01/T1w.nii'                                                                       
    '/Users/tomcat/BIDSTEST/Project_E_sameName/sub-02/T1w.nii'
    ... some directories
    '/Users/tomcat/BIDSTEST/Project_F_somedirs/mdir1/mdir2/sub-01_T1w.nii'                                                                           
    '/Users/tomcat/BIDSTEST/Project_F_somedirs/mdir1/mdir2/sub-02_T1w.nii'                                                                           
    ... previous derivatives 
    '/Users/tomcat/BIDSTEST/Project_G/mdir3/mdir4/derivatives/CAT12.9_2565/mri/p0sub-01_T1w.nii'                                            
    '/Users/tomcat/BIDSTEST/Project_G/mdir3/mdir4/derivatives/CAT12.9_2565/mri/p0sub-02_T1w.nii'      
    ... just some more directories 
    '/Users/tomcat/BIDSTEST/Project_H/mdir3/mdir4/mdir5/sub-01_T1w.nii'                                                                     
    '/Users/tomcat/BIDSTEST/Project_H/mdir3/mdir4/mdir5/sub-02_T1w.nii'  
    ... additional group dir
    '/Users/tomcat/BIDSTEST/Project_I_addGroup/group-01/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                           
    '/Users/tomcat/BIDSTEST/Project_I_addGroup/group-02/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                           
    ... doubled derivatives 
    '/Users/tomcat/BIDSTEST/Project_J_derivatives/derivatives/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                     
    '/Users/tomcat/BIDSTEST/Project_J_derivatives/derivatives/sub-02/ses-00/anat/sub-02_ses-01_T1w.nii'       
    ... multiple subjects
    '/Users/tomcat/BIDSTEST/Project_K_derivatives/derivatives/group-01/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii'                                     
    '/Users/tomcat/BIDSTEST/Project_K_derivatives/derivatives/group-01/sub-02/ses-00/anat/sub-02_ses-01_run-01_T1w.nii'                                    
    ... spaces
    '/Users/tom cat/BIDSTEST/Project L IXI2026/group-01/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii'                                     
    '/Users/tom cat/BIDSTEST/Project L IXI2026/group-01/sub-02/ses-00/anat/sub-02_ses-01_run-01_T1w.nii'                                    
    ... main path with/without BIDS dirs
    '/sub-01/ses-00/anat/sub-01_ses-01_run-01_T1w.nii'                                     
    '/sub-02/ses-00/anat/sub-02_ses-01_run-01_T1w.nii'   
    '/sub-01_ses-01_run-01_T1w.nii'                                     
    '/sub-02_ses-01_run-01_T1w.nii'
    ...
    '/mymri/derivatives/BIDS/sub-01/ses-00/anat/surf/lh.thickness.sub-01_ses-01_run-01_T1w'                                     
    '/mymri/derivatives/BIDS/sub-01/ses-00/anat/surf/lh.central.sub-01_ses-01_run-01_T1w.gii'   
    ...
    '/mymri/derivatives/noBIDS/surf/lh.thickness.T1w'                                     
    '/mymri/derivatives/noBIDS/surf/lh.central.T1w.gii'                                       
  };

  % keep it shorter
  files(2:2:end) = []; 
 
  % adopt for windows
  if ispc
    files  = cat_io_strrep(files, {'/Users/','/'},{'C:\Users\','\'});
    % replace the starting / by 
    cfiles = char(files); 
    corr   = cfiles(:,1)=='\'; 
    files(corr)  = cellstr( strcat( 'C:' , files(corr) ) );
  end

  if ~isnumeric(testcase) && strcmp(testcase,'testfiles')
    out = files; 
    return; 
  end
 
  % run it
  clear job
  job.extopts.resdircase = testcase;
  job.extopts.verbBIDS   = 0; 
  job.extopts.mkBIDSdir  = 0; 
  switch job.extopts.resdircase
    case 0, job.output.BIDS.BIDSno              = struct();
    case 1, job.output.BIDS.BIDSyes.BIDSfolder  = fullfile('../derivatives','CATdefbids');
    case 2, job.output.BIDS.BIDSrel.BIDSfolder  = fullfile('../derivatives','CATrelbids');
    case 3, job.output.BIDS.relative.BIDSfolder = fullfile('../derivatives','CATreldirs');
  end
    
  % show result
  BIDS = cat_io_BIDS(files,job); 
  out  = cat_io_BIDS(BIDS,'surfdir','prefix','lh.central.','suffix','.topofix','ext','');

end
