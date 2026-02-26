function out = cat_io_BIDS(files, job, subfield, varargin)
%Function to analyse files and setup related BIDS directories and files. 
%
%  out = cat_io_BIDS(files, job, action, varargin)
%
%  files          .. input files  or  call of test dataset 
%  job            .. CAT preprocessing job variable
%   .output       .. CAT GUI structure that is always prefered to extopts 
%   .extopts      .. CAT extopts struture partially defined in cat_defaults.m
%    .BIDS_folder .. "../derivatives/CAT#" 
%    .subfolder   .. 0 - none, 1 - use catsubfolder (default)                     
%    .resdircase  .. 0 - noBIDS
%                    1 - BIDSdir only if BIDS (default)
%                        no relative paths, i.e., ingnore "../"
%                    2 - BIDSdir allways
%                        with relative paths, i.e., use leading "../"
%                    3 - BIDSdir as pure relative dir independed of BIDS  
%    .createdirs  .. 0 - no
%                    1 - yes, dependent on job.output, e.g. ceate surf 
%                        if surfaces are written
%                    2 - yes, allways
%    .verbBIDS    .. be verbose and print out (0-no,1-yes)
%    .mkBIDSdir   .. create (BIDS) result directories (0-no, default; 1-yes)
%  subfield       .. field selctor to get only this field in the output structur
%  varargin       .. input for spm_file to add and maninpulate the filename 
%
%  out            .. full BIDS structure or BIDS.(subfield) element 
%   .BIDSdir      .. BIDS resultdirectory (/deriatives/CAT)
%   .files        .. cleaned up input files
%   .isBIDS       .. one if input file is BIDS
%   .BIDSsub      .. subject/file name
%   .devdir       .. main data directory (parent of sub in case of BIDS)
%                    with derivatives!
%   .rdevdir      .. relative directory (e.g. sub*/ses*/anat*)
%   .logdir       .. common directory for log files
%   .mdevdir      .. main data directory (parent of sub in case of BIDS)
%                    without derivatives
%
% Examples: 
%  1) Call testfunction for a set of BIDS / nonBIDS files:
%      cat_io_BIDS( 'test' , 
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


% TODO: 
% - implement/check directory creation 


  % == check input ==
  if nargin < 1, help cat_io_BIDS; end
  if ~exist('subfield','var'), subfield = ''; end
  if ~exist('job','var'), job = struct(); end

  % prepare/complete job variable
  job = updateJob(job);


  % == handle char/cell input and tests ==
  if ischar(files) && size(files,1) > 1
    % Creation of a data record out(:) for each input. 
    out = cat_io_BIDSchar(files, job, subfield); 
    out(:).main_logdir = out(1).main_logdir;
    return;
  elseif ischar(files) || isnumeric(files)
    % call testfunction 
    out = testfunction(files,job,subfield); 
    return
  else
    files = cellstr(files); 
  end


  % == main relative (BIDS) result directory ==
  % setup variables for subfolders (e.g, mri/surf) and 
  % the relative BIDSfolder_rel (e.g., ../derivatives/CAT0815)
  BIDS.main_ressubdir  = job.extopts.bids_folder;
  BIDS.main_resdircase = job.extopts.resdircase;
  BIDS.main_catfolders = job.extopts.subfolders;
  BIDS.main_logdircase = job.extopts.logdircase;
  BIDS.main_mkBIDSdir  = job.extopts.mkBIDSdir; 
  
  % update by GUI output variable
  if isfield(job,'output') && isfield(job.output,'BIDS')
    if isfield(job.output.BIDS,'BIDSno') 
      BIDS.main_ressubdir  = '';
      BIDS.main_resdircase = 0; 
    elseif isfield(job.output.BIDS,'BIDSyes') 
      BIDS.main_ressubdir  = job.output.BIDS.BIDSyes.BIDSfolder;
      BIDS.main_resdircase = 1; 
    elseif isfield(job.output.BIDS,'BIDSrel')
      BIDS.main_ressubdir  = job.output.BIDS.BIDSrel.BIDSfolder;
      BIDS.main_resdircase = 2;
    elseif isfield(job.output.BIDS,'relative')  
      BIDS.main_ressubdir  = job.output.BIDS.relative.BIDSfolder;
      BIDS.main_resdircase = 3;
    end
  end

  % In the default BIDS case, no-relative path is allowed!
  if BIDS.main_resdircase==0
    if contains( BIDS.main_ressubdir , ['..' filesep] )
      cat_io_cprintf('warn','Warning:  Relative BIDS folder are not supported for "autobids" (resdircase=1)! \n')
    end
    BIDS.main_ressubdir = strrep(BIDS.main_ressubdir,['..' filesep],''); 
  end

  % count the number of wished recursive directories
  % that we will have to add after the BIDS result folder together with 
  % the BIDS subject subfolders
  BIDS.main_adddirs = max([0,strfind(BIDS.main_ressubdir,['..' filesep])]);

  % basic catfolder definition
  catfolder = {'label','mri','surf','report','err'}; 
  
  % remove SPM image number
  BIDS.files = cellstr(files);
  BIDS.files = spm_file(BIDS.files, 'number', ''); 


  % == file specific evaluation == 
  for fi = 1:numel(files) 
    sdirs = strsplit(files{fi},filesep); 
    

    % Detect BIDS directories by looking for key words used in to define
    % the subject, session, and the anatomic directory.
    [BIDS.SUB{fi,1}, BIDS.SES{fi,1}, BIDS.ANA{fi,1}, BIDS.RUN{fi,1}, BIDS.MOD{fi,1}, ...
     BIDS.isSUB(fi,1), BIDS.isSES(fi,1), BIDS.isANA(fi,1), BIDS.isRUN(fi,1), BIDS.isMOD(fi,1), ...
     BIDS.BIDSsubdirs{fi,1}, BIDS.BIDSsubdirnumber(fi,1)] = getBIDSsubdirs(sdirs);
  
    % define if output can or has to be BIDS
    BIDS.isBIDS(fi,1) = BIDS.isSUB(fi) & ~isempty(BIDS.main_ressubdir) & BIDS.main_resdircase ~= 3; 

    % avoid multiple derivatives and remove the leading part
    if  BIDS.isBIDS(fi)
      sdirs = clearDoubleDerivatives(sdirs, BIDS.main_ressubdir, BIDS.BIDSsubdirnumber(fi));
    end

    % detect extra directories
    subhome = max(1,numel(sdirs) - BIDS.BIDSsubdirnumber(fi)); 
    BIDS.BIDSrawpath{fi,1} = [ repmat( filesep, isempty(sdirs{1}), isempty(sdirs{1})) ...
      fullfile( sdirs{ 1:max( 1, subhome - BIDS.main_adddirs - 1 ) } )]; 
    if BIDS.main_adddirs>1
      BIDS.adddirs{fi,1} = fullfile( sdirs{ max( 2, subhome - BIDS.main_adddirs) : max(2,subhome - 1) } ); 
    else
      BIDS.adddirs{fi,1} = '';
    end

    
    % define main resultdir and subdirs
    BIDS.resdir{fi,1}     = fullfile( strrep( BIDS.main_ressubdir , ...
      ['..' filesep], ''), BIDS.adddirs{fi} ); 
    BIDS.resdirpath{fi,1} = fullfile( BIDS.BIDSrawpath{fi} , ...
      strrep( BIDS.main_ressubdir ,['..' filesep], '') , BIDS.adddirs{fi} ); 
    for sfi = 1:numel(catfolder)
      % just the dirname
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'dir']){fi,1} = catfolder{sfi}; 
      else
        BIDS.([catfolder{sfi} 'dir']){fi,1} = ''; 
      end
    end      
    for sfi = 1:numel(catfolder)
      % full path
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'path']){fi,1} = fullfile( BIDS.BIDSrawpath{fi} , ...
          BIDS.resdir{fi,1} ,  BIDS.BIDSsubdirs{fi}, catfolder{sfi} ); 
      else
        BIDS.([catfolder{sfi} 'path']){fi,1} = fullfile( BIDS.BIDSrawpath{fi} , ...
          BIDS.resdir{fi,1} ,  BIDS.BIDSsubdirs{fi} ); 
      end
    end      
  end


  % == directory setup ==
  prepare_logdir(BIDS); 
  %prepare_resdir(BIDS); 


  % == specify output ==
  if isempty( subfield )
    % full BIDS structure
    out = BIDS; 
  else
    % only the relevant action subfield 
    if isfield( BIDS , subfield )
       out = BIDS.(subfield); 

       if nargin > 3
         [~,ff,ee] = fileparts(files); 
         out = spm_file( fullfile(out,strcat(ff,ee)) , varargin{:} ); 
       end
    else
      error('cat_io_BIDS:unkownField','Unknown field "%s". \n', subfield)
    end
  end

  % == debugging ==
  if job.extopts.verbBIDS 
    fprintf('\nVerbose cat_io_BIDS:\n')
    if iscell(out),   for fi = 1:size(out,1), cat_io_cprintf('blue','%4d) %s\n', fi, out{fi}); end; end
    if ischar(out),   for fi = 1:size(out,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out{fi}); end; end
    if isstruct(out), for fi = 1:size(out.resdirpath,1), cat_io_cprintf('blue','%4d)  %s\n', fi, out.resdirpath{fi}); end; end
    fprintf('\n'); 
  end
end
%==========================================================================
function main_logpath = prepare_logdir(BIDS)
%prepare_logdir.

  switch BIDS.main_logdircase
    case 1
      % one main directory for all cases
      if numel(sfiles) > 1
        [~,S] = spm_str_manip(files,'C'); 
        main_logpath = fullfile(S.s, BIDS.main_ressubdir, 'log');
      else
        main_logpath = fullfile(BIDS.BIDSrawpath{fi}, BIDS.main_ressubdir, 'log');
      end
    case 2
      % write into current directory
      if strcmp(pwd,fullfile(spm('dir'),'toolbox','CAT'))
        % in case of the CAT directory we file this in a dataspecific subdirectory
        main_logpath = fullfile(pwd, 'logs', char(datetime('now','Format','yyyyMMdd')) );
      else
        main_logpath = fullfile(pwd, BIDS.main_ressubdir,'log');
      end
    case 3
      % write into first subject directory
      main_logpath = fullfile(BIDS.BIDSrawpath{fi}, BIDS.main_ressubdir, 'log');
  end
  

  % avoid SPM dirs and write into a specific CAT subdir
  if cat_io_contains( main_logpath , spm('dir') )
    logdirfailed = main_logpath; 
    main_logpath = fullfile(BIDS.resdirpath{1}, 'log'); % #########
    cat_io_cprintf('warning',['cat_io_BIDS:  Creation of log directory "%s" failed. \n' ...
      'Write into first subject directory: "%s"\n'],logdirfailed,main_logpath)
  end
  

  % == create log and main result directory == 
  if ~exist(main_logpath,'dir')  &&  BIDS.main_mkBIDSdir 
    try
      mkdir(main_logpath); 
    catch
      cat_io_cprintf('warning','cat_run: creation of log directory "%s" failed. \n',main_logpath)
    end
  end

end
%==========================================================================
function out = cat_io_BIDSchar(files, job, subfield)
%cat_io_BIDSchar. Recursive call to create a data record out(f).
  for fi = 1:size(files,1)
    if isempty(subfield) % struct
      out(fi) = cat_io_BIDS({files(fi,:)},job,subfield); %#ok<*AGROW>
    else
      out1 = cat_io_BIDS({files(fi,:)},job,subfield); 
      if iscell(out1) || ischar(out1)
        out{fi,1} = out1; 
      elseif isscalar(out1)
        out(fi,1) = out1; 
      end
    end  
  end
end
%==========================================================================
function job = updateJob(job)
%updateJob. Check all varialbes. 

  def.extopts             = cat_get_defaults('extopts');
  def.extopts.bids_folder = '';  % main (BIDS) result folder
  def.extopts.resdircase  = 1; 
  def.extopts.logdircase  = 2; 
  def.extopts.mkBIDSdir   = 0; 
  def.extopts.verbBIDS    = 0; 

  job = cat_io_checkinopt(job,def); 
  
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
  BIDSsubpath, BIDSsubpathdepth] = getBIDSsubdirs(sdirs)
%getBIDSsubdirs. Get information about possible subject, session, and anat directories.  
% Detect BIDS directories by looking for key words used in to define
% the subject, session, and the anatomic directory.
%
%   [SUB, SES, ANA, RUN, isSUB, isSES, isANA, isRUN, BIDSsubpath, BIDSsubpathdepth] = getBIDSsubdirs(sdirs)
%

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
  isMOD = cat_io_contains(lower(sdirs{end}),{'_t1w','_t2w','_pdw','_flair'});

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
  if isRUN
    fparts = strsplit(sdirs{end},'_'); 
    subid  = find( cat_io_contains(lower(sdirs),{'sub-','sub_','ses-','ses_','sess-','sess_'})==1,1,'last');
    runid  = find( cat_io_contains(lower(sdirs{min(numel(sdirs),subid+1):end}),'run-')==1,1,'last');
    RUN    = fparts{runid};
  end

  % weighting
  if isMOD
    fparts = strsplit(sdirs{end},'_'); 
    modid  = cat_io_contains(lower(fparts),{'t1w','t2w','pdw','flair'});
    RUN    = strcat(fparts{modid});
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
function out = testfunction( testcase, job , subfield)
%%

  if ~exist('testcase','var')
    testcase = 1;
  else
    if isnumeric(testcase)
      testcase = max(0,min(3,testcase));
    else
      testcase = 1; 
    end
  end
  if ~exist('job','var') || ~isfield(job,'extopts'), job.extopts.resdircase = testcase; end
  if ~exist('subfield','var') || isempty(subfield), subfield = 'surfpath'; end

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

  % run it
  job.extopts.verbBIDS  = 1; 
  switch job.extopts.resdircase
    case 0, job.output.BIDS.BIDSno              = struct();
    case 1, job.output.BIDS.BIDSyes.BIDSfolder  = fullfile('../derivatives','CATdefbids');
    case 2, job.output.BIDS.BIDSrel.BIDSfolder  = fullfile('../derivatives','CATrelbids');
    case 3, job.output.BIDS.relative.BIDSfolder = fullfile('../derivatives','CATreldirs');
  end
    
  % show result
  out = cat_io_BIDS(files,job,subfield,'prefix','lh.central.','suffix','.topofix','ext','');

end
