function out = cat_io_BIDS(files, job, subfield, varargin)
%Function to analyse files and setup related BIDS directories and files. 
%
%  out = cat_io_BIDS(files, job, action, varargin)
%
%  files          .. input files  or  call of test dataset 
%
%  job            .. CAT preprocessing job variable
%   .output       .. CAT GUI structure that is always prefered to extopts 
%   .extopts      .. CAT extopts struture partially defined in cat_defaults.m
%    .BIDS_folder .. "../derivatives/CAT#" 
%    .subfolder   .. 0 - none, 1 - use catsubfolder (default)   
%    .logdircase  .. 1 - main common directory of all inputs
%                    2 - current working directory (default)
%                    3 - first input directory 
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
%
%  subfield       .. field selctor to get only this field in the output structur
%
%  varargin       .. input for spm_file to add and maninpulate the filename 
%
%
%  out            .. full BIDS structure or BIDS.(subfield) element 
%                    Variables with:
%                     - main_*:    are the same in all cases
%                     - *dir(s)*:  describe diretories 
%                     - *path*:    describe full paths 
%   .main_ressubdir   ..
%   .main_BIDSfolder  .. see job.extopts.bids_folder (in cat_default)
%   .main_logdircase  .. see job.extopts.logdircase
%   .main_mkBIDSdir   .. see job.extopts.mkBIDSdir
%   .main_npredirs    .. additional directories defined by the number of 
%                        relative directoires "../" in job.extopts.BIDS_folder
%
%   .files        .. cleaned up input files
%
%   .SUB          .. BIDS subject directory or filename if noBIDS
%   .SES          .. BIDS session directory or empty
%   .ANA          .. BIDS anat directory or empty
%   .RUN          .. BIDS run filename specifier for rescans
%   .MOD          .. MRI modality (T1w,T2w,PDw,FLAIR)
%
%   .isBIDS       .. definition if the input is or the setup requires BIDS
%   .is[SUB|SES|ANA|RUN|MOD] .. availability of a BIDS entry 
%
%   .postdirs     .. BIDS subdirectories, eg. sub*/ses*/anat
%   .predirs      .. addition relative directories ahead of the main 
%                    BIDS directory (see out.main_npredirs)
%   .npostdirs    .. number of BIDS subdirectories, eg., 3
%
%   .rawBIDSpath  .. path of the raw BIDS folder 
%                    (e.g, the parent of sub in case of BIDS)
%   .resBIDSpath  .. path of the new BIDS derivative folder (output)
% 
%   .main_logdir  .. common directory for log files
%    
%   .[label|mri|surf|report|err]path .. name of specific result directories
%   .[label|mri|surf|report|err]path .. path to specific result directories
%
%
% Examples: 
%  1) Call testfunction for a set of BIDS / nonBIDS files:
%      files = cat_io_BIDS( 'testfiles' ) % get testfiles
%
%  2) Get full BIDS structure for the first 3 test files
%      cat_io_BIDS( files(1:3) )
%
%  3) Create the surface filename for the fifth file: 
%      cat_io_BIDS( files(5), '', 'surfpath', 'prefix', 'lh.central.', 'suffix', '.topofix','ext','')
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
    out(:).main_logpath = out(1).main_logpath;
    return;
  elseif ischar(files)
    if strcmp(files,'files')
      out = testfunction(files,job,subfield); 
      return
    else
      files = cellstr(files); 
    end
  elseif isnumeric(files)
    % call testfunction 
    out = testfunction(files,job,subfield); 
    return
  else
    files = cellstr(files); 
  end


  % == main relative (BIDS) result directory ==
  % setup variables for subfolders (e.g, mri/surf) and 
  % the relative BIDSfolder_rel (e.g., ../derivatives/CAT0815)
  BIDS.main_BIDSfolder = job.extopts.bids_folder;
  BIDS.main_resdircase = job.extopts.resdircase;
  BIDS.main_catfolders = job.extopts.subfolders;
  BIDS.main_logdircase = job.extopts.logdircase;
  BIDS.main_mkBIDSdir  = job.extopts.mkBIDSdir; 
  
  % update by GUI output variable
  if isfield(job,'output') && isfield(job.output,'BIDS')
    if isfield(job.output.BIDS,'BIDSno') 
      BIDS.main_BIDSfolder = '';
      BIDS.main_resdircase = 0; 
    elseif isfield(job.output.BIDS,'BIDSyes') 
      BIDS.main_BIDSfolder = job.output.BIDS.BIDSyes.BIDSfolder;
      BIDS.main_resdircase = 1; 
    elseif isfield(job.output.BIDS,'BIDSrel')
      BIDS.main_BIDSfolder = job.output.BIDS.BIDSrel.BIDSfolder;
      BIDS.main_resdircase = 2;
    elseif isfield(job.output.BIDS,'relative')  
      BIDS.main_BIDSfolder = job.output.BIDS.relative.BIDSfolder;
      BIDS.main_resdircase = 3;
    end
  end

  % In the default BIDS case, no-relative path is allowed!
  if BIDS.main_resdircase==0
    if contains( BIDS.main_BIDSfolder , ['..' filesep] )
      cat_io_cprintf('warn','Warning:  Relative BIDS folder are not supported for "autobids" (resdircase=1)! \n')
    end
    BIDS.main_BIDSfolder = strrep(BIDS.main_BIDSfolder,['..' filesep],''); 
  end

  % count the number of wished recursive directories
  % that we will have to add after the BIDS result folder together with 
  % the BIDS subject subfolders
  BIDS.main_npredirs = max([0,strfind(BIDS.main_BIDSfolder,['..' filesep])]);

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
     BIDS.postdirs{fi,1}, BIDS.npostdirs(fi,1)] = getBIDSpostdirs(sdirs);
  
    % define if output can or has to be BIDS
    BIDS.isBIDS(fi,1) = BIDS.isSUB(fi) & ~isempty(BIDS.main_BIDSfolder) & BIDS.main_resdircase ~= 3; 

    % avoid multiple derivatives and remove the leading part
    if BIDS.isBIDS(fi)
      sdirs = clearDoubleDerivatives(sdirs, BIDS.main_BIDSfolder, BIDS.npostdirs(fi));
    end

    % detect extra directories
    subhome = max(1,numel(sdirs) - BIDS.npostdirs(fi)); 
    BIDS.rawBIDSpath{fi,1} = [ repmat( filesep, isempty(sdirs{1}), isempty(sdirs{1})) ...
      fullfile( sdirs{ 1:max( 1, subhome - BIDS.main_npredirs - 1 ) } )]; 
    if BIDS.main_npredirs>1
      BIDS.predirs{fi,1} = fullfile( sdirs{ max( 2, subhome - BIDS.main_npredirs) : max(2,subhome - 1) } ); 
    else
      BIDS.predirs{fi,1} = '';
    end

    
    % define main resultdir and subdirs
    BIDS.resdir{fi,1} = fullfile( strrep( BIDS.main_BIDSfolder , ...
      ['..' filesep], ''), BIDS.predirs{fi} ); 
    BIDS.resBIDSpath{fi,1} = fullfile( BIDS.rawBIDSpath{fi} , ...
      strrep( BIDS.main_BIDSfolder ,['..' filesep], '') , BIDS.predirs{fi} ); 
    for sfi = 1:numel(catfolder)
      % just the dirname
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'dir']){fi,1} = catfolder{sfi}; 
      else
        BIDS.([catfolder{sfi} 'dir']){fi,1} = ''; 
      end
    end  
   
    BIDS.mainpath{fi,1} = fullfile( BIDS.rawBIDSpath{fi} , BIDS.resdir{fi,1} , BIDS.postdirs{fi});  
    for sfi = 1:numel(catfolder)
      % full path
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'path']){fi,1} = fullfile( BIDS.rawBIDSpath{fi} , ...
          BIDS.resdir{fi,1} , BIDS.postdirs{fi}, catfolder{sfi} ); 
      else
        BIDS.([catfolder{sfi} 'path']){fi,1} = fullfile( BIDS.rawBIDSpath{fi} , ...
          BIDS.resdir{fi,1} , BIDS.postdirs{fi} ); 
      end
    end      
  end


  % == final directory setup & creation ==
  BIDS.main_logpath = prepare_logpath(BIDS); 
  create_resultdirs(BIDS,job); 


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
function main_logpath = prepare_logpath(BIDS)
%prepare_logdir.

  switch BIDS.main_logdircase
    case 1
      % one main directory for all cases
      if numel(sfiles) > 1
        [~,S] = spm_str_manip(files,'C'); 
        main_logpath = fullfile(S.s, BIDS.main_BIDSfolder, 'log');
      else
        main_logpath = fullfile(BIDS.resBIDSpath{fi}, 'log');
      end
    case 2
      % write into current directory
      main_logpath = fullfile(pwd, 'log');
    case 3
      % write into first subject directory
      main_logpath = fullfile(BIDS.resBIDSpath{1}, 'log');
  end
        
  if cat_io_contains(main_logpath, spm('dir'))
    % in case of the SPM directory we file this in a dataspecific subdirectory
    main_logpath = fullfile(spm('dir'),'toolbox','CAT', 'logs', char(datetime('now','Format','yyyyMMdd')) );
    
    % give some feedback if directories are created (only then otherwise the warning comes with every call of this function)
    if BIDS.main_mkBIDSdir 
      cat_io_cprintf('warn','cat_io_BIDS:  Creation of log directory in SPM/CAT path redirected to: \n  %s\n', main_logpath)
    end
  end  

  % == create log and main result directory == 
  if ~exist(main_logpath,'dir')  &&  BIDS.main_mkBIDSdir 
    try
      mkdir(main_logpath); 
    catch
      main_logpath0 = main_logpath;
      main_logpath  = fullfile(BIDS.rawBIDSpath{1}, BIDS.main_ressubdir, 'log');
      cat_io_cprintf('warning','cat_io_BIDS: Creation of log directory "%s" failed. \nUse path of first subject "%s".\n', ...
        main_logpath0, main_logpath)
    end
  end

end
%==========================================================================
function create_resultdirs(BIDS,job)
  % for longitudinal we need update

  if ~BIDS.main_mkBIDSdir; return; end

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

  % create data specific outputs 
  for fi = 1:numel(BIDS.files)
    if ~exist(BIDS.mainpath{fi},'dir')
      mkdir(BIDS.mainpath{fi}); 
    end
  
    if ~exist(BIDS.labelpath{fi},'dir') && ( writeVatlases || writeSatlases )
      mkdir(BIDS.labelpath{fi}); 
    end

    if ~exist(BIDS.mripath{fi},'dir') && writeVols
      mkdir(BIDS.mripath{fi}); 
    end

    if ~exist(BIDS.surfpath{fi},'dir') && ...
      isfield(job,'output') && isfield(job.output,'surface') && job.output.surface
      mkdir(BIDS.surfpath{fi}); 
    end

    if ~exist(BIDS.reportpath{fi},'dir')
      mkdir(BIDS.reportpath{fi}); 
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

  if isempty(job) && ~isstruct(job), clear job; job = struct(); end

  def.extopts             = cat_get_defaults('extopts');
 %def.extopts.bids_folder = cat_get_defaults('extopts.bids_folder');  % already defined
  def.extopts.resdircase  = 1; % 0-noBIDS, 1-onlyBIDS, 2-allways
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
  BIDSsubpath, BIDSsubpathdepth] = getBIDSpostdirs(sdirs)
%getpostdirs. Get information about possible subject, session, and anat directories.  
% Detect BIDS directories by looking for key words used in to define
% the subject, session, and the anatomic directory.
%
%   [SUB, SES, ANA, RUN, isSUB, isSES, isANA, isRUN, BIDSsubpath, BIDSsubpathdepth] = getpostdirs(sdirs)
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
%testfunction. Create a set of test files. 

  if ~exist('testcase','var')
    testcase = 1;
  else
    if isnumeric(testcase)
      testcase = max(0,min(3,testcase));
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

  if ~isnumeric(testcase) && strcmp(testcase,'testfiles')
    out = files; 
    return; 
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
