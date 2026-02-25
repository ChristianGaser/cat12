function out = cat_io_BIDS(files,job,action)
%
%  out = cat_io_BIDS(files,job,action)
%
%  sfiles   .. input files 
%  action   .. field selctor to get only this field in the output structur
%  job      .. CAT preprocessing job variable
%   .output.BIDS         .. CAT GUI structure that is always prefered
%   .extopts.BIDS_folder .. 
%   .extopts.resdircase  .. 0 - BIDSdir only if BIDS 
%                           1 - BIDSdir allways, 
%                           2 - BIDSdir as pure relative dir also in case of BIDS  
%
%  out          .. output file structure 
%   .BIDSdir    .. BIDS resultdirectory (/deriatives/CAT)
%   .files      .. cleaned up input files
%   .isBIDS     .. one if input file is BIDS
%  BIDSsub  .. subject/file name
%  devdir   .. main data directory (parent of sub in case of BIDS)
%              with derivatives!
%  rdevdir  .. relative directory (e.g. sub*/ses*/anat*)
%  logdir   .. common directory for log files
%  mdevdir  .. main data directory (parent of sub in case of BIDS)
%              without derivatives
%

  % == check input ==
  if ~exist('action','var'), action = ''; end
  if ~exist('job','var'), job = struct(); end


  % prepare/complete job variable
  job = updateJob(job);

  
  if ischar(files) && size(files,1) > 1
    % Creation of a data record out(:) for each input. 
    out = cat_io_BIDSchar(files, job, action); 
    return; 
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
  
  % update by GUI output variable
  if isfield(job,'output') && isfield(job.output,'BIDS')
    if isfield(job.output.BIDS,'BIDSyes') 
      BIDS.main_ressubdir  = job.output.BIDS.BIDSyes.BIDSfolder;
      BIDS.main_resdircase = 0; 
    elseif isfield(job.output.BIDS,'BIDSrel')
      BIDS.main_ressubdir  = job.output.BIDS.BIDSrel.BIDSfolder;
      BIDS.main_resdircase = 1;
    elseif isfield(job.output.BIDS,'relative')  
      BIDS.main_ressubdir  = job.output.BIDS.relative.BIDSfolder;
      BIDS.main_resdircase = 2;
    end
  end

  % In the default BIDS case, no-relative path is allowed!
  if BIDS.main_resdircase==0
    if contains( BIDS.main_ressubdir , ['..' filesep] )
      cat_io_cprintf('warn','Warning:  Relative BIDS folder are not supported for autobids! \n')
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
    [BIDS.SUB{fi,1}, BIDS.SES{fi,1}, BIDS.ANA{fi,1}, ...
     BIDS.isSUB(fi,1), BIDS.isSES(fi,1), BIDS.isANA(fi,1), ...
     BIDS.BIDSsubdirs{fi,1}, BIDS.BIDSsubdirnumber(fi,1)] = getBIDSsubdirs(sdirs);
  
    % avoid multiple derivatives and remove the leading part
    sdirs = clearDoubleDerivatives(sdirs,BIDS.main_ressubdir);

    % detect extra directories
    subhome            = max(1,numel(sdirs) - BIDS.BIDSsubdirnumber(fi)); 
    BIDS.BIDSrawpath{fi,1} = [filesep*isempty(sdirs{1}) fullfile( sdirs{ 1:max( 1, subhome - BIDS.main_adddirs - 1 ) } )]; 
    if BIDS.main_adddirs
      BIDS.adddirs{fi,1} = fullfile( sdirs{ max( 2, subhome - BIDS.main_adddirs) : max(2,subhome - 1) } ); 
    else
      BIDS.adddirs{fi,1} = '';
    end

    % define if output can or has to be BIDS
    BIDS.isBIDS(fi,1) = BIDS.isSUB(fi) & ~isempty(BIDS.main_ressubdir) & BIDS.main_resdircase ~= 2; 

    % define main resultdir and subdirs
    BIDS.resdir{fi,1}     = fullfile( strrep( BIDS.main_ressubdir ,['..' filesep], ''), BIDS.adddirs{fi} ); 
    BIDS.resdirpath{fi,1} = fullfile( BIDS.BIDSrawpath{fi} , strrep( BIDS.main_ressubdir ,['..' filesep], '') , BIDS.adddirs{fi} ); 
    for sfi = 1:numel(catfolder)
      % just the dirname
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'dir']) = catfolder{sfi}; 
      else
        BIDS.([catfolder{sfi} 'dir']) = ''; 
      end
    end      
    for sfi = 1:numel(catfolder)
      % full path
      if BIDS.main_catfolders && ~BIDS.isBIDS(fi)
        BIDS.([catfolder{sfi} 'path'])   = fullfile( BIDS.BIDSrawpath{fi} , BIDS.resdir{fi,1} ,  BIDS.BIDSsubdirs{fi}, catfolder{sfi} ); 
      else
        BIDS.([catfolder{sfi} 'path'])   = fullfile( BIDS.BIDSrawpath{fi} , BIDS.resdir{fi,1} ,  BIDS.BIDSsubdirs{fi} ); 
      end
    end      
  end


  % == log-directory setup ==
  switch BIDS.main_logdircase
    case 1
      % one main directory for all cases
      if numel(sfiles) > 1
        [~,S] = spm_str_manip(files,'C'); 
        BIDS.main_logpath = fullfile(S.s, BIDS.main_ressubdir, 'log');
      else
        BIDS.main_logpath = fullfile(BIDS.BIDSrawpath{fi}, BIDS.main_ressubdir, 'log');
      end
    case 2
      % write into current directory
      if strcmp(pwd,fullfile(spm('dir'),'toolbox','CAT'))
        % in case of the CAT directory we file this in a dataspecific subdirectory
        BIDS.main_logpath = fullfile(pwd, 'logs', char(datetime('now','Format','yyyyMMdd')) );
      else
        BIDS.main_logpath = fullfile(pwd, BIDS.main_ressubdir,'log');
      end
    case 3
      % write into first subject directory
      BIDS.main_logpath = fullfile(BIDS.BIDSrawpath{fi}, BIDS.main_ressubdir, 'log');
  end


  % == create log and main result directory == 
%### update  
  if 0
    if ~exist(logdir,'dir')
      try
        mkdir(logdir); 
      catch
        if cat_io_contains( logdir , spm('dir') )
          logdirfailed = logdir; 
          logdir = fullfile(pdir{1}, BIDSdir, 'log');
          cat_io_cprintf('warning',['cat_run: Creation of log directory "%s" failed. \n' ...
            'Write into first subject directory: "%s"'],logdirfailed,logdir)
        else
          cat_io_cprintf('warning','cat_run: Creation of log directory "%s" failed. \n',logdir)
        end
      end
    end

    for sfi = 1:numel(sfiles)
      if ~exist(devdir{sfi},'dir'); mkdir(devdir{sfi}); end
    end

    % for the subdirs we should maybe check if we need them. 

  end


  % == specify output ==
  if isempty( action )
    out = BIDS; 
  else
    if isfield( BIDS , action )
       out = BIDS.(action); 
    else
      error('cat_io_BIDS:unkownField', ...
        'Unknown field "%s". \n', action)
    end
  end
end
%==========================================================================
function out = cat_io_BIDSchar(files, job, action)
%cat_io_BIDSchar. Recursive call to create a data record out(f).
  for fi = 1:size(files,1)
    if isempty(action) % struct
      out(fi) = cat_io_BIDS({files(fi,:)},job,action); %#ok<*AGROW>
    else
      out1 = cat_io_BIDS({files(fi,:)},job,action); 
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

  if isempty(job) || ~isfield(job,'extopts')
    job.extopts = cat_get_defaults('extopts');
  end
  
  if ~isfield(job.extopts,'subfolders') % CAT subfolders
    job.extopts.subfolders = cat_get_defaults('extopts.subfolders');
  end

  if ~isfield(job.extopts,'bids_folder') % main (BIDS) result folder
    job.extopts.bids_folder = '';
  end

  if ~isfield(job.extopts,'resdircase') 
    job.extopts.resdircase = 0; 
  end

  if ~isfield(job.extopts,'logdircase') 
    job.extopts.logdircase = 2; 
  end
end
%==========================================================================
function sdirs = clearDoubleDerivatives(sdirs,ressubfolder)
%clearDoubleDerivatives. Remove double derivatives directories. 
% sdirs = clearDoubleDerivatives(sdirs)

  dev = find(strcmp('derivatives',sdirs));
  if ~isempty(ressubfolder) && ~isempty(dev)
    % ../derivatives/d1/derivatives/d2/.. >> ../derivatives/d2/..
    sdirs(dev(1):dev(end)-1) = []; 
  end
end
%==========================================================================
function [SUB, SES, ANA, isSUB, isSES, isANA, BIDSsubpath, BIDSsubpathdepth] = getBIDSsubdirs(sdirs)
%getBIDSsubdirs. Get information about possible subject, session, and anat directories.  
% Detect BIDS directories by looking for key words used in to define
% the subject, session, and the anatomic directory.
%
%   [SUB, SES, ANA, isSUB, isSES, isANA, BIDSsubpath, BIDSsubpathdepth] = getBIDSsubdirs(sdirs)
%

  SUB = spm_str_manip(sdirs{end},'r');  % use filename in case of non-BIDS (check for uniqueness?)
  SES = ''; 
  ANA = ''; 
  
  isANA = 0; 
  isSUB = 0; 
  isSES = 0; 
  
  if numel(sdirs) > 1
  % get the subpath with sub-*[/ses-*][/anat]
    isANA = cat_io_contains(sdirs{end-1}(1:min(4,numel(sdirs{end-1}))),'anat');
    isSES = any(cat_io_contains(sdirs{end-isANA-1}(1:min(4,numel(sdirs{end-1-isANA}))),{'ses-','ses_'}));
    isSUB = any(cat_io_contains(sdirs{end-isANA-1-isSES}(1:min(4,numel(sdirs{end-1-isANA-isSES}))),{'sub-','sub_'})); 
   
    if isANA, ANA = sdirs{end-1}; end
    if isSES, SES = sdirs{end-1-isANA}; end
    if isSUB, SUB = sdirs{end-1-isANA-isSES}; end

    % very tolerant sub defintion
    if ~isSES && ~isSUB
      subcan = cat_io_contains(sdirs(1:end-1),{'sub-','sub_'});
      for si = sort(find(subcan),'des')
        isSUB = (numel(sdirs) - si) * any(cat_io_contains(sdirs{si}(1:min(4,numel(sdirs{si}))),{'sub-','sub_'})); 
        if isSUB>0, break; end
      end
      
    end
  end
 
  BIDSsubpathdepth = isANA + isSES + isSUB; 
  if BIDSsubpathdepth > 0
    BIDSsubpath = fullfile( sdirs{ end - BIDSsubpathdepth : end - 1} );
  else
    BIDSsubpath = ''; 
  end
end 
%==========================================================================
function testfunction
%%
  files = {
    '/Users/tomcat/MRData/BIDStest/ProjectA_fulllong/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii'                                                    
    '/Users/tomcat/MRData/BIDStest/ProjectA_fulllong/sub-01/ses-02/anat/sub-01_ses-02_T1w.nii'                                                    
    '/Users/tomcat/MRData/BIDStest/ProjectA_fulllong/sub-02/ses-01/anat/sub-02_ses-01_T1w.nii'                                                    
    '/Users/tomcat/MRData/BIDStest/ProjectB_noDirs/sub-01_T1.nii'                                                                                 
    '/Users/tomcat/MRData/BIDStest/ProjectB_noDirs/sub-02_T1.nii'                                                                                 
    '/Users/tomcat/MRData/BIDStest/ProjectC_sameLongName_noAnat/sub-01/BL/sub-01_T1w.nii'                                                         
    '/Users/tomcat/MRData/BIDStest/ProjectC_sameLongName_noAnat/sub-01/FU/sub-01_T1w.nii'                                                         
    '/Users/tomcat/MRData/BIDStest/ProjectC_sameLongName_noAnat/sub-02/BL/sub-02_T1w.nii'                                                         
    '/Users/tomcat/MRData/BIDStest/ProjectC_sameLongName_noAnat/sub-02/FU/sub-02_T1w.nii'                                                         
    '/Users/tomcat/MRData/BIDStest/ProjectD_sameName_multiAnat/sub-01/ses-00/anat/T1w.nii'                                                        
    '/Users/tomcat/MRData/BIDStest/ProjectD_sameName_multiAnat/sub-01/ses-00/anat2/T1w.nii'                                                       
    '/Users/tomcat/MRData/BIDStest/ProjectD_sameName_multiAnat/sub-02/ses-00/anat/T1w.nii'                                                        
    '/Users/tomcat/MRData/BIDStest/ProjectD_sameName_multiAnat/sub-02/ses-00/anat2/T1w.nii'                                                       
    '/Users/tomcat/MRData/BIDStest/ProjectE_sameName/sub-01/ses-01/T1w.nii'                                                                       
    '/Users/tomcat/MRsdirsData/BIDStest/ProjectE_sameName/sub-02/ses-02/T1w.nii'                                                                       
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir1/mdir2/sub-01_T1w.nii'                                                                           
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir1/mdir2/sub-02_T1w.nii'                                                                           
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/mwp1sub-01_T1w.nii'                                          
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/mwp1sub-02_T1w.nii'                                          
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/mwp2sub-01_T1w.nii'                                          
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/mwp2sub-02_T1w.nii'                                          
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/p0sub-01_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/p0sub-02_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/wmsub-01_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/wmsub-02_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/y_sub-01_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/derivatives/CAT12.9_2565/mri/y_sub-02_T1w.nii'                                            
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/mdir5/sub-01_T1w.nii'                                                                     
    '/Users/tomcat/MRData/BIDStest/ProjectF/mdir3/mdir4/mdir5/sub-02_T1w.nii'                                                                     
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-01/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                           
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-01/sub-02/ses-00/anat/sub-02_ses-01_T1w.nii'                                           
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                           
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/derivatives/CAT12.9_2565/sub-01/ses-00/anat/mwp1sub-01_ses-01_T1w.nii'
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/derivatives/CAT12.9_2565/sub-01/ses-00/anat/mwp2sub-01_ses-01_T1w.nii'
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/derivatives/CAT12.9_2565/sub-01/ses-00/anat/p0sub-01_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/derivatives/CAT12.9_2565/sub-01/ses-00/anat/wmsub-01_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-01/ses-00/derivatives/CAT12.9_2565/sub-01/ses-00/anat/y_sub-01_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/anat/sub-02_ses-01_T1w.nii'                                           
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/derivatives/CAT12.9_2565/sub-02/ses-00/anat/mwp1sub-02_ses-01_T1w.nii'
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/derivatives/CAT12.9_2565/sub-02/ses-00/anat/mwp2sub-02_ses-01_T1w.nii'
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/derivatives/CAT12.9_2565/sub-02/ses-00/anat/p0sub-02_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/derivatives/CAT12.9_2565/sub-02/ses-00/anat/wmsub-02_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectG_addGroup/group-02/sub-02/ses-00/derivatives/CAT12.9_2565/sub-02/ses-00/anat/y_sub-02_ses-01_T1w.nii'  
    '/Users/tomcat/MRData/BIDStest/ProjectH_derivatives/derivatives/sub-01/ses-00/anat/sub-01_ses-01_T1w.nii'                                     
    '/Users/tomcat/MRData/BIDStest/ProjectH_derivatives/derivatives/sub-02/ses-00/anat/sub-02_ses-01_T1w.nii'                                    
  };

  %%
  clc, clear job; 
  %job.output.BIDS.BIDSyes.BIDSfolder = fullfile('../derivatives','CATnogroup');
  job.output.BIDS.BIDSrel.BIDSfolder    = fullfile('../derivatives','CATgroup');
  %job.output.BIDS.relative.BIDSfolder    = fullfile('../derivatives','CATrel');
  out = cat_io_BIDS(files(1),job)

end
