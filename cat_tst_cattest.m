function cat_tst_cattest(job)
%  _____________________________________________________________________
%  CAT test script. 
% 
%   cat_tst_cattest(job)
%   
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id: cat_run_job.m 1013 2016-09-22 11:49:13Z dahnke $
 
%  _____________________________________________________________________
%  The idea is to create a common test script that can be run by any user.
%  Different handling of user mode, i.e., the default user does not need 
%  to test expert functions? 
%   - Better to test all, because default users maybe switch to experts!
%   - Better to avoid the full expert output!
%   > GUI switch?
%
%  Development framework / versions:
%  ---------------------------------------------------------------------
%  V1:    processing of important CAT GUI functions with default parameter
%         for one subject
%         1) VBM/SBM preprocessing
%         2) Volume mapping from template to surface space
%            Volume mapping from surface to template space
%         3) Surface parameter estimation 
%            Volume to surface mapping (individual surface)
%            Volume to surface mapping (template surface)
%            Resample and smooth of surface data
%            S
%         4) Mapping tools 
%         ...
%  V2:    GUI with:
%           [job.resdir]   > cattestdir is ok 
%           job.userlevel  [ default | expert | developer ] 
%                          > maybe create a hard copy for the default user at least for preprocessing?
%           job.datalevel  [ basic | multisubject | long | non-human | full ]
%           job.paralevel  [ basic | enhanced | ... ]
%                          ... focus on cat_main paramter ...
%  V3:    processing of primate data 
%  V4:    processing with different important parameters
%         (e.g., NCstr, LASstr, gcutstr, cleanupstr, WHMC(str), vox)
%  V5:    statistical functions >> required futher data 
%          - we can add some further low res images
%          - we can link some sources that allow fast and simple access
%            and have good quality (IXI, OASIS)
%  V6:    Volume vs. Surface Smoothing
%  V*:    bad paramter input as internal test function to test error
%         handling
%
%  Planed data extensions:
%  ---------------------------------------------------------------------
%  * human example datasets (~ 1.0 MB / subject) 
%    - anatomy:   very old, very young (Berlin?), WMH (OASIS31?), tumor
%    ? quality:   strong bias, noisy, artefacts
%    ? AD:        10x10 subjects
%    ? aging:     10 subjects
%  * Adding of primate example datasets (~ 0.5 - 1.5 MB / subject)
%    1-? larger (kekla,molek,kenge/lana/laz,lorel)
%    1-? smaller (cleo, ...)
%  * Adding of T2/PD exampled datasets
%  * Adding of other modalities 
%  * Adding of low quality datasets?
%    - extreme noise 
%    - motion artifacts
%  * Adding of highres data? 
%    - CG, RD, Tohoku, ... with resolution test?
%  _____________________________________________________________________


%#ok<*ASGLU,*NASGU,*AGROW>
  clear

  if ~exist('job','var'), job=struct(); end
  [cv,rv] = cat_version;

  % defaults
  def.userlevel = 1; % [ default | expert | developer ] 
  def.datalevel = 1; % [ basic | multisubject | long | non-human | full ]
  def.paralevel = 1; % [ basic | enhanced | ... ]
  def.resdir    = fullfile(spm('dir'),'toolbox','cat12','cattest',[cv 'R' rv]);
  def.batchdir  = fullfile(spm('dir'),'toolbox','cat12','batches','cattest');
  % human test data
  def.data_human = {
    fullfile(spm('dir'),'canonical','single_subj_T1.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_4397-tfl.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_BUSS_2002_1YO_t1.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_OAS1_0031_MR1_mpr_n4_anon_sbj_111.nii'); 
    };
  def.data_human = {};
  % primate test data
  def.data_greaterapes = {
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_chimpanzee_kenge.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_chimpanzee_laz.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_orangutan_minyak.nii');
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_gorilla_kekla.nii');
    };
  def.data_oldworldmonkeys = {
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_baboon_F3S12s20130924_120712.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_gibbon_cleo.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_gibbon_gibbon4.nii');
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_mangabey_fso.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_rhesus_caretF99.nii');
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_primate_rhesus_research.nii');
    };
  job = cat_io_checkinopt(job,def);

  job.para = {
    ... scipt, testlevel, variable, further values
    'cat12_SBM_101_segment' 2 'spm.tools.cat.estwrite.extopts.APP'          {0}; 
    'cat12_SBM_101_segment' 1 'spm.tools.cat.estwrite.extopts.sanlm'        {0 2}; 
    'cat12_SBM_101_segment' 1 'spm.tools.cat.estwrite.extopts.NCstr'        {0 0.5 1}; 
    'cat12_SBM_101_segment' 1 'spm.tools.cat.estwrite.extopts.LASstr'       {0 1}; 
    'cat12_SBM_101_segment' 2 'spm.tools.cat.estwrite.extopts.gcutstr'      {0 1}; 
    'cat12_SBM_101_segment' 2 'spm.tools.cat.estwrite.extopts.cleanupstr'   {0 1}; 
    'cat12_SBM_101_segment' 3 'spm.tools.cat.estwrite.extopts.BVCstr'       {0.5 1}; 
    'cat12_SBM_101_segment' 3 'spm.tools.cat.estwrite.extopts.WMHCstr'      {0 1}; 
    'cat12_SBM_101_segment' 3 'spm.tools.cat.estwrite.extopts.WMHC'         {0 2}; 
    'cat12_SBM_101_segment' 2 'spm.tools.cat.estwrite.extopts.vox'          {1}; 
    'cat12_SBM_101_segment' 2 'spm.tools.cat.estwrite.extopts.pbtres'       {0.25 1};   
  };
    
  % add batch dir
  addpath(job.batchdir); 
  
  % check input files
  species = {'data_human','data_greaterapes','data_oldworldmonkeys'}; 
  for si = 1:numel(species)
    for j = numel(job.(species{si})):-1:1
      if ~exist(job.(species{si}){j},'file'), 
        job.(species{si}){j} = []; 
        %warning('ERROR:cat_tst_single:noExistingData','The input data "%s" does not exist!\n',job.(species{si}){j}); 
      end
    end  
  end
  
  % create output directory and copy files
  % the different species required another preprocessing, but the following
  % routines are identical!
  if ~exist(job.resdir,'dir'), mkdir(job.resdir); end 
  for si = 1:numel(species)
    eval(sprintf('files_%s = {};',species{si}(6:end)));  
    for j = numel(job.(species{si})):-1:1
      copyfile(job.(species{si}){j},job.resdir); 
      [pp,ff,ee] = fileparts(job.(species{si}){j});
      eval(sprintf('files_%s{j,1} = fullfile(job.resdir,[ff ee]);',species{si}(6:end)));  
    end
  end
  files = [files_human; files_greaterapes; files_oldworldmonkeys]; 
  if isempty(files), return; end
  
 
  % subdirs
  if cat_get_defaults('extopts.subfolders')
    mridir    = 'mri';
    surfdir   = 'surf';
    roidir    = 'label';
    reportdir = 'report';
  else
    mridir    = '';
    surfdir   = '';
    roidir    = '';
    reportdir = '';
  end  
  
  % get batches 
  batches = cat_vol_findfiles(job.batchdir,'*.m'); 

  % load batches and create maun batch
  mainbatch   = {};
  batchname   = {}; 
  for bi = 1:numel(batches)
    matlabbatch = {};
    [pp,ffbi]  = fileparts(batches{bi}); 
    
    % load batch
    eval(ffbi);
    
    if ~isempty(matlabbatch)
      for mbi = 1:numel(matlabbatch)
        mainbatch{end+1,1} = matlabbatch{mbi}; 
        batchname{end+1,1} = ffbi;
        batchname{end,2}   = mbi; 
        batchname{end,3}   = size(batchname,1);  
      end
    end
  end
  
  
  
  %% process mainbatch
  spm_jobman('initcfg');
  for mbi = 2 %numel(mainbatch)
    %%
    try
      spm_jobman('run',mainbatch(mbi));  
    catch e
      [mbi batchname(mbi,:)],
      rethrow(e);
    end
      
  end
end

