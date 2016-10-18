function cat_tst_cattest(job)
%  CAT test script. 
%  _____________________________________________________________________
%  
%  The idea is to create a common test script that can be run by any user.
%  It include the processing routines of important CAT GUI and background
%  functions with default (and modified) parameters for one (non)human 
%  subject (or multiple) subjects.
% 
%   cat_tst_cattest(job)
%   
%  Data:
%  1) Human:
%   * Adult:            single_T1subj  (Collins brain)
%                       Tohoku,CG,...? 
%   * Old:              OASIS31        (with WMHs)
%                       IXI?           (without WMHs)
%   * Infant:           Berlin? 
%                       IXI?
%   * Tumor:            tb09
%                       tp01
%  2) Primates (~1.5 mm resolution):
%   * Greater Apes:     chimpanzee_laz
%                       gorilla_kekla
%   * Lesser Apes:      gibbon_cleo (use oldworld monkey template)
%   * Oldworld Monkey:  rhesus_caretF99 (atlas)
%                       baboon_F3S12 / mangabey_fso
%
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id$



%  _____________________________________________________________________
%
%  Development framework / versions:
%  ---------------------------------------------------------------------
%  Further scripts:
%    *) Volume mapping from template to surface space
%       Volume mapping from surface to template space
%    *) Volume to surface mapping (individual surface)
%       Volume to surface mapping (template surface)
%    *) Projection tools 
%    *) Longitudinal Processing >> Data
%
%  Further functionality:
%    *) GUI with:
%           [job.resdir]   > cattestdir is ok ... not in general
%           job.userlevel  [ default | expert | developer ] 
%                          > maybe create a hard copy for the default user at least for preprocessing?
%           job.datalevel  [ basic | multisubject | long | non-human | full ]
%           job.paralevel  [ basic | enhanced | ... ] >> userlevel
%    *) processing with different important parameters
%         (e.g., NCstr, LASstr, gcutstr, cleanupstr, WHMC(str), vox)
%  	 *) Help Skripts, e.g., Volume vs. Surface Smoothing
%
%  Large extensions
%  V2:    statistical functions >> required futher data 
%          - we can add some further low res images
%          - we can link some sources that allow fast and simple access
%            and have good quality (IXI, OASIS)
%  V3:    bad paramter input as internal test function to test error
%         handling
%
%  Planed data extensions:
%  ---------------------------------------------------------------------
%  * human example datasets (~?1.0 MB / subject) 
%    - anatomy:   very old, very young (Berlin?), WMH (OASIS31?), tumor
%    ? quality:   strong bias, noisy, artefacts
%    ? AD:        10x10 subjects
%    ? aging:     10 subjects
%  * Adding of primate example datasets (~?0.5 - 1.5 MB / subject)
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
  def.computer  = computer; 
  def.userlevel = ''; % [ default | expert | developer ] 
  def.datalevel = ''; % [ basic | human | long | ape | monkey | full ]
  def.resdir    = fullfile(spm('dir'),'toolbox','cat12','cattest',[cv 'R' rv]);
  def.batchdir  = fullfile(spm('dir'),'toolbox','cat12','batches','cattest');

  % testdata definition
  % --------------------------------------------------------------------
  % single run human test data
  def.data_human = {
    fullfile(spm('dir'),'canonical','single_subj_T1.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_4397-tfl.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_BUSS_2002_1YO_t1.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x2_OAS1_0031_MR1_mpr_n4_anon_sbj_111.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x1_human_tumor_tb09.nii'); 
    fullfile(spm('dir'),'toolbox','cat12','data','uint8_lowresR2x2x1_human_tumor_tp01.nii'); 
    };
  % longitudinal test data 
  def.data_human_long = {
    };
  % dataset for tests of statistical functions  
  def.data_human_group = {
  };
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
  % test parameter
  def.para = {
    ... scipt, testlevel, variable, further values
    ''                      1  0 ''                                           ''           {0}; % default
    'cat12_SBM_101_segment' 1  2 'spm.tools.cat.estwrite.extopts.APP'         'APP'        {0}; 
    'cat12_SBM_101_segment' 1  1 'spm.tools.cat.estwrite.extopts.sanlm'       'sanlm'      {0 2}; 
    'cat12_SBM_101_segment' 1  1 'spm.tools.cat.estwrite.extopts.NCstr'       'NCstr'      {0 0.5 1}; 
    'cat12_SBM_101_segment' 1  1 'spm.tools.cat.estwrite.extopts.LASstr'      'LASstr'     {0 1}; 
    'cat12_SBM_101_segment' 1  2 'spm.tools.cat.estwrite.extopts.gcutstr'     'gcutstr'    {0 1}; 
    'cat12_SBM_101_segment' 1  2 'spm.tools.cat.estwrite.extopts.cleanupstr'  'cleanupstr' {0 1}; 
    'cat12_SBM_101_segment' 1  3 'spm.tools.cat.estwrite.extopts.BVCstr'      'BVCstr'     {0.5 1}; 
    'cat12_SBM_101_segment' 1  3 'spm.tools.cat.estwrite.extopts.WMHCstr'     'WMHCcstr'   {0 1}; 
    'cat12_SBM_101_segment' 1  3 'spm.tools.cat.estwrite.extopts.WMHC'        'WMHCstr'    {0 2}; 
    'cat12_SBM_101_segment' 1  2 'spm.tools.cat.estwrite.extopts.vox'         'vox'        {1}; 
    'cat12_SBM_101_segment' 1  2 'spm.tools.cat.estwrite.extopts.pbtres'      'pbtres'     {0.25 1};   
  };
  job = cat_io_checkinopt(job,def);
  
  
  % choose datalevel
  % --------------------------------------------------------------------
  if isempty(job.datalevel)
    job.datalevel = char(spm_input('Datalevel',1,'min|basic|human|prim|all', ...
      {'minimal','basic','human','primates','all'},1));
  end
  switch job.datalevel
    case 'minimal'
      job.data_human            = job.data_human(1);
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {}; 
    case 'basic'
      job.data_human            = job.data_human(1); 
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_greaterapes      = job.data_greaterapes(1); 
      job.data_oldworldmonkeys  = job.data_oldworldmonkeys(1); 
    case 'human'
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {};
    case 'long'
      job.data_human            = {};
      job.data_human_group      = {};
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {};
    case 'group'
      job.data_human            = {};
      job.data_human_long       = {};
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {};
    case 'ape'
      job.data_human            = {};
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_oldworldmonkeys  = {};
    case 'monkey'
      job.data_human            = {};
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_greaterapes      = {};
    case 'primates'
      job.data_human            = {};
      job.data_human_long       = {};
      job.data_human_group      = {};
  end
  
  % check input files
  species = {'data_human','data_human_long','data_human_group','data_greaterapes','data_oldworldmonkeys'}; 
  for si = 1:numel(species)
    for j = numel(job.(species{si})):-1:1
      if ~exist(job.(species{si}){j},'file'), 
        job.(species{si}){j} = []; 
        %warning('ERROR:cat_tst_single:noExistingData','The input data "%s" does not exist!\n',job.(species{si}){j}); 
      end
    end  
  end
  
  
  
  
  % choose userlevel
  % --------------------------------------------------------------------
  userlevels = {'default','expert','developer'};
  if isempty(job.userlevel)
    job.userlevel = spm_input('Userlevel',1,'Default|Expert|Developer',[0,1,2],1);
  end
  for pi = size(job.para,1):-1:1
    if job.para{pi,3}>job.userlevel, job.para(pi,:) = []; end
  end
  exp = job.userlevel; 
  if cat_get_defaults('extopts.expertgui')~=job.userlevel
    cat12(userlevels{job.userlevel+1});
  end
      
    
  for pi=1:size(job.para)
    for ppi=1:numel(job.para{pi,6})
      % create output directory and copy files
      % the different species required another preprocessing, but the following
      % routines are identical!
      if isempty(job.para{pi,5})
        subresdir = [job.resdir '_' job.computer];
      else
        subresdir = [job.resdir job.para{pi,5} num2str(job.para{pi,6}{ppi}) '_' job.computer];
      end  
      if ~exist(subresdir,'dir'), mkdir(subresdir); end 

      for si = 1:numel(species)
        eval(sprintf('files_%s = {};',species{si}(6:end)));  
        for j = numel(job.(species{si})):-1:1
          copyfile(job.(species{si}){j},subresdir); 
          [pp,ff,ee] = fileparts(job.(species{si}){j});
          eval(sprintf('files_%s{j,1} = fullfile(subresdir,[ff ee]);',species{si}(6:end)));  
        end
      end
      files = [files_human; files_greaterapes; files_oldworldmonkeys]; 
      if isempty(files), return; end


      % cat subdirs
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

      % add batch dir & get batches 
      addpath(job.batchdir); 
      batches = cat_vol_findfiles(job.batchdir,'*.m'); 

      % load batches and create main batch
      mainbatch{pi}{ppi} = {};
      batchname         = {}; 
      for bi = 1:numel(batches)
        matlabbatch = {};
        [pp,ffbi]  = fileparts(batches{bi}); 

        % load batch
        eval(ffbi);

        if strcmp(job.para{pi,1},ffbi)
          if isnumeric(job.para{pi,6}{ppi})
            eval(sprintf('matlabbatch{%d}.%s = %f;',job.para{pi,2},job.para{pi,4},job.para{pi,6}{ppi}));
          else
            error('bad parameter ... coding required')
          end
        end
        if ~isempty(matlabbatch)
          for mbi = 1:numel(matlabbatch)
            mainbatch{pi}{ppi}{end+1,1} = matlabbatch{mbi};
            
            batchname{end+1,1} = ffbi;
            batchname{end,2}   = mbi; 
            batchname{end,3}   = size(batchname,1);  
          end
        end
      end
    end
  end
  
  
  %% test compiled (and other) functions
  compile(0,1,job.userlevel+1);
  
  
  
  
  %% process main batch
  spm_jobman('initcfg');
  perror = {}; 
  for pi=1:size(job.para)
    for ppi=1:numel(job.para{pi,6})
      perror{pi}{ppi} = 2*ones(numel(mainbatch{pi}{ppi}),1);
      for mbi = 1:numel(mainbatch{pi}{ppi})
        try 
          spm_jobman('run',mainbatch{pi}{ppi}(mbi));  
          perror{pi}{ppi}(mbi)=0;
        catch
          perror{pi}{ppi}(mbi)=1;
        end
      end
    end
  end

  
  
  
  %% status report
  fprintf('CAT Test:\n')
  fprintf('  userlevel:  %s\n',userlevels{job.userlevel+1});
  fprintf('  datalevel:  %s\n',job.datalevel);
  for fi=1:numel(files)
    fprintf('    %s\n',files{fi});
  end
  fprintf('  paralevel: %s\n',job.paralevel);
  for fi=1:numel(files)
    fprintf('    %s\n',files{fi});
  end
  fprintf('\n');
  for pi=1:size(job.para)
    for ppi=1:numel(job.para{pi,6})
      fprintf('\n%45s: %s\n',sprintf('Script %s %0.2f',job.para{pi,5},job.para{pi,6}{ppi}),'Status');
      for mbi = 1:numel(perror{pi}{ppi})
        fprintf('%2d)%40s%02d: ',batchname{mbi,3},batchname{mbi,1},batchname{mbi,2}); 
        if perror{pi}{ppi}(mbi)==2
          cat_io_cprintf([1.0 0.5 0.0],'unprocessed\n'); 
        elseif perror{pi}{ppi}(mbi)==1
          cat_io_cprintf([0.8 0.0 0.0],'FAILED\n'); 
        else
          cat_io_cprintf([0.0 0.6 0.0],'OK\n'); 
        end    
      end
    end
  end
end

