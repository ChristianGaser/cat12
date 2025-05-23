function [mainbatch,perror] = cat_tst_cattest(job)
%  CAT test script. 
%  _____________________________________________________________________
%  
%  The idea is to create a common test script that can be run by any user.
%  It include the processing routines of important CAT GUI and background
%  functions with default (and modified) parameters for one (non)human 
%  subject (or multiple) subjects.
% 
%   [mainbatch,perror] = cat_tst_cattest(job)
%
%   job
%    .datalevel  .. test diffenent datasets
%                     working - {'basic'}
%                     planned - {'human','long','group','nonhuman'}
%    .paralevel  .. test different parameter settings
%                     1 - default, 2 - userlevel, 3 - full
%    .userlevel  .. test specific userlevel
%                     0 - default, 1 - expert, 2 - developer
%    .debug      .. using debuging mode
%
%  Examples:
%    cat_tst_cattest(struct('datalevel','basic','paralevel',1,'userlevel',0))
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$

 

%  _____________________________________________________________________
%
%  Development framework / versions:
%  ---------------------------------------------------------------------
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
%         (e.g. NCstr, LASstr, gcutstr, cleanupstr, WHMC(str), vox)
%  	 *) Help Skripts, e.g. Volume vs. Surface Smoothing
%
%  Large extensions
%  V2:    statistical functions >> required further data 
%          - we can add some further low res images
%          - we can link some sources that allow fast and simple access
%            and have good quality (IXI, OASIS)
%  V3:    bad parameter input as internal test function to test error
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
%  # delete data
%  # evaluation concept? > GT? > BWP?
%  _____________________________________________________________________


%#ok<*ASGLU,*NASGU,*AGROW>
  clc
  %clear
  spm_clf('Interactive'); 

  if nargin==0, help cat_tst_cattest; return; end 

  if ~exist('job','var'), job=struct(); end
  [cv,rv] = cat_version;
  
  % defaults
  def.computer  = computer; 
  def.resdir    = fullfile(fileparts(mfilename('fullpath')),'cattest',[cv 'R' rv]);
  def.batchdir  = fullfile(fileparts(mfilename('fullpath')),'batches','cattest');
  def.expert    = cat_get_defaults('extopts.expertgui');

  % testdata definition
  % --------------------------------------------------------------------
  % single run human test data
  def.data_human = {
    fullfile(spm('dir'),'canonical','single_subj_T1.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_4397-tfl.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_BUSS_2002_1YO_t1.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_OAS1_0031_MR1_mpr_n4_anon_sbj_111.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x1_human_tumor_tb09.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x1_human_tumor_tp01.nii'); 
    };
  % longitudinal test data 
  def.data_human_long = {
    };
  % dataset for tests of statistical functions  
  def.data_human_group = {
  };
  % primate test data
  def.data_greaterapes = {
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_chimpanzee_kenge.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_chimpanzee_laz.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_orangutan_minyak.nii');
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_gorilla_kekla.nii');
    };
  def.data_oldworldmonkeys = {
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_baboon_F3S12s20130924_120712.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_gibbon_cleo.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_gibbon_gibbon4.nii');
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_mangabey_fso.nii'); 
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_rhesus_caretF99.nii');
    fullfile(fileparts(mfilename('fullpath')),'data','uint8_lowresR2x2x2_primate_rhesus_research.nii');
    };
  % test parameter (especially for CAT preprocessing)
  def.para = {
    ... 'scriptname', [batchnumber], [testlevel], 'variable', {further values};  % comments
    ''                       1  0 ''                                             ''           {0};                       % default
    ...
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.opts.tpm'                           'tpm' ...       
      {{fullfile(spm('dir'),'tpm','TPM.nii')} ...  
       {fullfile(cat_get_defaults('extopts.pth_templates'),'TPM_Age11.5.nii')} ...                                     % children
      };     
    ... SPM bias
    ... SPM acc
    ... == parameter in default GUI ==
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.vox'                        'vox'        {2 1};            % def=1.5;
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.segmentation.APP'           'APP'        {0 1};     % def=1070;
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.segmentation.LASstr'        'LASstr'     {0 0.01 1.0};     % def=0.5;
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.segmentation.gcutstr'       'gcutstr'    {0 0.5 -1};       % def=2;
    ...'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.registration.darteltpm'  'template'   { ... 
    ...'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.registration.darteltpm'  'template'   { ... 
    ...  {cat_get_defaults('extopts.pth_templates'),'Template_0_IXI555_MNI152_GS.nii')}};   % Shooting template
    'cat12_101_MAIN_segment' 1  1 'spm.tools.cat.estwrite.extopts.registration.regstr'        'regstr'     {0.5 4};          % def=0;  
    ... == parameter in expert GUI ==
    'cat12_101_MAIN_segment' 1  2 'spm.tools.cat.estwrite.extopts.segmentation.cleanupstr'    'cleanupstr' {0.0 0.01 1.0};   % def=0.5;
    %{
    'cat12_101_MAIN_segment' 1  2 'spm.tools.cat.estwrite.extopts.segmentation.NCstr'         'NCstr'      {0 1 2 4};        % def=inf; 
    ... the restype and resval variable were used in the GUI to create the restypes structure that does not exist in the cat_defaults file
    'cat12_101_MAIN_segment' 1  2 'spm.tools.cat.estwrite.extopts.registration.restypes'      'restypes'   { ...
        struct('native',{})        ... native resolution - no changes
        struct('fixed' ,[1.0 0.1]) ... interpolation to 1 mm
        struct('fixed' ,[0.8 0.1]) ... interpolation to 0.8 mm 
        struct('best'  ,[0.5 0.1]) ... interpolation to best resolution but not more than 0.5 mm 
      };       
    ... == parameter in developer GUI ==
    'cat12_101_MAIN_segment' 1  2 'spm.tools.cat.estwrite.extopts.surface.pbtres'             'pbtres'     {1.00 0.25};      % def=0.5;  
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.segmentation.mrf'           'mrf'        {0.00 0.30};      % def=1; auto  
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.segmentation.BVCstr'        'BVCstr'     {0.5 0.01 1.0};   % def=0;
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.segmentation.WMHCstr'       'WMHCstr'    {0.0 0.01 1.0};   % def=0.5;
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.segmentation.WMHC'          'WMHC'       {0 2 3};          % def=1;
    ...'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.tca'           'tca'        {0 2};            % def=1;
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.admin.subfolders'           'subfolders' {0};              % def=1;
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.admin.verb'                 'verb'       {0 1 3};          % def=2;
    ... admin options
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.admin.expertimental'        'expertimental'     {1};              % def=0;
    'cat12_101_MAIN_segment' 1  3 'spm.tools.cat.estwrite.extopts.admin.ignoreErrors'         'ignoreErrors'      {1};              % def=0;
    ...
    %}
    ... 
    ... S1173
    %{
    'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.opts.tpm'                           'tpm' ...       
      {{fullfile(spm('dir'),'tpm','TPM.nii')} ...  
       {fullfile(cat_get_defaults('extopts.pth_templates'),'TPM_Age11.5.nii')} ...                                     % children
      };     
    ... SPM bias
    ... SPM acc
    ... == parameter in default GUI ==
    'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.vox'                        'vox'        {2 1};            % def=1.5;
'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.APP'                        'APP'        {4}; %0 1 2 3 4};     % def=1070; ERROR;
    'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.LASstr'                     'LASstr'     {0 0.01 1.0};     % def=0.5; OK;
    'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.gcutstr'                    'gcutstr'    {0 0.5 -1};       % def=2; OK;
    ...'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.registration.darteltpm'  'template'   { ... 
    ...'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.registration.darteltpm'  'template'   { ... 
    ...  {fullfile(cat_get_defaults('extopts.pth_templates'),'Template_0_IXI555_MNI152_GS.nii')}};   % Shooting template
    'cat12_105_MAIN_segment_1173plus' 1  1 'spm.tools.cat.estwrite1173plus.extopts.registration.regstr'        'regstr'     {0.5 4};          % def=0; OK;  
    ... == parameter in expert GUI ==
    %}
    %{
    'cat12_105_MAIN_segment_1173plus' 1  2 'spm.tools.cat.estwrite1173plus.extopts.segmentation.cleanupstr'    'cleanupstr' {0.0 0.01 1.0};   % def=0.5;
    'cat12_105_MAIN_segment_1173plus' 1  2 'spm.tools.cat.estwrite1173plus.extopts.segmentation.NCstr'         'NCstr'      {0 1 2 4};        % def=inf; 
    ... the restype and resval variable were used in the GUI to create the restypes structure that does not exist in the cat_defaults file
    'cat12_105_MAIN_segment_1173plus' 1  2 'spm.tools.cat.estwrite1173plus.extopts.registration.restypes'      'restypes'   { ...
        struct('native',{})        ... native resolution - no changes
        struct('fixed' ,[1.0 0.1]) ... interpolation to 1 mm
        struct('fixed' ,[0.8 0.1]) ... interpolation to 0.8 mm 
        struct('best'  ,[0.5 0.1]) ... interpolation to best resolution but not more than 0.5 mm 
      };       
    ... == parameter in developer GUI ==
    'cat12_105_MAIN_segment_1173plus' 1  2 'spm.tools.cat.estwrite1173plus.extopts.surface.pbtres'             'pbtres'     {1.00 0.25};      % def=0.5;  
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.segmentation.mrf'           'mrf'        {0.00 0.30};      % def=1; auto  
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.segmentation.BVCstr'        'BVCstr'     {0.5 0.01 1.0};   % def=0;
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.segmentation.WMHCstr'       'WMHCstr'    {0.0 0.01 1.0};   % def=0.5;
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.segmentation.WMHC'          'WMHC'       {0 2 3};          % def=1;
    ... admin options
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.admin.subfolders'           'subfolders' {0};              % def=1;
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.admin.verb'                 'verb'       {0 1 3};          % def=2;
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.admin.expertimental'        'expertimental'     {1};              % def=0;
    'cat12_105_MAIN_segment_1173plus' 1  3 'spm.tools.cat.estwrite1173plus.extopts.admin.ignoreErrors'         'ignoreErrors'      {1};              % def=0;
    %}
    }; 
  def.userlevel   = ''; 
  def.datalevel   = ''; 
  def.paralevel   = ''; 
  def.segmentonly = 1; % it is much faster, if we avoid surface reconstruction in most cases
  def.fixres      = 2; 
  job = cat_io_checkinopt(job,def);
  
  %{ 
  fprintf([
    'The datalevel controls the number and kind (e.g. animals) of test cases.\n' ... 
    'The userlevel controls the basic setting of the test parameters,\n' ... 
    'whereas the parameterlevel controls the definition of test parameters.']),
  %}
  
  % choose datalevel
  % --------------------------------------------------------------------
  % if there are multipe subjects an GUI for default user maybe useful
  if isempty(job.datalevel)
    job.datalevel = char(spm_input('Datalevel',1,'basic|human|animal|all', ...
      {'basic','human','primates','all'},1));
  end
  switch job.datalevel
    case 'basic'
      job.data_human            = job.data_human(1);
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {}; 
    case 'basica'
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
    case 'nonhuman'
      job.data_human            = {};
      job.data_human_long       = {};
      job.data_human_group      = {};
    otherwise % basic
      job.data_human            = job.data_human(1);
      job.data_human_long       = {};
      job.data_human_group      = {};
      job.data_greaterapes      = {};
      job.data_oldworldmonkeys  = {}; 
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
  
  
  
  % choose userlevel (only lower or equal level than the actual level)
  % --------------------------------------------------------------------
  if isempty(job.userlevel)
    userlevelsstr = 'Default|Expert|Developer'; 
    userlevels    = eval(sprintf('{''%s''}',strrep(lower(userlevelsstr),'|',''','''))); 
    job.userlevel = spm_input('Userlevel','+1',userlevelsstr,0:numel(userlevels)-1,1);
  else
    userlevels = {'Default','Expert','Developer'};   
  end
  exp = job.userlevel; 
  if cat_get_defaults('extopts.expertgui')~=job.userlevel
    cat12(userlevels{job.userlevel+1});
  end
  
  
  % choose parameter test level
  % --------------------------------------------------------------------
  % basic       = only defaults
  % userlevel   = only parameter defined for this userlevel without defaults
  % full        = all parameters defined for this and lower userlevels
  paralevelstr = 'Defaults|Userlevel|Full';
  if isempty(job.paralevel)
    paralevels   = eval(sprintf('{''%s''}',strrep(lower(paralevelstr),'|',''','''))); 
    job.paralevel = spm_input('Parameter test level','+1',paralevelstr,0:numel(paralevels)-1,1);
  end
  switch job.paralevel
    case 0 % just defaults
      for pi = size(job.para,1):-1:1
        if job.para{pi,3}>0, job.para(pi,:) = []; end
      end
    case 1 
      for pi = size(job.para,1):-1:1
        if (job.para{pi,3})==job.userlevel, job.para(pi,:) = []; end
      end
      job.para(1,:) = []; 
    case 2 
      for pi = size(job.para,1):-1:1
        if (job.para{pi,3}-1)>job.userlevel, job.para(pi,:) = []; end
      end
      job.para(1,:) = []; 
  end
  % 
  fprintf('Testparameter:\n'); 
  for pi=1:numel(job.para(:,5)), fprintf('  %s:%s\n',job.para{pi,1},job.para{pi,5}); end
  fprintf('\n');
  
  
  
  %%  
  for pi=1:size(job.para)
    for ppi=1:numel(job.para{pi,6})
      % create output directory and copy files
      % the different species required another preprocessing, but the following
      % routines are identical!
      if isempty(job.para{pi,5})
        subresdir = [job.resdir '_' job.computer];
      else
        if isnumeric(job.para{pi,6}{ppi})
          subresdir = [job.resdir '_' job.computer '_' job.para{pi,5} num2str(job.para{pi,6}{ppi}(1,:))];
        else
          subresdir = [job.resdir '_' job.computer '_' job.para{pi,5} num2str(ppi)];
        end
      end  
      if ~exist(subresdir,'dir'), mkdir(subresdir); end 
      cd(subresdir);

      
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

      clearvars cat cat1173 cat1173plus
      cat_get_defaults; cat_get_defaults1173; cat_get_defaults1173plus;
      
      files_human1173     = spm_file(files_human,'prefix','S1173_');
      files_human1173plus = spm_file(files_human,'prefix','S1173plus_');
      for fi = 1:numel(files_human1173)
        copyfile(files_human{fi},files_human1173{fi}); 
        copyfile(files_human{fi},files_human1173plus{fi}); 
      end
      
      % cat subdirs
      [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders;

      % add batch dir & get batches 
      addpath(job.batchdir); 
      batches = cat_vol_findfiles(job.batchdir,'*.m'); 

      % load batches and create main batch
      mainbatch{pi}{ppi} = {};
      batchname{pi}{ppi} = {}; 
      RBM = 1; 
      SBM = 1;     
      
      for bi = 1:numel(batches)
        matlabbatch = {};
        [pp,ffbi]  = fileparts(batches{bi}); 

        % load batch
        try
          eval(ffbi);
        catch
          fprintf('eval batch %d failed: %s\n',bi,ffbi);
          continue
        end
        if strcmp(job.para{pi,1},ffbi)
          if isnumeric(job.para{pi,6}{ppi})
            eval(sprintf('matlabbatch{%d}.%s = [%s];',job.para{pi,2},job.para{pi,4},sprintf('%f ',job.para{pi,6}{ppi})));
          elseif ischar(job.para{pi,6}{ppi}) 
            eval(sprintf('matlabbatch{%d}.%s = ''%s'';',job.para{pi,2},job.para{pi,4},job.para{pi,6}{ppi}));
          elseif iscellstr(job.para{pi,6}{ppi}) 
            eval(sprintf('matlabbatch{%d}.%s = {''%s''};',job.para{pi,2},job.para{pi,4},job.para{pi,6}{ppi}{1}));
          elseif isstruct(job.para{pi,6}{ppi})
            eval(sprintf('matlabbatch{%d}.%s = job.para{pi,6}{ppi};',job.para{pi,2},job.para{pi,4}));
          else
            error('bad parameter ... coding required')
          end
        end
        
        % use only 2 mm processing with 
        if job.fixres 
          for mbi=numel(matlabbatch):-1:1; 
            pp = getpp(matlabbatch,mbi);
%%
             if iscell(matlabbatch{mbi})
               % SPM surf pipeline does not contrain further test 
               for smbi=1:numel(matlabbatch{mbi})
                 if ~isempty(pp) && isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp),'extopts') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp),'output') %&& ...
                   %isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).output,'surface') && ...
                   %~matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).output.surface
                   if isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts,'segmentation') && ...
                      isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.segmentation,'vox') 
                     matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.segmentation.vox      = job.fixres;
                     matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.segmentation.restypes = struct('fixed',[job.fixres 0]); 
                   elseif isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts,'vox') 
                     matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.vox      = job.fixres;
                    matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.restypes = struct('fixed',[job.fixres 0]); 
                   end
                 end
               end
             else
               if ~isempty(pp) && isfield(matlabbatch{mbi}.spm.tools.cat.(pp),'extopts') && ...
                  isfield(matlabbatch{mbi}.spm.tools.cat.(pp),'output') %&& ...
                  %isfield(matlabbatch{mbi}.spm.tools.cat.(pp).output,'surface') && ...
                  %~matlabbatch{mbi}.spm.tools.cat.(pp).output.surface
                 if isfield(matlabbatch{mbi}.spm.tools.cat.(pp).extopts,'segmentation') && ...
                    isfield(matlabbatch{mbi}.spm.tools.cat.(pp).extopts.segmentation,'vox')
                   matlabbatch{mbi}.spm.tools.cat.(pp).extopts.segmentation.vox      = job.fixres;
                   matlabbatch{mbi}.spm.tools.cat.(pp).extopts.segmentation.restypes = struct('fixed',[job.fixres 0]); 
                 elseif isfield(matlabbatch{mbi}.spm.tools.cat.(pp).extopts,'vox')
                   matlabbatch{mbi}.spm.tools.cat.(pp).extopts.vox      = job.fixres;
                   matlabbatch{mbi}.spm.tools.cat.(pp).extopts.restypes = struct('fixed',[job.fixres 0]); 
                 end
               end
             end
          end
        end
           
       
        % set/unset surface processing
        if strcmp(job.para{pi,5},'pbtres') || isempty(job.para{pi,5})
          for mbi=numel(matlabbatch):-1:1
            pp = getpp(matlabbatch,mbi);
            
            % no surface-tools
            if iscell(matlabbatch{mbi})
              % SPM surf pipeline does not contrain further test 
              for smbi=1:numel(matlabbatch{mbi})
                % deactivate surface processing
                if ~isempty(pp) && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp),'output') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).output,'surface') 
                  if job.userlevel > 0
                    switch job.para{pi,5}
                      case 'pbtres',  matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).extopts.pbtres = job.para{pi,6}{ppi}; 
                    end
                  end
                  matlabbatch{mbi}{smbi}.spm.tools.cat.(pp).output.surface = 1; 
                  SBM = 1; 
                end
              end
           else
              % deactivate surface processing
              if ~isempty(pp) && ...
                 isfield(matlabbatch{mbi}.spm.tools.cat.(pp),'output') && ...
                 isfield(matlabbatch{mbi}.spm.tools.cat.(pp).output,'surface')
                if job.userlevel > 0
                  switch job.para{pi,5}
                    case 'pbtres',  matlabbatch{mbi}.spm.tools.cat.(pp).extopts.pbtres = job.para{pi,6}{ppi};
                  end
                end
                matlabbatch{mbi}.spm.tools.cat.(pp).output.surface = 1; 
                SBM = 1; 
              end
            end
          end
        elseif job.segmentonly
        % to avoid time intensive surface processing ... 
          %mainbatch{pi}{ppi}{1}.spm.tools.cat.(pp).output.surface = 0; 
          for mbi=numel(matlabbatch):-1:1; 
            % no surface-tools
            pp = getpp(matlabbatch,mbi);
            if iscell(matlabbatch{mbi})
              %{
              % SPM surf pipeline does not contrain further test ...
              for smbi=1:numel(matlabbatch{mbi})
                % deactivate surface processing
                if isfield(matlabbatch{mbi}{smbi},'spm') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm,'tools') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools,'cat') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat,'estwrite_spm') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp)_spm,'output') && ...
                   isfield(matlabbatch{mbi}{smbi}.spm.tools.cat.(pp)_spm.output,'surface') 
                  matlabbatch{mbi}{smbi}.spm.tools.cat.(pp)_spm.output.surface = 0; 
                  SBM = 0; 
                end
              end
              %}
            else
              % deactivate surface processing
              if ~isempty(pp) && ...
                 isfield(matlabbatch{mbi}.spm.tools.cat.(pp),'output') && ...
                 isfield(matlabbatch{mbi}.spm.tools.cat.(pp).output,'surface') 
                matlabbatch{mbi}.spm.tools.cat.(pp).output.surface = 0; 
                SBM = 0; 
              end
            end
          end
        end
        
        
        % no RBM tools if not ROI processing 
        for mbi=numel(matlabbatch):-1:1
          if ~isempty(pp) && ...
             isfield(matlabbatch{mbi}.spm.tools.cat.(pp),'output') && ...
             isfield(matlabbatch{mbi}.spm.tools.cat.(pp).output,'ROI') && ...
             matlabbatch{mbi}.spm.tools.cat.(pp).output.ROI==0;
            RBM = 0; 
          end
        end
        
        % continue if batch required surfaces/ROIs and not surface/ROIs are available
        if ~SBM && ~isempty(strfind(ffbi,'SBM')); continue; end
        if ~RBM && ~isempty(strfind(ffbi,'RBM')); continue; end
        if ~isempty(job.para{pi,1}) && ~strcmp(job.para{pi,1},ffbi); continue; end
        
        if ~isempty(matlabbatch)
          for mbi = 1:numel(matlabbatch)
            mainbatch{pi}{ppi}{end+1,1} = matlabbatch{mbi};
            
            batchname{pi}{ppi}{end+1,1} = ffbi;
            batchname{pi}{ppi}{end,2}   = mbi; 
            batchname{pi}{ppi}{end,3}   = size(batchname,1);  
            batchname{pi}{ppi}{end,4}   = matlabbatch{mbi};
          end
        end
      end
    end
  end
  
  
  %% test compiled (and other) functions
  compile(0,1,job.userlevel+1);
  spm_jobman('initcfg'); 
  
  
  
  %% process main batch
  col    = [0.0 0.2 0.6]; 
  perror = {}; 
  cat_io_cprintf('silentreset');
  for pi=1:size(job.para,1)
    for ppi=1:numel(job.para{pi,6})
      %%
      cat_io_cprintf(col,'\n\n\n------------------------------------------------------------------------\n');
      if isempty(job.para{pi,5}), cat_io_cprintf(col,sprintf('% 3d/%d) defaults ',pi,size(job.para,1)));
      else cat_io_cprintf(col,sprintf('% 3d/%d) %s = ',pi,size(job.para,1),job.para{pi,5}));
      end
      if ~isempty(job.para{pi,5})  
        if      isnumeric(job.para{pi,6}{ppi}),   cat_io_cprintf(col,sprintf('%0.2f ', job.para{pi,6}{ppi})); 
        elseif  ischar(job.para{pi,6}{ppi}),      cat_io_cprintf(col,sprintf('%s'    , job.para{pi,6}{ppi})); 
        elseif  iscell(job.para{pi,6}{ppi}),      cat_io_cprintf(col,sprintf('%s'    , job.para{pi,6}{ppi}{1}));  
        end
        cat_io_cprintf(col,sprintf(' (%d of %d)',ppi,numel(job.para{pi,6})));
      end
      cat_io_cprintf(col,'\n------------------------------------------------------------------------\n'); 
          
      %%
      perror{pi}{ppi} = 2*ones(numel(mainbatch{pi}{ppi}),1);
      for mbi = 1:numel(mainbatch{pi}{ppi})
        cat_io_cprintf([0 0.5 0],sprintf('\n  %s %d ',...
          batchname{pi}{ppi}{mbi,1}(7:end),batchname{pi}{ppi}{mbi,2}));
        try 
          if job.debug
            % -- debuging code 20161118--
            try
              fprintf('\nOutput Surface: %d', mainbatch{pi}{ppi}{mbi}.spm.tools.cat.(pp).output.surface); 
            end
            try
              fprintf('\nOutput ROI:     %d', mainbatch{pi}{ppi}{mbi}.spm.tools.cat.(pp).output.ROI);
            end
            if iscell(mainbatch{pi}{ppi}{mbi})
              for ci=1:numel(mainbatch{pi}{ppi}{mbi})
                try
                  fprintf('\nOutput Surface: %d', mainbatch{pi}{ppi}{mbi}{ci}.spm.tools.cat.(pp).spm.output.surface); 
                end
              end
              for ci=1:numel(mainbatch{pi}{ppi}{mbi})
                try
                  fprintf('\nOutput ROI:     %d', mainbatch{pi}{ppi}{mbi}{ci}.spm.tools.cat.(pp).spm.output.ROI);
                end
              end
            end           
            perror{pi}{ppi}(mbi)=0;
          else
            spm_jobman('run',mainbatch{pi}{ppi}(mbi));  
            perror{pi}{ppi}(mbi)=0;
          end
        catch e
          cat_io_cprintf([0.8 0 0],'ERROR\n\n');
          disp(e)
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
      if isnumeric(job.para{pi,6}{ppi}),  parastr = num2str(job.para{pi,6}{ppi});
      elseif ischar(job.para{pi,6}{ppi}), parastr = job.para{pi,6}{ppi};
      else                                parastr = '';
      end
      if ppi==1, fprintf('\n%42s: %s\n','Parameter','Status'); end
      if ppi==1 || ppi==numel(job.para{pi,6})
        fprintf('--------------------------------------------------\n');
      end
      fprintf('%42s: ',sprintf('%s = %s',job.para{pi,5},parastr));
      for mbi = 1:numel(perror{pi}{ppi})
        %fprintf('%40s%02d: ',batchname{pi}{ppi}{mbi,1},batchname{pi}{ppi}{mbi,2}); 
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
function pp = getpp(matlabbatch,mbi) 
  if isfield(matlabbatch{mbi},'spm') && ...
    isfield(matlabbatch{mbi}.spm,'tools') && ...
    isfield(matlabbatch{mbi}.spm.tools,'cat')
    if isfield(matlabbatch{mbi}.spm.tools.cat,'estwrite')
      pp = 'estwrite'; 
    elseif isfield(matlabbatch{mbi}.spm.tools.cat,'estwrite1173')
      pp = 'estwrite1173';
    elseif isfield(matlabbatch{mbi}.spm.tools.cat,'estwrite1173plus')
      pp = 'estwrite1173plus';
    elseif isfield(matlabbatch{mbi}.spm.tools.cat,'estwrite_spm')
      pp = 'estwrite_spm';
    else 
      pp = ''; 
    end
  else 
    pp = '';
  end
end
