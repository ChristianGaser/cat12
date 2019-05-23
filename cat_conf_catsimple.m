function [catsimple,catsimplelong] = cat_conf_catsimple(expert)
% Configuration file for simplyfied preprocessing jobs
% _________________________________________________________________________
% This is a simpliefied batch for the whole CAT preprocessing.  It should
% extract various cortical information suitable for standard VBM, SBM, and 
% RBM analysis.  There should be no detailed settings for any preprocessing
% subfuction (e.g., skull-stripping, bias-corrections, registration, ...), 
% their processing quality or their output including only "save" measures
% such as GM volume, thickness, or curvature. 
% The DEP output should only contain usefull sets of the most important 
% statistical files with smoohted and unsmoothed data. 
% _________________________________________________________________________
% Robert Dahnke
% $Id$


% _________________________________________________________________________
% Internal comments:
%
% * Update smoothing default values
%
% * Add further modalities mapping (version 2)
%
% * Subfunction in some other conf batch ? 
%   > Maybe tbx_cfg_cat?
%
% * Futher references?
%   > No the link to the CAT paper is enought
%   > Update cite!
%   > simple cites by high numbers (use char):
%     0-10:  176, 185, 178, 179, 8308-8313
%     ,():   183, 8317 ,8318
% _________________________________________________________________________


%% Data
%  ------------------------------------------------------------------------


  % files
  data              = cfg_files;
  data.tag          = 'data';
  data.name         = 'Volumes';
  data.filter       = 'image';
  data.ufilter      = '.*';
  data.num          = [1 Inf];
  data.help         = {'Select one high resolution T1 image for each subject. '};

  
  % files long with two different selection schemes
  % - timepoints-subjects
  % - subjects-timepoints
  timepoint         = data; 
  timepoint.tag     = 'timepoints';
  timepoint.name    = 'Timepoint';
  timepoints        = cfg_repeat;
  timepoints.tag    = 'timepoints';
  timepoints.name   = 'Timepoints';
  timepoints.values = {timepoint};
  timepoints.num    = [2 Inf];
  timepoints.help   = {'Specify the same number and order of subjects for each timepoint. '};

  subjlong          = data; 
  subjlong.num      = [2 Inf];
  subjlong.tag      = 'subjects';
  subjlong.name     = 'Subject';
  subjlong.help     = {'Select all longitudinal high resolution T1 images of one subject. '};
  subjects          = cfg_repeat;
  subjects.tag      = 'subjects';
  subjects.name     = 'Subjects';
  subjects.values   = {subjlong};
  subjects.num      = [1 Inf];
  subjects.help     = {'Specify the same number and order of subjects for each timepoint. '};

  datalong          = cfg_choice;
  datalong.tag      = 'datalong';
  datalong.name     = 'Longitudinal data';
  datalong.values   = {timepoints subjects};
  datalong.val      = {timepoints};
  datalong.help     = {
   ['Select mode of longitudinal data selection for timepoints or subjects. ' ...
    'In case of timepoints you can create multiple timepoints where each timepoints has to contain the same number and order of subjects. ' ...
    'If you have a variing number of timepoints you have to use the subject mode where you have to define the files of each subject. ']
  }; 


  if expert>1
    % additional modalities
    % ---------------------------------------------------------------------
    % This is just the basic concept to support handling of further 
    % modalities in future! The different setting may require to use 
    % multiple mod-type fields for task-bask / resting-state fMRI and 
    % GM- / WM-focused sMRI data to make allow the users to select the
    % best fitting case. However a super mod that support more options
    % is maybe also useful (at least for experts). However, do not forget
    % that this has to be as simplest as possible!
    % For most modalities we furst have to develop a general analyse
    % scheme.
    % 
    % mod
    %  - [rs|tb]-frmi with non-linear coreg
    %    - name
    %    - source
    %    - images
    %  - smri with coreg
    %    - name
    %    - masking
    %    - images
    %  - smri without coreg
    %    - name
    %    - masking
    %    - images
    % ---------------------------------------------------------------------
    mname             = cfg_entry;
    mname.tag         = 'name';
    mname.name        = 'Name';
    mname.strtype     = 's';
    mname.num         = [0 20];
    mname.val         = {'MRI'};
    mname.help        = {
     ['Name identifier of this modality use for volumes (e.g., "[s#]mwMRI*.nii") and ' ...
      'surface data (e.g., "[s#mm.mesh.]MRI.*[.gii]") and ROIs (e.g. "MRI").']
      ''
    };

    fname             = mname;
    fname.val         = {'fMRI'};
    fname.help        = strrep(fname.help,'MRI',fname.val{1});

    sname             = mname;
    sname.val         = {'sMRI'};
    sname.help        = strrep(fname.help,'MRI',fname.val{1});

    mdata             = data; 
    mdata.tag         = 'data';
    mdata.name        = 'Data';
    mdata.help        = {'Specify the same number and order of subjects for each additional image modality'};


    % -- masking --
    maskth            = cfg_entry;
    maskth.tag        = 'maskth'; 
    maskth.name       = 'Threshold for masking';
    maskth.strtype    = 'r';
    maskth.num        = [1 1];
    maskth.val        = {0.5};
    maskth.help       = {'Percentual level for tissue masking, where 0.80 means that the value has to belong to the tissue glass in 80% of the subjects.';''}; 

    none              = cfg_branch;
    none.tag          = 'none';
    none.name         = 'No masking';
    none.help         = {'Use no tissue for masking of VBM data and volume ROIs.'};

    gmmask            = cfg_branch;
    gmmask.tag        = 'gm';
    gmmask.name       = 'Masking by GM tissue';
    gmmask.val        = {maskth};
    gmmask.help       = {'Use GM tissue for masking of VBM data and volume ROIs.'};

    wmmask            = cfg_branch;
    wmmask.tag        = 'wm';
    wmmask.name       = 'Masking by WM tissue';
    wmmask.val        = {maskth};
    wmmask.help       = {'Use WM tissue for masking of VBM data and volume ROIs.'};

    masking           = cfg_choice;
    masking.tag       = 'masking'; 
    masking.name      = 'Volumetric group masking';
    masking.values    = {none,gmmask,wmmask};
    masking.val       = {none};
    masking.help      = {'Use group masking with a specific threshold.'}; 


    % -- spatial normalization --
    nonlin            = cfg_menu;
    nonlin.tag        = 'reg'; 
    nonlin.name       = 'Regularisation of spatial normalization';
    nonlin.labels     = {'low','high','very high'};
    nonlin.values     = {1,10,100};
    nonlin.val        = {1};
    nonlin.help       = {'Apply spatial normalization to reduce non-linear warping of this modality. ';''}; 

    nonlin2           = cfg_menu;
    nonlin2.tag       = 'reg'; 
    nonlin2.name      = 'Spatial Normalization';
    nonlin2.labels    = {'none','very low','low','normal'};
    nonlin2.values    = {inf,100,10,1};
    nonlin2.val       = {1};
    nonlin2.help      = {'Apply spatial normalization to reduce non-linear warping of this modality. ';''}; 

    avg               = cfg_branch;
    avg.tag           = 'avg';
    avg.name          = 'Average image';
    avg.help          = {'Use the average of the 4D input dataset for coregistration.'};

    first             = cfg_branch;
    first.tag         = 'first';
    first.name        = 'First image';
    first.help        = {'Use first image of the 4D input dataset for coregistration.'};

    source            = cfg_files;
    source.tag        = 'source';
    source.name       = 'External source images';
    source.filter     = 'image';
    source.ufilter    = '.*';
    source.num        = [0 Inf];
    source.val        = {''};
    source.preview    = @(f) spm_image('Display',char(f));
    source.help       = {'Select images that are jiggled about to best match the reference (e.g. mean EPI, B0 images). '};

    stype             = cfg_choice;
    stype.tag         = 'stype'; 
    stype.name        = 'Sourcetype';
    stype.values      = {avg,first,source};
    stype.val         = {first};
    stype.help        = {
     ['Select type of source image. Use "First" if the first image is suited as source and should not be used further in dataanlysis, ' ...
      'select "External source image" to specify other image or "Average" if no special source is given and the 4D average of the input data should be used. ']}; 


    % -- modality types --
    mod               = cfg_branch;
    mod.tag           = 'mod';
    mod.name          = 'Modality with coregistration';
    mod.val           = {mname masking nonlin2 stype mdata}; 
    mod.help          = {'Select data files, type of masking, spatial normalization, type of source image, and name the modality.'};

    moda              = cfg_branch;
    moda.tag          = 'moda';
    moda.name         = 'Modality without coregistration';
    moda.val          = {mname masking mdata}; 
    moda.help         = {'Select data files, type of masking, and name the modality.'};

    % main
    mods              = cfg_repeat;
    mods.tag          = 'mods';
    mods.name         = 'Additional modalities for surface projection';
    mods.values       = {mod,moda};
    mods.val          = {};
    mods.num          = [0 Inf];
    mods.help         = {
     ['Specify the same number and order of subjects for each additional image modality ' ...
      '(e.g., 3D sMRI or 4D fMRI data) that should be projected to the brain surface. ']
    };


    % additional longidudinal modalities
    % - modality-timepoints-subjects
    % - modality-subjects-timepoints
    mdatalong         = datalong; 
    mdatalong.tag     = 'data';
    mdatalong.name    = 'Modality';
    mdatalong.help    = {'Specify the same number and order of subjects for each additional image modality'};

    longmod           = cfg_branch;
    longmod.tag       = 'mod';
    longmod.name      = 'Modality';
    longmod.val       = {mname masking mdatalong}; 
    longmod.help      = {'Select timpoint/subject files and name the modality.'};

    longmod           = cfg_branch;
    longmod.tag       = 'modc';
    longmod.name      = 'sMRI modality without coregistration';

    longmodc          = longmod;
    longmodc.tag      = 'modc';
    longmodc.name     = 'sMRI modality with coregistration';

    longfmri          = longmod; 
    longfmri.tag      = 'modf';
    longmod.val       = {fname masking mdatalong}; 
    longfmri.name     = 'fMRI modality with coregistration';


    longmods          = cfg_repeat;
    longmods.tag      = 'mods';
    longmods.name     = 'Additional modalities for surface projection';
    longmods.values   = {longfmri,longmodc,longmod};
    longmods.num      = [0 Inf];
    longmods.help     = {
     ['Specify the same number and order of subjects for each additional image modality ' ...
      '(e.g., 3D T2 images or 4D fMRI data) that should be projected to the brain surface. ' ...
      'If no individual surfaces were created the data is normalized and a Template brain ' ...
      'surface is used for extraction that is less accurate. ' ]
    };
  end
  


%% Parameter
%  ------------------------------------------------------------------------

  % CAT preprocessing version 
  catversion          = cfg_menu;
  catversion.tag      = 'catversion';
  catversion.name     = 'CAT preprocessing version';
  catversion.labels   = {'CAT12.1 (2017/09)','CAT12.6 (2019/03)'};
  catversion.values   = {'estwrite1173','estwrite1445'};
  if expert
    catversion.labels = [catversion.labels(1) {'CAT12.3 (2018/12)'}  catversion.labels(2:end)];
    catversion.values = [catversion.values(1) {'estwrite1173plus'}   catversion.values(2:end)];
  end
  catversion.val      = {'estwrite1445'};
  if expert>1
    catversion.labels = [catversion.labels {'CAT12.# (current)'}];
    catversion.values = [catversion.values {'estwrite'}];
    catversion.val    = {'estwrite'};
  end
  catversion.help     = {[
    'To expand previously process datasets select the same version of CAT preprocessing used before. ' ...
    'Do not mix different versions!' ...
    ];''};
  
  
  % tpm:
  % the TPM may further support primate preprocessing in future
  tpm               = cfg_menu;
  tpm.tag           = 'tpm';
  tpm.name          = 'Tissue Probability Map';
  tpm.labels        = {'Children','Adults'};
  tpm.values        = {'children','adults'};
  tpm.val           = {'adults'};
  tpm.help          = {[
    'CAT uses the tissue probability map (TPM) only for the initial SPM ' ...
    'segmentation using a prior independent AMAP appraoch.  Although even ' ...
    'the standard TPM of SPM give accurate and robust results in general, ' ...
    'it is recommended to use a children specific TPM with smaller skull. '] ...
    ''
    ... further information about the SPM TPM?
    ... further information about the child TPM?
    };

  
  surface           = cfg_menu;
  surface.tag       = 'surface';
  surface.name      = 'Surface-based processing';
  surface.labels    = {'No','Yes'};
  surface.values    = {0 1};
  surface.def       = @(val)cat_get_defaults('output.surface', val{:});
  surface.help      = {
    'Use surface-based preprocessing to estimate cortical surface, area, volume, and folding. '
    'Please note that surface reconstruction additionally requires about 20-60 min of computation time.'
    ''
  };


  [vROI,sROI]         = cat_conf_ROI(expert);    
  
  
  ignoreErrors        = cfg_menu;
  ignoreErrors.tag    = 'ignoreErrors';
  ignoreErrors.name   = 'Ignore errors';
  ignoreErrors.labels = {'No','Yes'};
  ignoreErrors.values = {0 1};
  ignoreErrors.def    = @(val)cat_get_defaults('extopts.ignoreErrors', val{:});
  ignoreErrors.help   = {
    'Catch preprocessing errors and move on with the next subject'
    ''
  };


  % nproc - unse only a menu for simpler access? - could maybe a problem if numcore miss some cores
  cores             = feature('numcores'); 
  if 0 % enter value
    nproc           = cfg_entry;
    nproc.strtype   = 'w';
    nproc.val       = {cores};
    nproc.num       = [1 1];
  else
    % choose by menu
    nproc           = cfg_menu;
    nproc.labels    = {'no multi-threading',sprintf('quarter of available threads (%d processes)',...
      floor(cores/4)),sprintf('half of available threads (%d processes)',floor(cores/2)),...
      sprintf('all threads-1 (%d processes)',cores-1),sprintf('all available threads (%d processes)',cores)};
    nproc.values    = {0 floor(cores/4) floor(cores/2) cores-1 cores};
    if cores<=4, nproc.labels{2} = [nproc.labels{2}(1:end-3) ')']; end
    if cores<4,  nproc.labels(2:3) = []; nproc.values(2:3) = []; end
    nproc.val       = {cores};
  end
  nproc.tag         = 'nproc';
  nproc.name        = 'Split job into separate processes';
  nproc.help        = { 
   ['In order to use multi-threading, the CAT12 segmentation job with multiple ' ...
    'subjects can be split into separate processes that run in the background.  ' ...
    'If you do not want to run processes in the background then set this value to 0. ']
   ['Keep in mind that each process needs a CPU core and about 2GB of RAM, ' ...
    'which should be considered to choose the appropriate number of processes.']
    ''
    };

  
  % debugging mode for developer
  debug               = cfg_menu;
  debug.tag           = 'debug';
  debug.name          = 'Debugging';
  debug.labels        = {'No','Yes'};
  debug.values        = {0,1};
  debug.val           = {1};
  debug.help          = {'Use low resolution setting for fast tests of the whole pipeline.';''};

  
  % main 
  catsimple         = cfg_exbranch;
  catsimple.tag     = 'cat_simple';
  catsimple.name    = 'CAT12 Simple Preprocessing'; 
  catsimple.val     = {data catversion tpm surface};
  if expert
    catsimple.val   = [catsimple.val {vROI sROI ignoreErrors}];
  end
  if expert > 1 % further mods do not work right now!
    catsimple.val   = [catsimple.val(1) {mods} catsimple.val(2:end)];
  end
  if cores > 1 % use multithreading only if availabe
    catsimple.val   = [catsimple.val {nproc}];
  end
  if expert>1 % add final debuging option
    catsimple.val   = [catsimple.val {debug}];
  end
  catsimple.prog    = @cat_simple;
  catsimple.vout    = @(job) vout_catsimple(job);
  catsimple.help    = { 
   ['This batch is a fully standardized cross-sectional CAT preprocessing that prepare your data ' ...
    'for voxel- (VBM), surface- (SBM) and region-based morphometry analysis (RBM). ' ...
    'It classifies the GM and WM brain tissue (segmentation) and maps them to the template space (spacial registration), ' ...
    'where it is smoothed with 8 mm (see MRI subdirectory). ' ... % ### UPDATE SMOOTHING ###
    'Next the central cortical surface is reconstructed and cortical measures were estimated (such as thickness, area, volume, and gyrification), ' ...
    'registered to the template surface (spherical registration) and smoothed 12 mm (see surf subdirectory). ' ... % ### UPDATE SMOOTHING ###
    'For region-of-interest (ROI) analysis the volumetric Neuromorphometrics and surface-based Desikan atlas are applied (see label subdirectory) ' ...
    'Moreover, statistical nuisance parameter such as total intracranial volume (TIV) and total surface area (TSA) were estimated (see report directory).']
    ''
    'Main reference:'
    '  CAT Toolbox paper' % ### UPDATE PAPER ###
    ''
   };

 
  % main long
  catsimplelong       = cfg_exbranch;
  catsimplelong.tag   = 'cat_simple_long';
  catsimplelong.name  = 'CAT12 Simple Longitudinal Preprocessing';
  catsimplelong.val   = {datalong catversion tpm surface}; 
  if expert
    catsimplelong.val = [catsimplelong.val {vROI sROI ignoreErrors}];
  end
  if expert > 1
    catsimplelong.val = [catsimplelong.val(1) {longmods} catsimplelong.val(2:end)];
  end
  if feature('numcores') > 1 % use multithreading only if availabe
    catsimplelong.val = [catsimplelong.val {nproc}];
  end
  if expert>1
    catsimplelong.val = [catsimplelong.val {debug}];
  end
  catsimplelong.prog  = @cat_simple;
  catsimplelong.vout  = @(job) vout_catsimple(job);
  catsimplelong.help  = strrep(catsimple.help,'cross-sectional','longidudinal'); 
  catsimplelong.help  = [
    {[catsimplelong.help{1} ' It requires the same subjects in the same order! ']} 
    catsimplelong.help(2:end)
    ];
   
return
function dep = vout_catsimple(job)
% _________________________________________________________________________
% Do we need the full output here! For instance to remove files?
% No, because batcher should not use the simple script. 
% Hence, only statistic relevant thing should be outputed!
% This already leads to 3
%
%   s#gmv*            ... for GM VBM
%
%   s#thickness.*     ... for cortical SBM  
%   s#gyrification.*  ... for folding SBM
%
%     s#volume.*
%     s#area.*
%     s#myelination.*
%   
%     s#intvol#
%
%   roi-data?
%
% _________________________________________________________________________

% ### UPDATE SMOOTHING ###
  vsmooth = [4 8 12];
  ssmooth = [6 12 24];
  
  % volume data
  for si = 1:numel(vsmooth)
    if ~exist('cdep','var'), dep = cfg_dep; else, dep(end+1) = cfg_dep; end %#ok<AGROW>
    dep(end).sname      = sprintf('%dmm smoothed modulated GMV',vsmooth(si));
    dep(end).src_output = substruct('.',sprintf('s%dmwp1',vsmooth(si)),'()',{':'});
    dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  
  % surface data
  if job.surface
    measures = {'thickness','curvature'}; %,'area','volume', 'myelination'
    if isfield(job,'mod')
      fname    =  fieldnames( [job.mod(:).name] );
      for fi=1:numel(fname)
        fname{fi} = genvarname(fname{fi});
        eval(fname{fi});
      end
      measures = [measures,fname];
    end
    for mi = 1:numel(measures)
      for si = 1:numel(ssmooth)
        dep(end+1)          = cfg_dep; %#ok<AGROW>
        dep(end).sname      = sprintf('%dmm smoothed %s',ssmooth(si),measures{mi});
        dep(end).src_output = substruct('.',sprintf('s%d%s',ssmooth(si),measures{mi}),'()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
      end
    end
  end
  
  % further images 
  % ... use job information
  %if ( isfield(job,'modality') && job.modality>1 ) && ...
  %   ( isfield(job,'timepoint') && isfield(job.timepoint,'modality') && job.timepoint.modality>1 )
  %end
  
  % ROI data
  dep(end+1)          = cfg_dep; 
  dep(end).sname      = '';
  dep(end).src_output = substruct('.','catroi','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});

  
return
