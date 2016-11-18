function extopts = cat_conf_extopts(expert,spm)
% Configuration file for extended CAT options
%
% Christian Gaser
% $Id$
%#ok<*AGROW>

if ~exist('expert','var')
  expert = 0; % switch to de/activate further GUI options
end
if ~exist('spm','var')
  spm = 0; % SPM segmentation input
end

%_______________________________________________________________________
% options for output
%-----------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size for normalized images';
vox.strtype = 'r';
vox.num     = [1 1];
vox.def     = @(val)cat_get_defaults('extopts.vox', val{:});
vox.help    = {
    'The (isotropic) voxel sizes of any spatially normalised written images. A non-finite value will be replaced by the average voxel size of the tissue probability maps used by the segmentation.'
''
};

%---------------------------------------------------------------------

pbtres         = cfg_entry;
pbtres.tag     = 'pbtres';
pbtres.name    = 'Voxel size for thickness estimation';
pbtres.strtype = 'r';
pbtres.num     = [1 1];
pbtres.def     = @(val)cat_get_defaults('extopts.pbtres', val{:});
pbtres.help    = {
    'Internal isotropic resolution for thickness estimation in mm.'
''
};

%---------------------------------------------------------------------

bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.strtype = 'r';
bb.num     = [2 3];
bb.def     = @(val)cat_get_defaults('extopts.bb', val{:});
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'
''
};

%------------------------------------------------------------------------
% special expert and developer options 
%------------------------------------------------------------------------

ignoreErrors        = cfg_menu;
ignoreErrors.tag    = 'ignoreErrors';
ignoreErrors.name   = 'Ignore errors';
ignoreErrors.labels = {'No','Yes'};
ignoreErrors.values = {0 1};
ignoreErrors.def    = @(val)cat_get_defaults('extopts.ignoreErrors', val{:});
ignoreErrors.help   = {
  'Catch preprocessing errors and go on with the next subject'
};

verb         = cfg_menu;
verb.tag     = 'verb';
verb.name    = 'Verbose processing level';
verb.labels  = {'none','default','details'};
verb.values  = {0 1 2};
verb.def     = @(val)cat_get_defaults('extopts.verb', val{:});
verb.help    = {
  'Verbose processing.'
};

debug         = cfg_menu;
debug.tag     = 'debug';
debug.name    = 'Debuging level';
debug.labels  = {'none','light','details'};
debug.values  = {0 1 2};
debug.def     = @(val)cat_get_defaults('extopts.debug', val{:});
debug.help    = {
  'Debuging processing.'
};


%---------------------------------------------------------------------
% Resolution
%---------------------------------------------------------------------

resnative        = cfg_branch;
resnative.tag    = 'native';
resnative.name   = 'Native resolution ';
resnative.help   = {
    'Preprocessing with native resolution.'
    'In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). '
    ''
    'Examples:'
    '  native resolution       internal resolution '
    '   0.95 0.95 1.05     >     0.95 0.95 1.05'
    '   0.45 0.45 1.70     >     0.45 0.45 1.50 (if voxel size for normalized images is 1.5mm)'
    '' 
  }; 

resbest        = cfg_entry;
resbest.tag    = 'best';
resbest.name   = 'Best native resolution';
resbest.def    = @(val)cat_get_defaults('extopts.resval', val{:});
resbest.num    = [1 2];
resbest.help   = {
    'Preprocessing with the best (minimal) voxel dimension of the native image.'
    'The first parameters defines the lowest spatial resolution for every dimension, while the second is used to avoid tiny interpolations for almost correct resolutions.'
    'In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). '
    ''
    'Examples:'
    '  Parameters    native resolution       internal resolution'
    '  [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.00 1.00'
    '  [1.00 0.10]    0.45 0.45 1.50     >     0.45 0.45 1.00'
    '  [0.75 0.10]    0.45 0.45 1.50     >     0.45 0.45 0.75'  
    '  [0.75 0.10]    0.45 0.45 0.80     >     0.45 0.45 0.80'  
    '  [0.00 0.10]    0.45 0.45 1.50     >     0.45 0.45 0.45'  
    ''
  }; 

resfixed        = cfg_entry;
resfixed.tag    = 'fixed';
resfixed.name   = 'Fixed resolution';
resfixed.def    = @(val)cat_get_defaults('extopts.resval', val{:});
resfixed.num    = [1 2];
resfixed.help   = {
    'This options prefers an isotropic voxel size that is controled by the first parameters.  '
    'The second parameter is used to avoid tiny interpolations for almost correct resolutions. ' 
    'In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). '
    ''
    'Examples: '
    '  Parameters     native resolution       internal resolution'
    '  [1.00 0.10]     0.45 0.45 1.70     >     1.00 1.00 1.00'
    '  [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.02]     0.95 1.05 1.25     >     1.00 1.00 1.00'
    '  [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [0.75 0.10]     0.75 0.95 1.25     >     0.75 0.75 0.75'
  }; 


restype        = cfg_choice;
restype.tag    = 'restypes';
restype.name   = 'Internal resampling for preprocessing';
switch cat_get_defaults('extopts.restype')
  case 'native', restype.val = {resnative};
  case 'best',   restype.val = {resbest};
  case 'fixed',  restype.val = {resfixed};
end
restype.values = {resnative resbest resfixed};
restype.help   = {
    'There are 3 major ways to control the internal spatial resolution ''native'', ''best'', and ''fixed''. In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). The minimum spatial resolution is 0.5mm. '
    ''
    'We commend to use ''best'' option to ensure optimal quality for preprocessing. ' 
}; 


%------------------------------------------------------------------------
% Cleanup
%------------------------------------------------------------------------

cleanupstr         = cfg_entry;
cleanupstr.tag     = 'cleanupstr';
cleanupstr.name    = 'Strength of Final Clean Up';
cleanupstr.strtype = 'r';
cleanupstr.num     = [1 1];
cleanupstr.def     = @(val)cat_get_defaults('extopts.cleanupstr', val{:});
cleanupstr.help    = {
  'Strength of tissue clean up after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. The default of 0.5 was successfully tested on a variety of scans. Use smaller values (>=0) for small changes and higher values (<=1) for stronger corrections. '
  ''
  'The strength changes multiple internal parameters: '
  ' 1) Size of the correction area'
  ' 2) Smoothing parameters to controll the opening processes to remove thin structures '
  ''
  'If parts of brain tissue were missing than decrease the strength.  If to many meninges are vissible than increase the strength. '
''
};

%------------------------------------------------------------------------
% Skull-stripping
%------------------------------------------------------------------------

gcutstr         = cfg_entry;
gcutstr.tag     = 'gcutstr';
gcutstr.name    = 'Strength of Skull-Stripping';
gcutstr.strtype = 'r';
gcutstr.num     = [1 1];
gcutstr.def     = @(val)cat_get_defaults('extopts.gcutstr', val{:});
gcutstr.help    = {
  'Strength of skull-stripping before AMAP segmentation, with 0 for a more liberal and wider brain masks and 1 for a more aggressive skull-stripping. The default of 0.5 was ' 
  'successfully tested on a variety of scans. '
  ''
  'The strength changes multiple internal parameters: '
  ' 1) Intensity thresholds to deal with blood-vessels and meninges '
  ' 2) Distance and growing parameters for the graph-cut/region-growing '
  ' 3) Closing parameters that fill the sulci'
  ' 4) Smoothing parameters that allow sharper or wider results '
  ''
  'If parts of the brain were missing in the brain mask than decrease the strength.  If the brain mask of your images contains parts of the head, than increase the strength. '
''
};

%------------------------------------------------------------------------
% Noise correction
%------------------------------------------------------------------------
sanlm        = cfg_menu;
sanlm.tag    = 'sanlm';
sanlm.name   = 'Use SANLM de-noising filter';
sanlm.labels = {'No denoising','SANLM denoising','ISARNLM denoising'};
sanlm.values = {0 1 2};
sanlm.def    = @(val)cat_get_defaults('extopts.sanlm', val{:});
sanlm.help   = {
    'This function applies an spatial adaptive non local means (SANLM) or the iterative spatial resolution adaptive non local means (ISARNLM) denoising filter to the data. Using the ISARNLM filter is only required for high resolution data with parallel image artifacts or strong noise. Both filter will remove noise while preserving edges. Further modification of the strength of the noise correction is possible by the NCstr parameter. '
    'The following options are available: '
    '  * No noise correction '
    '  * SANLM '
    '  * ISARNLM ' 
};

NCstr         = cfg_entry;
NCstr.tag     = 'NCstr';
NCstr.name    = 'Strength of Noise Corrections';
NCstr.strtype = 'r';
NCstr.num     = [1 1];
NCstr.def     = @(val)cat_get_defaults('extopts.NCstr', val{:});
NCstr.help    = {
  'Strength of the SANLM noise correction. The default "inf" uses an adaptive noise correction and was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction!'
''
};

%------------------------------------------------------------------------
% Blood Vessel Correction
%------------------------------------------------------------------------

BVCstr        = cfg_entry;
BVCstr.tag    = 'BVCstr';
BVCstr.name    = 'Strength of Blood Vessel Corrections';
BVCstr.strtype = 'r';
BVCstr.num     = [1 1];
BVCstr.def     = @(val)cat_get_defaults('extopts.BVCstr', val{:});
BVCstr.help    = {
  'Strength of the Blood Vessel Correction. The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction!'
  ''
};

%------------------------------------------------------------------------
% Local Adapative Segmentation
%------------------------------------------------------------------------

LAS        = cfg_menu;
LAS.tag    = 'LAS';
LAS.name   = 'Local Adaptive Segmentation';
LAS.labels = {'No','Yes'};
LAS.values = {0 1};
LAS.def    = @(val)cat_get_defaults('extopts.LAS', val{:});
LAS.help   = {
  'Correction of local intensity changes with medium/low spatial frequencies. This will affect mostly subcortical GM areas. This function will also utilize the inhomogeneity correction. '
  ''
  'See also ...'
  ''
};

LASstr         = cfg_entry;
LASstr.tag     = 'LASstr';
LASstr.name    = 'Strength of Local Adaptive Segmentation';
LASstr.strtype = 'r';
LASstr.num     = [1 1];
LASstr.def     = @(val)cat_get_defaults('extopts.LASstr', val{:});
LASstr.help    = {
  'Strength of the modification by the Local Adaptive Segmentation (LAS). The default 0.5 was successfully tested on a large variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger corrections. The value 0 will deactive LAS.'
  ''
};

%------------------------------------------------------------------------
% WM Hyperintensities:
%------------------------------------------------------------------------

wmhc        = cfg_menu;
wmhc.tag    = 'WMHC';
wmhc.name   = 'WM Hyperintensity Correction (WMHC)';
wmhc.labels = { ...
  'no correction' ...
  'WMHC - correction of WM only for spatial normalization' ... 
  'WMHC - correction of WM segmentations' ...
  'WMHC - correction as separate class' ...
};
wmhc.values = {0 1 2 3};
wmhc.def    = @(val)cat_get_defaults('extopts.WMHC', val{:});
wmhc.help   = {
  'In aging or diseases WM intensity be strongly reduces in T1 or increased in T2/PD images. These so called WM hyperintensies (WMHs) can lead to preprocessing errors. Large GM areas next to the ventricle can cause normalization problems. Therefore, a temporary correction for the normalization is meaningfull, if WMHs were expected. As far as these changes are an important marker, CAT allows different ways to handel WMHs. '
  ''
  ' 0) No Correction. '
  '     - Take care of large WMHs that might cause normalization problems. '
  '     - Consider that GM in unexpected regions represent WMCs.  '
  ' 1) Temporary correction for spatial normalization. '
  '     - Consider that GM in unexpected regions represent WMCs.  '
  ' 2) Correction of WM segmentations (like SPM). ' 
  ' 3) Correction as separate class. '
  ''
  'See also ...'
''
};

WMHCstr         = cfg_entry;
WMHCstr.tag     = 'WMHCstr';
WMHCstr.name    = 'Strength of WMH Correction';
WMHCstr.strtype = 'r';
WMHCstr.num     = [1 1];
WMHCstr.def     = @(val)cat_get_defaults('extopts.WMHCstr', val{:});
WMHCstr.help    = {
  'Strength of the modification of the WM Hyperintensity Correction (WMHC). The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger corrections. The value 0 will deactive WMHC.'
''
};

%------------------------------------------------------------------------

print        = cfg_menu;
print.tag    = 'print';
print.name   = 'Display and print results';
print.labels = {'No','Yes'};
print.values = {0 1};
print.def    = @(val)cat_get_defaults('extopts.print', val{:});
print.help   = {
'The result of the segmentation can be displayed and printed to a ps-file. This is often helpful to check whether registration and segmentation were successful. Furthermore, some additional information about the used versions and parameters are printed.'
''
};

%------------------------------------------------------------------------

darteltpm         = cfg_files;
darteltpm.tag     = 'darteltpm';
darteltpm.name    = 'Spatial normalization Template';
darteltpm.filter  = 'image';
darteltpm.ufilter = 'Template_1'; % the string Template_1 is SPM default
darteltpm.def     = @(val)cat_get_defaults('extopts.darteltpm', val{:});
darteltpm.num     = [1 1];
darteltpm.help    = {
  'Selected Dartel/Shooting template must be in multi-volume nifti format and should contain GM and WM segmentations. The template of the first iteration (indicated by "Template_1" for Dartel and "Template_0" for Shooting) must be selected, but the templates of all iterations must be existing.'
  ''
  'Please note that the use of an own DARTEL template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ.'
  ''
};

%------------------------------------------------------------------------

% for other species
cat12atlas         = cfg_files;
cat12atlas.tag     = 'cat12atlas';
cat12atlas.name    = 'CAT12 ROI atlas';
cat12atlas.filter  = 'image';
cat12atlas.ufilter = '_1';
cat12atlas.def     = @(val)cat_get_defaults('extopts.cat12atlas', val{:});
cat12atlas.num     = [1 1];
cat12atlas.help    = {
  'CAT12 atlas file to handle major regions.'
};

%------------------------------------------------------------------------

% for other species
brainmask         = cfg_files;
brainmask.tag     = 'brainmask';
brainmask.name    = 'Brainmask';
brainmask.filter  = 'image';
brainmask.ufilter = '_1';
brainmask.def     = @(val)cat_get_defaults('extopts.brainmask', val{:});
brainmask.num     = [1 1];
brainmask.help    = {
  'Initial brainmask.'
};

%------------------------------------------------------------------------

% for other species
T1         = cfg_files;
T1.tag     = 'T1';
T1.name    = 'T1';
T1.filter  = 'image';
T1.ufilter = '_1';
T1.def     = @(val)cat_get_defaults('extopts.T1', val{:});
T1.num     = [1 1];
T1.help    = {
  'Affine registration template.'
};

%------------------------------------------------------------------------
appfull          = cfg_entry;
appfull.tag      = 'APP';
appfull.name     = 'Affine Preprocessing (APP) code';
appfull.strtype  = 'n';
appfull.num      = [1 1];
appfull.def      = @(val)cat_get_defaults('extopts.APP', val{:});
appfull.help     = { ...
  'Affine alignment and SPM preprocessing can fail in untypical subjects with deviating anatomy (other species/neonates) or in images with strong signal inhomogeneities or untypical intensities like in synthetic images). An initial bias correction, a head/brain masking/extraction can help to reduce such problems. ' ...
  '' ... 
  'Before changing this parameter please check image orientation.' ...
  'Strong and heavy APP uses only fine affine registration that typically reqired a good intial orientation.' ...
  '' ...
  'Enter a number with 3 Digits with the first digit for biascorr, digit 2 for masking, and digit 3 for affreg with:' ...
  '' ...
  ' biascorr: ' ...
  '  0 = none '...
  '  1 = initial, only for affine registration (default) '...
  '  2 = initial, full use '...
  '  3 = fine, only for affine registration '...
  '  4 = fine, full use '...
  '' ...
  ' masking: ' ...
  '  0 = none '...
  '  1 = head, msk only (default) '...
  '  2 = head, apply mal use '...
  '  3 = brain, msk only '...
  '  4 = brain, apply mal use '...
  '' ...
  ' affreg:' ...
  '  0 = none '...
  '  1 = use affine preregistration (default) '...
  '' ...
  '' ...
  'Standard APP values for default/expert user menu:' ...
  '  none:      0 = 000 - old default' ...
  '  light:     1 = 111' ...
  '  medium:    2 = 211' ...
  '  strong:    3 = 210' ...
  '  heavy:     4 = 430' ...
  '  nonhuman:  5 = 440' ...
  '' ...
};

applight        = cfg_menu;
applight.tag    = 'APP';
applight.name   = 'Affine Preprocessing (APP)';
applight.labels = { ...
  'none'  ...
  'light' ... 
  'heavy' ... 
};
applight.values = {0 1 4};
applight.def    = @(val)cat_get_defaults('extopts.APP', val{:});
applight.help   = { ...
  'Affine alignment and SPM preprocessing can fail in untypical subjects with deviating anatomy (other species/neonates) or in images with strong signal inhomogeneities or untypical intensities like in synthetic images). An initial bias correction, a head/brain masking/extraction can help to reduce such problems. ' ...
  '' ... 
  'Before changing this paramter please check image orientation.' ...
  'Strong and heavy APP uses only fine affine registration that typically reqired a good intial orientation.' ...
  '' ...
};

app        = cfg_menu;
app.tag    = 'APP';
app.name   = 'Affine Preprocessing (APP)';
app.labels = { ...
  'none' ... 
  'light' ... 
  'medium' ... 
  'strong' ... 
  'heavy' ... 
};
app.values = {0 1 2 3 4};
app.def    = @(val)cat_get_defaults('extopts.APP', val{:});
app.help   = { ...
  'Affine alignment and SPM preprocessing can fail in untypical subjects with deviating anatomy (other species/neonates) or in images with strong signal inhomogeneities or untypical intensities like in synthetic images). An initial bias correction, a head/brain masking/extraction can help to reduce such problems. ' ...
  '' ... 
  'Before changing this paramter please check image orientation.' ...
  'Strong and heavy APP uses only fine affine registration that typically reqired a good intial orientation.' ...
  '' ...
  'APP is still in development and these options may change. ' ...
  '' ...
  '  APP: biascorr, masking, affreg' ...
  '  -------------------------------' ...
  '  none: none, none, yes' ...
  '  light: just for affreg, head msk, yes' ...
  '  medium: permanent (init), head msk, yes' ...
  '  strong: permanent (fine), brain msk, no'  ...
  '  heavy: permanent (fine), brain msk, no' ...
  '' ...
};

%------------------------------------------------------------------------

lazy         = cfg_menu;
lazy.tag     = 'lazy';
lazy.name    = 'Lazy processing';
lazy.labels  = {'yes','no'};
lazy.values  = {1,0};
lazy.val     = {0};
lazy.help    = {
  'Do not process data if the result exist. '
};

extopts       = cfg_branch;
extopts.tag   = 'extopts';
extopts.name  = 'Extended options for CAT12 segmentation';
if ~spm
  if expert>=2 % experimental expert options
    extopts.val   = {lazy,appfull,sanlm,NCstr,LASstr,gcutstr,cleanupstr,BVCstr,WMHCstr,wmhc,...
                     darteltpm,cat12atlas,brainmask,T1,...
                     restype,vox,pbtres,ignoreErrors,debug,verb}; 
  elseif expert==1 % working expert options
    extopts.val   = {app,sanlm,NCstr,LASstr,gcutstr,cleanupstr,WMHCstr,wmhc,darteltpm,restype,vox,ignoreErrors}; 
  else
    extopts.val   = {applight,LASstr,gcutstr,cleanupstr,darteltpm,vox}; 
  end
else
  if expert>=2 % experimental expert options
    extopts.val   = {lazy,darteltpm,cat12atlas,brainmask,T1,vox,pbtres,ignoreErrors,debug,verb}; 
  elseif expert==1 % working expert options
    extopts.val   = {darteltpm,vox,ignoreErrors}; 
  else
    extopts.val   = {darteltpm,vox}; 
  end 
end
extopts.help  = {'Using the extended options you can adjust special parameters or the strength of different corrections ("0" means no correction and "0.5" is the default value that works best for a large variety of data).'};
