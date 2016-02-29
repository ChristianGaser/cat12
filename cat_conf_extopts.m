function extopts = cat_conf_extopts(expert)
% Configuration file for extended CAT options
%
% Christian Gaser
% $Id$
%#ok<*AGROW>

if ~exist('expert','var')
  expert = 0; % switch to de/activate further GUI options
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
  'Strengh of tissue clean up after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. The default of 0.5 was successfully tested on a variety of scans. Use smaller values (>=0) for small changes and higher values (<=1) for stronger corrections. '
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
  'Strengh of skull-stripping before AMAP segmentation, with 0 for a more liberal and wider brain masks and 1 for a more aggressive skull-stripping. The default of 0.5 was ' 
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
  'Strengh of the SANLM, ORNLM and MRF noise correction. The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction!'
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
  'Strengh of the Blood Vessel Correction. The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction!'
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
  'Strengh of the modification by the Local Adaptive Segmentation (LAS). The default 0.5 was successfully tested on a large variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger corrections. The value 0 will deactive LAS.'
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
  'Strengh of the modification of the WM Hyperintensity Correction (WMHC). The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger corrections. The value 0 will deactive WMHC.'
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
darteltpm.ufilter = '_1';
darteltpm.def     = @(val)cat_get_defaults('extopts.darteltpm', val{:});
darteltpm.num     = [1 1];
darteltpm.help    = {
   'Selected Dartel/Shooting template must be in multi-volume nifti format and should contain GM and WM segmentations. The template of the first iteration (indicated by "_1") must be selected, but the templates of all iterations must be present.'
''
};

%------------------------------------------------------------------------

app        = cfg_menu;
app.tag    = 'APP';
app.name   = 'Affine Preprocessing (APP)';
app.labels = { ...
  'none (only SPM preprocessing)' ... the old default 
  'APP1 (fast bias correction only for affine registration)' ... just to test it ... the second BC should be much more exact and stable
  'APP2 (fast bias correction)' ... just to test it ... the second BC should be much more exact and stable
  'APP3 (fine bias correction)' ... just to test it ... the hard skull-stripping can lead to problems if brain tissue is missing, but actual it seams to work very well
  'APP4 (fine bias correction, skull-stripping, no initial registration)' ... I expect that this is only important for animals and can maybe controlled by the species parameter
};
app.values = {0 1 2 3 4};
app.def    = @(val)cat_get_defaults('extopts.APP', val{:});
app.help   = {
  'Affine alignment and SPM preprocessing can fail in untypical subjects with deviating anatomy (other species/neonates) or in images  with strong signal inhomogeneities. ' ...
  'An initial rough bias correction and the extraction of the brain can reduce problems (APP = Affine PreProcessing). ' ...
  'Because the first affine registration that is required for the brain mask can also fail in non-humans the option (APP3) without initial registation is available that requires excact AD-PC alignment by the user (i.e. by SPM Display). ' ... 
  ''
};

%------------------------------------------------------------------------

extopts       = cfg_branch;
extopts.tag   = 'extopts';
extopts.name  = 'Extended options for CAT12 segmentation';
if expert>2 % experimental expert options
  extopts.val   = {app,sanlm,NCstr,LASstr,gcutstr,cleanupstr,BVCstr,WMHCstr,wmhc,darteltpm,restype,vox,pbtres,ignoreErrors,debug,verb,subfolders}; 
elseif expert==1 % working expert options
  extopts.val   = {sanlm,NCstr,LASstr,gcutstr,cleanupstr,WMHCstr,wmhc,darteltpm,restype,vox,ignoreErrors}; 
else
  extopts.val   = {NCstr,LASstr,gcutstr,cleanupstr,darteltpm,vox}; 
end
extopts.help  = {'Using the extended options you can adjust special parameter or the strength of different corrections ("0" means no correction and "0.5" is the default value that works best for a large variety of data).'};
