function extopts = cg_vbm_extopts(expert)
% Configuration file for extended VBM options
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
vox.def     = @(val)cg_vbm_get_defaults('extopts.vox', val{:});
vox.help    = {
    'The (isotropic) voxel sizes of any spatially normalised written images. A non-finite value will be replaced by the average voxel size of the tissue probability maps used by the segmentation.'
''
};


%---------------------------------------------------------------------

bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.strtype = 'r';
bb.num     = [2 3];
bb.def     = @(val)cg_vbm_get_defaults('extopts.bb', val{:});
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'
''
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
resbest.def    = @(val)cg_vbm_get_defaults('extopts.resval', val{:});
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
resfixed.def    = @(val)cg_vbm_get_defaults('extopts.resval', val{:});
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
restype.tag    = 'restype';
restype.name   = 'Internal resampling for preprocessing';
switch cg_vbm_get_defaults('extopts.restype')
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
cleanupstr.def     = @(val)cg_vbm_get_defaults('extopts.cleanupstr', val{:});
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
gcutstr.def     = @(val)cg_vbm_get_defaults('extopts.gcutstr', val{:});
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
sanlm.labels = {'No denoising','SANLM denoising','SANLM denoising (multi-threaded)',...
    'SANLM denoising + ORNLM','SANLM denoising (multi-threaded) + ORNLM','ORNLM denoising'
};
sanlm.values = {0 1 2 3 4 5};
sanlm.def    = @(val)cg_vbm_get_defaults('extopts.sanlm', val{:});
sanlm.help   = {
    'This function applies an spatial adaptive non local means (SANLM) denoising filter to the data. This filter will remove noise while preserving edges. The smoothing filter size is automatically estimated based on the local variance in the image. Optionally, a second ORNLM filter is  used after inhomogeneity correction and intensity scaling in order to remove increased noise after local scaling. Using an option with ORNLM filter allows further modification of the strength of the noise correction. '
    'The following options are available: '
    '  0) No noise correction '
    '  1) SANLM '
    '  2) SANLM (multi-threaded) '
    '  3) SANLM + ORNLM 4) SANLM (multi-threaded) + ORNLM'
    '  5) Only ORNLM (with temporarily SANLM filter, but only one final noise correction of the original data)' 
};

NCstr         = cfg_entry;
NCstr.tag     = 'NCstr';
NCstr.name    = 'Strength of Noise Corrections';
NCstr.strtype = 'r';
NCstr.num     = [1 1];
NCstr.def     = @(val)cg_vbm_get_defaults('extopts.NCstr', val{:});
NCstr.help    = {
  'Strengh of the SANLM, ORNLM and MRF noise correction. The default 0.5 was successfully tested on a variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction!'
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
LAS.def    = @(val)cg_vbm_get_defaults('extopts.LAS', val{:});
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
LASstr.def     = @(val)cg_vbm_get_defaults('extopts.LASstr', val{:});
LASstr.help    = {
  'Strengh of the modification by the Local Adaptive Segmentation (LAS). The default 0.5 was successfully tested on a large variety of scans. Use smaller values (>0) for small changes and higher values (<=1) for stronger corrections. The value 0 will deactive LAS.'
  ''
};

%------------------------------------------------------------------------
% WM Hyperintensities:
%------------------------------------------------------------------------

WMHC        = cfg_menu;
WMHC.tag    = 'WMHC';
WMHC.name   = 'WM Hyperintensity Correction (WMHC)';
WMHC.labels = {'no correction','WMHC - correction to WM like SPM','WMHC - correction to a separate class'};
WMHC.values = {0 1 2 3};
WMHC.def    = @(val)cg_vbm_get_defaults('extopts.WMHC', val{:});
WMHC.help   = {
  'In aging or diseases WM intensity be strongly reduces in T1 or increased in T2/PD images. These so called WM hyperintensies (WMHs) can lead to preprocessing errors. Large GM areas next to the ventricle can cause normalization problems. Therefore, a temporary correction for the normalization is meaningfull, if WMHs were expected. As far as these changes are an important marker, VBM allow different ways to handel WMHs. '
  ''
  ' 0) No Correction (like VBM8). '
  '     - Take care of large WMHs with normalization problems. '
  '     - Consider that GM in unexpected regions represent WMCs.  '
  ' 1) Temporary Correction for normalization. '
  '     - Consider that GM in unexpected regions represent WMCs.  '
  ' 2) Correction to WM (like SPM). ' 
  ' 3) Correction to separate class. '
  ''
  'See also ...'
''
};

WMHCstr         = cfg_entry;
WMHCstr.tag     = 'WMHCstr';
WMHCstr.name    = 'Strength of WMH Correction';
WMHCstr.strtype = 'r';
WMHCstr.num     = [1 1];
WMHCstr.def     = @(val)cg_vbm_get_defaults('extopts.WMHCstr', val{:});
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
print.def    = @(val)cg_vbm_get_defaults('extopts.print', val{:});
print.help   = {
'The result of the segmentation can be displayed and printed to a ps-file. This is often helpful to check whether registration and segmentation were successful. Furthermore, some additional information about the used versions and parameters are printed.'
''
};

%------------------------------------------------------------------------

darteltpm         = cfg_files;
darteltpm.tag     = 'darteltpm';
darteltpm.name    = 'Dartel Template';
darteltpm.filter  = 'image';
darteltpm.ufilter = '_1_';
darteltpm.def     = @(val)cg_vbm_get_defaults('extopts.darteltpm', val{:});
darteltpm.num     = [1 1];
darteltpm.help    = {
   'Selected Dartel template must be in multi-volume nifti format and should contain GM and WM segmentations. Here the template of the first iteration (indicated by _1_) should be selected, but the templates of all 6 iterations should be present.'
''
};

%------------------------------------------------------------------------

extopts       = cfg_branch;
extopts.tag   = 'extopts';
extopts.name  = 'Extended options for VBM12 segmentation';
if expert
  extopts.val   = {sanlm,NCstr,LASstr,gcutstr,cleanupstr,darteltpm,restype,vox,print}; 
else
  extopts.val   = {NCstr,LASstr,gcutstr,cleanupstr,darteltpm,vox,print}; 
end
extopts.help  = {'Using the extended options you can adjust the strength of different corrections ("0" means no correction and "0.5" is the default value that works best for a large variety of data).'};
