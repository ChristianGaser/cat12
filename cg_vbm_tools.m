function tools = cg_vbm_tools
% wrapper for calling VBM utilities
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

%_______________________________________________________________________

data         = cfg_files; 
data.tag     = 'data';
data.name    = 'Volumes';
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
data.help    = {
'Select raw data (e.g. T1 images) for processing. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not yet implemented for this method.'};

%------------------------------------------------------------------------

data_T2x         = cfg_files;
data_T2x.tag     = 'data_T2x';
data_T2x.name    = 'Volumes';
data_T2x.filter  = 'image';
data_T2x.ufilter = '^spmT.*\.[in][im][gi]$';
data_T2x.num     = [1 Inf];
data_T2x.help    = {'Select spmT-images to transform or convert.'};

sel        = cfg_menu;
sel.name   = 'Convert t value to';
sel.tag    = 'sel';
sel.labels = {'p','-log(p)','correlation coefficient cc','effect size d','apply thresholds without conversion'};
sel.values = {1,2,3,4,5};
sel.val    = {2};
sel.help   = {'Select conversion of t-value'};

thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'r';
thresh.num     = [1 1];
thresh.val     = {0.05};

thresh2         = cfg_entry;
thresh2.tag     = 'thresh2';
thresh2.name    = 'Threshold';
thresh2.help    = {''};
thresh2.strtype = 'r';
thresh2.num     = [1 1];
thresh2.val     = {0.001};

kthresh         = cfg_entry;
kthresh.tag     = 'kthresh';
kthresh.name    = 'Extent (voxels)';
kthresh.help    = {'Enter the extent threshold in voxels'};
kthresh.strtype = 'r';
kthresh.val     = {0};
kthresh.num     = [1 1];

noniso        = cfg_menu;
noniso.name   = 'Correct for non-isotropic smoothness';
noniso.tag    = 'noniso';
noniso.labels = {'yes','no'};
noniso.values = {1,0};
noniso.val    = {1};
noniso.help  = {'Correct for non-isotropic smoothness for cluster extent thresholds.'};

none         = cfg_const;
none.tag     = 'none';
none.name    = 'None';
none.val     = {1};
none.help    = {'No threshold'};

k         = cfg_branch;
k.tag     = 'k';
k.name    = 'k-value';
k.val     = {kthresh, noniso };
k.help    = {''};

fwe         = cfg_branch;
fwe.tag     = 'fwe';
fwe.name    = 'FWE';
fwe.val     = {thresh };
fwe.help    = {''};

fdr         = cfg_branch;
fdr.tag     = 'fdr';
fdr.name    = 'FDR';
fdr.val     = {thresh };
fdr.help    = {''};

fwe2         = cfg_branch;
fwe2.tag     = 'fwe2';
fwe2.name    = 'FWE';
fwe2.val     = {thresh, noniso };
fwe2.help    = {''};

uncorr         = cfg_branch;
uncorr.tag     = 'uncorr';
uncorr.name    = 'unocrrected';
uncorr.val     = {thresh2 };
uncorr.help    = {''};

uncorr2         = cfg_branch;
uncorr2.tag     = 'uncorr2';
uncorr2.name    = 'unocrrected';
uncorr2.val     = {thresh2, noniso };
uncorr2.help    = {''};

En         = cfg_branch;
En.tag     = 'En';
En.name    = 'Expected voxels per cluster';
En.val     = {noniso };
En.help    = {''};

inverse        = cfg_menu;
inverse.name   = 'Show also inverse effects (e.g. neg. values)';
inverse.tag    = 'inverse';
inverse.labels = {'yes','no'};
inverse.values = {1,0};
inverse.val    = {0};
inverse.help   = {'Show also inverse effects (e.g. neg. values). This is not valid if you convert to (log) p-values.'};

threshdesc        = cfg_choice;
threshdesc.name   = 'Threshold type peak-level';
threshdesc.tag    = 'threshdesc';
threshdesc.values = {none uncorr fdr fwe};
threshdesc.val    = {uncorr};
threshdesc.help   = {'Select method for voxel threshold'};

cluster        = cfg_choice;
cluster.name   = 'Cluster extent threshold';
cluster.tag    = 'cluster';
cluster.values = {none k En uncorr2 fwe2};
cluster.val    = {none};
cluster.help   = {'Select method for extent threshold'};

conversion         = cfg_branch;
conversion.tag     = 'conversion';
conversion.name    = 'Conversion';
conversion.val     = {sel threshdesc inverse cluster};
conversion.help    = {''};

T2x      = cfg_exbranch;
T2x.tag  = 'T2x';
T2x.name = 'Threshold and transform spmT-maps';
T2x.val  = {data_T2x,conversion};
T2x.prog = @cg_spmT2x;
T2x.help = {
          'This function transforms t-maps to P, -log(P), r or d-maps.'
          'The following formulas are used:'
          '--------------------------------'
          'correlation coefficient:'
          '          sign(t)'
          'r = ------------------'
          '           df'
          '    sqrt(------ + 1)'
          '          t*t'
          'effect-size'
          '           2r'
          'd = ----------------'
          '    sqrt(1-sqr(r))'
          'p-value'
          'p = 1-spm_Tcdf'
          'log p-value'
          '-log10(1-P) = -log(1-spm_Tcdf)'
          'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
          'Examples:'
          'p-value  -log10(1-P)'
          '0.1      1'
          '0.05     1.3'
          '0.01     2'
          '0.001    3'
          '0.0001   4'
          'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
          'Naming convention of the transformed files:'
          '   Type_Contrast_Pheight_Pextent_K_Neg'
          '   Type:      P    - p-value'
          '              logP - log p-value'
          '              R    - correlation coefficient'
          '              D    - effect size'
          '              T    - t-value'
          '   Contrast:  name used in the contrast manager with replaced none valid'
          '              strings'
          '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
          '              pFWE - p-value with FWE correction in %'
          '              pFDR - p-value with FDR correction in %'
          '   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")'
          '              pkFWE - extent p-value with FWE correction in %'
          '   K:         extent threshold in voxels'
          '   Neg:       image also shows thresholded inverse effects (e.g. neg. '
          '              values) '
}';

%------------------------------------------------------------------------

data_F2x         = cfg_files;
data_F2x.tag     = 'data_F2x';
data_F2x.name    = 'Volumes';
data_F2x.filter  = 'image';
data_F2x.ufilter = '^spmF.*\.[in][im][gi]$';
data_F2x.num     = [1 Inf];
data_F2x.help    = {'Select spmF-images to select.'};

sel        = cfg_menu;
sel.name   = 'Convert F value to';
sel.tag    = 'sel';
sel.labels = {'p','-log(p)','coefficient of determination R^2'};
sel.values = {1,2,3};
sel.val    = {2};
sel.help   = {'Select conversion of F-value'};

none         = cfg_const;
none.tag     = 'none';
none.name    = 'None';
none.val     = {1};
none.help    = {'No threshold'};

cluster        = cfg_choice;
cluster.name   = 'Cluster extent threshold';
cluster.tag    = 'cluster';
cluster.values = {none k};
cluster.val    = {none};
cluster.help  = {'Select method for extent threshold'};

conversion         = cfg_branch;
conversion.tag     = 'conversion';
conversion.name    = 'Conversion';
conversion.val     = {sel threshdesc cluster};
conversion.help    = {''};

F2x      = cfg_exbranch;
F2x.tag  = 'F2x';
F2x.name = 'Threshold and transform spmF-maps';
F2x.val  = {data_F2x,conversion};
F2x.prog = @cg_spmF2x;
F2x.help = {
          'This function transforms F-maps to P, -log(P), or R2-maps.'
          'The following formulas are used:'
          '--------------------------------'
          'coefficient of determination R2:'
          '          F*(n-1)'
          'R2 = ------------------'
          '        n-p + F*(n-1)'
          'p-value:'
          'p = 1-spm_Fcdf'
          'log p-value:'
          '-log10(1-P) = -log(1-spm_Fcdf)'
          'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
          'Examples:'
          'p-value  -log10(1-P)'
          '0.1      1'
          '0.05     1.3'
          '0.01     2'
          '0.001    3'
          '0.0001   4'
          'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
          'Naming convention of the transformed files:'
          '   Type_Contrast_Pheight_K'
          '   Type:      P    - p-value'
          '              logP - log p-value'
          '              R2   - coefficient of determination'
          '   Contrast:  name used in the contrast manager with replaced none valid'
          '              strings'
          '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
          '              pFWE - p-value with FWE correction in %'
          '              pFDR - p-value with FDR correction in %'
          '   K:         extent threshold in voxels'
}';

%------------------------------------------------------------------------

data.help = {
'Select all images. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images)'};

c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of nuisance values'};
c.strtype = 'r';
c.num     = [Inf 1];

nuisance       = cfg_branch;
nuisance.tag   = 'nuisance';
nuisance.name  = 'Nuisance';
nuisance.val   = {c};
nuisance.help  = {'Add a nuisance parameter to be removed from data'};

slice         = cfg_entry;
slice.tag     = 'slice';
slice.name    = 'Show slice (in mm)?';
slice.strtype = 'r';
slice.num     = [1 1];
slice.val     = {0};
slice.help    = {'Choose slice in mm.'};

gap         = cfg_entry;
gap.tag     = 'gap';
gap.name    = 'Separation';
gap.strtype = 'n';
gap.num     = [1 1];
gap.val     = {3};
gap.help    = {
  'To speed up calculations you can define that only every x voxel correlation is estimated. Smaller sampling distances gives slightly more accurate correlations, but will be much slower.'};

scale        = cfg_menu;
scale.tag    = 'scale';
scale.name   = 'Proportional scaling?';
scale.labels = {'no','yes'};
scale.values = {0 1};
scale.val    = {0};
scale.help   = {'This option should be only used if image intensity is not scaled (e.g. T1 images) or if images have to be scaled during statistical analysis (e.g. modulated images).'};

transform         = cfg_repeat;
transform.tag     = 'transform';
transform.name    = 'Nuisance';
transform.values  = {nuisance};
transform.num     = [0 Inf];
transform.help    = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to the calculation of the correlation.'};

check_cov      = cfg_exbranch;
check_cov.tag  = 'check_cov';
check_cov.name = 'Check homogeneity using sample correlation';
check_cov.val  = {data,gap,transform};
check_cov.prog = @cg_check_cov;
check_cov.help  = {
'If you have a reasonable sample size artefacts are easily overseen. In order to identify images with poor image quality or even artefacts you can use this function. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images). The idea of this tool is to check the correlation of all files across the sample.'
''
'The correlation is calculated between all images and the mean for each image is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if you have selected the images in the order of different sub-groups. Furthermore this is also useful for fMRI images which can be also used with this tool. The proportional scaling option should be only used if image intensity is not scaled (e.g. T1 images) or if images have to be scaled during statistical analysis (e.g. modulated images).'};

data_surf         = cfg_files;
data_surf.tag     = 'data_surf';
data_surf.name    = 'Surfaces';
data_surf.filter  = 'gifti';
data_surf.ufilter = '^[lr]h.central';
data_surf.num     = [1 Inf];
data_surf.help    = {'Select surfaces to extract values.'};

GI        = cfg_menu;
GI.name   = 'Gyrification index';
GI.tag    = 'GI';
GI.labels = {'none','yes'};
GI.values = {0,1};
GI.val    = {1};
GI.help   = {'Extract gyrification index (GI) based on absolute mean curvature. The method is described in Luders et al. NeuroImage, 29: 1224-1230, 2006.'};

FD        = cfg_menu;
FD.name   = 'Cortical complexity (fractal dimension)';
FD.tag    = 'FD';
FD.labels = {'none','yes'};
FD.values = {0,1};
FD.val    = {0};
FD.help   = {'Extract Cortical complexity (fractal dimension) based on absolute mean curvature. The method is described in Yotter et al. Neuroimage, 56(3): 961-973, 2011.'
             ''
             'Warning: Estimation of cortical complexity is very slow!'
''
};

SD        = cfg_menu;
SD.name   = 'Sulcus depth';
SD.tag    = 'SD';
SD.labels = {'none','yes'};
SD.values = {0,1};
SD.val    = {1};
SD.help   = {'Extract log10-transformed sulcus depth based on the euclidian distance between the central surface and its convex hull.'
             ''
             'Log-transformation is used to render the data more normally distributed.'
''
};

SA        = cfg_menu;
SA.name   = 'Surface area';
SA.tag    = 'SA';
SA.labels = {'none','yes'};
SA.values = {0,1};
SA.val    = {1};
SA.help   = {'Extract log10-transformed local surface area using re-parameterized tetrahedral surface. The method is described in Winkler et al. NeuroImage, 61: 1428–1443, 2012.',
             ''
             'Log-transformation is used to render the data more normally distributed.'
''
};

surfextract      = cfg_exbranch;
surfextract.tag  = 'surfextract';
surfextract.name = 'Extract surface parameters';
surfextract.val  = {data_surf,GI,FD,SD};
surfextract.prog = @vbm_surf_parameters;
surfextract.help = {'Using this option several surface parameters can be extracted that can be further analyzed.'};

data_surf         = cfg_files;
data_surf.tag     = 'data_surf';
data_surf.name    = 'Surfaces parameters';
data_surf.filter  = 'any';
data_surf.ufilter = '^[lr]h.[tgfl][hyro][irag][cias]';
data_surf.num     = [1 Inf];
data_surf.help    = {'Select surfaces parameter files for resampling to template space.'};

fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'Smoothing filter size in fwhm';
fwhm.strtype = 'r';
fwhm.num     = [1 1];
fwhm.val     = {15};
fwhm.help    = {
'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm.'};

surfresamp      = cfg_exbranch;
surfresamp.tag  = 'surfresamp';
surfresamp.name = 'Resample and smooth surface parameters';
surfresamp.val  = {data_surf,fwhm};
surfresamp.prog = @vbm_surf_resamp;
surfresamp.help = {
'In order to analyze surface parameters all data have to be resampled into template space and the resampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};

data_fs         = cfg_files;
data_fs.tag     = 'data_fs';
data_fs.name    = 'Freesurfer subject directories';
data_fs.filter  = 'dir';
data_fs.ufilter = '.*';
data_fs.num     = [1 Inf];
data_fs.help    = {'Select subject folders of freesurfer data to resample thickness data.'};

outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
outdir.help    = {'Select a directory where files are written.'};

surfresamp_fs      = cfg_exbranch;
surfresamp_fs.tag  = 'surfresamp_fs';
surfresamp_fs.name = 'Resample and smooth existing freesurfer thickness data';
surfresamp_fs.val  = {data_fs,fwhm,outdir};
surfresamp_fs.prog = @vbm_surf_resamp_freesurfer;
surfresamp_fs.help = {
'If you have existing freesurfer thickness data this function can be used to resample these data, smooth the resampled data, and convert freesurfer data to gifti format.'};

data_surf         = cfg_files;
data_surf.tag     = 'data_surf';
data_surf.name    = 'Surfaces';
data_surf.filter  = 'gifti';
data_surf.ufilter = 'resampled';
data_surf.num     = [1 Inf];
data_surf.help    = {'Select surfaces parameter files.'};

check_mesh_cov      = cfg_exbranch;
check_mesh_cov.tag  = 'check_mesh_cov';
check_mesh_cov.name = 'Check homogeneity using sample correlation';
check_mesh_cov.val  = {data_surf,transform};
check_mesh_cov.prog = @cg_check_cov;
check_mesh_cov.help = {
'If you have a reasonable sample size artefacts are easily overseen. In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces have to be resampled to the template space (e.g. normalized images). The idea of this tool is to check the correlation of all files across the sample.'
''
'The correlation is calculated between all images and the mean for each image is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this surface from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if you have selected the images in the order of different sub-groups.'
''
};

%------------------------------------------------------------------------

data.help = {
'Select images for quality assurance.'};

qa        = cfg_exbranch;
qa.tag    = 'qa';
qa.name   = 'VBM Quality Assurance';
qa.val    = {data};
qa.prog   = @vbm_tst_qa;
qa.vfiles = @vfiles_qa;
qa.help   = {'VBM Quality Assurance of T1 images. '};

%------------------------------------------------------------------------

data.help = {
'Select all images. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images)'};

showslice      = cfg_exbranch;
showslice.tag  = 'showslice';
showslice.name = 'Display one slice for all images';
showslice.val  = {data,scale,slice};
showslice.prog = @cg_showslice_all;
showslice.help = {'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality for a large number of files in a circumscribed region (slice).'};

%------------------------------------------------------------------------

data.help = {
'Select images for filtering'};

rician         = cfg_menu;
rician.tag     = 'rician';
rician.name    = 'Rician noise';
rician.labels  = {'Yes' 'No'};
rician.values  = {1 0};
rician.val     = {0};
rician.help    = {'MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected.'
''
'Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts (e.g. cortex can be affected by very high values in the scalp or in blood vessels.'
''
};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'sanlm_'};
prefix.help    = {'Specify the string to be prepended to the filenames of the smoothed image file(s). Default prefix is ''samlm_''.'};

sanlm        = cfg_exbranch;
sanlm.tag    = 'sanlm';
sanlm.name   = 'Spatially adaptive non-local means denoising filter';
sanlm.val    = {data prefix rician};
sanlm.prog   = @vbm_vol_sanlm;
sanlm.vfiles = @vfiles_sanlm;
sanlm.help   = {
'This function applies an spatial adaptive non-local means denoising filter to the data. This filter will remove noise while preserving edges. The filter strength is automatically estimated based on the standard deviation of the noise. '
'',
'This filter is internally used in the segmentation procedure anyway. Thus, it is not neccessary (and not recommended) to apply the filter before segmentation.'
''
};

%------------------------------------------------------------------------
calcvol_files         = cfg_files;
calcvol_files.tag     = 'data';
calcvol_files.name    = 'Volumes';
calcvol_files.filter  = '*';
calcvol_files.ufilter = 'seg.*\.txt$';
calcvol_files.num     = [1 Inf];
calcvol_files.help    = {
'Select all *_seg.txt files containing raw volumes, which were saved by VBM12 toolbox.'};

calcvol_name         = cfg_entry;
calcvol_name.tag     = 'calcvol_name';
calcvol_name.name    = 'Output file';
calcvol_name.strtype = 's';
calcvol_name.num     = [1 Inf];
calcvol_name.val     = {'raw_volumes.txt'};
calcvol_name.help    = {
'The output file is written to current working directory unless a valid full pathname is given'};

calcvol       = cfg_exbranch;
calcvol.tag   = 'calcvol';
calcvol.name  = 'Read raw volumes (GM/WM/CSF/Total)';
calcvol.val   = {calcvol_files,calcvol_name};
calcvol.prog  = @execute_calcvol;
calcvol.help  = {
'This function reads raw volumes for GM/WM/CSF/Total and saves values in a txt-file. These values can be read with the matlab command: vol = spm_load. The values for GM/WM/CSF/TOTAL are now saved in vol(:,1) vol(:,2) vol(:,3) and vol(:,4).'
''
'You can use these variables either as nuisance in an AnCova model or as user-specified globals with the "global calculation" option. Depending on your hypothesis and/or your data you can just use gray matter ("gm") or calculate the sum of gray/white matter with "gm+wm". The use of raw volumes as nuisance or globals is only recommended for modulated data. These data are corrected for size changes due to spatial  normalization and are thought to be in raw (un-normalized) space. In contrast, un-modulated data are yet corrected for differences in size due to spatial normalization to a reference brain and there is no need to correct for these differences again.'
''
};

%------------------------------------------------------------------------

field         = cfg_files;
field.tag     = 'field';
field.name    = 'Deformation Field';
field.filter  = 'image';
field.ufilter = '.*y_.*\.nii$';
field.num     = [1 Inf];
field.help    = {
'Deformations can be thought of as vector fields. These can be represented by three-volume images.'};

field1         = cfg_files;
field1.tag     = 'field1';
field1.name    = 'Deformation Field';
field1.filter  = 'image';
field1.ufilter = '.*y_.*\.nii$';
field1.num     = [1 1];
field1.help    = {
'Deformations can be thought of as vector fields. These can be represented by three-volume images.'};

images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are deformation fields, such that each flow field warps one image.'};
images1.filter = 'image';
images1.ufilter = '.*';
images1.num     = [1 Inf];

images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'The flow field deformations can be applied to multiple images. At this point, you are choosing how many images each flow field should be applied to.'};
images.values  = {images1 };
images.num     = [1 Inf];

interp        = cfg_menu;
interp.name   = 'Interpolation';
interp.tag    = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-spline',...
'3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {0,1,2,3,4,5,6,7};
interp.def    = @(val)cg_vbm_get_defaults('defs.interp',val{:});
interp.help   = {
'The method by which the images are sampled when being written in a different space.'
'    Nearest Neighbour:     - Fastest, but not normally recommended.'
'    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
'    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially       with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '
}';

modulate        = cfg_menu;
modulate.tag    = 'modulate';
modulate.name   = 'Modulate image (preserve volume)';
modulate.labels = {'yes','no'};
modulate.values = {1 0};
modulate.val    = {0};
modulate.help   = {
'"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). The SPM default is to adjust spatially normalised grey matter (or other tissue class) by using both terms and the resulting modulated images are preserved for the total amount of grey matter. Thus, modulated images reflect the grey matter volumes before spatial normalisation. However, the user is often interested in removing the confound of different brain sizes and there are many ways to apply this correction.'};

defs        = cfg_exbranch;
defs.tag    = 'defs';
defs.name   = 'Apply Deformations (Many images)';
defs.val    = {field1,images1,interp,modulate};
defs.prog   = @cg_vbm_defs;
defs.vfiles = @vfiles_defs;
defs.help   = {'This is a utility for applying a deformation field of one subject to many images.'};

defs2        = cfg_exbranch;
defs2.tag    = 'defs2';
defs2.name   = 'Apply Deformations (Many subjects)';
defs2.val    = {field,images,interp,modulate};
defs2.prog   = @cg_vbm_defs;
defs2.vfiles = @vfiles_defs2;
defs2.help   = {'This is a utility for applying deformation fields of many subjects to images.'};

%------------------------------------------------------------------------
realign = cg_cfg_realign;
bias    = cg_vbm_bias;
long    = cg_vbm_longitudinal_multi;
%------------------------------------------------------------------------

tools = cfg_choice;
tools.name   = 'Tools';
tools.tag    = 'tools';
tools.values = {showslice,check_cov,qa,calcvol,T2x,F2x,sanlm,bias,realign,long,defs,defs2,surfextract,surfresamp,surfresamp_fs,check_mesh_cov};

return

%_______________________________________________________________________

function vf = vfiles_defs(job)

PU = job.field1;
PI = job.images;

vf = cell(numel(PI),1);
for i=1:numel(PU),
    [pth,nam] = spm_fileparts(PU{i});
    for m=1:numel(PI),
        [pth1,nam,ext,num] = spm_fileparts(PI{m});
        if job.modulate,
            fname = fullfile(pth,['mw' nam ext]);
        else
            fname = fullfile(pth,['w' nam ext]);
        end;
        vf{m} = fname;
    end
end

return;
%_______________________________________________________________________

function vf = vfiles_defs2(job)

PU = job.field;
PI = job.images;

vf = cell(numel(PU),numel(PI));
for i=1:numel(PU),
    [pth,nam] = spm_fileparts(PU{i});
    for m=1:numel(PI),
        [pth1,nam,ext,num] = spm_fileparts(PI{m}{i});
        if job.modulate,
            fname = fullfile(pth,['mw' nam ext]);
        else
            fname = fullfile(pth,['w' nam ext]);
        end;
        vf{i,m} = fname;
    end
end

return;
%_______________________________________________________________________

function vf = vfiles_sanlm(job)
vf = {};

s  = strvcat(job.data);
for i=1:size(s,1),
    [pth,nam,ext,num] = spm_fileparts(s(i,:));
    vf = {vf{:}, fullfile(pth,[job.prefix,nam,ext,num])};
end;
return;
%_______________________________________________________________________
function vf = vfiles_qa(job)
vf = {};

s  = strvcat(job.data);
for i=1:size(s,1),
    [pth,nam,ext,num] = spm_fileparts(s(i,:));
    vf = {vf{:}, fullfile(pth,['p0',nam,'ext'])};
end;
return;
%_______________________________________________________________________

%------------------------------------------------------------------------
function execute_calcvol(p)
%
% calculate raw volumes all tissue classes
%

fprintf('%35s\t%5s\t%5s\t%5s\t%5s\n','Name','GM','WM','CSF','Total');
fid = fopen(p.calcvol_name,'w');

if fid < 0
	error('No write access: check file permissions or disk space.');
end

for i=1:length(p.data)
	tmp = load(deblank(p.data{i}));
    [pth,nam]     = spm_fileparts(p.data{i});
	if numel(tmp)==3
		fprintf(fid,'%3.2f\t%3.2f\t%3.2f\t%3.2f\n',tmp(1),tmp(2),...
			tmp(3),sum(tmp));
		fprintf('%35s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n',nam(1:end-4),tmp(1),tmp(2),...
			tmp(3),sum(tmp));
	else
		error(['Wrong format in ' p.data{i}]);
	end
end
if fclose(fid)==0
	fprintf('\nValues saved in %s.\n',p.calcvol_name);
end

return
%------------------------------------------------------------------------
