%-----------------------------------------------------------------------
% Job saved on 05-Aug-2021 10:59:09 by cfg_util (rev $Rev$)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% Batch for preprocess high-resolution ex-vivo PD data with SPM and CAT 
% to create volume and surfaces data for structural analyses. The batch 
% presented only a first attempt and to process the data presented in 
% (Edlow et al., 2019). 
%
%   Edlow, B.L., Mareyam, A., Horn, A., Polimeni, J.R., Witzel, T., 
%   Tisdall, M.D., Augustinack, J.C., Stockmann, J.P., Diamond, B.R., 
%   Stevens, A., Tirrell, L.S., Folkerth, R.D., Wald, L.L., Fischl, B.R.,
%   van der Kouwe, A., 2019. 
%   7 Tesla MRI of the ex vivo human brain at 100 micron resolution.
%   Sci. Data 6, 1-10. <a>doi:10.1038/s41597-019-0254-8
%
%   Data: 
%     https://kottke.org/19/07/the-highest-resolution-mri-scan-of-a-human-brain	
% 
%
% Use: 
% * Start SPM/CAT with: 
%   spm 'fmri'; cat12('developer')
% * Specify your own directories below 
% * Update some variables such as "resolution", "spm_res_factor" and
%   "preview_surfaces" that are defined for a fast test run
%   (currently it is only reducing resolution with preview surface estimation)
% * Open the SPM batch mode via SPM GUI
% * Open this batch in the SPM batch mode 
% * Start the batch (if it is not ready or create errors then start 
%   checking the input files)
%
%
% Known limitations and problems (20210808):
% * I am currently limited to only 2 test datasets. 
%
% * There were multiple problems due to the segmentation, even with reduced
%   classes.  I tried two models: (1) a 3-class model with GM, WM, and  
%   CSF/BG (background) and (2) a 4-class model with GM, WM, CSF, and BG.
%   Both are not really perfect but the simpler model fit a bit better.
%
% * The 4-class model can also be used to skull-stripping the image. This 
%   can help in some way (e.g. bias correction) but there are some strange 
%   side effect and in the stripped background the TPM probability appeared 
%   (i.e. there was soft GM ribbon around the old brainmask). 
%
% * I tried different number of Gaussians per tissue but it did not work in 
%   the way I see by eye in the images. Besides the SPM/CAT defaults [1 1 2]
%   for [GM,WM,CSF/BG] I got good results for [2,1,3].
%
% * Intensity normalization/limitation helped overall to remove some crazy 
%   outliers and I added it multiple times.
%
% * Moreover, I observed that increased segmentation sampling resolution 
%   (samp<3) result in larger problems especially in images with higher
%   resolution (the background was often miss-classified as GM or WM). 
%   I also tried high numbers of Gaussians per tissue class without success. 
%
% In conclusion, I think the samp parameter of the unified segmentation has
% the strongest effect and I was not able to find a way to use higher 
% sampling (samp<3) here that generally improve segmentation quality. 
% In addition, the changed gaussians, as well as the intensity limitation 
% improved the result
%
%
% Methods: 
% (1) Trimming
% (2) Denoising
% (3) Downsampling
% (4) Data range limitation 
% (5) Unified segmentation 
% (6) CAT12 preprocessing for SPM segmentation 
% (7) CAT12 preprocessing
%
% Further directions/optimizations: 
% * The high-res PD data may allow an enhanced cortical model with two or 
%   three GM MRI layer with very bright layer (lamina 1-2?), a medium dark 
%   layer (lamina 3,5 and 6?) and a dark layer (lamina 4). On the other 
%   side, a higher number of Gaussians is quite similar and also not really 
%   working, i.e., high affords but low benefit . 
% * Moreover, some WM regions could be defined as separate class, as well 
%   as some damaged areas (lesions?). However, this would require to create 
%   a standard GM,WM,CSF/background input model. 
%
% * Moreover, the surface reconstruction has to use the segmentation and 
%   not the original data as far as this pipeline (cat_surf_create[2]) is
%   optimized for T1 data. 
%   > a simulated T1 dataset (inverted image) maybe works better
%
%-----------------------------------------------------------------------
% Robert Dahnke 2020/09 (robert.dahnke@uni-jena.de)
%-----------------------------------------------------------------------
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
%-----------------------------------------------------------------------
% $Id$


% Todo: 
% - extend to strongly trimmed images or better set up CAT processing
% - add non-developer version

% set random generators for SPM (not really working)
if exist('rng','file') == 2, rng('default'); rng(0); else, rand('state',0); randn('state',0); end
warning('OFF','MATLAB:RandStream:ActivatingLegacyGenerators');   

% clear old stuff
clear matlabbatch; 


% (0) Define the path to your files and directories here!
data  = {
  {'/Volumes/WD4TBE/MRData/202106 UHR PD/Acquired_FA25_downsampled_200um.nii,1'}
  }; 
TPM                = fullfile(spm('dir'),'TPM','TPM.nii');
template_dir       = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym'); 
resolution         = 0.5;  %0.8;  % define here the resolution you want (0 - original resolution)
                           % for fast test use 0.8 mm 
samp               = 3.0;  % 3 mm is ok,
tol                = 1e-32;  % super accurate to correct even strong inhomogeneities but it will take a while
%spm_res_factor     = 1.5/resolution;    % use x-times lower sampling resolution (depending on the resolution variable) for SPM for faster tests 
                           % for a resolution of 0.8 and spm_res_factor of 2 SPM will use 1.6 mm (default is 3 mm for human data and)                            
nproc              = 0;    % number of parallel processes of the CAT preprocessing
preview_surfaces   = 5;    % (0 - none, 1 - yes, 5 - fast preview surfaces for visual checks but not for statistical analyses
                           % use 5 for fast reviews to optimize the pipeline 
reduce_mesh        = 5;    % more accurate but maybe with fatal Matlab crash (just try to rerun) otherwise use 1 (optimal by volume reduction)
%ngaus              = [1 1 2 3];  % [GM WM CSF[/BG | BG]] ... [2 1 3 4] seems also to work
ngaus              = [2 1 3 4];  % [GM WM CSF[/BG | BG]] ... [2 1 3 4] seems also to work
BGmodel            = 2;    % 0 - 4-class, 1 - 4-class-skull-stripped, 2 - 3-class(CSFBG)
BC                 = 0; 
lazy               = 1;    % do not reprocess data if the input was not modified (default = 0)


% (0) Display parameter ...
%{
fprintf('Ultra-high resolution preprocessing batch: ')
for fi = 1:numel(data)
  [pp,ff,ee] = spm_fileparts(data{fi});
  file = fullfile(pp,[ff ee]); 
  if ~exist(file,'file')
    cat_io_printf('err',sprintf('ERROR: Image %d/%d - "%s" does not exist!\n',fi,numel(data),)
  end
end
%}

% (1) Trimm the image and remove empty space around the brain to save 
%     memory and increase processing time. There are many parameters but 
%     no changes are required. 
% ######### 
% TODO: 
%  * It would be better to use mm rather than voxel for the limitation
%  * It could be helpful to guaranty some boundary (e.g. for resampling)
% ######### 
mi=1; 
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.image_selector.subjectimages  = data;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.prefix                        = 'trimmed_';
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.mask                          = 1;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.suffix                        = '';
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.intlim1                       = 90;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.pth                           = 0.4;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.avg                           = 0;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.open                          = 2;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.addvox                        = 50; % default 2; 0.2 mm > 1 cm
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.spm_type                      = 0;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.intlim                        = 100;
matlabbatch{mi}.spm.tools.cat.tools.datatrimming.lazy                          = lazy;




% (2) De-noise the image but also add some basic noise in the background 
%     that is important for SPM to fit Gaussian curves there. Noise is only
%     added in regions without noise!
mi = mi + 1; misanlm = mi; 
matlabbatch{mi}.spm.tools.cat.tools.sanlm.data(1)                              = cfg_dep('Image data trimming: first images of all subjects', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','image_selector', '.','firstimages'));
matlabbatch{mi}.spm.tools.cat.tools.sanlm.spm_type                             = 16;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.prefix                               = 'sanlm_';
matlabbatch{mi}.spm.tools.cat.tools.sanlm.suffix                               = '';
matlabbatch{mi}.spm.tools.cat.tools.sanlm.intlim                               = 100;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.addnoise                             = 0.5; 
    % We have to guaranty some principle noise in noise-free regions.
    % However, I am not sure if this is optimal in this special case.
matlabbatch{mi}.spm.tools.cat.tools.sanlm.rician                               = 0;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.replaceNANandINF                     = 1;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.NCstr               = -Inf;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.iter                = 0;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.iterm               = 0;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.outlier             = 1;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.relativeIntensityAdaption     = 1;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.relativeIntensityAdaptionTH   = 2;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.relativeFilterStengthLimit    = 1;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.resolutionDependency          = 0;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.resolutionDependencyRange     = [1 2.5];
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.red                           = 0;
matlabbatch{mi}.spm.tools.cat.tools.sanlm.nlmfilter.expert.lazy                          = lazy;




% (3) Create an new CSF/background class to have a suitable SPM unified
%     segmentation model. We prepare a 3 class model with GM, WM, and a
%     CSF/background class (1 - (GM + WM)). The result is used as 3rd 
%     tissue class in the unified segmentation below.
%     
% 3 class model 
mi = mi + 1; miBG = mi; 
pp = spm_fileparts(data{1}{1}); % we write into the default output directory
if BGmodel == 2 
  matlabbatch{mi}.spm.util.imcalc.input        = {
                                                  [TPM ',1']; 
                                                  [TPM ',2']; 
                                                 };
  matlabbatch{mi}.spm.util.imcalc.output       = 'TPM3of3_BG';
  matlabbatch{mi}.spm.util.imcalc.expression   = '1 - (i1 + i2)';
else
  matlabbatch{mi}.spm.util.imcalc.input        = {
                                                [TPM ',1']; 
                                                [TPM ',2']; 
                                                [TPM ',3']; 
                                               };
  matlabbatch{mi}.spm.util.imcalc.output       = 'TPM4of4_BG';
  matlabbatch{mi}.spm.util.imcalc.expression   = '1 - (i1 + i2 + i3)';
end
matlabbatch{mi}.spm.util.imcalc.outdir         = {pp};
matlabbatch{mi}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{mi}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{mi}.spm.util.imcalc.options.mask   = 0;
matlabbatch{mi}.spm.util.imcalc.options.interp = 1;
matlabbatch{mi}.spm.util.imcalc.options.dtype  = 4;




% (4) Resize the image for (faster) tests with lower resolution.
%     This is done by a special function that smooth the data to simulate
%     the partial volume effect (PVE) for denoising the data. Otherwise, 
%     only the sample points are used and noise is preserved. 
mi = mi + 1; miResize = mi; 
matlabbatch{mi}.spm.tools.cat.tools.resize.data(1)     = cfg_dep('Spatially adaptive non-local means (SANLM) denoising filter: SANLM Images', substruct('.','val', '{}',{misanlm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{mi}.spm.tools.cat.tools.resize.restype.res = resolution;
matlabbatch{mi}.spm.tools.cat.tools.resize.interp      = -2005;  
matlabbatch{mi}.spm.tools.cat.tools.resize.prefix      = 'auto'; % this will create some resolution specific values
matlabbatch{mi}.spm.tools.cat.tools.resize.outdir      = {''};
% (5) Normalize and limit image intensities.
%     The SPM segmentation seems to have some issues with strong outliers.  
mi = mi + 1; mihlim = mi;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.data(1)    = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{miResize}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'}));
matlabbatch{mi}.spm.tools.cat.tools.spmtype.ctype      = 16;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.prefix     = 'hlim_';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.suffix     = '';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.range      = 99.9999;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.intscale   = 1;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.lazy       = 1;



if BC
  % -----------------------------------------------------------------------
  % only bias correction with very high sampling (<= 1 mm) that was not
  % working quite well for this segmentation itself but helps to model 
  % the brain in the following steps.
  % -----------------------------------------------------------------------
  mi = mi + 1; segment0 = mi; 
  matlabbatch{mi}.spm.spatial.preproc.channel.vols(1)    = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
  matlabbatch{mi}.spm.spatial.preproc.channel.biasreg    = 0;
  matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm   = 30;     % default = 60 (for 1.5 Tesla and precorrected 3.0 Tesla (T1) data in humans)
  matlabbatch{mi}.spm.spatial.preproc.channel.write      = [0 1];  % we also write the bias corrected image
  matlabbatch{mi}.spm.spatial.preproc.tissue(1).tpm      = {[TPM ',1']};
  matlabbatch{mi}.spm.spatial.preproc.tissue(1).ngaus    = ngaus(1);      
  matlabbatch{mi}.spm.spatial.preproc.tissue(1).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(1).warped   = [0 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(2).tpm      = {[TPM ',2']};
  matlabbatch{mi}.spm.spatial.preproc.tissue(2).ngaus    = ngaus(2);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(2).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(2).warped   = [0 0];
  if BGmodel == 2 
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM3of3_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
  else
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm      = {[TPM ',3']};
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(4).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM4of4_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{mi}.spm.spatial.preproc.tissue(4).ngaus    = ngaus(4);    
    matlabbatch{mi}.spm.spatial.preproc.tissue(4).native   = [0 0];
    matlabbatch{mi}.spm.spatial.preproc.tissue(4).warped   = [0 0];
  end
  matlabbatch{mi}.spm.spatial.preproc.warp.mrf           = 0;
  matlabbatch{mi}.spm.spatial.preproc.warp.cleanup       = 0;
  matlabbatch{mi}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
  matlabbatch{mi}.spm.spatial.preproc.warp.affreg        = 'mni';
  matlabbatch{mi}.spm.spatial.preproc.warp.fwhm          = 0; 
  matlabbatch{mi}.spm.spatial.preproc.warp.samp          = 0.5; % there are problems for high-resolution combined with high sampling ... no idea why
  matlabbatch{mi}.spm.spatial.preproc.warp.write         = [0 0];
  matlabbatch{mi}.spm.spatial.preproc.warp.vox           = nan;
  matlabbatch{mi}.spm.spatial.preproc.warp.bb            = [nan nan nan; nan nan nan]; 
  % Create a label map p0 for fast visual checks.
  mi = mi + 1; 
  [pp,ff,ee] = spm_fileparts(data{1}{1}); % we write into the default output directory
  matlabbatch{mi}.spm.util.imcalc.input(1)               = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
  matlabbatch{mi}.spm.util.imcalc.input(2)               = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
  matlabbatch{mi}.spm.util.imcalc.input(3)               = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
  matlabbatch{mi}.spm.util.imcalc.output                 = spm_file([ff ee],'prefix',sprintf('p0masked_r%04.0f_sanlm_trimmed_',resolution*1000));
  matlabbatch{mi}.spm.util.imcalc.outdir                 = {pp};
  matlabbatch{mi}.spm.util.imcalc.expression             = 'i1*2 + i2*3 + i3';
  matlabbatch{mi}.spm.util.imcalc.var                    = struct('name', {}, 'value', {});
  matlabbatch{mi}.spm.util.imcalc.options.dmtx           = 0;
  matlabbatch{mi}.spm.util.imcalc.options.mask           = 0;
  matlabbatch{mi}.spm.util.imcalc.options.interp         = 1;
  matlabbatch{mi}.spm.util.imcalc.options.dtype          = 4;
  % Intensity normalization and limitiation.
  mi = mi + 1; mihlim0 = mi;
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.data(1)    = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.ctype      = 16;
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.prefix     = 'hlimbc_';
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.suffix     = '';
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.range      = 99.999;
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.intscale   = 1;
  matlabbatch{mi}.spm.tools.cat.tools.spmtype.lazy       = lazy;
  % Delete c1-c4 that we don't need anymore (we will do another SPM segmentation) 
  mi = mi + 1; 
  matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
  matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
  matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{segment0}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
  matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
  % -----------------------------------------------------------------------
end




% (6) Unified segmentation with 3 or 4 classes with a strong bias correction. 
%     First a pre-correction that helps to build a more accurate second level model. 
%     Check the p0 images to see the difference (especially the WM in the occipital lobe.
mi = mi + 1; segment1 = mi; 
if BC
  matlabbatch{mi}.spm.spatial.preproc.channel.vols(1)  = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim0}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
 else
  matlabbatch{mi}.spm.spatial.preproc.channel.vols(1)  = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
end
matlabbatch{mi}.spm.spatial.preproc.channel.biasreg    = 0;
matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm   = 30;     % default = 60 (for 1.5 Tesla and precorrected 3.0 Tesla (T1) data in humans)
matlabbatch{mi}.spm.spatial.preproc.channel.write      = [0 1];  % we also write the bias corrected image
matlabbatch{mi}.spm.spatial.preproc.tissue(1).tpm      = {[TPM ',1']};
matlabbatch{mi}.spm.spatial.preproc.tissue(1).ngaus    = ngaus(1);      
matlabbatch{mi}.spm.spatial.preproc.tissue(1).native   = [1 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(1).warped   = [0 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(2).tpm      = {[TPM ',2']};
matlabbatch{mi}.spm.spatial.preproc.tissue(2).ngaus    = ngaus(2);    
matlabbatch{mi}.spm.spatial.preproc.tissue(2).native   = [1 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(2).warped   = [0 0];
if BGmodel == 2 
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM3of3_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
else
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm      = {[TPM ',3']};
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM4of4_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).ngaus    = ngaus(4);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).warped   = [0 0];
end
matlabbatch{mi}.spm.spatial.preproc.warp.mrf           = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.cleanup       = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
matlabbatch{mi}.spm.spatial.preproc.warp.affreg        = 'mni';
matlabbatch{mi}.spm.spatial.preproc.warp.fwhm          = 0; 
matlabbatch{mi}.spm.spatial.preproc.warp.samp          = samp; % there are problems for high-resolution combined with high sampling ... no idea why
matlabbatch{mi}.spm.spatial.preproc.warp.write         = [0 0];
matlabbatch{mi}.spm.spatial.preproc.warp.vox           = nan;
matlabbatch{mi}.spm.spatial.preproc.warp.bb            = [nan nan nan; nan nan nan]; 
% (7) Create a label map p0 for fast visual checks.
mi = mi + 1; 
[pp,ff,ee] = spm_fileparts(data{1}{1}); % we write into the default output directory
matlabbatch{mi}.spm.util.imcalc.input(1)       = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.input(2)       = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.input(3)       = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.output         = spm_file([ff ee],'prefix',sprintf('p0masked_r%04.0f_sanlm_trimmed_',resolution*1000));
matlabbatch{mi}.spm.util.imcalc.outdir         = {pp};
matlabbatch{mi}.spm.util.imcalc.expression     = 'i1*2 + i2*3 + i3';
matlabbatch{mi}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{mi}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{mi}.spm.util.imcalc.options.mask   = 0;
matlabbatch{mi}.spm.util.imcalc.options.interp = 1;
matlabbatch{mi}.spm.util.imcalc.options.dtype  = 4;
% (8) Intensity normalization and limitiation.
mi = mi + 1; mihlim2 = mi;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.data(1)    = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{mi}.spm.tools.cat.tools.spmtype.ctype      = 16;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.prefix     = 'hlim_';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.suffix     = '';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.range      = 99.9999;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.intscale   = 1;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.lazy       = lazy;
% (9) Create a brain masked image for CAT (really required)
if BGmodel == 1  
  mi = mi + 1; miMSK = mi; 
  [pp,ff,ee] = spm_fileparts(data{1}{1}); % we write into the default output directory
  matlabbatch{mi}.spm.util.imcalc.input(1)       = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
  matlabbatch{mi}.spm.util.imcalc.input(2)       = cfg_dep('Segment: c4 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{4}, '.','c', '()',{':'}));
  matlabbatch{mi}.spm.util.imcalc.output         = spm_file([ff ee],'prefix',sprintf('masked_r%04.0f_sanlm_trimmed_',resolution*1000));
  matlabbatch{mi}.spm.util.imcalc.outdir         = {pp};
  matlabbatch{mi}.spm.util.imcalc.expression     = 'i1 .* (i2 < 0.8)'; % X(:,:,:,1) .* smooth3(X(:,:,:,2) < 0.5)';
  matlabbatch{mi}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
  matlabbatch{mi}.spm.util.imcalc.options.dmtx   = 0;
  matlabbatch{mi}.spm.util.imcalc.options.mask   = 0;
  matlabbatch{mi}.spm.util.imcalc.options.interp = 1;
  matlabbatch{mi}.spm.util.imcalc.options.dtype  = 16;
end
% (10) Delete c1-c4 that we don't need anymore (we will do another SPM segmentation) 
mi = mi + 1; 
matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
if BGmodel < 2 
  matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Segment: c4 Images', substruct('.','val', '{}',{segment1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{4}, '.','c', '()',{':'}));
end
matlabbatch{mi}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
% ############
% * a morph-ops batch would be helpful
% * a gradient / divergence batch maybe too
% * cleanup module 
% * bias correction
% * intensity normalization with SPM mat or segments?
% ###########




% (11) Final SPM segmentation
mi = mi + 1; segmentf = mi;
if BGmodel == 1 % skull-stripped
  matlabbatch{mi}.spm.spatial.preproc.channel.vols(1)    = cfg_dep(sprintf('Image Calculator: ImCalc Computed Image: %s',matlabbatch{miMSK}.spm.util.imcalc.output), substruct('.','val', '{}',{miMSK}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
else
  matlabbatch{mi}.spm.spatial.preproc.channel.vols(1)    = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
end
matlabbatch{mi}.spm.spatial.preproc.channel.biasreg    = 0;
matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm   = 30;     % default = 60 (for 1.5 Tesla and precorrected 3.0 Tesla (T1) data in humans)
matlabbatch{mi}.spm.spatial.preproc.channel.write      = [0 1];  % we also write the bias corrected image
matlabbatch{mi}.spm.spatial.preproc.tissue(1).tpm      = {[TPM ',1']};
matlabbatch{mi}.spm.spatial.preproc.tissue(1).ngaus    = ngaus(1);     
matlabbatch{mi}.spm.spatial.preproc.tissue(1).native   = [1 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(1).warped   = [0 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(2).tpm      = {[TPM ',2']};
matlabbatch{mi}.spm.spatial.preproc.tissue(2).ngaus    = ngaus(2);     
matlabbatch{mi}.spm.spatial.preproc.tissue(2).native   = [1 0];
matlabbatch{mi}.spm.spatial.preproc.tissue(2).warped   = [0 0];
if BGmodel == 2
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM3of3_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
else
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).tpm(1)   = {[TPM ',3']};
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).ngaus    = ngaus(3);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(3).warped   = [0 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).tpm(1)   = cfg_dep('Image Calculator: ImCalc Computed Image: TPM4of4_BG', substruct('.','val', '{}',{miBG}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).ngaus    = ngaus(4);    
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).native   = [1 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(4).warped   = [0 0];
end
matlabbatch{mi}.spm.spatial.preproc.warp.mrf           = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.cleanup       = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
matlabbatch{mi}.spm.spatial.preproc.warp.affreg        = 'mni';
matlabbatch{mi}.spm.spatial.preproc.warp.fwhm          = 0;
matlabbatch{mi}.spm.spatial.preproc.warp.samp          = samp; % max(0.5,min(2,resolution * spm_res_factor * 2));
matlabbatch{mi}.spm.spatial.preproc.warp.write         = [0 0];
matlabbatch{mi}.spm.spatial.preproc.warp.vox           = nan;
matlabbatch{mi}.spm.spatial.preproc.warp.bb            = [nan nan nan; nan nan nan]; 
if 0 %BGmodel == 1 
  % apply skull-stripping!
  for ci = 1:3
    mi = mi + 1; 
    matlabbatch{mi}.spm.util.imcalc.input(1)       = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
    matlabbatch{mi}.spm.util.imcalc.input(2)       = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
    matlabbatch{mi}.spm.util.imcalc.output         = spm_file([ff ee],'prefix',sprintf('p0mmasked_r%04.0f_sanlm_trimmed_',resolution*1000));
    matlabbatch{mi}.spm.util.imcalc.outdir         = {pp};
    matlabbatch{mi}.spm.util.imcalc.expression     = 'i1 .* i2';
    matlabbatch{mi}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{mi}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{mi}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{mi}.spm.util.imcalc.options.interp = 1;
    matlabbatch{mi}.spm.util.imcalc.options.dtype  = 4;
  end
end
  
% (12) Create a label map p0* for visual checks. 
mi = mi + 1; 
matlabbatch{mi}.spm.util.imcalc.input(1)       = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.input(2)       = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.input(3)       = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.util.imcalc.output         = spm_file([ff ee],'prefix',sprintf('p0mmasked_r%04.0f_sanlm_trimmed_',resolution*1000));
matlabbatch{mi}.spm.util.imcalc.outdir         = {pp};
matlabbatch{mi}.spm.util.imcalc.expression     = 'i1*2 + i2*3 + i3';
matlabbatch{mi}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
matlabbatch{mi}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{mi}.spm.util.imcalc.options.mask   = 0;
matlabbatch{mi}.spm.util.imcalc.options.interp = 1;
matlabbatch{mi}.spm.util.imcalc.options.dtype  = 4;
% normalize intensities
mi = mi + 1; mihlim3 = mi;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.data(1)    = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.ctype      = 16;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.prefix     = 'hlim_';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.suffix     = '';
matlabbatch{mi}.spm.tools.cat.tools.spmtype.range      = 99.9;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.intscale   = 1;
matlabbatch{mi}.spm.tools.cat.tools.spmtype.lazy       = lazy;



% (6) CAT12 preprocessing with given SPM segmentation without CAT voxel-
%     based preprocessing (especially the AMAP segmentation) to be more 
%     flexible with different image modalities (CAT prefers T1).
mi = mi + 1; 
matlabbatch{mi}.spm.tools.cat.estwrite_spm.data(1)                          = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{segmentf}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.tools.cat.estwrite_spm.nproc                            = nproc;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.T1          = {fullfile(template_dir,'Template_T1.nii')};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.brainmask   = {fullfile(template_dir,'brainmask.nii')};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.cat12atlas  = {fullfile(template_dir,'cat.nii')};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.darteltpm   = {fullfile(template_dir,'Template_1.nii')};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.shootingtpm = {fullfile(template_dir,'Template_0_GS.nii')};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.regstr      = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.bb          = 12;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.vox         = 15; %0.5; 
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.vox                      = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.bb                       = 12;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.pbtres           = 0.8; %min(0.5,max(0.3,resolution));
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.pbtmethod        = 'pbt2x';
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.SRP              = 22;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.reduce_mesh      = reduce_mesh;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.vdist            = 2;  
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.scale_cortex     = 0.7;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.add_parahipp     = 0.1; 
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.close_parahipp   = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.experimental       = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.new_release        = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.lazy               = lazy;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.ignoreErrors       = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.verb               = 2;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.print              = 2;
% ############ Use this part to define the output of both CAT preprocessings ############## 
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.BIDS.BIDSno               = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.surface                   = preview_surfaces;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.lpba40             = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.cobra              = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.hammers            = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.thalamus           = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.ibsr               = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.aal3               = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.mori               = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.anatomy3           = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.julichbrain        = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.atlases.ownatlas  = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.warped                 = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.mod                    = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.dartel                 = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.warped                 = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.mod                    = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.dartel                 = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.warped                = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.mod                   = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.dartel                = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.native              = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.warped              = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.dartel              = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.labelnative               = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.jacobianwarped            = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.warps                     = [1 0];




% (7) Full CAT preprocessing with AMAP segmentation. 
%     However, I am not really satisfied with the results and the GM seems
%     to be underestimated. 
%%{
mi = mi + 1; 
matlabbatch{mi}.spm.tools.cat.estwrite.data(1)                              = cfg_dep('Image data type converter: Converted Images', substruct('.','val', '{}',{mihlim3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
matlabbatch{mi}.spm.tools.cat.estwrite.data_wmh                             = {''};
matlabbatch{mi}.spm.tools.cat.estwrite.nproc                                = nproc;
matlabbatch{mi}.spm.tools.cat.estwrite.useprior                             = '';
matlabbatch{mi}.spm.tools.cat.estwrite.opts.tpm                             = {TPM};
matlabbatch{mi}.spm.tools.cat.estwrite.opts.affreg                          = 'none';
matlabbatch{mi}.spm.tools.cat.estwrite.opts.ngaus                           = [ngaus 4 4]; % default parameters are ok here 
matlabbatch{mi}.spm.tools.cat.estwrite.opts.warpreg                         = [0 0.001 0.5 0.05 0.2];
matlabbatch{mi}.spm.tools.cat.estwrite.opts.bias.spm.biasfwhm               = 30;          % strong bias correction
matlabbatch{mi}.spm.tools.cat.estwrite.opts.bias.spm.biasreg                = 0.00001;     % strong bias correction
matlabbatch{mi}.spm.tools.cat.estwrite.opts.acc.spm.samp                    = samp; %max(0.3,min(1.5,resolution * spm_res_factor));  % higher resolution for SPM preprocessing
matlabbatch{mi}.spm.tools.cat.estwrite.opts.acc.spm.tol                     = tol;        % higher accuracy for SPM preprocessing
matlabbatch{mi}.spm.tools.cat.estwrite.opts.redspmres                       = 0;           
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.restypes.best   = [0.5 0.3];
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.setCOM          = 1;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.APP             = 1070;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.affmod          = 0;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.NCstr           = -Inf;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.LASstr          = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr       = 0;           % not working in PD ?
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr         = -2;          % post-mortem skull-stripping model
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr      = 1.0;         % strong cleanup (default = 0.5)
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr          = 1.0;         % strong correction
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.WMHC            = 0;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.SLC             = 0;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.segmentation.mrf             = 1;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.registration                 = matlabbatch{mi-1}.spm.tools.cat.estwrite_spm.extopts.registration;
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.surface                      = matlabbatch{mi-1}.spm.tools.cat.estwrite_spm.extopts.surface; 
matlabbatch{mi}.spm.tools.cat.estwrite.extopts.admin                        = matlabbatch{mi-1}.spm.tools.cat.estwrite_spm.extopts.admin;
matlabbatch{mi}.spm.tools.cat.estwrite.output                               = matlabbatch{mi-1}.spm.tools.cat.estwrite_spm.output; 

%}