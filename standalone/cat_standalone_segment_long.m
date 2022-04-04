% Batch file for CAT12 longitudinal segmentation for SPM12/CAT12 standalone installation
%
%_______________________________________________________________________
% $Id$

% first undefined data field, that will be dynamically replaced by cat_standalone.sh
% The different definitions of the subjects-field are necessary to be compatible 
% with CAT12 longitudinal batch (using "{}") and cat_standalone where the
% UNDEFINED field is necessary. The clear command prevents error due to different
% datatypes and the comented out part for cat_standalone will be removed in the
% shell script and the last definition of subjects is finally used. Looks weird,
% but only works in that way.
matlabbatch{1}.spm.tools.cat.long.datalong.subjects = {};
clear matlabbatch
%matlabbatch{1}.spm.tools.cat.long.datalong.subjects = '<UNDEFINED>';

% Entry for choosing longitudinal model
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
% (0) large changes with brain/head growth (i.e. developmental effects)
% (1) small changes (i.e. plasticity/learning effects)
% (2) large changes (i.e. aging effects)
% (3) save results for both models 1 and 2
%matlabbatch{1}.spm.tools.cat.long.longmodel = '<UNDEFINED>';

% Entry for choosing TPM
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 2nd parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.tools.cat.long.opts.tpm = '<UNDEFINED>';

% use priors for longitudinal data
matlabbatch{1}.spm.tools.cat.long.enablepriors = 1;

% Remove comments and edit entry if you would like to change the Dartel/Shooting approach
% Otherwise the default value from cat_defaults.m is used.
% entry for choosing shooting approach
%matlabbatch{1}.spm.tools.cat.long.extopts.registration.regmethod.shooting.shootingtpm = {fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_0_GS.nii')};
% entry for choosing dartel approach
%matlabbatch{1}.spm.tools.cat.long.extopts.registration.regmethod.dartel.darteltpm = {fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_1_Dartel.nii')};
% Strength of Shooting registration: 0 - Dartel, eps (fast), 0.5 (default) to 1 (accurate) optimized Shooting, 4 - default Shooting; default 0.5
%matlabbatch{1}.spm.tools.cat.long.extopts.registration.regmethod.shooting.regstr = 0.5;

% additional bounding box
matlabbatch{1}.spm.tools.cat.long.extopts.registration.bb = 12;

% Affine regularisation (SPM12 default = mni) - '';'mni';'eastern';'subj';'none';'rigid'
matlabbatch{1}.spm.tools.cat.long.opts.affreg = 'mni';

% Strength of the bias correction that controls the biasreg and biasfwhm parameter (CAT only!)
% 0 - use SPM parameter; eps - ultralight, 0.25 - light, 0.5 - medium, 0.75 - strong, and 1 - heavy corrections
% job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
% job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*job.opts.biasstr ));  
matlabbatch{1}.spm.tools.cat.long.opts.biasstr = 0.5;

matlabbatch{1}.spm.tools.cat.long.opts.accstr = 0.5;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.spm_kamap = 0;

% surface options
matlabbatch{1}.spm.tools.cat.long.extopts.surface.pbtres = 0.5;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{1}.spm.tools.cat.long.extopts.surface.SRP = 22;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.vdist = 1.33333333333333;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.scale_cortex = 0.7;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.add_parahipp = 0.1;
matlabbatch{1}.spm.tools.cat.long.extopts.surface.close_parahipp = 0;

matlabbatch{1}.spm.tools.cat.long.extopts.admin.experimental = 0;
matlabbatch{1}.spm.tools.cat.long.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.long.extopts.admin.lazy = 0;
matlabbatch{1}.spm.tools.cat.long.extopts.admin.ignoreErrors = 1;
matlabbatch{1}.spm.tools.cat.long.extopts.admin.verb = 2;
matlabbatch{1}.spm.tools.cat.long.extopts.admin.print = 2;

matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.NCstr = -Inf;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.BVCstr = 0.5;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.WMHC = 2;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.SLC = 0;
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.mrf = 1;

% Affine PreProcessing (APP) with rough bias correction and brain extraction for special anatomies (nonhuman/neonates) 
% 0 - none; 1070 - default; [1 - light; 2 - full; 1144 - update of 1070, 5 - animal (no affreg)]
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.APP = 1070;

% Strength of the local adaptation: 0 to 1; default 0.5
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.LASstr = 0.5;

% Strength of skull-stripping: 0 - SPM approach; eps to 1  - gcut; 2 - new APRG approach; -1 - no skull-stripping (already skull-stripped); default = 2
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.gcutstr = 2;

% voxel size for normalized data (EXPERIMENTAL: inf - use Template values)
matlabbatch{1}.spm.tools.cat.long.extopts.registration.vox = 1.5;

% resolution handling: 'native','fixed','best', 'optimal'
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.restypes.optimal = [1 0.3];

% use center-of-mass approach for estimating origin
matlabbatch{1}.spm.tools.cat.long.extopts.segmentation.setCOM = 1;

% surface and thickness creation:   0 - no (default), 1 - lh+rh, 2 - lh+rh+cerebellum, 
% 3 - lh, 4 - rh, 5 - lh+rh (fast, no registration, only for quick quality check and not for analysis),
% 6 - lh+rh+cerebellum (fast, no registration, only for quick quality check and not for analysis)
% 9 - thickness only (for ROI analysis, experimental!)
% +10 to estimate WM and CSF width/depth/thickness (experimental!)
matlabbatch{1}.spm.tools.cat.long.output.surface = 1;

% BIDS output
matlabbatch{1}.spm.tools.cat.long.output.BIDS.BIDSno = 1;

% define here volume atlases
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.aal3 = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.long.ROImenu.atlases.ownatlas = {''};

% create and use longitudinal TPM to get more stable segmentations
matlabbatch{1}.spm.tools.cat.long.longTPM = 1;

% apply modulation
matlabbatch{1}.spm.tools.cat.long.modulate = 1;

% save dartel export
matlabbatch{1}.spm.tools.cat.long.dartel = 0;

% delete temporary files
matlabbatch{1}.spm.tools.cat.long.delete_temp = 1;
