% batch file of CAT12 simple processing for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.tools.cat.cat_simple.data = '<UNDEFINED>';         

% template for initial affine registration/segmentation; 'adult','children'
matlabbatch{1}.spm.tools.cat.cat_simple.tpm = 'adults';               

matlabbatch{1}.spm.tools.cat.cat_simple.registration.regmethod.shooting.shootingtpm = {fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_0_GS.nii')};

% voxel size and bounding box
matlabbatch{1}.spm.tools.cat.cat_simple.registration.vox = 1.5;
matlabbatch{1}.spm.tools.cat.cat_simple.registration.bb = 12;

% smoothing filter size for volumes
matlabbatch{1}.spm.tools.cat.cat_simple.fwhm_vol = 6;                 

% define here volume atlases
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.hammers = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.ibsr = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.aal3 = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.mori = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.thalamus = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.anatomy3 = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.julichbrain = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.ownatlas = {''};

% catch errors: 0 - stop with error (default); 1 - catch preprocessing errors (requires MATLAB 2008 or higher); 
matlabbatch{1}.spm.tools.cat.cat_simple.ignoreErrors = 1;

% define here surface atlases
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Desikan = 1;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.HCP = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Destrieux = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Schaefer2018_100P_17N = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Schaefer2018_200P_17N = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Schaefer2018_400P_17N = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Schaefer2018_600P_17N = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.ownatlas = {''};

% smoothing filter size for cortical thickness (fwhm_surf1) and gyrification (fwhm_surf2)
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.fwhm_surf1 = 12;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.fwhm_surf2 = 20;

% in order to skip surface processing remove this comment and add comments
% to all line with parameter ".surface.yes"
% matlabbatch{1}.spm.tools.cat.cat_simple.surface.no = 1;

% disable parallel processing
