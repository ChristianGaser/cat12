% batch file of CAT12 simple processing for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id: cat_batch_simple_standalone.m 1510 2019-10-16 10:12:29Z gaser $

% data field, will be dynamically replaced by cat_batch_standalone.sh
matlabbatch{1}.spm.tools.cat.cat_simple.data = '<UNDEFINED>';         

% CAT12 version, see simple processing batch for more options
% other catversions are 'estwrite1445'
matlabbatch{1}.spm.tools.cat.cat_simple.catversion = 'estwrite';      

% template for initial affine registration/segmentation; 'adult','children'
matlabbatch{1}.spm.tools.cat.cat_simple.tpm = 'adults';               

% smoothing filter size for volumes
matlabbatch{1}.spm.tools.cat.cat_simple.fwhm_vol = 8;                 

% define here volume atlases
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.lpba40 = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.cobra = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.ROImenu.atlases.hammers = 0;

% define here surface atlases
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Desikan = 1;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.HCP = 0;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.sROImenu.satlases.Destrieux = 0;

% smoothing filter size for cortical thickness (fwhm_surf1) and gyrification (fwhm_surf2)
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.fwhm_surf1 = 12;
matlabbatch{1}.spm.tools.cat.cat_simple.surface.yes.fwhm_surf2 = 20;

% in order to skip surface processing remove this comment and add comments
% to all line with parameter ".surface.yes"
% matlabbatch{1}.spm.tools.cat.cat_simple.surface.no = 1;

