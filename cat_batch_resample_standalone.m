% batch file of CAT12 Resample & Smooth for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id: cat_batch_resample_standalone.m 1510 2019-10-16 10:12:29Z gaser $

% data field, will be dynamically replaced by cat_batch_standalone.sh
matlabbatch{1}.spm.tools.cat.stools.surfresamp.data_surf = {'<UNDEFINED>'};

% merge hemispheres and use 32k mesh from HCP
matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = 1;

% smoothing filter size surface values
% use 12-15mm for cortical thickness and 20-25mm for folding measures
matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = 12;
