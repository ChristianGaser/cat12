% Batch file for CAT12 Resample & Smooth for SPM12/CAT12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf = '<UNDEFINED>';

% Entry for choosing smoothing filter size surface values
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
% use 12-15mm for cortical thickness and 20-25mm for folding measures
%matlabbatch{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = '<UNDEFINED>';

% Entry for using 32k mesh from HCP or 164k mesh from Freesurfer
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 2nd parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.tools.cat.stools.surfresamp.mesh32k = '<UNDEFINED>';

% merge hemispheres?
matlabbatch{1}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;

% set this to 1 for skipping processing if already processed data exist
matlabbatch{1}.spm.tools.cat.stools.surfresamp.lazy = 0;

% disable parallel processing
matlabbatch{1}.spm.tools.cat.stools.surfresamp.nproc = 0;
