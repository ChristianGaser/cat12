% Batch file for importing DICOM data for SPM12/CAT12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.util.import.dicom.data = '<UNDEFINED>';

% Entry for choosing directory structure
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.util.import.dicom.root = '<UNDEFINED>';

% Entry for choosing output directory
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 2nd parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.util.import.dicom.outdir = '<UNDEFINED>';

% protocol name filter
matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';

% output image format
matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';

% export metadata
matlabbatch{1}.spm.util.import.dicom.convopts.meta = 1;

% use IDEDims in filename
matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

