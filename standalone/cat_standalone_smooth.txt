% batch file of volume smoothing for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.spatial.smooth.data = '<UNDEFINED>';

% Entry for choosing smoothing filter size
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.spatial.smooth.fwhm = '<UNDEFINED>';

% Entry for prepending string for smoothed file (e.g. 's6')
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.spatial.smooth.prefix = '<UNDEFINED>';

matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;

