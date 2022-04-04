% Batch file for getting ROI values for SPM12/CAT12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.tools.cat.tools.calcroi.roi_xml = '<UNDEFINED>';

% Entry for output-file string
% Remove comments and edit entry if you would like to change the parameter.
% Otherwise the default value from cat_defaults.m is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.tools.cat.tools.calcroi.calcroi_name = '<UNDEFINED>';

matlabbatch{1}.spm.tools.cat.tools.calcroi.point = '.';
matlabbatch{1}.spm.tools.cat.tools.calcroi.outdir = {''};
