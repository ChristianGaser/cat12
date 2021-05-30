% batch file of TFCE statistics for SPM12 standalone installation
%
%_______________________________________________________________________
% $Id$

% data field, that will be dynamically replaced by cat_standalone.sh
matlabbatch{1}.spm.tools.tfce_estimate.data = '<UNDEFINED>';

% Entry for contrast number
% Remove comments and edit entry if you would like to change the contrast number.
% Otherwise the default value is used.
% Or use 1st parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.tools.tfce_estimate.conspec.contrasts = '<UNDEFINED>';

% Entry for number of permutations
% Remove comments and edit entry if you would like to change the number of permutations.
% Otherwise the default value is used.
% Or use 2nd parameter field, that will be dynamically replaced by cat_standalone.sh
%matlabbatch{1}.spm.tools.tfce_estimate.conspec.n_perm = '<UNDEFINED>';

matlabbatch{1}.spm.tools.tfce_estimate.nproc = 0;
matlabbatch{1}.spm.tools.tfce_estimate.mask = '';
matlabbatch{1}.spm.tools.tfce_estimate.conspec.titlestr = '';
matlabbatch{1}.spm.tools.tfce_estimate.nuisance_method = 2;
matlabbatch{1}.spm.tools.tfce_estimate.tbss = 0;
matlabbatch{1}.spm.tools.tfce_estimate.E_weight = 0.5;
matlabbatch{1}.spm.tools.tfce_estimate.singlethreaded = 0;

