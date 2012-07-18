function cg_spm12_batch(batchname)
% wrapper for using spm12 batch mode (see cg_vbm_batch.sh)
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cg_spm12_batch(batchname)\n');
	exit
end

spm_get_defaults;
global defaults

if ~exist(batchname,'file')
	fprintf('Batchfile %s not found\n',batchname);
	exit
end

eval(batchname)

if ~exist('matlabbatch','var')
	fprintf('Batchfile %s did not returned variable matlabbatch.\n', batchname);
	exit
end

warning off
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

exit
