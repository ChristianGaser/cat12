function cg_spm8_batch(batchname)
% wrapper for using spm8 batch mode (see cg_vbm8_batch.sh)
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cg_spm8_batch(batchname)\n');
	exit
end

spm_get_defaultsg
lobal defaults

if ~exist(batchname,'file')
	fprintf('Batchfile %s not found\n',batchname);
	exit
end

eval(batchname)

if ~exist('matlabbatch','var')
	fprintf('Batchfile %s did not returned variable matlabbatch.\n', batchname);
	exit
end

spm_jobman('initcfg');
spm_jobman('run_nogui',matlabbatch)

exit
