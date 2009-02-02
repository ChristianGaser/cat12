function cg_vbm8_batch(pattern,writeonly)
% wrapper for using batch mode (see cg_vbm8_batch.sh)
%_______________________________________________________________________
% $ID$

if nargin < 1
	fprintf('Syntax: cg_vbm8_batch(pattern)\n');
	return
end

if nargin < 2
	writeonly = 0;
end

addpath(fullfile(spm('dir'),'toolbox','vbm8'));
global defaults
spm_defaults
cg_vbm8_defaults

warning off
% extract folder
folder = fileparts(pattern);
d = dir(pattern);
n = length(d);

if n == 0, error(sprintf('No file %s found in %s.\n',pattern,folder)); end

for i=1:n
	name = fullfile(folder,d(i).name);
	if writeonly
		matlabbatch{1}.spm.tools.vbm8.write = defaults.vbm8;
		matlabbatch{1}.spm.tools.vbm8.write.data{i} = name;
	else
		matlabbatch{1}.spm.tools.vbm8.estwrite = defaults.vbm8;
		matlabbatch{1}.spm.tools.vbm8.estwrite.data{i} = name;
	end
end

spm_jobman('initcfg');
spm_jobman('run_nogui',matlabbatch)

exit
