function cg_vbm8_batch(namefile,writeonly)
% wrapper for using batch mode (see cg_vbm8_batch.sh)
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cg_vbm8_batch(namefile)\n');
	return
end

if nargin < 2
	writeonly = 0;
end

addpath(fullfile(spm('dir'),'toolbox','vbm8'));
spm_defaults
cg_vbm8_defaults
global defaults
spm_jobman('initcfg');

names = textread(namefile,'%s');
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end

if writeonly
	matlabbatch{1}.spm.tools.vbm8.write = defaults.vbm8;
else
	matlabbatch{1}.spm.tools.vbm8.estwrite = defaults.vbm8;
end

for i=1:n
	if writeonly
		matlabbatch{1}.spm.tools.vbm8.write.data{i} = names{i};
	else
		matlabbatch{1}.spm.tools.vbm8.estwrite.data{i} = names{i};
	end
end

spm_jobman('run_nogui',matlabbatch)
spm_unlink(namefile)

exit
