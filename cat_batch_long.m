function cat_batch_long(namefile)
% wrapper for using spm8 batch mode (see cat_batch_cat12.sh)
%
% namefile      - array of file names
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cat_batch_long(namefile)\n');
	exit
end

names = textread(namefile,'%s');
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end

spm_get_defaults;
cat_get_defaults;

global defaults cat matlabbatch

for i=1:n
  matlabbatch{1}.spm.tools.cat.tools.long.subj.mov{i} = names{i};
end

% always deselect print option
matlabbatch{1}.spm.tools.cat.tools.long.extopts.print = 0;

warning off
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

spm_unlink(char(namefile))

exit
