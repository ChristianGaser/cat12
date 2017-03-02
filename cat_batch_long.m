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

fid = fopen(namefile,'r');
names = textscan(fid,'%s');
names = names{:};
fclose(fid);

n = length(names);

if n == 0, error('No file found in %s.\n',namefile); end

spm_get_defaults;
cat_get_defaults;

global defaults cat12 matlabbatch

matlabbatch{1}.spm.tools.cat.tools.long.subj.mov = cell(n,1);
for i=1:n
  matlabbatch{1}.spm.tools.cat.tools.long.subj.mov{i} = names{i};
end

matlabbatch{1}.spm.tools.cat.tools.long.modulate = 1;
matlabbatch{1}.spm.tools.cat.tools.long.warps = 0;

% always deselect print option
matlabbatch{1}.spm.tools.cat.tools.long.extopts.print = 0;

warning off
try
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
catch %#ok<CTCH> % catch with lasterror is necessary for old matlab versions
  caterr = lasterror;  %#ok<LERR>
  sprintf('\n%s\nCAT Preprocessing error: %s:\n%s\n', repmat('-',1,72),caterr.identifier,caterr.message,repmat('-',1,72));
  for si=1:numel(caterr.stack), cat_io_cprintf('err',sprintf('%5d - %s\n',caterr.stack(si).line,caterr.stack(si).name)); end;
  cat_io_cprintf('err',sprintf('%s\\n',repmat('-',1,72)));  
  error('Batch failed.');
end

spm_unlink(char(namefile))

exit
