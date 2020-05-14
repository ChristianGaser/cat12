function cat_batch_long(namefile,output_surface,cat_defaults)
% wrapper for using batch mode (see cat_batch_long.sh)
%
% namefile       - array of file names
% output_surface - enable surface estimation
% cat_defaults   - use this default file instead of cat_defaults.m
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cat_batch_long(namefile)\n');
	exit
end

if nargin < 2
  output_surface = 0;
else
  % string argument has to be converted 
  if isstr(output_surface)
    output_surface = str2num(output_surface);
  end
end

fid = fopen(namefile,'r');
names = textscan(fid,'%s');
names = names{:};
fclose(fid);

n = length(names);

if n == 0, error('No file found in %s.\n',namefile); end

global defaults cat matlabbatch

spm_get_defaults;

if nargin < 3
    cat_get_defaults;
else
    if isempty(cat_defaults)
        cat_get_defaults;
    else
        fprintf('Use defaults in %s.\n',cat_defaults);
        [pp, name] = spm_fileparts(cat_defaults);
        clear cat_defaults
        oldpath = pwd;
        cd(pp)
        eval(name);
        cd(oldpath)
    end
end

matlabbatch{1}.spm.tools.cat.long.subj.mov = cell(n,1);
matlabbatch{1}.spm.tools.cat.long.nproc = 0;

for i=1:n
  matlabbatch{1}.spm.tools.cat.long.subj.mov{i} = names{i};
end

matlabbatch{1}.spm.tools.cat.long.modulate = 1;

if output_surface == 1
  matlabbatch{1}.spm.tools.cat.long.output.surface = 1;
end

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
