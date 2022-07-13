function cat_batch_cat(namefile,cat_defaults)
% wrapper for using batch mode (see cat_batch_cat.sh)
%
% namefile      - array of file names or text file with file names
% cat_defaults  - use this default file instead of cat_defaults.m
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

 %#ok<*TRYNC>
 
if nargin < 1
  fprintf('Syntax: cat_batch_cat(namefile,cat_defaults)\n');
  return
end

addpath(fileparts(which(mfilename)))

[t,pid] = system('echo $$');
fprintf('cat_batch_cat: \n  PID = %s\n\n',pid);

global defaults cat matlabbatch %#ok<NUSED>

spm_get_defaults;

if nargin < 2
    cat_get_defaults;
else
    if isempty(cat_defaults)
        cat_get_defaults;
    else
        fprintf('Use defaults in %s.\n',cat_defaults);
        [pp, name] = spm_fileparts(cat_defaults);
        clear cat_defaults
        oldpath = pwd;
        if ~isempty(pp), cd(pp); end
        eval(name);
        cd(oldpath)
    end
end

if ~iscell(namefile)
  [pth,nam,ext] = spm_fileparts(namefile);
end

% check whether namefile is a cell of filenames, a nifti filename,
% or a text file with filenames
if iscell(namefile)
  names0 = namefile;
  is_filelist = 1;
elseif strcmp(ext,'.nii') | strcmp(ext,'.img')
  names0 = cellstr(namefile);
  is_filelist = 1;
else % or use list of names in text file
  fid = fopen(namefile,'r');
  names0 = textscan(fid,'%s');
  names0 = names0{:};
  fclose(fid);
  is_filelist = 0;
end

n = length(names0);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end %#ok<SPERR>

i = 1;
while i <= n
  % if no .nii or .img was found assume that the filenames contains spaces and is therefore divided into
  % different cells
  if isempty(strfind(names0{i},'.nii')) && isempty(strfind(names0{i},'.img')) && i<length(names0)
    names{i,1} = [deblank(names0{i}) ' ' deblank(names0{i+1})];
    i = i+1;
  else
    names{i,1} = deblank(names0{i});
  end
  i = i+1;
end

matlabbatch{1}.spm.tools.cat.estwrite = cat;
matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr(names);

% remove fields to suppress warnings
matlabbatch{1}.spm.tools.cat.estwrite = rmfield(matlabbatch{1}.spm.tools.cat.estwrite,'extopts');
matlabbatch{1}.spm.tools.cat.estwrite = rmfield(matlabbatch{1}.spm.tools.cat.estwrite,'output');
matlabbatch{1}.spm.tools.cat.estwrite = rmfield(matlabbatch{1}.spm.tools.cat.estwrite,'opts');

% deselect multi-threading for batch
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;

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

% delete text file with filenames
if ~is_filelist, spm_unlink(char(namefile)); end

warning off
exit
