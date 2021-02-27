function cat_batch_cat(namefile,cat_defaults)
% wrapper for using batch mode (see cat_batch_cat.sh)
%
% namefile      - array of file names or text file with file names
% cat_defaults  - use this default file instead of cat_defaults.m
%
%_______________________________________________________________________
% $Id$

 %#ok<*TRYNC>
 
if nargin < 1
  fprintf('Syntax: cat_batch_cat(namefile,cat_defaults)\n');
  return
end

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

% remove extopts fields
extopts_fields = char('NCstr','BVCstr','regstr','WMHC','WMHCstr','mrf','INV','restype','resval','species','darteltpm','shootingtpm',...
            'cat12atlas','brainmask','T1','pbtres','close_parahipp','scale_cortex','add_parahipp','colormap','verb','ignoreErrors',...
            'expertgui','subfolders','experimental','atlas','LAB','print','cleanupstr','SLC','spm_kamap','fontsize','satlas',...
            'send_info','pbtlas','thick_measure','thick_limit','collcorr','nproc','gifti_dat','reduce_mesh','vdist','setCOM',...
            'shootingsurf','pth_templates','shootingT1','report','vox');
for i=1:size(extopts_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.extopts = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.extopts,deblank(extopts_fields(i,:)));
  end
end

% remove output fields
output_fields = char('atlas','te','pc','WMH','ROI','TPMC','label','CSF','WM','GM','las','bias','ct','SL','jacobian','atlases','pp');
for i=1:size(output_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.output = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output,deblank(output_fields(i,:)));
  end
end

% remove opts fields
opts_fields = char('ngaus','warpreg','biasreg','biasfwhm','samp','redspmres','tol','accstr','biasstr');
for i=1:size(opts_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.opts = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.opts,deblank(opts_fields(i,:)));
  end
end

tmp_fields  = char('mod','native','warped','dartel');
tmp_map = char('GM','WM','CSF','bias','las','WMH','label','jacobian','ct','SL');
for i=1:size(tmp_map,1)
  for j=1:size(tmp_fields,1)  
    if isfield(matlabbatch{1}.spm.tools.cat.estwrite.output,(deblank(tmp_map(i,:))))
      if isfield(matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))),deblank(tmp_fields(j,:)))
        matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))) = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))),deblank(tmp_fields(j,:)));
      end
    end
  end
end

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
