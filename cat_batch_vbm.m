function cat_batch_vbm(namefile,cat_defaults)
% wrapper for using batch mode (see cat_batch_cat12.sh)
%
% namefile      - array of file names
% cat_defaults  - use this default file instead of cat_defaults.m
%
%_______________________________________________________________________
% $Id$

 %#ok<*TRYNC>
 
if nargin < 1
	fprintf('Syntax: cat_batch_vbm(namefile)\n');
	return
end

[t,pid] = system('echo $$');
fprintf('cat_batch_vbm: \n  PID = %s\n\n',pid);

spm_get_defaults;

if nargin < 2
    cat_get_defaults;
else
    if isempty(cat_defaults)
        cat_get_defaults;
    else
        fprintf('Use defaults in %s.\n',cat_defaults);
        [pp, name] = spm_fileparts(cat_defaults);
        oldpath = pwd;
        cd(pp)
        eval(name);
        cd(oldpath)
    end
end
global defaults cat12 matlabbatch %#ok<NUSED>

fid = fopen(namefile,'r');
names = textscan(fid,'%s');
names = names{:};
fclose(fid);
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end %#ok<SPERR>

matlabbatch{1}.spm.tools.cat.estwrite = cat12;
matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr(names);

tmp_fields = char('darteltpm','gcutstr','cleanupstr','mrf','NCstr','BVCstr','LASstr','restype','resval','species',...
              'WMHC','WMHCstr','pbtres','INV','colormap','atlas','print','debug','verb','ignoreErrors',...
              'QAcleanup','QAcleanupth','LAB','vox','bb','cat12atlas','sanlm','expertgui','brainmask','T1','APP','subfolders');
              
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.extopts = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.extopts,deblank(tmp_fields(i,:)));
  end
end

tmp_fields = char('atlas','te','pc','WMH','ROI','TPMC');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.output = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output,deblank(tmp_fields(i,:)));
  end
end


tmp_fields = char('opts','bias','realign','defs','nproc');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite = rmfield(matlabbatch{1}.spm.tools.cat.estwrite,deblank(tmp_fields(i,:)));
  end
end

try 
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.GM,'mod');
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.WM,'mod');
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.CSF,'mod');
end

try 
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.GM,'native');
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.GM,'warped');
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.WM,'native');
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.WM,'warped');
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.bias,'native');
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias  = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.bias,'dartel');
  matlabbatch{1}.spm.tools.cat.estwrite.output = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output,'CSF');
  matlabbatch{1}.spm.tools.cat.estwrite.output = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output,'label');
  matlabbatch{1}.spm.tools.cat.estwrite = rmfield(matlabbatch{1}.spm.tools.cat.estwrite,'estwrite');
end

% deselect multithreading for batch
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

spm_unlink(char(namefile))

exit
