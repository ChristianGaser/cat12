function cg_vbm_batch(namefile,vbm_defaults)
% wrapper for using batch mode (see cg_vbm_batch.sh)
%
% namefile      - array of file names
% vbm_defaults  - use this default file instead of cg_vbm_defaults.m
%
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cg_vbm_batch(namefile)\n');
	return
end

[t,pid]=system('echo $$');
fprintf('cg_vbm_batch: \n  PID = %s\n\n',pid);

spm_get_defaults;

if nargin < 3
    cg_vbm_get_defaults;
else
    if isempty(vbm_defaults)
        cg_vbm_get_defaults;
    else
        fprintf('Use defaults in %s.\n',vbm_defaults);
        [pp, name] = spm_fileparts(vbm_defaults);
        oldpath = pwd;
        cd(pp)
        eval(name);
        cd(oldpath)
    end
end
global defaults vbm matlabbatch

% always deselect print option
vbm.extopts.print = 0;

names = textread(namefile,'%s');
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end

matlabbatch{1}.spm.tools.vbm.estwrite = vbm;

for i=1:n
	matlabbatch{1}.spm.tools.vbm.estwrite.data{i} = names{i};
end

tmp_fields = char('darteltpm','gcutstr','cleanupstr','mrf','NCstr','BVCstr','LASstr','restype','resval','species',...
              'WMHC','WMHCstr','pbtres','INV','colormap','atlas','print','debug','verb','ignoreErrors',...
              'QAcleanup','QAcleanupth','LAB','vox','bb','vbm12atlas','sanlm','gui','brainmask','T1');
              
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.vbm.estwrite.extopts = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.extopts,deblank(tmp_fields(i,:)));
  end
end

tmp_fields = char('atlas','te','pc','WMH');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.vbm.estwrite.output = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output,deblank(tmp_fields(i,:)));
  end
end


tmp_fields = char('opts','bias','realign','defs');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.vbm.estwrite = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite,deblank(tmp_fields(i,:)));
  end
end

try
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM  = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.GM,'mod');
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM  = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.WM,'mod');
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF,'mod');
end

%save(fullfile(spm('dir'),'tmp.vbm.batch.mat'));

try
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
catch
  %vbmerr = lasterror; 
  %sprintf('\n%s\nVBM Preprocessing error: %s:\n%s\n', repmat('-',1,72),vbmerr.identifier,vbmerr.message,repmat('-',1,72));
  %for si=1:numel(vbmerr.stack), vbm_io_cprintf('err',sprintf('%5d - %s\n',vbmerr.stack(si).line,vbmerr.stack(si).name)); end;
  %vbm_io_cprintf('err',sprintf('%s\\n',repmat('-',1,72))); exit; end; 
  error('Batch failed.');
end

spm_unlink(char(namefile))

exit
