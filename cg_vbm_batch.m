function cg_vbm_batch(namefile,writeonly,vbm_defaults)
% wrapper for using batch mode (see cg_vbm_batch.sh)
%
% namefile      - array of file names
% writeonly     - if "1" do not estimate segmentations
% vbm_defaults - use this default file instead of cg_vbm_defaults.m
%
%_______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cg_vbm_batch(namefile)\n');
	return
end

if nargin < 2
	writeonly = 0;
end

spm_get_defaults;

if nargin < 3
    cg_vbm_get_defaults;
else
    if isempty(vbm_defaults)
        cg_vbm_get_defaults;
    else
        fprintf('Use defaults in %s.\n',vbm_defaults);
        [path, name] = spm_fileparts(vbm_defaults);
        oldpath = pwd;
        eval(['cd ' path])
        eval(name);
        eval(['cd ' oldpath])
    end
end
global defaults vbm matlabbatch

% always deselect print option
vbm.extopts.print = 0;

names = textread(namefile,'%s');
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end

if writeonly
	matlabbatch{1}.spm.tools.vbm.write = vbm;
else
	matlabbatch{1}.spm.tools.vbm.estwrite = vbm;
end

for i=1:n
	if writeonly
		matlabbatch{1}.spm.tools.vbm.write.data{i} = names{i};
	else
		matlabbatch{1}.spm.tools.vbm.estwrite.data{i} = names{i};
	end
end

tmp_fields = char('darteltpm','finalmask','gcut','kmeans','mrf','bias_fwhm','BVC','pbtres','mask','INV','colormap');
if writeonly
    matlabbatch{1}.spm.tools.vbm.write.extopts = rmfield(matlabbatch{1}.spm.tools.vbm.write.extopts,tmp_fields);
else
    matlabbatch{1}.spm.tools.vbm.estwrite.extopts = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.extopts,tmp_fields);
end

tmp_fields = char('l1','ml','pc','te','surf');
try
  if writeonly
    matlabbatch{1}.spm.tools.vbm.write.output = rmfield(matlabbatch{1}.spm.tools.vbm.write.output,tmp_fields);
  else
    matlabbatch{1}.spm.tools.vbm.estwrite.output = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output,tmp_fields);
  end
end

tmp_fields = char('opts','bias','realign','defs');
try
  if writeonly
    matlabbatch{1}.spm.tools.vbm.write = rmfield(matlabbatch{1}.spm.tools.vbm.write,tmp_fields);
  else
    matlabbatch{1}.spm.tools.vbm.estwrite = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite,tmp_fields);
  end
end

try
  if writeonly
    matlabbatch{1}.spm.tools.vbm.write.output.GM  = rmfield(matlabbatch{1}.spm.tools.vbm.write.output.GM,'mod');
    matlabbatch{1}.spm.tools.vbm.write.output.WM  = rmfield(matlabbatch{1}.spm.tools.vbm.write.output.WM,'mod');
    matlabbatch{1}.spm.tools.vbm.write.output.CSF = rmfield(matlabbatch{1}.spm.tools.vbm.write.output.CSF,'mod');
    matlabbatch{1}.spm.tools.vbm.write.output.mgT = rmfield(matlabbatch{1}.spm.tools.vbm.write.output.mgT,'mod');
  else
    matlabbatch{1}.spm.tools.vbm.estwrite.output.GM  = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.GM,'mod');
    matlabbatch{1}.spm.tools.vbm.estwrite.output.WM  = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.WM,'mod');
    matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF,'mod');
    matlabbatch{1}.spm.tools.vbm.estwrite.output.mgT = rmfield(matlabbatch{1}.spm.tools.vbm.estwrite.output.mgT,'mod');
  end
end

try
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
catch
  error('Batch failed.');
end

spm_unlink(char(namefile))

exit
