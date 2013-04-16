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
    cg_vbm_defaults;
else
    if isempty(vbm_defaults)
        cg_vbm_defaults;
    else
        fprintf('Use defaults in %s.\n',vbm_defaults);
        [path, name] = spm_fileparts(vbm_defaults);
        oldpath = pwd;
        eval(['cd ' path])
        eval(name);
        eval(['cd ' oldpath])
    end
end
global defaults vbm

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

try
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
end

spm_unlink(char(namefile))

exit
