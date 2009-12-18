function cg_vbm8_update
% check for new updates
%
% FORMAT cg_vbm8_update
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the VBM8 installation 
% currently in the MATLAB path and will display the outcome.
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_vbm8_update.m 183 2009-10-07 21:10:26Z gaser $

rev = '$Rev: 183 $';

r = 0;

% get current release number
A = ver;
for i=1:length(A)
  if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
    r = str2double(A(i).Version);
  end
end

url = 'http://dbm.neuro.uni-jena.de/vbm8/';

% get new release number
[s,sts] = urlread(url);
if ~sts, error('Cannot access the SBM server.'); end
n = regexp(s,'vbm8_r(\d.*?)\.zip','tokens','once');
if isempty(n)
  fprintf('There are no new releases available yet.\n');
  return;
else
  n = str2double(n{1});
end

if n > r
  fprintf('A new version of VBM8 is available on: %s\n',url);
  fprintf('Your version: %d - New version: %d\n',r,n);

  d = fullfile(spm('Dir'),'toolbox','vbm8'); 
  overwrite = spm_input('Update',1,'m','Download zip-file only|Overwrite old VBM8 installation',[0 1],2);
  if overwrite
    try
      s = unzip([url sprintf('vbm8_r%d.zip',n)], d);
      fprintf('%d files have been updated.\n',numel(s));
    catch
      fprintf('Update failed: check file permissions. Download only\n');
      web([url sprintf('vbm8_r%d.zip',n)],'-browser');
    end
  else
    web([url sprintf('vbm8_r%d.zip',n)],'-browser');
  end
else
  fprintf('You already have the newest version %d\n',r);
end