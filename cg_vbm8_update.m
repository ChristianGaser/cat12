function cg_vbm8_update(update)
% check for new updates
%
% FORMAT cg_vbm8_update(update)
% update - allow installation of update
% 
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the VBM8 installation 
% currently in the MATLAB path and will display the outcome.
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_vbm8_update.m 183 2009-10-07 21:10:26Z gaser $

rev = '$Rev: 183 $';

if nargin == 0
  update = 0;
end

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
if ~sts
  fprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\n',url); 
  return
end

n = regexp(s,'vbm8_r(\d.*?)\.zip','tokens');
if isempty(n)
  fprintf('There are no new releases available yet.\n');
  return;
else
  % get largest release number
  rnew = [];
  for i=1:length(n)
    rnew = [rnew str2double(n{i})];
  end
  rnew = max(rnew);
end

if rnew > r
  fprintf('A new version of VBM8 is available on: %s\n',url);
  fprintf('Your version: %d - New version: %d\n',r,rnew);
  if ~update
    fprintf('In order to update use Toolbox|VBM8|Check for updates\n',r,rnew);
  end

  if update
    d = fullfile(spm('Dir'),'toolbox'); 
    overwrite = spm_input('Update',1,'m','Download zip-file only|Overwrite old VBM8 installation',[0 1],2);
    if overwrite
      try
        % list mex-files and delete these files to prevent that old
        % compiled files are used
        mexfiles = dir(fullfile(d,'vbm8','*.mex*'));
        for i=1:length(mexfiles)
          name = fullfile(d,'vbm8',mexfiles(i).name);
          spm_unlink(name);
          disp(name)
        end
        s = unzip([url sprintf('vbm8_r%d.zip',rnew)], d);
        fprintf('%d files have been updated.\nSPM should be restarted.\n',numel(s));
        restart = spm_input('Restart SPM',1,'m','no|yes',[0 1],2);
        if restart
          eval(['spm fmri']);
        end
      catch
        fprintf('Update failed: check file permissions. Download zip-file only.\n');
        web([url sprintf('vbm8_r%d.zip',rnew)],'-browser');
        fprintf('Unzip file to %s\n',d);
      end
    else
      web([url sprintf('vbm8_r%d.zip',rnew)],'-browser');
      fprintf('Unzip file to %s\n',d);
    end
  end
else
  if ~update
    fprintf('You already have the newest version %d.\n',r);
  end
end