function varargout = cat_update(update)
% check for new CAT updates
%
% FORMAT [sts, msg] = cat_update(update)
% sts    - status code:
%        NaN - server not accessible
%        Inf - no updates available
%        0   - CAT installation up-to-date
%        n   - new revision <n> is available for download
% msg    - string describing outcome, that would otherwise be displayed.
% update - allow installation of update
% 
% This function will connect to the SBM server, compare the
% version number of the updates with the one of the CAT12 installation 
% currently in the MATLAB path and will display the result.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

if isdeployed
  sts= Inf;
  msg = 'Update function is not working for compiled CAT12. Please check for a new compiled CAT12 version.';
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return;
end

% Github release url
url_github = 'https://api.github.com/repos/ChristianGaser/cat12/releases';

if ~nargin
  update = false;
else
  update = true;
end

% get current release
r = cat_version;
r = str2double(strrep(r,'CAT',''));

% get new release
try
  [jsonStr, sts] = urlread(url_github,'Timeout',2);
catch
  [jsonStr, sts] = urlread(url_github);
end

if ~sts
  msg = sprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\nYou can download your update at %s\n',url,url); 
  if ~nargout, fprintf(msg); else varargout = {NaN, msg}; end
  return
end

data = jsondecode(jsonStr);
% get largest release number
rnew = [];
for i = 1:length(data)
  rnew = [rnew str2double(data(i).tag_name)];
end
rnew = max(rnew);

if rnew > r
  sts = rnew;
  msg = sprintf('         A new version of CAT12 is available on:\n');
  msg = [msg sprintf('   %s\n',url_github)];
  msg = [msg sprintf('        (Your version: %g - New version: %g)\n',r,rnew)];
  if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
  sts = 0;
  msg = sprintf('Your version of CAT12 is up-to-date.');
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return
end

url = sprintf('https://github.com/ChristianGaser/cat12/releases/download/%g/cat%g.zip',rnew,rnew);

if update
  overwrite = spm_input(sprintf('Update to %g',rnew),1,'m','yes|no|download only',[1 -1 0],1);
  d0 = spm('Dir');
  d  = fileparts(fileparts(which('cat12')));
  
  if overwrite
    try
      % list mex-files and delete these files to prevent that old
      % compiled files are used
      mexfiles = dir(fullfile(fileparts(mfilename('fullpath')),'*.mex*'));
      for i=1:length(mexfiles)
        name = fullfile(fileparts(mfilename('fullpath')),mexfiles(i).name);
        spm_unlink(name);
      end
      
      % delete old html folder
      htmldir = fullfile(fileparts(mfilename('fullpath')),'html');
      if exist(htmldir,'dir'), rmdir(htmldir, 's'); end

      % delete old CAT12 manual
      pdffile = fullfile(fileparts(mfilename('fullpath')),'CAT12-Manual.pdf');
      spm_unlink(pdffile);

      % delete old atlas files
      atlasfiles = dir(fullfile(fileparts(mfilename('fullpath')),'atlases_surfaces','*.*'));
      for i=1:length(atlasfiles)
        name = fullfile(fileparts(mfilename('fullpath')),'atlases_surfaces',atlasfiles(i).name);
        spm_unlink(name);
      end

      % delete old atlas files with 32k meshes
      atlasfiles = dir(fullfile(fileparts(mfilename('fullpath')),'atlases_surfaces_32k','*.*'));
      for i=1:length(atlasfiles)
        name = fullfile(fileparts(mfilename('fullpath')),'atlases_surfaces_32k',atlasfiles(i).name);
        spm_unlink(name);
      end

      % delete old surface template files
      templatefiles = dir(fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','*.*'));
      for i=1:length(templatefiles)
        name = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces',templatefiles(i).name);
        spm_unlink(name);
      end

      % delete old surface template files with 32k meshes
      templatefiles = dir(fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k','*.*'));
      for i=1:length(templatefiles)
        name = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k',templatefiles(i).name);
        spm_unlink(name);
      end

      % delete old volume template files 
      templatefiles = dir(fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym','*.*'));
      for i=1:length(templatefiles)
        name = fullfile(fileparts(mfilename('fullpath')),'templates_volumes',templatefiles(i).name);
        spm_unlink(name);
      end

      templatefiles = dir(fullfile(cat_get_defaults('extopts.pth_templates'),'*.*'));
      for i=1:length(templatefiles)
        name = fullfile(fileparts(mfilename('fullpath')),'templates_volumes',templatefiles(i).name);
        spm_unlink(name);
      end

      lastwarn('');
      warning off
      delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath'); drawnow
      m = '          Download and install CAT12...\n';
      if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
            
      s = unzip(url, d);
      m = sprintf('         Success: %d files have been updated.\n',numel(s));
      if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
      addpath(d0);
      rehash;
      rehash toolboxcache;
      if exist('toolbox_path_cache','file'), toolbox_path_cache; end
      spm fmri; clear cat_version; spm_cat12
      warning on
    catch
      le = lasterror;
      switch le.identifier
          case 'MATLAB:checkfilename:urlwriteError'
              fprintf('          Update failed: cannot download update file.\n');
          otherwise
              fprintf('\n%s\n',le.message);
      end
    end
    
    % open version information if difference between release numbers 
    % is large enough
    if rnew > r+20
      web(fullfile(fileparts(mfilename('fullpath')),'doc','index.html#version'));
    end
    
    [warnmsg, msgid] = lastwarn;
    switch msgid
      case 'MATLAB:extractArchive:unableToCreate'
        fprintf('          Update failed: check folder permission.\n');
      case 'MATLAB:extractArchive:unableToOverwrite'
        fprintf('          Update failed: check file permissions.\n');
      otherwise
        fprintf('          Warning %s.\n',warnmsg);
    end  
  elseif overwrite == 0
    web(url,'-browser');
    fprintf('Unzip file to %s\n',d);
  end
end
