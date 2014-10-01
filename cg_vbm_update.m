function varargout = cg_vbm_update(update)
% check for new VBM updates
%
% FORMAT [sts, msg] = cg_vbm_update(update)
% sts    - status code:
%        NaN - SPM server not accessible
%        Inf - no updates available
%        0   - SPM installation up to date
%        n   - new revision <n> is available for download
% msg    - string describing outcome, that would otherwise be displayed.
% update - allow installation of update
% 
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the VBM12 installation 
% currently in the MATLAB path and will display the result.
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

url = 'http://dbm.neuro.uni-jena.de/vbm12/';

if ~nargin
    update = false;
else
    update = true;
end

r = 0;

% get current release number
A = ver;
for i=1:length(A)
  if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
    r = str2double(A(i).Version);
  end
end

% get new release numbers
[s,sts] = urlread(url);
if ~sts
  sts = NaN;
  msg = sprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\n. You can download your update at %s\n',url,url); 
  if ~nargout, error(msg); else varargout = {sts, msg}; end
  return
end

n = regexp(s,'vbm12_r(\d.*?)\.zip','tokens');
if isempty(n)
  sts= Inf;
  msg = 'There are no new releases available yet.';
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
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
    sts = n;
    msg = sprintf('         A new version of VBM12 is available on:\n');
    msg = [msg sprintf('   %s\n',url)];
    msg = [msg sprintf('        (Your version: %d - New version: %d)\n',r,rnew)];
    if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
    sts = 0;
    msg = sprintf('Your version of VBM12 is up to date.');
    if ~nargout, fprintf([blanks(9) msg '\n']);
    else varargout = {sts, msg}; end
    return
end

if update
    d = fullfile(spm('Dir'),'toolbox'); 
    overwrite = spm_input('Update',1,'yes|no',[1 0],1);
    if overwrite
      try
        lastwarn('');
        % list mex-files and delete these files to prevent that old
        % compiled files are used
        mexfiles = dir(fullfile(d,'vbm12','*.mex*'));
        for i=1:length(mexfiles)
          name = fullfile(d,'vbm12',mexfiles(i).name);
          spm_unlink(name);
        end
        m = '          Download and install VBM12...\n';
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
        s = unzip([url sprintf('vbm12_r%d.zip',rnew)], d);
        m = sprintf('         Success: %d files have been updated.\n',numel(s));
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
        rehash
        rehash toolboxcache;
        toolbox_path_cache
        eval(['spm fmri']);
      catch
        le = lasterror;
        switch le.identifier
            case 'MATLAB:checkfilename:urlwriteError'
                fprintf('          Update failed: cannot download update file.\n');
            otherwise
                fprintf('\n%s\n',le.message);
        end
      end
      [warnmsg, msgid] = lastwarn;
      switch msgid
        case ''
        case 'MATLAB:extractArchive:unableToCreate'
            fprintf('          Update failed: check folder permission.\n');
        case 'MATLAB:extractArchive:unableToOverwrite'
            fprintf('          Update failed: check file permissions.\n');
        otherwise
            fprintf('          Update failed: %s.\n',warnmsg);
      end
    end
end
