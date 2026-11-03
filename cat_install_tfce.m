function cat_install_tfce(install)
% Download the latest release of the TFCE toolbox from Github and install
% it in the SPM toolbox folder
%
% FORMAT cat_install_tfce(install)
% install - install without asking (otherwise the user is asked)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin == 0
  install = spm_input('Install TFCE toolbox',1,'yes|no',[1 0],1);
end

if ~install, return; end

url_api      = 'https://api.github.com/repos/ChristianGaser/tfce/releases/latest';
url_releases = 'https://github.com/ChristianGaser/tfce/releases';
d0 = spm('Dir');
d  = fullfile(d0,'toolbox');

% get the zip-file of the latest TFCE release from Github
try
  release = webread(url_api, weboptions('Timeout',10));
  assets  = release.assets;
  if ~iscell(assets), assets = num2cell(assets); end
  url = '';
  for i = 1:numel(assets)
    [~,~,ext] = fileparts(assets{i}.name);
    if strcmpi(ext,'.zip')
      url = assets{i}.browser_download_url;
      break
    end
  end
  if isempty(url)
    error('cat_install_tfce:noZip','Release %s contains no zip-file.',release.tag_name);
  end
catch err
  fprintf('          Cannot get TFCE release from Github: %s\n',err.message);
  fprintf('          Please download the zip-file from %s\n          and unzip it to %s.\n',url_releases,d);
  return
end

% check write permission before downloading (e.g. SPM installed by an admin)
testfile = fullfile(d,'cat_install_tfce.tmp');
fid = fopen(testfile,'w');
if fid < 0
  fprintf('          Installation failed: no write permission for %s.\n',d);
  fprintf('          Please download %s\n          and unzip it to %s.\n',url,d);
  return
end
fclose(fid); delete(testfile);

% the zip-file contains the folder TFCE
try
  fprintf('          Download and install TFCE %s...\n',release.tag_name);
  lastwarn('');
  s = unzip(url, d);
  warnmsg = lastwarn;
catch err
  fprintf('          Installation failed: %s\n',err.message);
  fprintf('          Please download %s\n          and unzip it to %s.\n',url,d);
  return
end

% unzip only warns if files cannot be written (the warning identifiers
% differ between MATLAB versions)
if ~isempty(warnmsg)
  fprintf('          Installation failed: %s\n',warnmsg);
  fprintf('          Please check the file permissions in %s.\n',fullfile(d,'TFCE'));
  return
end
fprintf('          Success: %d files have been downloaded.\n',numel(s));

% restart SPM and CAT to add TFCE to the SPM toolbox menu
try
  addpath(d0);
  rehash
  rehash toolboxcache;
  if exist('toolbox_path_cache','file'), toolbox_path_cache; end
  spm fmri; clear cat_version; spm_CAT
catch err
  fprintf('\n%s\n',err.message);
end
