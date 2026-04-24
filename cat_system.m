function [status,result] = cat_system(cmd,verb,trerr)
% ______________________________________________________________________
% CAT wrapper for system calls
% This is necessary because windows does not allow spaces in system
% calls. Thus, we have to cd into that folder and call the command
% from this folder.
%
% [status,result] = cat_system(cmd,verb,trerr)
% cmd            .. system call;
% verb           .. verbosity
% trerr          .. trough an error message (default), else just display error
% status, result .. system call outputs [status,result] = system('...');
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

if nargin == 0
  error('Argument is missing');
end
if nargin < 3 trerr = 1; end
if nargin < 2 verb = 0; end

CATDir = fullfile(fileparts(mfilename('fullpath')),'CAT');

% replace spaces in directory name
if ~ispc
  CATDir = strrep(CATDir,' ','\ ');
end

if ispc
  CATDir = [CATDir '.w32'];
elseif ismac
  [stat, output] = system('uname -v');
  % try to recognize new Apple arm64 processor
  if ~stat && ~isempty(strfind(output,'ARM64'))
    CATDir = [CATDir '.maca64'];
  else
    CATDir = [CATDir '.maci64'];
  end
elseif isunix
  CATDir = [CATDir '.glnx86'];
end  

% On Linux, MATLAB prepends its own library paths to LD_LIBRARY_PATH, which
% can cause CAT binaries to pick up MATLAB's older libstdc++.so.6 instead of
% the system one (missing GLIBCXX_3.4.29 and higher).  Unsetting
% LD_LIBRARY_PATH for the subprocess lets the OS use its ldconfig cache and
% find the correct system libraries, while leaving MATLAB's own environment
% untouched.
if isunix && ~ismac
  ldfix = 'env -u LD_LIBRARY_PATH ';
else
  ldfix = '';
end

if ispc
  olddir = pwd;
  cd(CATDir);
  [ST, RS] = system(cmd);
  cd(olddir);
else
  cmdfull = [ldfix fullfile(CATDir,cmd)];
  warning off % this is to prevent warnings by calling cat12 from the shell script
  [ST, RS] = system(cmdfull);
  warning on
end

% if not successful try it again after changing the file attributes to "+x"
if ST > 1
  if ispc
    try, fileattrib(fullfile(CATDir,'*'),'+x'); end
    olddir = pwd;
    cd(CATDir);
    [ST, RS] = system(cmd);
    cd(olddir);
  else
    try, fileattrib(fullfile(CATDir,'*'),'+x','a'); end
    cmdfull = [ldfix fullfile(CATDir,cmd)];
    warning off % this is to prevent warnings by calling cat12 from the shell script
    [ST, RS] = system(cmdfull);
    warning on
  end
end

if nargout > 0
  [status,result] = cat_check_system_output(ST,RS,verb,trerr);
else
  cat_check_system_output(ST,RS,verb,trerr);
end

% for mac we need to enable execution because of Apple Gatekeeper
if ismac && (ST == 137 || ST == 127)
  [fixStatus1, ~] = system(sprintf('xattr -dr com.apple.quarantine "%s"', strrep(CATDir,'\ ',' ')));
  [fixStatus2, ~] = system(sprintf('chmod -R a+x "%s"', strrep(CATDir,'\ ',' ')));
  if fixStatus1 == 0 && fixStatus2 == 0
    % retry the command after clearing quarantine
    warning off
    [ST, RS] = system(cmdfull);
    warning on
  end
end

if ST > 1 && ST~=139 % 139: data setup error
  if ispc
    [ST, RS] = system('systeminfo.exe');
  else
    [ST, RS] = system('uname -a');
  end
  str = sprintf('\nWARNING: Surface processing will not work because\n(1) File permissions are not correct (for unix use sudo chmod a+x) or\n(2) CAT binaries are not compatible to your system or\n(3) Antivirus software in Windows or Gatekeeper in macOS is blocking to execute binaries\nSystem: %s\n',RS);
  cat_io_cmd(str,'warning');
  helpdlg(str,'Error Using Surface Tools');
  
  % check Gatekeeper on macOS
  if ismac && ~isdeployed
    fprintf(2, '\n========================================================================\n');
    fprintf(2, 'CAT: Critical Permission Error\n');
    fprintf(2, 'macOS is blocking CAT binaries because they are quarantined.\n');
    fprintf(2, 'Please run this command in your Terminal to fix this:\n\n');
    fprintf(2, '     sudo xattr -dr com.apple.quarantine "%s"\n\n', strrep(CATDir,'\ ',' '));
    fprintf(2, 'Then restart MATLAB.\n');
    fprintf(2, '========================================================================\n\n');
  end
  fprintf('\n\nFor future support of your system please send this message to christian.gaser@uni-jena.de\n\n');
end

