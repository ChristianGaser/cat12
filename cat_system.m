function [status,result] = cat_system(cmd,verb,trerr)
% ______________________________________________________________________
% CAT12 wrapper for system calls
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
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
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

CATDir = fullfile(spm('dir'),'toolbox','cat12','CAT');

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

if ispc
  olddir = pwd;
  cd(CATDir);
  [ST, RS] = system(cmd);
  cd(olddir);
else
  cmd = fullfile(CATDir,cmd);
  warning off % this is to prevent warnings by calling cat12 from the shell script
  [ST, RS] = system(cmd);
  warning on
end

if nargout > 0
  [status,result] = cat_check_system_output(ST,RS,verb,trerr);
else
  cat_check_system_output(ST,RS,verb,trerr);
end

% for mac we need to enable execution because of Apple Gatekeeper
if ismac && ST == 137
  web('https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)#Troubleshooting');
  CATDir = fullfile(spm('dir'),'toolbox','cat12','CAT');
  cat_io_cmd(sprintf('\nThe following commands will be executed as administrator to allow execution of CAT12 binaries and mex-files.\n Please now type admin password to call sudo\n'),'warning');
  cat_io_cmd(sprintf('You can also break that command here and run the commands that are listed on the open website under Troubleshooting manually.\n'),'warning');
  cmd = ['sudo xattr -r -d com.apple.quarantine ' CATDir];
  system(cmd); fprintf([cmd '\n']);
  cmd = ['sudo find ' CATDir ' -name *.mexmac* -exec spctl --add {} \;'];
  system(cmd); fprintf([cmd '\n']);
  cmd = ['sudo chmod a+x ' CATDir '/CAT.mac*/CAT*'];
  system(cmd); fprintf([cmd '\n']);
  cmd = ['sudo find ' CATDir ' -name *.mexmac* -exec xattr -d com.apple.quarantine {} \;'];
  system(cmd); fprintf([cmd '\n']);
  ST = system(fullfile(CATDir,'CAT_3dVol2Surf'));
end

if ST > 1
  if ispc
    [ST, RS] = system('systeminfo.exe');
  else
    [ST, RS] = system('uname -a');
  end
  str = sprintf('\nWARNING: Surface processing will not work because\n(1) CAT binaries are not compatible to your system or\n(2) File permissions are not correct (for unix use chmod a+x) or\n(3) Antivirus software in Windwos or Gatekeeper in MAC OS is blocking to execute binaries\nSystem: %s\n',RS);
  cat_io_cmd(str,'warning');
  helpdlg(str,'Error Using Surface Tools');
  
  % check Gatekeeper on MAC OS
  if ismac
    [ST, RS] = system('spctl --status');
    if ~isempty(strfind(RS,'enabled'))
      str = 'Please disable Gatekeeper on MAC OS!';
      fprintf('\n\n%s\n',str);
      helpdlg(str,'Gatekeeper error');
      web('https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)#Troubleshooting');
    end
  end
  fprintf('\n\nFor future support of your system please send this message to christian.gaser@uni-jena.de\n\n');
end

