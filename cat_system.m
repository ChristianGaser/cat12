function [ST, RS] = cat_system(varargin)
% ______________________________________________________________________
% CAT12 wrapper for system calls
% This is necessary because windows does not allow spaces in system
% calls. Thus, we have to cd into that folder and call the command
% from this folder.
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

CATDir = fullfile(spm('dir'),'toolbox','cat12','CAT');

% replace spaces in directory name
if ~ispc
  CATDir      = strrep(CATDir,' ','\ ');
end

if ispc
  CATDir = [CATDir '.w32'];
elseif ismac
  CATDir = [CATDir '.maci64'];
elseif isunix
  CATDir = [CATDir '.glnx86'];
end  

if ispc
  olddir = pwd;
  cd(CATDir);
  [ST, RS] = system(varargin{1});
  cd(olddir);
else
  cmd = fullfile(CATDir,varargin{1});
  warning off % this is to prevent warnings by calling cat12 from the shell script
  [ST, RS] = system(cmd);
  warning on
end
