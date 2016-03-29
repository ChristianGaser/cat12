function [ST, RS] = cat_system(varargin)
% ______________________________________________________________________
% CAT12 wrapper for system calls
% This is necessary because windows does not allow spaces in system
% calls. Thus, we have to cd into that folder and call the command
% from this folder.
% ______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 0
  error('Argument is missing');
end

CATDir      = fullfile(spm('dir'),'toolbox','cat12','CAT');
if ispc
  CATDir = [CATDir '.w32'];
elseif ismac
  CATDir = [CATDir '.maci64'];
elseif isunix
  CATDir = [CATDir '.glnx86'];
end  

olddir = pwd;
cd(CATDir);
[ST, RS] = system(varargin{1});
cd(olddir);

