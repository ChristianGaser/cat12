function [ST, RS] = cat_system(varargin)
% ______________________________________________________________________
% CAT12 wrapper for system calls
% This is necessary because windows does not allow spaces in system
% calls. Thus, we have to cd into that folder and call the command
% from this folder.
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_system.m 892 2016-03-09 17:57:21Z dahnke $

rev = '$Rev: 892 $';

if nargin == 0
  error('Argument is missing');
end

if ispc
    olddir = pwd;
    [pth, nam] = fileparts(varargin{1}); 
    cd(pth);
    [ST, RS] = system(nam);
    cd(olddir);
else
    [ST, RS] = system(varargin{1});
end
