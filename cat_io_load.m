function x = cat_io_load(file)
% ______________________________________________________________________
% Read values from txt-files and replace empty lines/entries by NaNs 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~nargin
    [file,sts] = spm_select(1,'^.*\.txt$');
    if ~sts, x = []; return; end
end

if ~exist(file,'file')
    error('Unable to find file ''%s''',file);
end

try
  str = fileread(file);
  str = regexprep(str, '(?<=\r?\n)[ ]*(?=\r?\n)', 'NaN', 'emptymatch' );
  x = str2num(str);
catch
  fprintf('Error: Only txt-files without characters are allowed.\n');
  x = [];
end
