function T = cat_io_structEqual(S1,S2)
% _________________________________________________________________________
% Test if two structures are equal (independent of the order of the fields).
% Two structures are equal if they have the same size and fields and if
% all fields are equal:
%  - numeric and logical values by isequal (i.e., NaN is not equal to NaN)
%  - char by strcmp and strings by their text
%  - structures recursively (also arrays of structures)
%  - cellstr by their joined text (i.e., independent of the cell shape)
%  - other cells element-wise (with the same size)
% Values of different types are not equal (e.g., 1 and true or 1 and '1').
%
%   T = cat_io_structEqual(S1,S2)
%
%   S1, S2 .. structures
%   T      .. true if the structures are equal (false for non-structures)
%
% See also isequal, cat_io_updateStruct.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  T = isstruct(S1) && isstruct(S2) && isequal(size(S1),size(S2));
  if ~T, return; end

  FN1 = sort(fieldnames(S1));
  FN2 = sort(fieldnames(S2));
  if numel(FN1) ~= numel(FN2) || ~all(strcmp(FN1,FN2)), T = false; return; end

  for si = 1:numel(S1)
    for fni = 1:numel(FN1)
      T = valueEqual(S1(si).(FN1{fni}), S2(si).(FN1{fni}));
      if ~T, return; end
    end
  end
end
% =========================================================================
function T = valueEqual(v1,v2)
%valueEqual. Test if two values are equal (see cat_io_structEqual).
  if (isnumeric(v1) && isnumeric(v2)) || (islogical(v1) && islogical(v2))
    T = isequal(v1,v2);
  elseif ischar(v1) && ischar(v2)
    T = strcmp(v1,v2);
  elseif (ischar(v1) || isstring(v1)) && (ischar(v2) || isstring(v2))
    T = isequal(string(v1),string(v2));
  elseif isstruct(v1) && isstruct(v2)
    T = cat_io_structEqual(v1,v2);
  elseif iscellstr(v1) && iscellstr(v2)
    T = strcmp(char(join(v1)),char(join(v2)));
  elseif iscell(v1) && iscell(v2)
    T = isequal(size(v1),size(v2));
    for ci = 1:numel(v1)
      if ~T, return; end
      T = valueEqual(v1{ci},v2{ci});
    end
  else
    T = false;
  end
end
