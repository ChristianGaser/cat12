function [FN,C] = cat_io_struct2cell(S,separator)
% cat_io_struct2cell. Recursive use of struct2cell with fieldnames.
%
%  [FN,C] = cat_io_struct2cell(S,separator)
% 
%  FN        .. list of full fieldname 
%  C         .. data entry 
%  S         .. struct
%  separator .. seperator of fields in FN (default = '.')
%
% Example: 
%  S      = struct('a',1,'b',struct('b1',2,'b2',3))
%  [FN,C] = cat_io_struct2cell(S)
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

  if ~exist('separator','var')
    separator = '.';
  end

  if numel(S) > 1
    error( ...
      sprintf('%s:nelment',mfilename),...
      'Only 1 element per (sub)struct supportet yet!');
  end

  FN = fieldnames(S); 
  C  = struct2cell(S); 
  for ci = numel(C):-1:1
    if isstruct(C{ci})
      % recursive call
      [FNi,Ci] = cat_io_struct2cell(C{ci}); 
      
      % extend fieldnames
      for FNii = 1:numel(FNi)
        FNi{FNii} = sprintf('%s%s%s',FN{ci},separator,FNi{FNii});
      end

      % update fieldnames and cell
      FN = [FN(1:ci-1); FNi; FN(ci+1:end)];
      C  = [C(1:ci-1);  Ci;  C(ci+1:end)]; 
    end
  end

end