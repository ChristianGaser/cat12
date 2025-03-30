function varargout = cat_check_system_output(status,result,debugON,trerr)
%_______________________________________________________________________
% cat_check_system_output check of system commands and returned result 
%
% cat_check_system_output(status,result,debugON,trerr)
%
% status, result .. system call outputs [status,result] = system('...');
% debugON        .. display result
% trerr          .. trough an error message (default), else just display 
%                   error
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('debugON','var'), debugON=0; end
  if ~exist('trerr','var'), trerr=1; end
  if nargout>0, varargout{1} = false; varargout{2} = result; end 
  
  % replace special characters
  result = genstrarray(result);
  
  if status > 1 || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'))
    if nargout>0, varargout{1} = true; varargout{2} = result; end
    if trerr
      try
        error('CAT:system_error',sprintf('CAT System_error: %s',result)); 
      catch
        fprintf('CAT System_error: %s',sprintf(result)); 
      end
    else
      cat_io_cprintf('warn','CAT:system_error:%s',sprintf(result)); 
    end
  end
  if nargin > 2
    if debugON>0 && ~strcmp(result,'')
      fprintf('%s',sprintf(result)); 
    end
  end
end

function str = genstrarray(stritem)
% generate a string of properly quoted strings 

  str = strrep(stritem, '''', '''''');
  if ~any(str  == char(0)) &&  ~any(str  == char(9)) && ~any(str  == char(10)) && ~strcmp(str,'')
    str  = sprintf('''%s''', str ); 
  else
    % first, quote sprintf special chars % and \
    % second, replace special characters by sprintf equivalents
    replacements = {'%', '%%'; ...
        '\', '\\'; ...
        char(0), '\0'; ...
        char(9), '\t'; ...
        char(10), '\n';...
        '\S', ''};
    for cr = 1:size(replacements, 1)
        str  = strrep(str , replacements{cr,:});
    end
  end
end