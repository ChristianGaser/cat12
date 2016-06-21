function varargout = cat_check_system_output(status,result,debugON,trerr)
%_______________________________________________________________________
% cat_check_system_output check of system commands and returned result 
%
% cat_check_system_output(status,result,debugON,trerr)
%
% status, result .. system call outputs [status,result] = system('...');
% debugON        .. dipslay result
% trerr          .. trough an error message (default), else just display 
%                   error
%_______________________________________________________________________
% Christian Gaser
% $Id$

  if ~exist('debugON','var'), debugON=0; end
  if ~exist('trerr','var'), trerr=1; end
  if nargout>0, varargout{1} = false; end 
  
  % replace special characters
  result = genstrarray(result);
  
  if status || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'))
   if nargout>0, varargout{1} = true; end
    if trerr
      error('CAT:system_error',result); 
    else
      cat_io_cprintf('err','CAT:system_error\n');
    end
  end
  if nargin > 2
    if debugON, disp(result); end
  end
end

function str = genstrarray(stritem)
% generate a string of properly quoted strings 

  str = strrep(stritem, '''', '''''');
  if ~any(str  == char(0)) &&  ~any(str  == char(9)) && ~any(str  == char(10))
    str  = sprintf('''%s''', str );
  else
    % first, quote sprintf special chars % and \
    % second, replace special characters by sprintf equivalents
    replacements = {'%', '%%'; ...
        '\', '\\'; ...
        char(0), '\0'; ...
        char(9), '\t'; ...
        char(10), '\n'};
    for cr = 1:size(replacements, 1)
        str  = strrep(str , replacements{cr,:});
    end
  end
end