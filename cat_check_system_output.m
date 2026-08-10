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
  
  % Detect error messages in the raw output, ie. before genstrarray replaces
  % the line breaks and quotes the string.
  % The CAT binaries are not consistent in the notation of their error
  % messages and most of them use "Error: ..." rather than "ERROR: ...",
  % eg. CAT_VolMarchingCubes prints "Error: Label image must have the same
  % dimensions as the input image." and exits with a status of 1, that is
  % not counted as failure here. Such errors passed the case-sensitive test
  % unnoticed and the calling function silently continued with the (old or
  % missing) output file. Hence, the test is case-insensitive but anchored
  % at the begin of a line to avoid false positives by messages that only
  % mention an error value, eg. "root mean square error" or "errors : %d"
  % (the lookahead excludes the plural and is used rather than a word
  % boundary because MATLAB does not support \b and Octave does not \>).
  haserror = status > 1 || ...
     ~isempty(regexpi(result,'(^|[\r\n])\s*error(?![a-zA-Z])','once')) || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'));

  % replace special characters
  result = genstrarray(result);

  if haserror
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