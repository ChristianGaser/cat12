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

  if status==1 || ...
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
