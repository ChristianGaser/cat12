function vbm_check_system_output(status,result,debugON)
%vbm_check_system_output check of system commands and returned result 

%_______________________________________________________________________
% Christian Gaser
% $Id$


  if status==1 || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'))
    error('VBM:system_error',result); 
  end
  if nargin > 2
    if debugON, disp(result); end
  end
end
