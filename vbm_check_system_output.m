function vbm_check_system_output(status,result,debugON)
%vbm_check_system_output check of system commands and returned result 

%_______________________________________________________________________
% Christian Gaser
% $Id: vbm_check_system_output.m 584 2014-03-12 16:59:22Z gaser $


  if status==1 || ...
     ~isempty(strfind(result,'ERROR')) || ...
     ~isempty(strfind(result,'Segmentation fault'))
    error('VBM:system_error',result); 
  end
  if nargin > 2
    if debugON, disp(result); end
  end
end
