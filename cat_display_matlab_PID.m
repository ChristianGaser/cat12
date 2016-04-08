function varagout = cat_display_matlab_PID
% cat_display_matlab_PID (in development!)
% ______________________________________________________________________
%
% Display exor return PID of this (Linux/Mac) or the last started 
% (Windows) MATLAB instace.
% ______________________________________________________________________
%
%   Robert Dahnke - robert.dahnke@uni-jena.de
%   Center of Neuroimaging 
%   Department of Psychiatry and Psychotherapy 
%   University Hostpital Jena
% ______________________________________________________________________
% $Id$

  % get PID
  if ispc
    pid = []; pidi = 0; 
    while isempty(pid) && pidi<100 
      s = pidi; 
      [t,pids] = system(sprintf(['tasklist /FO TABLE /NH /FI "imagename eq MATLAB.exe" ' ...
        '/FI "CPUTime lt 00:%02d:%02d" /FI "MemUsage lt 500000"'],round(s/60),mod(s,60)));
      pid = textscan(pids,'MATLAB.exe %d Console'); 
      if ~isempty(pid); pid = pid{1}; end
      pidi = pidi+1;
    end
  else
    [t,pid] = system('echo $$'); 
  end
  pids = pid; 
  pid = str2double(pid);
  
  % display PID
  if nargout==0
    if isnumeric(pid) && ~isempty(pid)
      fprintf('CAT parallel processing with MATLAB PID: %d\n',pid);
    else
      fprintf('CAT parallel processing with MATLAB PID: unknown %s PID %d\n',computer,pids);
    end
  else 
    if isnumeric(pid) && ~isempty(pid)
      varagout{1} = pid; 
    else
      varagout{1} = pids; 
    end  
  end
  
end
