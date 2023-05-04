tfce = tfceMex_pthread(t, dh, E, H, calc_neg, single_threaded)
% Apply threshold-free cluster enhancement (TFCE)
% FORMAT tfce = tfceMex(t, dh, E, H, calc_neg, single_threaded)
% Estimate TFCE
% t         - T map 
% dh        - step size (e.g. dh = max(abs(t))/100)
% E         - TFCE parameter for extent
% H         - TFCE parameter for height
% calc_neg  - also calc neg. TFCE values (default)
% single_threaded - use single thread only
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev: 125 $';

disp('Compiling tfceMex_pthread.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);

if strcmpi(spm_check_version,'octave')
  mexcmd = 'mkoctfile --mex';
  mexflag=' -O -DOCTAVE ';
else
  mexcmd = 'mex';
  mexflag=' -O -largeArrayDims COPTIMFLAGS=''-O3 -fwrapv -DNDEBUG'' CFLAGS=''$CFLAGS -pthread -Wall -ansi -pedantic -Wextra'' ';
end

eval([mexcmd ' ' mexflag ' tfceMex_pthread.c'])

cd(p_path);

tfce = tfceMex_pthread(t, dh, E, H, calc_neg, single_threaded);

return
