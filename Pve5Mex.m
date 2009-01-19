function Pve5Mex(src, prob, label, mean, BG);
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling Pve5Mex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O Pve5Mex.c Pve5.c 
cd(p_path);

Pve5Mex(src, prob, label, mean);

return
