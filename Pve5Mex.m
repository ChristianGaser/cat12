function Pve5Mex(src, prob, label, mean, BG);
%
% Christian Gaser
% $Id: Pve5Mex.m 8 2008-12-10 20:21:27Z gaser $

rev = '$Rev:$';

disp('Compiling Pve5Mex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O Pve5Mex.c Pve5.c 
cd(p_path);

Pve5Mex(src, prob, label, mean, BG);

return
