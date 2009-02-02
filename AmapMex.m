function [prob, mean] = AmapMex(src, label, nc, niters, sub, pve)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling AmapMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O AmapMex.c Amap.c MrfPrior.c Pve6.c
cd(p_path);

[prob, mean] = AmapMex(src, label, nc, niters, sub, pve);

return
