function [prob, mean] = PveAmapMex(src, label, mask, vx)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling PveAmapMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O PveAmapMex.c PveAmap.c Amap.c MrfPrior.c Pve5.c Kmeans.c optimizer3d.c diffeo3d.c SplineSmooth.cc -lEBTKS
cd(p_path);

[prob, mean] = PveAmapMex(src, label, mask, vx);

return
