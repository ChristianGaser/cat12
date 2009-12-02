function [prob, mean] = AmapMex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize)
% FORMAT [prob, mean] = AmapMex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling AmapMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O AmapMex.c KmeansProper.c Amap.c MrfPrior.c Pve.c 
cd(p_path);

try 
[prob, mean] = AmapMex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize);
end

return
