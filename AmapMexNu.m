function [prob, mean] = AmapMexNu(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize)
% FORMAT [prob, mean] = AmapMexNu(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize)
%
% Christian Gaser
% $Id: AmapMexNu.m 215 2009-11-15 23:01:31Z gaser $

rev = '$Rev: 215 $';

disp('Compiling AmapMexNu.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
ext = mexext;
eval(['mex -O AmapMexNu.c Kmeans.c Amap.c MrfPrior.c Pve.c SplineSmooth.cc -I. ' fullfile('.',ext(4:end),'libEBTKS.a')])
cd(p_path);

try 
[prob, mean] = AmapMexNu(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize);
end

return
