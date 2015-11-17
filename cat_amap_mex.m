function [prob, mean] = cat_amap_mex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
% FORMAT [prob, mean] = cat_amap_mex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling cat_amap_mex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cat_amap_mex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
cd(p_path);

[prob, mean] = cat_amap_mex(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm);

return
