function [prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
% FORMAT [prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

disp('Compiling cat_amap.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cat_amap.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
cd(p_path);

[prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm);

return
