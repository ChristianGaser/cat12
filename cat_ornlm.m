function out = cat_ornlm(in, v, f, h)
% FORMAT out = cat_ornlm(in, v, f, h)
% 
% Optimized Blockwise Non Local Means Denoising Filter
%
% v - size of search volume (M in paper)
% f - size of neighborhood (d in paper)
% h - smoothing parameter
%
%                          Details on ONLM filter                        
% ***************************************************************************
%  The ONLM filter is described in:                                       
%                                                                         
%  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     
%  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic
%  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, 
%  April 2008                                                             
% ***************************************************************************
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

rev = '$Rev$';

disp('Compiling cat_ornlm.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cat_ornlm.c ornlm_float.c 
cd(p_path);

out = cat_ornlm(in, v, f, h);

return
