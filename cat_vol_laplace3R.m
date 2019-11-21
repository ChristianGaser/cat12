%cat_vol_laplace3R Volumetric Laplace filter with Dirichlet boundary. 
%  Filter SEG within the intensity range of low and high until the changes
%  are below TH. 
% 
%  L = cat_vol_laplace3(SEG,R,TH)
%
%  SEG  .. 3D single input matrix
%  R    .. 3D boolean volume to describe the filter area
%  TH   .. threshold to controll the number of interations
%          maximum change of an element after interation
%
%  See also cat_vol_laplace3, compile.
%  ________________________________________________________________________
%  Robert Dahnke 2009/02
%  $Id: cat_vol_laplace3.m 1519 2019-11-19 10:48:29Z gaser $
