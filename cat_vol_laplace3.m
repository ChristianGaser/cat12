%cat_vol_laplace3 Volumetric Laplace filter with Dirichlet boundary. 
%  Filter SEG within the intensity range of low and high until the changes
%  are below TH. 
% 
%  L = cat_vol_laplace3(SEG,low,high,TH)
%
%  SEG  .. 3D single input matrix
%  low  .. low boundary threshold
%  high .. high boundary threshold
%  TH   .. threshold to controll the number of interations
%          maximum change of an element after interation
%
%  See also cat_vol_laplace3R, compile.
%  ________________________________________________________________________
%  Robert Dahnke 2009/01
%  $Id: cat_vol_laplace3.m 1519 2019-11-19 10:48:29Z gaser $
