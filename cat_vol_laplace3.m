%cat_vol_laplace3 Volumetric Laplace filter with Dirichlet boundary. 
%  Filter SEG within the intensity range of low and high until the changes
%  are below TH. 
% 
%  L = cat_vol_laplace3(SEG,low,high,TH)
%
%  SEG  .. 3D single input matrix
%  low  .. low boundary threshold
%  high .. high boundary threshold
%  TH   .. threshold to control the number of iterations
%          maximum change of an element after iteration
%
% Example: 
%   A = zeros(50,50,3,'single'); A(10:end-9,10:end-9,2)=0.5; 
%   A(20:end-19,20:end-19,2)=1;
%   C = cat_vol_laplace3(A,0,1,0.001); ds('d2smns','',1,A,C,2); 
% 
%  See also cat_vol_laplace3R, compile.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
