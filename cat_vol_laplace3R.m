%cat_vol_laplace3R Volumetric Laplace filter with Dirichlet boundary. 
%  Filter SEG within the intensity range of low and high until the changes
%  are below TH. 
% 
%  L = cat_vol_laplace3R(SEG,R,TH)
%
%  SEG  .. 3D single input matrix
%  R    .. 3D boolean volume to describe the filter area
%  TH   .. threshold to control the number of iterations
%          maximum change of an element after iteration
%
% Example: 
%   A = zeros(50,50,3,'single'); A(10:end-9,10:end-9,2)=0.5; 
%   A(20:end-19,20:end-19,2)=1;
%   B = A==0.5; 
%   C = cat_vol_laplace3R(A,B,0.001); ds('d2smns','',1,A,C,2); 
%
%  See also cat_vol_laplace3, compile.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
