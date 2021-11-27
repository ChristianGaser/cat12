%cat_vol_eidist Eikonal distance.
%  Estimates eikonal distance to an object (negative values) with the speed 
%  given by the positive values
% 
%   D = cat_vol_eidist(O,F,vx_vol) 
% 
%   D (3D single)       .. distance map
%   O (3D single)       .. object map
%   F (3D single)       .. speed map
%   vx_vol (1x3 double) .. voxel size
%
%  Example: 
%   1) distance from two points with a simple mask 
%    A = zeros(50,50,3,'single'); A(15,25,2)=-1; A(35,25,2)=-1; 
%    B = nan(size(A),'single'); B(5:end-5,5:end-5,:) = rand(size(A)-[9 9 0]);
%    D = cat_vol_eikonal3(A,B);  ds('d2smns','',1,A+B,D,2);
%
%   2) not working mask definition
%    B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
%    B(10:end-10,10:end-10,:) = false; 
%    D = cat_vol_eikonal3(A,B); ds('d2smns','',1,A+B,D,2);
%
%  See also cat_vol_eidist, compile.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
