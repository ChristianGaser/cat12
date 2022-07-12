%cat_vbdist Voxel-based euclidean distance calculation.
% Calculates the euclidean distance without PVE to an object in P with a 
% boundary of 0.5 for all voxels within a given mask M that should define
% a convex hull with direct connection between object and estimation 
% voxels.
% 
%  [D,I,L] = vbdist(P[,M])
%  
%  P (single)  input image with zero for non elements 
%  M (logical) mask to limit the distance calculation roughly, e.g., to 
%              the brain defined by a convex hull
%              WARNING: Voxels have to see objects within the mask!
%  D (single)  distance image
%  L (uint8)   label map
%  I (uint32)  index of nearest point
%
% Examples: 
%  % (1) distance from two points with a simple mask 
%   A = zeros(50,50,3,'single'); A(15,25,2)=1; A(35,25,2)=1; 
%   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
%   D = cat_vbdist(A,B); 
%   ds('d2sm','',1,A/3+B/3+1/3,D/20,2);
%
%  % (2) not working mask definition
%   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
%    B(10:end-10,10:end-10,:) = false; 
%   D = cat_vbdist(A,B); 
%   ds('d2sm','',1,A/3+B/3+1/3,D/20,2);
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
