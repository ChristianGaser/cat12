%cat_vol_median3c Median filter for label maps.
%  Median Filter for 3D single image D. Bi is used to mask voxels for the 
%  filter process, whereas Bn is used to mask voxels that are used as 
%  neighbors in the filter process. 
% 
%  This is a special subversion to filter label maps!
% 
%  M = cat_vol_median3c(D[,Bi,Bn])
% 
%  D  (single)  .. 3D matrix for filter process 
%  Bi (logical) .. 3D matrix that mark voxel that should be filtered
%  Bn (logical) .. 3D matrix that mark voxel that are used as neighbors 
% 
%  Examples: 
%   1)
%     A = round(smooth3(rand(50,50,3,'single')*3));
%     B = false(size(A)); B(5:end-4,5:end-4,:)=true; 
%     C = cat_vol_median3c(A,B); C = cat_vol_median3c(C,B); 
%     ds('d2smns','',1,A+B,C,2);
%
%  See also cat_vol_median3, compile.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id: 2558 2024-02-28 $
