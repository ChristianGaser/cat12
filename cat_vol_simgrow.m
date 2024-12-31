%cat_vol_simgrow Volumetric region-growing.
%  
%  [SLAB,DIST] = cat_vol_simgrow(ALAB,SEG[,d,dims,dd]);
%  
%  SLAB (3D  single) .. output label map  
%  DIST (3D  single) .. distance map from region-growing
%  ALAB (3D  single) .. input label map
%  SEG  (3D  single) .. input tissue map
%  d    (1x1 double) .. growing treshhold parameter (max local gradient)
%                       in SEG
%  dims (1x3 double) .. voxel dimensions (default [1,1,1])
%  dd   (1x2 double) .. general growing limits in SEG 
%
% Examples:
%  1) 
%    A = zeros(50,50,3,'single'); A(:,1:25,:)=0.25; A(:,25:end,:)=0.75; 
%    A = A + (rand(size(A),'single')-0.5)*0.05; 
%    B = zeros(50,50,3,'single'); B(15:35,15:20,:)=1; B(15:35,30:35,:)=2; 
%    [C,D] = cat_vol_simgrow(B,A,1); ds('d2smns','',1,A,C,2);
%
%  See also cat_vol_downcut.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
