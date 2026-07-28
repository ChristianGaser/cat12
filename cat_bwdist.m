%cat_bwdist Euclidean distance to a binary object, distance only.
% Toolbox free replacement for bwdist of the Image Processing Toolbox.
% Calculates the euclidean distance to an object in P, where the object is
% given by all voxels with P>0.5. In contrast to cat_vbdist it returns only
% the distance and no index or label map, which allows a separable algorithm
% that is about ten times faster and exact.
%
%  D = cat_bwdist(P[,vx_vol])
%
%  P      (single) input image, the object is P>0.5
%  vx_vol (double) 1x3 voxel size in mm, default [1 1 1]
%  D      (single) distance image in mm
%
% Use this function whenever only the distance is needed, e.g. to threshold
% it for a morphological operation. Use cat_vbdist if the index of the
% nearest object voxel (I), the label map (L), or the mask argument (M) that
% restricts the distance calculation are required.
%
% Differences to bwdist of the Image Processing Toolbox:
%   - the object is P>0.5 and not every nonzero voxel, which follows the
%     convention of cat_vbdist and makes this a drop in replacement for it
%   - only the distance is returned, there is no second output IDX
%   - only the euclidean metric, no cityblock/chessboard/quasi-euclidean
%   - the voxel size can be given, so the distance is in mm and not in
%     voxels, and anisotropic voxels are handled correctly
%
% The algorithm is the separable distance transform of Felzenszwalb &
% Huttenlocher (2012): the squared distance is obtained with three one
% dimensional passes, one per dimension, each in linear time. cat_vbdist
% instead propagates the vector to the nearest object voxel, which is exact
% for isotropic voxels but approximate for anisotropic ones and accumulates
% a small error with increasing distance to large objects.
%
% Examples:
%  % (1) distance from two points
%   A = zeros(50,50,3,'single'); A(15,25,2)=1; A(35,25,2)=1;
%   D = cat_bwdist(A);
%
%  % (2) anisotropic voxel size, distance in mm
%   D = cat_bwdist(A,[0.5 0.5 1.5]);
%
%  % (3) dilation of a binary mask by 5mm
%   Ydil = cat_bwdist(single(Ymask),vx_vol) <= 5;
%
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$

%-This is merely the help file for the compiled routine
error('cat_bwdist.c not compiled - see compile.m')
