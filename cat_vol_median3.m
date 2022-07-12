%cat_vol_median3 Median filter for label maps.
%  Median Filter for a 3D single image D. Bi is used to mask voxels for the 
%  filter process, whereas Bn is used to mask voxels that are used as 
%  neighbors in the filter process. Both mask can changed by intensity 
%  threshold (Bi_low,Bi_high,Bn_low,Bn_high) for D. Local NAN and INF
%  values can be replaced, if a neighbor has a non NAN/INF value and is
%  within the defined maskes and boundaries. 
% 
%   M = cat_vol_median3(D[, Bi, Bn, sf, Bi_low, Bi_high, Bn_low, Bn_high, 
%        filterNaNandINF])
% 
%   D      (single)  .. 3d matrix for filter process 
%   Bi     (logical) .. 3d matrix that marks voxels that should be filtered
%   Bn     (logical) .. 3d matrix that marks voxels that are used to filter
%   sf     (double)  .. threshold that is used to filter the result
%                          sf=0: no filter
%                          sf<0: only smaller changes
%                          sf>0: only bigger changes
%   Bi_low  (double) .. low  threshold in D for filtering (add to Bi)
%   Bi_high (double) .. high threshold in D for filtering (add to Bi)
%   Bn_low  (double) .. low  threshold in D for neighbors (add to Bn)
%   Bn_high (double) .. high threshold in D for neighbors (add to Bn)
%   filterNaNandINF (double ) .. replace NaN or Inf by the median of non
%                       NaN/INF voxels (default=0)
%
%  Used slower quicksort for median calculation, because the faster median 
%  of the median estimation leaded to incorrect results. 
% 
%  Example: 
%   A is the image that should be filter and that may contain NaN and Inf
%   values, whereas B defines the regions that should be filtered and spend
%   values. 
%
%     A = randn(50,50,3,'single');
%     B = false(size(A)); B(5:end-4,5:end-4,:)=true; 
%     N = rand(size(A),'single'); 
%     A(N>0.9 & N<1.0) = NaN; A(N<0.1 & N>0) = -inf; A(N<0.05 & N>0) = inf; 
%
%   1) simple cases without limits
%     C = cat_vol_median3(A,B); ds('d2smns','',1,A+B,C,2);
%
%   2) simple case without limits bud with NaN that are replaced by default
%     C = cat_vol_median3(A,B,B,0,-inf,inf,-inf,inf,1); ds('d2smns','',1,A+B,C,2);
%
%   3) Replace only small changes in C1, eg. to filter within tissue classes.
%      Replace only large outlier in C2, eg. to remove outlier like salt &
%      pepper noise. In both cases NANs/INFs were replaced.    
%     C1 = cat_vol_median3(A,B,B, -1.0 ,-inf,inf,-inf,inf,1 ); 
%     C2 = cat_vol_median3(A,B,B,  1.0 ,-inf,inf,-inf,inf,1 ); 
%     ds('d2smns','',1,C1,C2,2); 
%
%  See also cat_vol_median3c, compile.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
