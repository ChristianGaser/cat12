%cat_vol_localstat Local mean, minimum, maximum, SD, and peak estimation.
* Estimates specific functions in a volume V within a mask region M. 
% For each voxel v of M, the values of the neigbors of v that belong to M 
% and are within a distance smaller than nb where used (i.e. masked voxels
% within a sphere of radius nb). 
%
% If V contains NaNs, -INFs or INFs these values are ignored and added to 
% the mask M.  Masked voxels in the output volume were defined depending 
% on the mskval  variable, i.e. as zeros (0, default), input (1), NANs (2),
% -INF (3), or INF(4).  
%
% The function was designed to extract the inhomogeneity in noisy data in 
% a well-known area, such a tissue class in structural images. In general  
% the mean estimated in a local neighborhood nb or after some interations
% iter. However, tissues are often affected by partial volume effects near
% the tissue boundaries and minimum/maximum can be quite helpful to reduce
% such effects. For noise estimation the local variance/standard deviation 
% is also quite useful. 
% Besides the mean value the local peak of the histogram can also work 
%
%  S = cat_vol_localstat(V,M[,nb,stat,iter,filter0,verb])
% 
%  V    (single)    input volume
%  M    (logical)   mask volume
%  nb   (double)    neigbhour distance (1 .. 10)
%  stat (double)    1-mean, 2-min, 3-max, 4-std  
%                   5-peak1, 6-peak2, 7-peak3    (experimental)
%                   8-median                     
%                   9-hist                       (experimental)
% iter             number of iterations          (default=1)
% filter0 (double) originally values <=0 were ignored (default=0)
% mskval  (double) setting of masked voxels
%                   (0-zeros,1-input,2-NAN,3--INF,4-INF)
% verb (double)    verbose output for debugging
%
%
% Examples: 
%  Here are some simple samples to outline the subfunctions. The mask area
%  is defined by NaN. The simulated data of A is between -1 and 1 and B is 
%  a locial mask.
%
%  == input variables ==
%  A  = rand(20,20,3,'single') - 1;           
%    for i=1:size(A,2), A(:,i,:) = A(:,i,:) + (( i/size(A,2) ) - 0.5); end
%  B  = smooth3(smooth3(rand(size(A))))>0.5;
%
%   
%  == function calls ==
%  (1) MEAN:     values around 0
%      C = cat_vol_localstat(A,B,2,1,2,2); ds('d2smns','',1,A,C,2); 
%
%  (2) MINIMUM:  values trending torwards -1
%      C = cat_vol_localstat(A,B,2,2,2,2); ds('d2smns','',1,A,C,2); 
%
%  (3) MAXIMUM:  values trending torwards 1
%      C = cat_vol_localstat(A,B,2,3,2,2); ds('d2smns','',1,A,C,2); 
%
%  (4) STANDARD DEVIATION:  values about 0.5 
%      C = cat_vol_localstat(A,B,2,4,1,2); ds('d2smns','',1,A,C,2); 
%
%  (8) MEDIAN:   values around 0
%      C = cat_vol_localstat(AL,B,2,8,1,2); ds('d2smns','',1,A,C,2); 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%  $Id$
