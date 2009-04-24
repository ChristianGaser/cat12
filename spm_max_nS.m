function [N,Z,M,A,R] = spm_max_nS(X,L,VRPV)
%
% function spm_max_nS.m
% Finds sizes, maxima, locations, and resel counts for local excursion sets.
% This is an extension of spm_max function, with an additional functionality
% to calculate cluster sizes in terms of resels. 
%
% FORMAT [N Z M A R] = spm_max_nS(X,L,VRPV)
%  X     - values of 3-D field
%  L     - locations [x y x]' {in voxels}
%  VRPV  - pointer for the RPV image
%  N     - size of region {in voxels)
%  Z     - Z values of maxima
%  M     - location of maxima {in voxels}
%  A     - region number
%  R     - resel count for the region
%________________________________________________________________________________
%
% spm_max_nS is the non-stationary version of spm_max. This function calls
% spm_max to obtain cluster sizes, locations, and local maxima. This function
% also calculates cluster sizes in terms of resels by summing RPV 
% (resel-per-voxel) values in the RPV image provided by the user. If the RPV
% value is not available for all the voxels in a cluster, then the RPV value
% is imputed by the average RPV value of that cluster.
%
%________________________________________________________________________________
% spm_max_nS.m  0.75b
% Written by Satoru Hayasaka   October 23, 2006
%

%-Running the spm_max
[N Z M A] = spm_max(X,L);


%-Reading in the RPV image
VRPV = spm_vol(VRPV.fname);
tmpx = L(1,:);
tmpy = L(2,:);
tmpz = L(3,:);
XRPV = spm_sample_vol(VRPV,tmpx,tmpy,tmpz,1);


%-First, calculate average RPV for all the ST voxels
dAll = find(~isnan(XRPV));
if dAll
     mRPV = mean(XRPV(dAll));
else
     mRPV = 0;
end


%-Clustering
A2     = spm_clusters(L);


%-Checking number of clusters
numCl  = max(A);
numCl2 = max(A2);
if numCl ~= numCl2
     error('Number of clusters does not much in spm_max and spm_clusters!');
end


%-Calculating resels for clusters
R      = zeros(size(A));

for icl = 1:max(A2)
    d     = find(A2==icl);
    dXYZ  = L(:,d);
    dRPV  = XRPV(d);
    isRPV = find(~isnan(dRPV));
    if length(isRPV) == 0
       clresel = mRPV*length(d);
       %clresel = length(d);
    else
       clresel = (sum(dRPV(isRPV))/length(isRPV))*length(d);
       %clresel = length(d);
    end

    %-finding the corresponding cluster from the spm_max results
    jcl   = 0;
    bcl   = 0;
    while bcl == 0
        jcl  = jcl + 1;
        jvox = min(find(A==jcl));
        jXYZ = M(:,jvox);
        bcl  = ismember(jXYZ',dXYZ','rows');
    end

    jclind    = find(A==jcl);
    R(jclind) = clresel;
end;

return;

