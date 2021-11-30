function [prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
% FORMAT [prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm)
%
%  WARNING: The function changes the input data src and label!
%
%  Examples: 
%   % create a labelmap with 3 large areas
%    L = zeros(50,50,9,'single'); 
%    L(round(size(L,1)*0.27), round(size(L,2)*0.21), round(size(L,3)/2)+1) = 1; 
%    L(round(size(L,1)*0.85), round(size(L,2)*0.45), round(size(L,3)/2)+1) = 2; 
%    L(round(size(L,1)*0.37), round(size(L,2)*0.76), round(size(L,3)/2)+1) = 3; 
%    B = false(size(L)); B(5:end-4,5:end-4,:) = true; 
%    [D,I] = cat_vbdist(L,B); L=L(I); L = uint8(round(L));
% 
%   % simulate real image with noise and bias
%    A = double(L);  
%    for i=1:size(A,2), A(:,i,:) = A(:,i,:) .* (( i/size(A,2) ) + 0.5) / 3; end
%    A = A + 0.1 * rand(size(L)); 
%    ds('d2smns','',1,A,L,round(size(L,3)/2)+1);
%
%   % use AMAP to segment the simulated image
%    Atmp = A + 0; % create copy!
%    [prob, mean] = cat_amap(Atmp,L, 3,16,32, 5,0, 0.1,ones(1,3),50,60);
%    S = single(prob(:,:,:,1))*1/3/255 + single(prob(:,:,:,2))*2/3/255 + single(prob(:,:,:,3))/255;
%    ds('d2sm','',1,A,S,round(size(L,3)/2)+1);
%
%   % another test
%    B = double(L)/3 + 0.03 * randn(size(L)); Btmp = B + 0; 
%    [prob, mean] = cat_amap(Btmp,L, 3,16,32, 5,0, 0.1,ones(1,3),50,60);
%    S = single(prob(:,:,:,1))*1/3/255 + single(prob(:,:,:,2))*2/3/255 + single(prob(:,:,:,3))/255;
%    ds('d2sm','',1,B,S,round(size(L,3)/2)+1);
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  %% THIS CODE IS NOT PROCESSED BECAUSE THIS FILE IS ONLY TO DISPLAY THE HELP
  rev = '$Rev$';

  disp('Compiling cat_amap.c')

  pth = fileparts(which(mfilename));
  p_path = pwd;
  cd(pth);
  mex -O cat_amap.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
  cd(p_path);
  
  [prob, mean] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm);

  
  if 0
    %% call test case
    % create a labelmap with 3 large areas
    L = zeros(50,50,9,'single'); 
    L(round(size(L,1)*0.27), round(size(L,2)*0.21), round(size(L,3)/2)+1) = 1; 
    L(round(size(L,1)*0.85), round(size(L,2)*0.45), round(size(L,3)/2)+1) = 2; 
    L(round(size(L,1)*0.37), round(size(L,2)*0.76), round(size(L,3)/2)+1) = 3; 
    B = false(size(L)); B(5:end-4,5:end-4,:) = true; 
    [D,I] = cat_vbdist(L,B); L=L(I); L = uint8(round(L));

    %% simulate real image with noise and bias
    A = double(L);  
    for i=1:size(A,2), A(:,i,:) = A(:,i,:) .* (( i/size(A,2) ) + 0.5) / 3; end
    A = A + 0.1 * rand(size(L)); 
    ds('d2smns','',1,A,L,round(size(L,3)/2)+1);

    % use AMAP to segment the simulated image
    Atmp = A + 0; Ltmp = L + 0; % create copies!
    [prob, mean] = cat_amap(Atmp,Ltmp, 3,16,32, 5,0, 0.1,ones(1,3),50,60);
    S = single(prob(:,:,:,1))*1/3/255 + single(prob(:,:,:,2))*2/3/255 + single(prob(:,:,:,3))/255;
    ds('d2sm','',1,A,S,round(size(L,3)/2)+1);

    %% another test
    B = double(L)/3 + 0.01 * randn(size(L)); Btmp = B + 0;  Ltmp = L + 0;
    amapres = evalc('[prob, mean] = cat_amap(Btmp,Ltmp, 3,16,32, 5,0, 0.1,ones(1,3),50,60)');
    S = single(prob(:,:,:,1))*1/3/255 + single(prob(:,:,:,2))*2/3/255 + single(prob(:,:,:,3))/255;
    ds('d2sm','',1,B,S,round(size(L,3)/2)+1);
    %rmse = sum( (S(:) - (single(L(:))/3)).^2 ).^(1/2)

  end

  
return
