function h = gaussian_noise_estimation(ima)
% FORMAT h = gaussian_noise_estimation(ima)
%
% Estimate gaussian noise in homogenious region of WM
% 
% ***************************************************************************
%  The noise estimation is described in:                                       
%                                                                         
%  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     
%  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic
%  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, 
%  April 2008                                                             
% ***************************************************************************
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

% convolve image using 6 neighbors
kernel = zeros(3,3,3);
kernel(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
kernel(:,:,2) = [0 1 0; 1 0 1; 0 1 0];
kernel(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

fima = (convn(ima,kernel,'same')/sum(kernel(:)));

% kmeans labeling into 3 classes
label = KmeansMex(fima, 3);

% use class with the highest value (e.g. WM)
ind = find(label==3);

% calculate pseudo residuals
residual = sqrt(6/7)*(ima(ind) - fima(ind));

% use sqrt of standard deviation of noise
h = sqrt(mean(residual.^2));