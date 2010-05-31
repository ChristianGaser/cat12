function h = cg_noise_estimation(ima)
% FORMAT h = cg_noise_estimation(ima)
%
% h    - Noise estimate (SD)
%
% based on RicianSTD.m from 
% Pierrick Coupe - pierrick.coupe@gmail.com                                  
% Jose V. Manjon - jmanjon@fis.upv.es                                        
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon  
%
% ***************************************************************************
%  The noise estimation is described in:                                       
%                                                                         
% Coupe P, Manjon JV, Gedamu E, Arnold D, Robles M, Collins DL.                                                  
% Robust Rician noise estimation for MR images.
% Med Image Anal. 2010 
% ***************************************************************************
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

s = size(ima);

% Estimation of the zeros pading size
p(1) = 2^(ceil(log2(s(1))));
p(2) = 2^(ceil(log2(s(2))));
p(3) = 2^(ceil(log2(s(3))));

% Zeros Pading
pad1 = zeros(p(1),p(2),p(3));
pad1(1:s(1),1:s(2),1:s(3)) = ima;

% Wavelet Transform
[af, sf] = farras;
w1 = dwt3D(pad1,1,af);

clear pad1

% Removing region corresponding to zeros pading
tmp  = w1{1}{7};
tmp2 = w1{2};
clear w1

tmp  = tmp(1:round((s(1)-1)/2),1:round((s(2)-1)/2),1:round((s(3)-1)/2));
tmp2 = tmp2(1:round((s(1)-1)/2),1:round((s(2)-1)/2),1:round((s(3)-1)/2));

% Detection of the object in the LLL subband
mu = kmeans3D(tmp2,2);
th = mean(mu);
map = tmp2 > th;

% Detection of the High gradient area in the LLL subband
[PX,PY,PZ] = gradient(tmp2);
clear tmp2

GR = sqrt(PX.^2 + PY.^2 + PZ.^2);
clear PX PY PZ

m = median(GR(map));
map2 = GR < m;
clear GR

% Map containing Object without strong edges
map = map & map2;

% Estimation of the magnitude noise STD in HHH subband
Nsig = median(abs(tmp(map)))/0.6745;
clear tmp

% Computation of SNR on object 
fima = convn(ima,ones(3,3,3),'same');
mu = kmeans3D(fima,2);
th = mean(mu);
map = find(fima > th);
SNR = mean(ima(map)) / Nsig;

% Iterative estimation of truth SNR based on Koay method
for un = 1:500
    SNR2 = sqrt(epsi(SNR)*(1 + mean(ima(map))^2  / Nsig^2 )-2);
    
    if abs(SNR-SNR2) < 0.000000001 
        break;
    end
    
    SNR = SNR2;
end
    
h = sqrt((Nsig^2 / epsi(SNR)));

PSNR = 20*log10(max(ima(:))/h);

return

function res = epsi(SNR)

% Pierrick Coupe - pierrick.coupe@gmail.com                                                                        
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe 

% Based on Koay estimation of truth SNR from Magnitude SNR.

if (SNR > 37) res = 1;

else
    
res = 2 + SNR^2 - pi/8 * exp(-(SNR^2/2))*((2+SNR^2)*besseli(0,(SNR^2/4)) + SNR^2*besseli(1,(SNR^2/4)))^2;

end
