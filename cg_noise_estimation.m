function h = cg_noise_estimation(ima)
% FORMAT h = cg_noise_estimation(ima)
%
% 
% ***************************************************************************
%  The noise estimation is described in:                                       
%                                                                         
%  S. Aja-Fern‡ndez, A. Trist‡n-Vega, C. Alberola-L—pez.     
%  Noise estimation in single- and multiple-coil magnetic resonance data 
%  based on statistical models (2009).
%  Magnetic Resonance Imaging, 27, 1397-1409.                                                             
% ***************************************************************************
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

% estimate local mean
k = 3;
localMean = (convn(single(ima),ones(k,k,k),'same')/k^3);

% find non-zero regions
ind = find(localMean>0);

% If image has non-zero background (we check that less than 5% of the image are zero) 
% we can assume Rayleigh PDF in the background and noise estimation can be based on
% mode of local mean (equation 11 in Aja-Fern‡ndez et al. 2009)
if length(ind)>0.95*prod(size(ima))
  h = sqrt(2/pi)*moda(localMean(ind),1000);
else % otherwise use mode of local variance (equation 15 in Aja-Fern‡ndez et al. 2009)
  localVar = (convn(single(ima).^2,ones(k,k,k),'same')/k^3) - localMean.^2;
  h = sqrt(moda(localVar(ind),1000));
end

return

function m = moda(u,N)
% MODA   Mode of a distribution
%
%    m=MODE(u,N) calculates the mode of the set of data "u" using the histogram.
%    To avoid outliers, for the calculation are only taken into account those
%    values less than mean+2sigma;
%
%    INPUT:
%
%	- u (set of data)
%       - N: Number of points for the histogram. If N=0 then 5000 points are
%            considered
%
%   Author: Santiago Aja Fernandez
%   LOCAL STATISTICS TOOLBOX
%
%   Modified: Feb 01 2008
%

if N==0
	N = 1000;
end

u = single(u(:));

M1 = mean(u);
V1 = std(u);
C2 = u(u<=(M1+2*V1));
whos
[h,x] = hist(C2,N);
[M,M2] = max(h);

m = x(M2);
