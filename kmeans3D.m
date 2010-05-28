function mu = kmeans3D(ima,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   kmeans image segmentation
%
%   Input:
%          ima: grey color image
%          k: Number of classes
%   Output:
%          mu: vector of class means 
%
%   Author: Jose Vicente Manjon Herrera
%    Email: jmanjon@fis.upv.es
%     Date: 27-08-2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_______________________________________________________________________
% Christian Gaser
% $Id$

scale = 1000;       % scale to that value to have enough values in the histogram

% check image
ima = double(ima);
ima = ima(:);       % vectorize ima
mx  = max(ima);
ima = scale*ima/mx; % scale to have enough values in the histogram
mi  = min(ima);     % deal with negative 
ima = ima-mi+1;     % and zero values
ima = round(ima);

s   = length(ima);

% create image histogram

m  = max(ima)+1;
h  = zeros(1,m);
hc = zeros(1,m);

for i = 1:s
  if(ima(i)>0) h(ima(i)) = h(ima(i))+1;end;
end

ind = find(h);
hl  = length(ind);

% initiate centroids

mu = (1:k)*m/(k+1);

% start process

while(true)
  
  oldmu = mu;
  
  % current classification  
  for i = 1:hl
      c = abs(ind(i)-mu);
      cc = find(c == min(c));
      hc(ind(i)) = cc(1);
  end
  
  % recalculation of means  
  for i = 1:k, 
      a = find(hc == i);
      mu(i) = sum(a.*h(a))/sum(h(a));
  end
  
  if(mu == oldmu) break; end;
  
end

mu = mu+mi-1;   % recover real range
mu = mx*mu/scale;

