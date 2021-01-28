function [mu,su,nu] = cat_stat_kmeans(y,k,s)
% K-means clustering
%_______________________________________________________________________
% FORMAT [mu,su,nu] = cat_stat_kmeans(y,k[,s])
% 
%  y  .. data 
%  k  .. Number of components
%  s  .. select peak (s>0 & s<=k) to have some side peaks just for
%        stabilization, s==0 select maximum peak
% 
%
%  mu .. vector of class means 
%  su .. vector of class std 
%  nu .. vector of class number of voxel
%
% modified version of
% spm_kmeans1.m 1143 2008-02-07 19:33:33Z spm $
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id$

if nargin<1, help cat_stat_kmeans; return; end
if nargin<2, k=1; end

k = max(1,k); 

dt = class(y); 
y = double(y); 
y = y(:)';
y(isnan(y))=[]; % remove NaNs 
if numel(y)<=0
  mu = nan(1,k);
  su = mu; 
  nu = mu; 
  return; 
end

N = length(y);

% Spread seeds evenly according to CDF
x = sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i = N; %(ones(1,N);
mu = x(seeds);
su = zeros(size(mu));
nu = ones(size(mu));

d = zeros(k,length(y));
for loops = 1:1000  

  for j=1:k
    d(j,:) = (y-mu(j)).^2;
  end
  [tmp,i] = min(d); clear tmp %#ok<ASGLU>
  if sum(i - last_i)==0 || isempty(y(i==j))
    % If assignment is unchanged
    break;
  else
   % Recompute centres
   for j=1:k
     mu(j) = mean(y(i==j));
   end
   last_i=i;
  end
end  

% Compute variances and mixing proportions
for j=1:k
  if isempty(y(i==j))
    su(j) = std(d(j,:));
    nu(j) = sum(std(d(j,:)))./numel(y(:));
  else
    su(j) = std(y(i==j));
    nu(j) = sum(i==j)./numel(y(:));
  end
end

if exist('s','var')
  if s>k 
    error('s has to be >=0 and <=k.'); 
  end 
  if s==0 
    [tmp,s] = max(nu); %#ok<ASGLU> % select maximum
  end
  mu = mu(s);
  su = su(s);
  nu = nu(s); 
end
  
if strcmp(dt,'double') %#ok<STISA>
  feval(dt,mu); 
  feval(dt,su); 
  feval(dt,nu); 
end
    
    