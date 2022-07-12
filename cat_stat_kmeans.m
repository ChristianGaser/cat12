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
%  nu .. vector of class percentage of values
%
% modified version of
% spm_kmeans1.m 1143 2008-02-07 19:33:33Z spm $
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
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
  else
    su(j) = std(y(i==j));
  end
end

bd = nan(k,2); 
for j=1:k
  % lower boundary 
  if j==1, bd(j,1) = -inf; end
  if j>1,  bd(j,1) = mean( [ mu(j-1) + su(j-1) , mu(j) - su(j) ] ); end
  % upper boudnary  
  if j<k,  bd(j,2) = mean( [ mu(j) + su(j) , mu(j+1) - su(j+1) ] ); end
  if j==k, bd(j,2) = inf; end
  % proportion of elements
  nu(j) = sum( y>bd(j,1) & y<=bd(j,2) )/numel(y);  
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
    
    