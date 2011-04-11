function mu = kmeans3D(y,k)
% K-means clustering
% FORMAT mu = kmeans3D (y,k)
% 
% y          data 
% k          Number of components
%
% mu         vector of class means 
%
% modified version of
% $Id$
%_______________________________________________________________________
% Christian Gaser
% $Id$

y=y(:)';
N=length(y);

% Spread seeds evenly according to CDF
[x,i]=sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i=ones(1,N);
mu=x(seeds);
for loops=1:100,  
 for j=1:k,
   d(j,:)=(y-mu(j)).^2;
 end
 [tmp,i]=min(d);
 if sum(i-last_i)==0
   % If assignment is unchanged
   break;
 else
   % Recompute centres
   for j=1:k,
     mu(j)=mean(y(i==j));
   end
   last_i=i;
 end
end  
