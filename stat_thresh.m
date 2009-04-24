function [peak_threshold, extent_threshold, peak_threshold_1, extent_threshold_1] = ...
   stat_thres(search_volume, num_voxels, fwhm, df, p_val_peak, ...
   cluster_threshold, p_val_extent, nconj, nvar, EC_file, expr)
%
% stat_thresh.m
% Modified version of STAT_THRESHOLD function by Keith Worsley. The original
% STAT_THRESHOLD function is part of the FMRISTAT packgage. The main use of the
% original version is to calculate peak and cluster thresholds, both corrected
% and uncorrected. Details on input and output arguments are found in the
% original STAT_THRESHOLD function available from Keith Worsley's web site
% at McGill University's Math & Stat department. 
%
% This stat_thresh.m function is a customized version to be called by a function
% spm_list_nS for non-stationarity correction of RFT-based cluster size test.
% The input and output of this function is therefore modified for producing cluster
% p-values (corrected) under non-stationarity. The modification includes:
%   -supressing the output from being displayed.
%   -the number of cluster p-values it can calculate has been increased to 500 clusters
%    (the default in the original version was 5).
%   -the p_val_extent is treated as extent, no matter how small it is.
%
% stat_thresh is called by spm_list_nS in the following format:
% [PEAK_P CLUSTER_P PEAK_P_1 CLUSTER_P_1] = 
%    stat_thresh(V_RESEL,NUM_VOX,1,[DF_ER;DF_RPV],ALPHA,CL_DEF_TH,CL_RESEL);
% PARAMETERS:
%    V_RESEL:      The seach volume in terms of resels. It is a 1x4 vector
%                  describing the topological characteristic of the search
%                  volume.
%    NUM_VOX:      The number of voxels in the search volume.
%    DF_ER:        Degrees of freedom of error
%    DF_RPV:       Degrees of freedom of RPV image estimation. Usually the same
%                  as the error df.
%    ALPHA:        The significance level of the peak (arbitrarily set to 0.05) 
%    CL_DEF_TH:    The cluster defining threshold. Can be entered in terms of
%                  a p-value (uncorrected) or a t-value.
%    CL_RESEL:     The cluster size in terms of resel
%
%    PEAK_P:       Peak p-value (FWE corrected). Not used for our purpose.
%    PEAK_P_1:     Peak p-value (uncorrected). Not used for our purpose.
%    CLUSTER_P:    Cluster p-value (FWE corrected)
%    CLUSTER_P_1:  Cluster p-value (uncorrected)
%    
%                           ----------------
%
% More etails on non-stationary cluster size test can be found in
%
% Worsley K J, Andermann M, Koulis T, MacDonald D and Evans A C
%   Detecting Changes in Nonisotropic Images
%   Human Brain Mapping 8: 98-101 (1999)
%
% Hayasaka S, Phan K L, Liberzon I, Worsley K J, and Nichols T E
%   Nonstationary cluster size inference with random-field and permutation methods
%   NeuroImage 22: 676-687 (2004)
%
%
%-----------------------------------------------------------------------------------
% Version 0.76b   Feb 19, 2007  by Satoru Hayasaka
%

%############################################################################
% COPYRIGHT:   Copyright 2003 K.J. Worsley 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              keith.worsley@mcgill.ca , www.math.mcgill.ca/keith
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

%############################################################################
% UPDATES:
%
%          Variable nvar is rounded so that it is recognized as an integer.
%          Feb 19, 2007 by Satoru Hayasaka
%############################################################################


% Defaults:
if nargin<1;  search_volume=[];  end
if nargin<2;  num_voxels=[];  end
if nargin<3;  fwhm=[];  end
if nargin<4;  df=[];  end
if nargin<5;  p_val_peak=[];  end
if nargin<6;  cluster_threshold=[];  end
if nargin<7;  p_val_extent=[];  end
if nargin<8;  nconj=[];  end
if nargin<9;  nvar=[];  end
if nargin<10;  EC_file=[];  end
if nargin<11;  expr=[];  end

if isempty(search_volume);  search_volume=1000000;  end
if isempty(num_voxels);  num_voxels=1000000;  end
if isempty(fwhm);  fwhm=0.0;  end
if isempty(df);  df=Inf;  end
if isempty(p_val_peak);  p_val_peak=0.05;  end
if isempty(cluster_threshold);  cluster_threshold=0.001;  end
if isempty(p_val_extent);  p_val_extent=0.05;  end
if isempty(nconj);  nconj=1;  end
if isempty(nvar);  nvar=1;  end

if size(fwhm,1)==1; fwhm(2,:)=fwhm; end
if size(fwhm,2)==1; scale=1; else scale=fwhm(1,2)/fwhm(1,1); fwhm=fwhm(:,1); end;
isscale=(scale>1); 

if length(num_voxels)==1; num_voxels(2,1)=1; end

if size(search_volume,2)==1
   radius=(search_volume/(4/3*pi)).^(1/3);
   search_volume=[ones(length(radius),1) 4*radius 2*pi*radius.^2 search_volume];
end
if size(search_volume,1)==1
   search_volume=[search_volume; [1 zeros(1,size(search_volume,2)-1)]];
end
lsv=size(search_volume,2);
fwhm_inv=all(fwhm>0)./(fwhm+any(fwhm<=0));
resels=search_volume.*repmat(fwhm_inv,1,lsv).^repmat(0:lsv-1,2,1);
invol=resels.*(4*log(2)).^(repmat(0:lsv-1,2,1)/2);
for k=1:2
   D(k,1)=max(find(invol(k,:)))-1;
end

% determines which method was used to estimate fwhm (see fmrilm or multistat): 
df_limit=4;

% max number of pvalues or thresholds to print:
% it can print out a ton of stuff! (the original default was 5)
nprint=500;

if length(df)==1; df=[df 0]; end
if size(df,1)==1; df=[df; Inf Inf]; end
if size(df,2)==1; df=[df [0; df(2,1)]]; end

% is_tstat=1 if it is a t statistic
is_tstat=(df(1,2)==0);
if is_tstat
   df1=1;
   df2=df(1,1);
else
   df1=df(1,1);
   df2=df(1,2);
end
if df2 >= 1000; df2=Inf; end
df0=df1+df2;

dfw1=df(2,1);
dfw2=df(2,2);
if dfw1 >= 1000; dfw1=Inf; end
if dfw2 >= 1000; dfw2=Inf; end

if length(nvar)==1; nvar(2,1)=df1; end
nvar = round(nvar); %-to make sure that nvar is integer!

if isscale & (D(2)>1 | nvar(1,1)>1 | df2<Inf)
   D
   nvar
   df2
   fprintf('Cannot do scale space.');
   return
end

Dlim=D+[scale>1; 0];
DD=Dlim+nvar-1;

% Values of the F statistic:
t=((1000:-1:1)'/100).^4;
% Find the upper tail probs cumulating the F density using Simpson's rule:
if df2==Inf
   u=df1*t;
   b=exp(-u/2-log(2*pi)/2+log(u)/4)*df1^(1/4)*4/100;
else  
   u=df1*t/df2;
   b=exp(-df0/2*log(1+u)+log(u)/4-betaln(1/2,(df0-1)/2))*(df1/df2)^(1/4)*4/100;
end
t=[t; 0];
b=[b; 0];
n=length(t);
sb=cumsum(b);
sb1=cumsum(b.*(-1).^(1:n)');
pt1=sb+sb1/3-b/3;
pt2=sb-sb1/3-b/3;
tau=zeros(n,DD(1)+1,DD(2)+1);
tau(1:2:n,1,1)=pt1(1:2:n);
tau(2:2:n,1,1)=pt2(2:2:n);
tau(n,1,1)=1;
tau(:,1,1)=min(tau(:,1,1),1);

% Find the EC densities:
u=df1*t;
for d=1:max(DD)
   for e=0:min(min(DD),d)
      s1=0;
      cons=-((d+e)/2+1)*log(pi)+gammaln(d)+gammaln(e+1);
      for k=0:(d-1+e)/2
         [i,j]=ndgrid(0:k,0:k);
         if df2==Inf
            q1=log(pi)/2-((d+e-1)/2+i+j)*log(2);
         else
            q1=(df0-1-d-e)*log(2)+gammaln((df0-d)/2+i)+gammaln((df0-e)/2+j) ...
               -gammalni(df0-d-e+i+j+k)-((d+e-1)/2-k)*log(df2);
         end
         q2=cons-gammalni(i+1)-gammalni(j+1)-gammalni(k-i-j+1) ...
            -gammalni(d-k-i+j)-gammalni(e-k-j+i+1);
         s2=sum(sum(exp(q1+q2)));
         if s2>0
            s1=s1+(-1)^k*u.^((d+e-1)/2-k)*s2;
         end
      end
      if df2==Inf
         s1=s1.*exp(-u/2);
      else
         s1=s1.*exp(-(df0-2)/2*log(1+u/df2));
      end
      if DD(1)>=DD(2)
         tau(:,d+1,e+1)=s1;
         if d<=min(DD)
            tau(:,e+1,d+1)=s1;
         end
      else
         tau(:,e+1,d+1)=s1;      
         if d<=min(DD)
            tau(:,d+1,e+1)=s1;
         end
      end
   end
end

% For multivariate statistics, add a sphere to the search region:
a=zeros(2,max(nvar));
for k=1:2
   j=(nvar(k)-1):-2:0;
   a(k,j+1)=exp(j*log(2)+j/2*log(pi) ...
      +gammaln((nvar(k)+1)/2)-gammaln((nvar(k)+1-j)/2)-gammaln(j+1));
end
rho=zeros(n,Dlim(1)+1,Dlim(2)+1);
for k=1:nvar(1)
   for l=1:nvar(2)
      rho=rho+a(1,k)*a(2,l)*tau(:,(0:Dlim(1))+k,(0:Dlim(2))+l);
   end
end

if is_tstat
   t=[sqrt(t(1:(n-1))); -flipdim(sqrt(t),1)];
   rho=[rho(1:(n-1),:,:); flipdim(rho,1)]/2;
   for i=0:D(1)
      for j=0:D(2)
         rho(n-1+(1:n),i+1,j+1)=-(-1)^(i+j)*rho(n-1+(1:n),i+1,j+1);
      end
   end
   rho(n-1+(1:n),1,1)=rho(n-1+(1:n),1,1)+1;
   n=2*n-1;
end

% For scale space:
if scale>1
   kappa=D(1)/2;
   tau=zeros(n,D(1)+1);
   for d=0:D(1)
      s1=0;
      for k=0:d/2
         s1=s1+(-1)^k/(1-2*k)*exp(gammaln(d+1)-gammaln(k+1)-gammaln(d-2*k+1) ...
            +(1/2-k)*log(kappa)-k*log(4*pi))*rho(:,d+2-2*k,1);
      end
      if d==0
         cons=log(scale);
      else
         cons=(1-1/scale^d)/d;
      end
      tau(:,d+1)=rho(:,d+1,1)*(1+1/scale^d)/2+s1*cons;
   end
   rho(:,1:(D(1)+1),1)=tau;
end

if D(2)==0
   d=D(1);
   if nconj>1
      % Conjunctions:
      b=gamma(((0:d)+1)/2)/gamma(1/2);
      for i=1:d+1
         rho(:,i,1)=rho(:,i,1)/b(i);
      end
      m1=zeros(n,d+1,d+1);
      for i=1:d+1
         j=i:d+1;
         m1(:,i,j)=rho(:,j-i+1,1);
      end
      for k=2:nconj
         for i=1:d+1
            for j=1:d+1
               m2(:,i,j)=sum(rho(:,1:d+2-i,1).*m1(:,i:d+1,j),2);
            end
         end
         m1=m2;
      end
      for i=1:d+1
         rho(:,i,1)=m1(:,1,i)*b(i);
      end
   end
   
   if ~isempty(EC_file)
      if d<3
         rho(:,(d+2):4,1)=zeros(n,4-d-2+1);
      end
      fid=fopen(EC_file,'w');
      % first 3 are dimension sizes as 4-byte integers:
      fwrite(fid,[n max(d+2,5) 1],'int');
      % next 6 are bounding box as 4-byte floats: 
      fwrite(fid,[0 0 0; 1 1 1],'float');
      % rest are the data as 4-byte floats:
      if ~isempty(expr)
         eval(expr);
      end
      fwrite(fid,t,'float');
      fwrite(fid,rho,'float');
      fclose(fid);
   end
end

if all(fwhm>0)
   pval_rf=zeros(n,1);
   for i=1:D(1)+1
      for j=1:D(2)+1
         pval_rf=pval_rf+invol(1,i)*invol(2,j)*rho(:,i,j);
      end
   end
else
   pval_rf=Inf;
end

% Bonferroni 
pt=rho(:,1,1);
pval_bon=abs(prod(num_voxels))*pt;

% Minimum of the two:
pval=min(pval_rf,pval_bon);

tlim=1;
if p_val_peak(1) <= tlim
   peak_threshold=minterp1(pval,t,p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold;
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak=interp1(t,pval,p_val_peak);
   peak_threshold=P_val_peak;
   if length(p_val_peak)<=nprint
      P_val_peak;
   end
end

if fwhm<=0 | any(num_voxels<0)
   peak_threshold_1=p_val_peak+NaN;
   extent_threshold=p_val_extent+NaN;
   extent_threshold_1=extent_threshold;
   return
end

% Cluster_threshold:
%###-Changed so that cluster_threshold is considered as cluster extent no matter what.
if cluster_threshold > eps
   tt=cluster_threshold;
else
   % cluster_threshold is treated as a probability:
   tt=minterp1(pt,t,cluster_threshold);
   Cluster_threshold=tt;
end

d=D(1);
rhoD=interp1(t,rho(:,d+1,1),tt);
p=interp1(t,pt,tt);

% Pre-selected peak:

pval=rho(:,d+1,1)./rhoD;

if p_val_peak(1) <= tlim 
   peak_threshold_1=minterp1(pval,t, p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold_1;
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak_1=interp1(t,pval,p_val_peak);
   peak_threshold_1=P_val_peak_1;
   if length(p_val_peak)<=nprint
      P_val_peak_1;
   end
end

if  D(1)==0 | nconj>1 | nvar(1)>1 | D(2)>0 | scale>1
    extent_threshold=p_val_extent+NaN;
    extent_threshold_1=extent_threshold;
    if length(p_val_extent)<=nprint
       extent_threshold;
       extent_threshold_1;
    end
    return
end

% Expected number of clusters:

%###-Change tlim to a small number so that p_val_extent is considered as cluster extent
tlim = eps;

EL=invol(1,d+1)*rhoD;
cons=gamma(d/2+1)*(4*log(2))^(d/2)/fwhm(1)^d*rhoD/p;

if df2==Inf & dfw1==Inf
   if p_val_extent(1) <= tlim 
      pS=-log(1-p_val_extent)/EL;
      extent_threshold=(-log(pS)).^(d/2)/cons;
      pS=-log(1-p_val_extent);
      extent_threshold_1=(-log(pS)).^(d/2)/cons;
      if length(p_val_extent)<=nprint
         extent_threshold;
         extent_threshold_1;
      end
   else
      % p_val_extent is now treated as a spatial extent:
      pS=exp(-(p_val_extent*cons).^(2/d));
      P_val_extent=1-exp(-pS*EL);
      extent_threshold=P_val_extent;
      P_val_extent_1=1-exp(-pS);
      extent_threshold_1=P_val_extent_1;
      if length(p_val_extent)<=nprint
         P_val_extent;
         P_val_extent_1;
      end
   end
else
   % Find dbn of S by taking logs then using fft for convolution:
   ny=2^12;
   a=d/2;
   b2=a*10*max(sqrt(2/(min(df1+df2,dfw1))),1);
   if df2<Inf
      b1=a*log((1-(1-0.000001)^(2/(df2-d)))*df2/2);
   else
      b1=a*log(-log(1-0.000001));
   end
   dy=(b2-b1)/ny;
   b1=round(b1/dy)*dy;
   y=((1:ny)'-1)*dy+b1;
   numrv=1+(d+1)*(df2<Inf)+d*(dfw1<Inf)+(dfw2<Inf);
   f=zeros(ny,numrv);
   mu=zeros(1,numrv);
   if df2<Inf
      % Density of log(Beta(1,(df2-d)/2)^(d/2)):
      yy=exp(y./a)/df2*2;  
      yy=yy.*(yy<1);
      f(:,1)=(1-yy).^((df2-d)/2-1).*((df2-d)/2).*yy/a;
      mu(1)=exp(gammaln(a+1)+gammaln((df2-d+2)/2)-gammaln((df2+2)/2)+a*log(df2/2));
   else
      % Density of log(exp(1)^(d/2)):
      yy=exp(y./a);   
      f(:,1)=exp(-yy).*yy/a;
      mu(1)=exp(gammaln(a+1));
   end
   
   nuv=[];
   aav=[];
   if df2<Inf
      nuv=[df1+df2-d  df2+2-(1:d)];
      aav=[a repmat(-1/2,1,d)]; 
   end
   if dfw1<Inf
      if dfw1>df_limit
         nuv=[nuv dfw1-dfw1/dfw2-(0:(d-1))];
      else
         nuv=[nuv repmat(dfw1-dfw1/dfw2,1,d)];
      end
      aav=[aav repmat(1/2,1,d)];
   end   
   if dfw2<Inf
      nuv=[nuv dfw2];
      aav=[aav -a];
   end   
   
   for i=1:(numrv-1)
      nu=nuv(i);
      aa=aav(i);
      yy=y/aa+log(nu);
      % Density of log((chi^2_nu/nu)^aa):
      f(:,i+1)=exp(nu/2*yy-exp(yy)/2-(nu/2)*log(2)-gammaln(nu/2))/abs(aa);
      mu(i+1)=exp(gammaln(nu/2+aa)-gammaln(nu/2)-aa*log(nu/2));
   end
   % Check: plot(y,f); sum(f*dy,1) should be 1
      
   omega=2*pi*((1:ny)'-1)/ny/dy;
   shift=complex(cos(-b1*omega),sin(-b1*omega))*dy;
   prodfft=prod(fft(f),2).*shift.^(numrv-1);
   % Density of Y=log(B^(d/2)*U^(d/2)/sqrt(det(Q))):
   ff=real(ifft(prodfft));
   % Check: plot(y,ff); sum(ff*dy) should be 1
   mu0=prod(mu);
   % Check: plot(y,ff.*exp(y)); sum(ff.*exp(y)*dy.*(y<10)) should equal mu0   
   
   alpha=p/rhoD/mu0*fwhm(1)^d/(4*log(2))^(d/2);
   
   % Integrate the density to get the p-value for one cluster: 
   pS=cumsum(ff(ny:-1:1))*dy;
   pS=pS(ny:-1:1);
   % The number of clusters is Poisson with mean EL:
   pSmax=1-exp(-pS*EL);
   
   if p_val_extent(1) <= tlim 
      yval=minterp1(-pSmax,y,-p_val_extent);
      % Spatial extent is alpha*exp(Y) -dy/2 correction for mid-point rule:
      extent_threshold=alpha*exp(yval-dy/2);
      % For a single cluster:
      yval=minterp1(-pS,y,-p_val_extent);
      extent_threshold_1=alpha*exp(yval-dy/2);
      if length(p_val_extent)<=nprint
         extent_threshold;
         extent_threshold_1;
      end
   else
      % p_val_extent is now treated as a spatial extent:
      P_val_extent=interp1(y,pSmax,log(p_val_extent/alpha)+dy/2);
      extent_threshold=P_val_extent;
      % For a single cluster:
      P_val_extent_1=interp1(y,pS,log(p_val_extent/alpha)+dy/2);
      extent_threshold_1=P_val_extent_1;
      if length(p_val_extent)<=nprint
         P_val_extent;
         P_val_extent_1;
      end
   end
   
end
   
return

function x=gammalni(n);
i=find(n>=0);
x=Inf+n;
if ~isempty(i)
   x(i)=gammaln(n(i));
end
return

function iy=minterp1(x,y,ix);
% interpolates only the monotonically increasing values of x at ix
n=length(x);
mx=x(1);
my=y(1);
xx=x(1);
for i=2:n
   if x(i)>xx
      xx=x(i);
      mx=[mx xx];
      my=[my y(i)];
   end
end
iy=interp1(mx,my,ix);
return
