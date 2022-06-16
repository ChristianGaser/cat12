function hAxes = cat_plot_scatter(X,Y, varargin)
% cat_plot_scatter creates a scatter plot coloured by density.
%
% cat_plot_scatter(X,Y) creates a scatterplot of X and Y at the locations
% specified by the vectors X and Y (which must be the same size), colored
% by the density of the points.
%
% cat_plot_scatter(...,'MARKER',M) allows you to set the marker for the
% scatter plot. Default is 's', square.
%
% cat_plot_scatter(...,'MSIZE',MS) allows you to set the marker size for the
% scatter plot. Default is 10.
%
% cat_plot_scatter(...,'FILLED',false) sets the markers in the scatter plot to be
% outline. The default is to use filled markers.
%
% cat_plot_scatter(...,'FIT_POLY',N) fit polynomial with degree N
% The default is 0.
%
% cat_plot_scatter(...,'CI',false) add plot of confidence interval for
% polynomial fit. The default is true.
%
% cat_plot_scatter(...,'CMAP',cmap) defines colormap
% The default is parula (if available) or jet.
%
% cat_plot_scatter(...,'fig',fig) defines figure handle
%
% cat_plot_scatter(...,'jitter',true) adds jitter on x-axsis for
% categorical x-variables
%
% cat_plot_scatter(...,'PLOTTYPE',TYPE) allows you to create other ways of
% plotting the scatter data. Options are 'image','contour','dscatter', and 'scatter'.
% These create surf, mesh and contour plots colored by density of the
% scatter data. The default is 'dscatter'.
%
% cat_plot_scatter(...,'BINS',[NX,NY]) allows you to set the number of bins used
% for the 2D histogram used to estimate the density. The default is to
% use the number of unique values in X and Y up to a maximum of 200.
%
% cat_plot_scatter(...,'SMOOTHING',LAMBDA) allows you to set the smoothing factor
% used by the density estimator. The default value is 20 which roughly
% means that the smoothing is over 20 bins around a given point.
%
% cat_plot_scatter(...,'LOGY',true) uses a log scale for the yaxis.
%
% Examples:
%
%   data randn(1000,2);
%   cat_plot_scatter(data(:,1),10.^(data(:,2)/256),'log',1)
%   % Add contours
%   hold on
%   cat_plot_scatter(data(:,1),10.^(data(:,2)/256),'log',1,'plottype','contour')
%   hold off
%   
% See also SCATTER.
%
% modified by Christian Gaser (christian.gaser@uni-jena.de) and
% original dscatter version was written by Paul H. C. Eilers
%
% Copyright 2003-2004 The MathWorks, Inc.
% $Revision$   $Date$
% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

lambda      = [];
nbins       = [];
plottype    = 'dscatter';
contourFlag = false;
msize       = 10;
marker      = 's';
logy        = false;
filled      = true;
fit_poly    = 0;
ci          = true;
jitter      = false;
if exist('parula')
  cmap        = 'parula';
else
  cmap        = 'jet';
end

if nargin > 2
  if rem(nargin,2) == 1
    error('IncorrectNumberOfArguments',...
      'Incorrect number of arguments to %s.',mfilename);
  end
  okargs = {'smoothing','bins','plottype','logy','contourFlag','marker','msize','filled','fit_poly','ci','cmap','fig','jitter'};
  for j=1:2:nargin-2
    pname = varargin{j};
    pval = varargin{j+1};
    k = strmatch(lower(pname), okargs); %#ok
    if isempty(k)
      error('UnknownParameterName',...
        'Unknown parameter name: %s.',pname);
    elseif length(k)>1
      error('AmbiguousParameterName',...
        'Ambiguous parameter name: %s.',pname);
    else
      switch(k)
        case 1  % smoothing factor
          if isnumeric(pval)
            lambda = pval;
          else
            error('InvalidScoringMatrix','Invalid smoothing parameter.');
          end
        case 2
          if isscalar(pval)
            nbins = [ pval pval];
          else
            nbins = pval;
          end
        case 3
          plottype = pval;
        case 4
          logy = pval;
          Y = log10(Y);
        case 5
          contourFlag = pval;
        case 6
          marker = pval;
        case 7
          msize = pval;
        case 8
          filled = pval;
        case 9
          fit_poly = pval;
        case 10
          ci = pval;
        case 11
          cmap = pval;
        case 12
          fig = pval;
        case 13
          jitter = pval;
      end
    end
  end
end

X0 = X; Y0 = Y;
ind = ~isnan(X) & ~isnan(Y);
X = X(ind); Y = Y(ind); 

minx = min(X(:));
maxx = max(X(:));
miny = min(Y(:));
maxy = max(Y(:));

if isempty(nbins)
  nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
end

if isempty(lambda)
  lambda = 20;
end

edges1 = linspace(minx, maxx, nbins(1)+1);
ctrs1 = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(miny, maxy, nbins(2)+1);
ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];
[n,p] = size(X);
bin = zeros(n,2);

% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
[dum,bin(:,2)] = histc(X,edges1);
[dum,bin(:,1)] = histc(Y,edges2);

% remove zero histogram entries that can't be used
ind = find(bin(:,1)==0 | bin(:,2)==0);
bin(ind,:) = [];
X(ind) = [];
Y(ind) = [];

H = accumarray(bin,1,nbins([2 1])) ./ n;
G = smooth1D(H,nbins(2)/lambda);
F = smooth1D(G',nbins(1)/lambda)';

if logy
  ctrs2 = 10.^ctrs2;
  Y = 10.^Y;
end

if exist('fig','var')
  figure(fig)
else
  figure
end

% check whether X has limited number of entries (and is rather categorical)
n_categories_X = numel(unique(X));
add_range = min([1,log10(n_categories_X)/10]);

% polynomial fit and confidence interval
if fit_poly
  clf
  [p,S] = polyfit(X,Y,fit_poly);
  
  % add some more range because of jittered data
  if n_categories_X < 30 && jitter
    xfit = linspace(minx-2*add_range,maxx+2*add_range,100);
  else
    xfit = linspace(minx,maxx,100);
  end
  yfit = polyval(p,xfit);

  if ci
    [Y2,DELTA] = cg_polyconf(p,xfit,S);
    plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);
    plot_variance(xfit,Y2+DELTA,Y2-DELTA,[0.75 0.75 0.75])
  end

  % estimate R^2
  R2 = 1 - (S.normr/norm(Y - mean(Y)))^2;
  [cc,pp] = corrcoef(X,Y);
  fprintf('R^2 = %g\tr = %g\tp = %g\n',R2,cc(1,2),pp(1,2))
  fprintf('Coefficients = %g\n',p)

  hold on
end

% check whether X has limited number of entries (and is rather categorical)
% and add jitter if defined
if n_categories_X < 30 && jitter
  X = X + add_range*randn(size(X));
end

okTypes = {'contour','image','dscatter','scatter'};
k = strmatch(lower(plottype), okTypes); %#ok

if isempty(k)
  error('cat_plot_scatter:UnknownPlotType',...
    'Unknown plot type: %s.',plottype);
elseif length(k)>1
  error('cat_plot_scatter:AmbiguousPlotType',...
    'Ambiguous plot type: %s.',plottype);
else
  switch(k)
    case 1 %'contour'
      [dummy, h] = contour(ctrs1,ctrs2,F);
    case 2 %'image'
      nc = 256;
      F = F./max(F(:));
      colormap(repmat(linspace(1,0,nc)',1,3));
      h = image(ctrs1,ctrs2,floor(nc.*F) + 1);
    case 3 %'dscatter'
      F = F./max(F(:));
      ind = sub2ind(size(F),bin(:,1),bin(:,2));
      col = F(ind);
      if filled
        h = scatter(X,Y,msize,col,marker,'filled');
      else
        h = scatter(X,Y,msize,col,marker);
      end
    case 4 %'scatter'
      if filled
        h = scatter(X,Y,msize,marker,'filled');
      else
        h = scatter(X,Y,msize,marker);
      end
  end
end

colormap(cmap)

if fit_poly
  pl = plot(xfit,yfit,'k');
  set(pl,'LineWidth',2)
  if minx > 0 && minx < 1e-4
    xl = xlim;
    xlim([0 xl(2)])
  end
  if miny > 0 && miny < 1e-4
    yl = ylim;
    ylim([0 yl(2)])
  end
    hold off
end

if logy
  set(gca,'yscale','log');
end

if nargout > 0
  hAxes = get(h,'parent');
end

%%%% This method is quicker for symmetric data.
% function Z = filter2D(Y,bw)
% z = -1:(1/bw):1;
% k = .75 * (1 - z.^2);
% k = k ./ sum(k);
% Z = filter2(k'*k,Y);
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;

function [y, delta] = cg_polyconf(p,x,S,varargin)
%POLYCONF Polynomial evaluation and confidence interval estimation.
%   Y = POLYCONF(P,X) returns the value of a polynomial P evaluated at X. P
%   is a vector of length N+1 whose elements are the coefficients of the
%   polynomial in descending powers.
%
%       Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
%
%   If X is a matrix or vector, the polynomial is evaluated at all points
%   in X.  See also POLYVALM for evaluation in a matrix sense.
%
%   [Y,DELTA] = POLYCONF(P,X,S) uses the optional output, S, created by
%   POLYFIT to generate 95% prediction intervals.  If the coefficients in P
%   are least squares estimates computed by POLYFIT, and the errors in the
%   data input to POLYFIT were independent, normal, with constant variance,
%   then there is a 95% probability that Y +/- DELTA will contain a future
%   observation at X.
%
%   [Y,DELTA] = POLYCONF(P,X,S,'NAME1',VALUE1,'NAME2',VALUE2,...) specifies
%   optional argument name/value pairs chosen from the following list.
%   Argument names are case insensitive and partial matches are allowed.
%
%      Name       Value
%     'alpha'     A value between 0 and 1 specifying a confidence level of
%                 100*(1-alpha)%.  Default is alpha=0.05 for 95% confidence.
%     'mu'        A two-element vector containing centering and scaling
%                 parameters as computed by polyfit.  With this option,
%                 polyconf uses (X-MU(1))/MU(2) in place of x.
%     'predopt'   Either 'observation' (the default) to compute intervals for
%                 predicting a new observation at X, or 'curve' to compute
%                 confidence intervals for the polynomial evaluated at X.
%     'simopt'    Either 'off' (the default) for non-simultaneous bounds,
%                 or 'on' for simultaneous bounds.
%
%   See also POLYFIT, POLYTOOL, POLYVAL, INVPRED, POLYVALM.

%   For backward compatibility we also accept the following:
%   [...] = POLYCONF(p,x,s,ALPHA)
%   [...] = POLYCONF(p,x,s,alpha,MU)

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision$  $Date$

error(nargchk(2,Inf,nargin,'struct'));

alpha = [];
mu = [];
doobs = true;   % predict observation rather than curve estimate
dosim = false;  % give non-simultaneous intervals
if nargin>3
    if ischar(varargin{1})
        % Syntax with parameter name/value pairs
        okargs =   {'alpha' 'mu' 'predopt' 'simopt'};
        defaults = {0.05    []   'obs'     'off'};
        [eid emsg alpha mu predopt simopt] = ...
                internal.stats.getargs(okargs,defaults,varargin{:});
        if ~isempty(eid)
            error(sprintf('stats:polyconf:%s',eid),emsg);
        end
        
        i = find(strncmpi(predopt,{'curve';'observation'},length(predopt)));
        if ~isscalar(i)
            error('stats:polyconf:BadPredOpt', ...
           'PREDOPT must be one of the strings ''curve'' or ''observation''.');
        end
        doobs = (i==2);
        
        i = find(strncmpi(simopt,{'on';'off'},length(simopt)));
        if ~isscalar(i)
            error('stats:polyconf:BadSimOpt', ...
           'SIMOPT must be one of the strings ''on'' or ''off''.');
        end
        dosim = (i==1);
    else
        % Old syntax
        alpha = varargin{1};
        if numel(varargin)>=2
            mu = varargin{2};
        end
    end
end
if nargout > 1
    if nargin < 3, S = []; end % this is an error; let polyval handle it
    if nargin < 4 || isempty(alpha)
        alpha = 0.05;
    elseif ~isscalar(alpha) || ~isnumeric(alpha) || ~isreal(alpha) ...
                            || alpha<=0          || alpha>=1
        error('stats:polyconf:BadAlpha',...
              'ALPHA must be a scalar between 0 and 1.');
    end
    if isempty(mu)
        [y,delta] = polyval(p,x,S);
    else
        [y,delta] = polyval(p,x,S,mu);
    end
    if doobs
        predvar = delta;                % variance for predicting observation
    else
        s = S.normr / sqrt(S.df);
        delta = delta/s;
        predvar = s*sqrt(delta.^2 - 1); % get uncertainty in curve estimation
    end
    if dosim
        k = length(p);
        crit = sqrt(k * finv(1-alpha,k,S.df)); % Scheffe simultaneous value
    else
        crit = tinv(1-alpha/2,S.df);           % non-simultaneous value
    end
    delta = crit * predvar;
else
    if isempty(mu)
        y = polyval(p,x);
    else
        y = polyval(p,x,[],mu);
    end
end

function x = tinv(p,v);
% TINV   Inverse of Student's T cumulative distribution function (cdf).
%   X=TINV(P,V) returns the inverse of Student's T cdf with V degrees 
%   of freedom, at the values in P.
%
%   The size of X is the common size of P and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
% This is an open source function that was assembled by Eric Maris using
% open source subfunctions found on the web.

if nargin < 2, 
    error('Requires two input arguments.'); 
end

[errorcode p v] = distchck(2,p,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
x=zeros(size(p));

k = find(v < 0  | v ~= round(v));
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end

k = find(v == 1);
if any(k)
  x(k) = tan(pi * (p(k) - 0.5));
end

% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
k0 = find(p == 0);
if any(k0)
    tmp   = Inf;
    x(k0) = -tmp(ones(size(k0)));
end
k1 = find(p ==1);
if any(k1)
    tmp   = Inf;
    x(k1) = tmp(ones(size(k1)));
end

k = find(p >= 0.5 & p < 1);
if any(k)
    z = betainv(2*(1-p(k)),v(k)/2,0.5);
    x(k) = sqrt(v(k) ./ z - v(k));
end

k = find(p < 0.5 & p > 0);
if any(k)
    z = betainv(2*(p(k)),v(k)/2,0.5);
    x(k) = -sqrt(v(k) ./ z - v(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION distchck
%%%%%%%%%%%%%%%%%%%%%%%%%

function [errorcode,varargout] = distchck(nparms,varargin)
%DISTCHCK Checks the argument list for the probability functions.

errorcode = 0;
varargout = varargin;

if nparms == 1
    return;
end

% Get size of each input, check for scalars, copy to output
isscalar = (cellfun('prodofsize',varargin) == 1);

% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(isscalar)), return; end

n = nparms;

for j=1:n
   sz{j} = size(varargin{j});
end
t = sz(~isscalar);
size1 = t{1};

% Scalars receive this size.  Other arrays must have the proper size.
for j=1:n
   sizej = sz{j};
   if (isscalar(j))
      t = zeros(size1);
      t(:) = varargin{j};
      varargout{j} = t;
   elseif (~isequal(sizej,size1))
      errorcode = 1;
      return;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION betainv
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = betainv(p,a,b);
%BETAINV Inverse of the beta cumulative distribution function (cdf).
%   X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%   parameters A and B at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINV uses Newton's method to converge to the solution.

%   Reference:
%      [1]     M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964.

%   B.A. Jones 1-12-93

if nargin < 3, 
    error('Requires three input arguments.'); 
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

%   Initialize x to zero.
x = zeros(size(p));

%   Return NaN if the arguments are outside their respective limits.
k = find(p < 0 | p > 1 | a <= 0 | b <= 0);
if any(k),
   tmp = NaN;
   x(k) = tmp(ones(size(k))); 
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0), 
    x(k0) = zeros(size(k0)); 
end

k1 = find(p==1);
if any(k1), 
    x(k1) = ones(size(k1)); 
end

% Newton's Method.
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;

k = find(p > 0 & p < 1 & a > 0 & b > 0);
pk = p(k);

%   Use the mean as a starting guess. 
xk = a(k) ./ (a(k) + b(k));


% Move starting values away from the boundaries.
if xk == 0,
    xk = sqrt(eps);
end
if xk == 1,
    xk = 1 - sqrt(eps);
end


h = ones(size(pk));
crit = sqrt(eps); 

% Break out of the iteration loop for the following:
%  1) The last update is very small (compared to x).
%  2) The last update is very small (compared to 100*eps).
%  3) There are more than 100 iterations. This should NEVER happen. 

while(any(abs(h) > crit * abs(xk)) & max(abs(h)) > crit    ...
                                 & count < count_limit), 
                                 
    count = count+1;    
    h = (betacdf(xk,a(k),b(k)) - pk) ./ betapdf(xk,a(k),b(k));
    xnew = xk - h;

% Make sure that the values stay inside the bounds.
% Initially, Newton's Method may take big steps.
    ksmall = find(xnew < 0);
    klarge = find(xnew > 1);
    if any(ksmall) | any(klarge)
        xnew(ksmall) = xk(ksmall) /10;
        xnew(klarge) = 1 - (1 - xk(klarge))/10;
    end

    xk = xnew;  
end

% Return the converged value(s).
x(k) = xk;

if count==count_limit, 
    fprintf('\nWarning: BETAINV did not converge.\n');
    str = 'The last step was:  ';
    outstr = sprintf([str,'%13.8f'],h);
    fprintf(outstr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION betapdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = betapdf(x,a,b)
%BETAPDF Beta probability density function.
%   Y = BETAPDF(X,A,B) returns the beta probability density 
%   function with parameters A and B at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.33.

if nargin < 3, 
   error('Requires three input arguments.');
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y = zeros(size(x));

% Return NaN for parameter values outside their respective limits.
k1 = find(a <= 0 | b <= 0 | x < 0 | x > 1);
if any(k1)
   tmp = NaN;
    y(k1) = tmp(ones(size(k1))); 
end

% Return Inf for x = 0 and a < 1 or x = 1 and b < 1.
% Required for non-IEEE machines.
k2 = find((x == 0 & a < 1) | (x == 1 & b < 1));
if any(k2)
   tmp = Inf;
    y(k2) = tmp(ones(size(k2))); 
end

% Return the beta density function for valid parameters.
k = find(~(a <= 0 | b <= 0 | x <= 0 | x >= 1));
if any(k)
    y(k) = x(k) .^ (a(k) - 1) .* (1 - x(k)) .^ (b(k) - 1) ./ beta(a(k),b(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION betacdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = betacdf(x,a,b);
%BETACDF Beta cumulative distribution function.
%   P = BETACDF(X,A,B) returns the beta cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINC does the computational work.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.5.

if nargin < 3, 
   error('Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
   error('Requires non-scalar arguments to match in size.');
end

% Initialize P to 0.
p = zeros(size(x));

k1 = find(a<=0 | b<=0);
if any(k1)
   tmp = NaN;
   p(k1) = tmp(ones(size(k1))); 
end

% If is X >= 1 the cdf of X is 1. 
k2 = find(x >= 1);
if any(k2)
   p(k2) = ones(size(k2));
end

k = find(x > 0 & x < 1 & a > 0 & b > 0);
if any(k)
   p(k) = betainc(x(k),a(k),b(k));
end

% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
p(k) = ones(size(k));
