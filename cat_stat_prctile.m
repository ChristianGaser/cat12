function y = cat_stat_prctile(x,p)
% ----------------------------------------------------------------------
% Percentiles of all values of x (ignoring NaNs) like prctile(x(:),p) but
% without the Statistics and Machine Learning Toolbox that is required by
% prctile, quantile and iqr before MATLAB R2022a. As in prctile, the
% sorted values are taken as the 100*(0.5:n-0.5)/n percentiles, values in
% between are linearly interpolated and values outside are set to the
% minimum or maximum.
%
% y = cat_stat_prctile(x,p)
%
%   x .. data array (all elements are used, NaNs are ignored)
%   p .. percentages between 0 and 100
%   y .. percentiles (same size as p, NaN if x has no values; single for
%        single x as in prctile and double otherwise)
%
% Example - interquartile range like iqr(x):
%   diff(cat_stat_prctile(x,[25 75]))
%
% See also cat_stat_nanmedian, cat_stat_nanstat1d.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  x = x(:);
  x = sort(x(~isnan(x)));
  n = numel(x);
  y = nan(numel(p),1);
  if n > 0
    % position of the percentiles in the sorted data limited to [1,n]
    r  = min(max(double(p(:))/100*n + 0.5, 1), n);
    lo = floor(r);
    w  = r - lo;

    % interpolate only between two values (avoids 0*inf for infinite values)
    y    = double(x(lo));
    i    = w > 0;
    y(i) = y(i) + w(i) .* (double(x(lo(i) + 1)) - y(i));
  end
  y = reshape(y,size(p));
  if isa(x,'single'), y = single(y); end
end
