function out = vbm_stat_nanmean(in, dim)
% ----------------------------------------------------------------------
% Average, not considering NaN values. Similar usage like mean() or 
% MATLAB nanmedian of the statistic toolbox.
%
% out nanmean(in,dim)
%
% Example:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = vbm_stat_nanmean(a,3); 
%   am = nanmean(a,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
% ----------------------------------------------------------------------
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
% ----------------------------------------------------------------------
% $Id$

  if nargin < 1
    help vbm_stat_nanmean;
    return;
  end;
  
  if nargin < 2
    if size(in,1) ~= 1
      dim = 1;
    elseif size(in,2) ~= 1
      dim = 2;
    else 
      dim = 3; 
    end;
  end;
  
  if isempty(in), out = nan; return; end
  
  % estimate mean
  tmpin = in;
  tmpin(isnan(in(:))) = 0;
  if sum(~isnan(in),dim)==0
    out = nan(size(sum(tmpin, dim)));
  else
    out = sum(tmpin, dim) ./ max(eps,sum(~isnan(in),dim));
  end
end