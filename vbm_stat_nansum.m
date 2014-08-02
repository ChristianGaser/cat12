function out = vbm_stat_nansum(in, dim)
% ----------------------------------------------------------------------
% Average, not considering NaN values. Similar usage like sum() or 
% MATLAB nansum of the statistic toolbox.
%
% out nansum(in,dim)
%
% Example:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = vbm_stat_nansum(a,3); 
%   am = nansum(a,3); % of the statistical toolbox ...
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
  
  if isempty(in), out = 0; return; end
  
  % estimate mean
  tmpin = in;
  tmpin(isnan(in(:))) = 0;
  out = sum(tmpin, dim);
end