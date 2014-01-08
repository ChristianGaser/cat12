function out = vbm_stat_nanstd(in, dim)
% ----------------------------------------------------------------------
% Standard deviation, not considering NaN values. Similar usage like 
% mean() or MATLAB nanmedian of the statistic toolbox.
%
% out nanmean(in,dim)
%
% Example:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = vbm_stat_nanstd(a,3); 
%   am = nanstd(a,0,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
% ----------------------------------------------------------------------
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
% ----------------------------------------------------------------------
% $Id$

  if nargin < 1
    help nanmean;
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
  
  % estimate mean
  tmpin = in;
  tmpin(isnan(in(:))) = 0;
  mn = vbm_stat_nanmean(in,dim);
 
  dm = size(in); dm(setdiff(1:numel(dm),dim)) = 1;
  tmpmn = repmat(mn,dm);
  tmpmn(isnan(in(:))) = 0; 
    
  % estimate std
  out = (sum( (tmpin-tmpmn).^2 , dim) ./ max(1,(size(in,dim) - sum(isnan(in),dim))-1)).^0.5;
  out(isnan(mn))=nan;
end