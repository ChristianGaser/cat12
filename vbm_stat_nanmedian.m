function out = vbm_stat_nanmedian(in, dim)
% ----------------------------------------------------------------------
% Median, not considering NaN values. Similar usage like median() or 
% MATLAB nanmedian of the statistic toolbox.
%
% out = vbm_stat_nanmedian(in,dim)
%
% Example:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = vbm_stat_nanmedian(a,3); 
%   am = nanmedian(a,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
% ----------------------------------------------------------------------
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
% ----------------------------------------------------------------------
% $Id$

  if nargin < 1
    help vbm_stat_nanmedian;
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
  
  sz  = size(in);

  if isempty(in), out = nan; return; end
  
  % reduce to 2d matrix
  pm = [dim:max(length(size(in)),dim) 1:dim-1];
  in = reshape(permute(in,pm),size(in,dim),prod(sz)/size(in,dim));

  in = sort(in,1);
  s  = size(in,1) - sum(isnan(in));
  
  % estimate median in loop
  out = zeros(size(s));
  for i = 1:length(s)
    if s(i)>0, out(i) = vbm_stat_nanmean([in(floor((s(i)+1)/2),i),in(ceil((s(i)+1)/2),i)]); else out(i)=nan; end
  end
  
  % estimate median as matrix ... doesn't work :/
  %out = nan(size(s)); si=1:numel(s);
  %out(s>0) = nanmean( [ in(([floor((s(s>0)+1)/2);si(s>0)])'); in(ceil((s(s>0)+1)/2),si(s>0)) ] );
  
  % correct for permutation
  sz(dim) = 1; out = ipermute(reshape(out,sz(pm)),pm);
end