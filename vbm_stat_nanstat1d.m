function x=vbm_stat_nanstat1d(x,action)
% replace nan* functions of the stat toolbox for the 1d case 
  x=x(:); x(isnan(x) | isinf(x))=[];
  if ~exist('action','var'), action='nanmean'; end

  switch lower(action)
    case {'mean'   'nanmean'},   x = mean(x);
    case {'median' 'nanmedian'}, x = median(x);
    case {'std'    'nanstd'},    x = std(x);
    case {'min'    'nanmin'},    x = min(x);
    case {'max'    'nanmax'},    x = max(x);
    case {'var'    'nanvar'},    x = var(x);
    case {'sum'    'nansum'},    x = sum(x);
    case {'prod'   'nanprod'},   x = prod(x);
    case {'cumsum' 'nancumsum'}, x = cumsum(x); 
  end
end

