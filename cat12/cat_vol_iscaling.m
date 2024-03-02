function clim = cat_vol_iscaling(cdata,plim)
% clim = cat_vol_iscaling(cdata,plim). Intensity scaling. 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

  cdata(isnan(cdata) | isinf(cdata))=[]; 
  ASD = min(0.02,max(eps,0.05*std(cdata))/max(abs(cdata))); 
  if ~exist('plim','var'), plim = [ASD 1-ASD]; end 

  bcdata  = [min(cdata) max(cdata)]; 
  if bcdata(1) == bcdata(2)
    clim = bcdata + [-eps eps];
  else
    range   = bcdata(1):diff(bcdata)/1000:bcdata(2);
    hst     = hist(cdata,range);
    clim(1) = range(max(1,find(cumsum(hst)/sum(hst)>plim(1),1,'first')));
    clim(2) = range(min([numel(range),find(cumsum(hst)/sum(hst)>plim(2),1,'first')]));
  end
end