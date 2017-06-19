function Ydiv = cat_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var'), vx_vol = repmat(1.5,1,3); end % no reduction
  Ym = single(Ym); 
  [D,I] = cat_vbdist(single(~isnan(Ym))); Ym = Ym(I); % replace nan
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,min(1.5,vx_vol*3),32); % don't forget small animals...
  [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = cat_vol_resize(Ydivr,'dereduceV',resT2); 
return