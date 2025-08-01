function Ydiv = cat_vol_div(Ym,vx_vol,vx_volr,norm,sm)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
%
%   Ydiv = cat_vol_div(Ym[,vx_vol,vx_volr])
%  
%   Ym        .. input image
%   vx_vol    .. voxel resolution
%   vx_volr   .. lower voxel resolution for faster processing 
%                and smoother results
%   norm      .. normalize gradients (default = 0)
%   sm        .. smooth (default = 1)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if nargin==0, help cat_vol_div; return; end
  if ~exist('vx_vol','var'), vx_vol = repmat(1.5,1,3); end % no reduction
  if isscalar(vx_vol), vx_vol = repmat(vx_vol,1,3); end
  if ~exist('vx_volr','var'), vx_volr = min(1.5,vx_vol*3); end
  if isscalar(vx_volr), vx_volr = repmat(vx_volr,1,3); end
  if ~exist('norm','var'), norm = 0; end
  if ~exist('sm','var'), sm = 1; end
  
  Ym   = single(Ym); 
  Ynan = isnan(Ym); 
  
  % replace nan by neighbor values
  [D,I] = cat_vbdist(single(~Ynan),cat_vol_morph(~Ynan,'d',2)); Ym(D<2) = Ym(I(D<2)); 
  clear D I
  
  % limit brain (nan or zero) and use low resolution for processing speed... % don't forget small animals...
  [Ym,BB]     = cat_vol_resize(Ym,'reduceBrain', vx_vol, max(vx_vol/vx_volr)*2, Ym~=0 & ~Ynan); 
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV', vx_vol, vx_volr, 16, 'cubic'); 
  clear Ym
  
  % gradients
  if norm == 1 
    [gx,gy,gz]  = cat_vol_gradient3(Ymr); clear Ymr
    gs = max(1e-1,(gx.^2 + gy.^2 + gz.^2).^0.5); 
    gx = gx ./ gs; gy = gy ./ gs; gz = gz ./ gs; 
  else
    [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); clear Ymr
  end
  
  % divergence function was too memory demanding for some systems
  [px,~,~] = cat_vol_gradient3(gx./resT2.vx_volr(1)); clear gx 
  [~,py,~] = cat_vol_gradient3(gy./resT2.vx_volr(2)); clear gy 
  [~,~,pz] = cat_vol_gradient3(gz./resT2.vx_volr(3)); clear gz 
  Ydivr = single(px) + single(py) + single(pz); clear px qy gz
  
  if sm, Ydivr = smooth3(Ydivr); end
  
  % original resolution
  Ydiv  = cat_vol_resize(Ydivr, 'dereduceV', resT2, 'cubic'); 
  Ydiv  = cat_vol_resize(Ydiv,  'dereduceBrain', BB); 
  
  Ydiv(Ynan) = 0;  % restore nan
return