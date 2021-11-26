function Yg = cat_vol_grad(Ym,vx_vol,method,repNaN)
% ----------------------------------------------------------------------
% Simple gradient map for edge description. Default method is the absolute
% sum of all vectors but this should be replaced in future by the length 
% of the gradient (RD202111). The function (temporary) replaces NaN to 
% reduce boundary effects.
%
%   Yg = cat_vol_grad(Ym,vx_vol,method)
%
%   Ym      .. input image
%   Yg      .. gradient map 
%   vx_vol  .. voxel size (default=[1 1 1]);
%   method  .. averaging method (0-abs sum, 1-sum, 2-grad legth, default=0)
%   repNaN  .. replace NaN (0-no, 1-permanent, 2-temporary, default=1)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('vx_vol','var'), vx_vol = ones(1,3); end
  if ~exist('method','var'), method = 1; end
  if ~exist('repNaN','var'), repNaN = 1; end
  if numel(vx_vol) == 1
    vx_vol(2:3) = vx_vol(1); 
  elseif numel(vx_vol) ~= 3
    error('cat_vol_grad:vx_vol','The size of the second input (vx_vol) has to be 1 or 3.\n'); 
  end
  
  Ym   = single(Ym); 
 
  % temporary replace NaNs by nearest values
  if repNaN
    Ynan  = isnan(Ym); 
    [D,I] = cat_vbdist(single(~Ynan),cat_vol_morph(~Ynan,'d',1)); Ym(D<2) = Ym(I(D<2)); % replace nan
    clear D I; 
  end
  
  [gx,gy,gz] = cat_vol_gradient3(single(Ym)); 
 
  % averaging 
  switch method
    case 0, Yg = gx./vx_vol(1) + gy./vx_vol(2) + gz./vx_vol(3);                 % simple sum
    case 1, Yg = abs(gx./vx_vol(1)) + abs(gy./vx_vol(2)) + abs(gz./vx_vol(3));  % absolute sum
    case 2, Yg = (gx.^2 + gy.^2 + gz.^2).^(0.5);                                % gradient length
  end
  
  % restore nan 
  if repNaN==1
    Yg(Ynan) = 0;
  elseif repNaN==2
    Yg(Ynan) = nan;
  end
return
