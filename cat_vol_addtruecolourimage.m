function cat_vol_addtruecolourimage(P,cmap)
% ______________________________________________________________________
% Function to show overlays of images either as multiple views or for
% up to 3 overlay images as RGB view
%
% Input:
% P    - char array of 2..12 filenames  
% cmap - colourmap for overlay (default jet)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 

  % transparency
  prop = 0.2;
  RGBstr = {'R','G','B'};
  
  if nargin < 1
    P = spm_select([2 12],'image','Select anatomical image and image(s) to overlay',...
        {fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_T1.nii')});
  elseif ~ischar(P)
    P = spm_select([2 12],'image','Select anatomical image and image(s) to overlay');
  end
    
  V = spm_vol(P);
  n = numel(V);
  rgb_overlay = 0;
  spm_orthviews('Reset');

  if n > 2 && n < 5
    rgb_overlay = spm_input('Overlay',1,'m','Single RGB overlay|Multiple overlays',[true false], 1);  
  end  
  
  if rgb_overlay
    cmap = [((1:64)/64)' zeros(64,2)];
  elseif ~exist('cmap','var')
    cmap = jet(64);
  end
  
  if nargin < 2
    if ~rgb_overlay
      cmap = spm_input('Colourmap',1,'e','jet');
    end
  end

  if rgb_overlay || n == 2
    spm_check_registration(deblank(P(1,:)));
  else
    spm_check_registration(repmat(deblank(P(1,:)),n-1,1));
  end

  
  for i=2:n
    if rgb_overlay
      if i==3, cmap = [zeros(64,1) ((1:64)/64)' zeros(64,1)]; end
      if i==4, cmap = [zeros(64,2) ((1:64)/64)']; end
      handle = 1;
    else
      handle = i - 1;
    end
    [mn,mx] = mn_mx_val(V(i));
    if (mn < 0) && (mx/mn > -100)
      mn = min([-mn mx]);
      mx = mn;
      mn = -mn;
    end
    
    spm_orthviews('addtruecolourimage',handle,deblank(P(i,:)),cmap,prop,mx,mn);
    
    [pp,nam] = spm_fileparts(deblank(P(i,:)));
    if ~rgb_overlay
      spm_orthviews('Caption',i-1,{nam},'FontSize',12,'FontWeight','Bold');
    else
      fprintf('%s: %s\n',RGBstr{i-1},nam)
    end
  end
  
  spm_orthviews('redraw');
end

%_______________________________________________________________________
function [mn,mx] = mn_mx_val(vol)
  mn = Inf;
  mx = -Inf;
  for i=1:vol.dim(3)
    tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
    imx = max(tmp(isfinite(tmp)));
    if ~isempty(imx),mx = max(mx,imx);end
    imn = min(tmp(isfinite(tmp)));
    if ~isempty(imn),mn = min(mn,imn);end
  end
end