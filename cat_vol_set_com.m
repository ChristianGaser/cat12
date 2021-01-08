function Affine = cat_vol_set_com(V)
% use center-of-mass (COM) to roughly correct for differences in the
% position between image and template
% ______________________________________________________________________
% FORMAT:  Affine = cat_vol_set_com(varargin)
%
% V      - mapped images or filenames 
% Affine - affine transformation to roughly correct origin 
% 
% If no output filed is defined the estimated transformation is applied 
% to the images
% ______________________________________________________________________
% Christian Gaser
% $Id$


if nargin == 1
  if isstruct(V)
    V = V;
  else
    P = char(V);
    V = spm_vol(P);
  end
else
  P = spm_select(Inf,'image','Select images to filter');
  V = spm_vol(P);
end
n = numel(V);

% pre-estimated COM of MNI template
com_reference = [0 -20 -15];

for i=1:n
  fprintf('Correct center-of-mass: %s',spm_str_manip(V(i).fname,'k41'));
  Affine = eye(4);
  if isfield(V(i),'dat')
    vol(:,:,:) = V(i).dat(:,:,:);
  else
    vol = spm_read_vols(V(i));
  end
  avg = mean(vol(:));
  avg = mean(vol(vol>avg));
  
  % don't use background values
  [x,y,z] = ind2sub(size(vol),find(vol>avg));
  com = V(i).mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
  com = com';

  M = spm_get_space(V(i).fname);
  Affine(1:3,4) = (com - com_reference)';
  if nargin < 1
    spm_get_space(V(i).fname,Affine\M);
    fprintf('\n');
  end
end
