function cat_vol_laterality_index(P)
% ______________________________________________________________________
% Calculation of laterality index for images in MNI152NLin2009cAsym space
% LI = (L-R)/(R+L)
%
% The result is indicated with a prepended 'LI_' in the dataname
% of the file.
% Please note that only the data of the left hemipshere is stored, since
% the values in the opposite hemisphere would be simply inverted and other-
% wise identical except for the sign.
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

fprintf('Warning: This only works with spatially registered images in MNI152NLin2009cAsym space (CAT12.8 or newer)\n');

if ~nargin
  P = spm_select(Inf,'image','Select images for LI estimation',{},pwd,'^T*');
end

V = spm_vol(P);
n = size(P,1);

sym_template = '/Users/gaser/Dropbox/GitHub/cat12-templates-atlases/templates_external/y_MNI152NLin2009cAsym_to_MNI152NLin2009cSym.nii';
Vsym = spm_vol(sym_template);
wvol = cat_vol_defs(struct('field1',{{sym_template}},'images',{{P}},'interp',5,'modulate',0));

for i = 1:n
  vol = wvol{1}{i};
  
  % estimate new dimensions for left/right hemisphere
  xdim  = size(vol,1);
  left_xdim = floor(xdim/2);

  % consider off x-dimensions
  if rem(xdim,2)
    right_xdim = ceil(xdim/2) + 1;
  else
    right_xdim = ceil(xdim/2);
  end
  
  % flip values
  left_data = vol(1:left_xdim,:,:);
  right_data  = flipud(vol(right_xdim:xdim,:,:)); % image should be flipped

  % estimate laterality index
	LI = (left_data-right_data)./(left_data+right_data+eps);
	  
  % rename dataname
  [pth,nm,xt] = spm_fileparts(deblank(P(i,:)));
  flipped_name = fullfile(pth, ['LI_' nm xt]);
  
  % we need origin and voxel size from template
  Vout = Vsym;
  Vout.fname = flipped_name;
  Vout.dim = size(LI);
  Vout.dt(1) = 16;
  Vout.pinfo(1) = 1;
  Vout.mat(1,:) = -Vout.mat(1,:);
  Vout.mat(1,4) = Vout.mat(1,4) - 1.5;

  Vout.descrip = 'Laterality index';
  
  spm_write_vol(Vout,left_data);
  fprintf('Save LI in %s\n',flipped_name);
end