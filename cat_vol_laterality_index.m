function cat_vol_laterality_index(P)
% ______________________________________________________________________
% Calculation of laterality index for images in MNI152NLin2009cAsym space
% LI = (L-R)/(R+L)
%
% We apply the deformation from MNI152NLin2009cAsym to the symmetrical template
% MNI152NLin2009cSym in order to obtain a symmetrical image, which then can be
% used to estimate laterality index.
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
  P = spm_select(Inf,'image','Select images for LI estimation',{},pwd);
end

Def_sym = fullfile(cat_get_defaults('extopts.pth_templates'),'y_MNI152NLin2009cAsym_to_MNI152NLin2009cSym.nii');
Vsym = spm_vol(Def_sym);

n = size(P,1);

cat_progress_bar('Init',n);

for i = 1:n
  name = deblank(P(i,:));
  
  % apply deformation to selected images
  wvol = cat_vol_defs(struct('field1',{{Def_sym}},'images',{{name}},'interp',5,'modulate',0));

  vol = wvol{1}{1};
  
  % estimate new dimensions for left/right hemisphere
  xdim  = size(vol,1);
  right_xdim = floor(xdim/2);

  % consider odd x-dimensions
  if rem(xdim,2)
    left_xdim = ceil(xdim/2) + 1;
  else
    left_xdim = ceil(xdim/2);
  end
  
  % flip values
  left_data  = vol(left_xdim:xdim,:,:); % image should be flipped
  right_data = flipud(vol(1:right_xdim,:,:));
  mx = max(vol(:));

  % estimate laterality index
	LI = (left_data-right_data)./(left_data+right_data+eps);
	
	% and exclude areas where intensity (i.e. of tissue) is quite low
	LI(abs(left_data+right_data)<0.01*mx) = 0;

  % rename dataname
  [pth,nm,xt] = spm_fileparts(deblank(P(i,:)));
  flipped_name = fullfile(pth, ['LI_' nm xt]);
  
  % we need origin and voxel size from template
  Vout = Vsym;
  Vout.fname = flipped_name;
  Vout.dim = size(LI);
  Vout.dt(1) = 16;
  Vout.pinfo(1) = 1;
  Vout.descrip = 'Laterality index (L-R)/(R+L)';

  % we have to correct origin
  vx_vol  = sqrt(sum(Vout.mat(1:3,1:3).^2));
  Vout.mat(1,4) = Vout.mat(1,4) - vx_vol(1)*left_xdim;
  Vout.private.mat(1,4) = Vout.private.mat(1,4) - vx_vol(1)*left_xdim;

  spm_write_vol(Vout,LI);
  fprintf('Save LI in %s\n',flipped_name);
  cat_progress_bar('Set',i);
end

cat_progress_bar('Clear');
