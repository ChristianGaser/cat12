function Atlas_correct_scaling
% Correct scaling factors for atlases to 1
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

csv_file = spm_select('FPList',cat_get_defaults('extopts.pth_templates'),'.csv');

for i = 1:size(csv_file,1)
  [pth,nam,ext] = spm_fileparts(deblank(csv_file(i,:)));
  atlas_file = fullfile(pth,[nam '.nii']);
  N = nifti(atlas_file);
  atlas = N.dat(:,:,:);
  max(atlas(:))
  N.dat.scl_slope = 1;
  N.dat.scl_inter = 0;
  create(N);
  N.dat(:,:,:) = round(atlas);
  fprintf('Save corrected atlas %s\n',atlas_file);
end