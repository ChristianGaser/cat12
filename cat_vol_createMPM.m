function cat_vol_create_MPM(Label, Deform, vox, thresholds)
% Create Maximum Probability Map (Label) of labels in native space and deformation 
% fields
%
% FORMAT cat_vol_create_MPM(Label,Deform)
% Label   - char array of filenames of labels
% Deform  - char array of filenames of deformation fields (leave empty if labels
%           are already normalized and no deformations are needed)
% vox     - voxel size (use NaNs to use voxel size of deformation fields)
%
% some of the subfunctions are modified versions from spm_deformations.m
%__________________________________________________________________________
% Christian Gaser
% $Id$

refine     = 1;   % always use refinement with slight smoothing and median filtering

if nargin < 1
  Label  = spm_select(Inf,'image','Select native label maps');
end

Vlabel = spm_vol(Label);
n_subjects = numel(Vlabel);

if nargin < 2
  Deform  = spm_select([0 n_subjects],'image','Select deformation fields (or press done if no deformations are needed)',{},pwd,'^y_');
end

% voxel size
if nargin < 3 & ~isempty(Deform)
  vox = spm_input('Voxel size','1','r',[NaN NaN NaN],[1,3]);
end

% thresholds for average probability to exclude non-brain areas
if nargin < 4
  thresholds = spm_input('Threshold(s)','2','r',[0.1:0.1:0.5]);
end

% check whether only one value was defined
if size(vox,1) ==1 & size(vox,2) == 1
  vox = [vox vox vox];
end

% transpose if necessary
if size(vox,1) > size(vox,2)
  vox = vox';
end

% find all unique values in structures
structures = spm_read_vols(Vlabel(1));
datarange  = 0:max(structures(:));
H = hist(structures(:),datarange);
ind = [1 find(H==0)];
datarange(ind) = [];

n_structures = numel(datarange);

if ~isempty(Deform)
  Vdeform = spm_vol(Deform);
  V = Vdeform;
  [Def,mat] = get_comp(Vdeform(1).fname,vox);
  sz = size(Def);
  V(1).dim(1:3) = sz(1:3);
  V(1).mat = mat;
else
  V = Vlabel;
end

% set data type w.r.t. maximum value
max_val = max(datarange);
if max_val < 2^8
  data_type   = 'uint8';
  fprintf('Set data type to uint8\n.')
elseif max_val < 2^16
  data_type   = 'uint16';
  fprintf('Set data type to uint16\n.')
else 
  data_type   = 'float32';
  fprintf('Set data type to float32\n.')
end

[tmp, name] = spm_str_manip(spm_str_manip(Label,'t'),'C');

watlas = zeros([V(1).dim(1:3) n_structures],'single');

for i=1:n_subjects
  fprintf('.');
  if ~isempty(Deform)
    [Def,mat] = get_comp(Vdeform(i).fname,vox);    
  end
  vol = spm_read_vols(Vlabel(i));
  for j=1:n_structures;
    if ~isempty(Deform)
      dat = apply_def(Def,double(round(vol)==datarange(j)),3,Vlabel(i).mat);
    else
      dat = double(round(vol)==datarange(j));
    end
    dat(isnan(dat)) = 0;
    watlas(:,:,:,j) = watlas(:,:,:,j) + single(dat);
  end
end
fprintf('\n');

% apply median filtering to each label and slight smoothing
if refine
  for j=1:n_structures;
    tmp = cat_vol_median3c(single(watlas(:,:,:,j)));
    spm_smooth(tmp,tmp,2);
    watlas(:,:,:,j) = tmp;
  end
end

[max_atlas, index_atlas_orig] = max(watlas,[],4);
avg_atlas = sum(watlas,4)/n_subjects;

for i=1:numel(thresholds)
  threshold = thresholds(i);
  
  index_atlas = index_atlas_orig;
  index_atlas(find(max_atlas<1 | isnan(max_atlas) | avg_atlas<threshold)) = 0;
  
  index_atlas0 = index_atlas;
  
  for j=1:n_structures;
    index_atlas(find(index_atlas0 == j)) = datarange(j);
  end
  
  % replace remaining holes with median value
  holes = index_atlas > 0;
  holes = (cat_vol_morph(holes,'c') - holes) > 0;
  tmp = cat_vol_median3c(single(index_atlas));
  index_atlas(holes) = double(tmp(holes));
  
  Vo = struct('fname',['MPM_th' sprintf('%0.2f',threshold) '_' name.s strrep(name.e,',1','')],...
              'dim',size(index_atlas),...
              'dt',[spm_type(data_type)  spm_platform('bigend')],...
              'pinfo',[1 0 352]',...
              'mat',V(1).mat,...
              'n',V(1).n,...
              'descrip',['n=' num2str(n_subjects)]);
  Vo = spm_create_vol(Vo);
  spm_write_vol(Vo, index_atlas);
  fprintf('%s saved.\n',Vo.fname);
end

%_______________________________________________________________________
function [Def,mat,vx,bb] = get_def(field)
% Load a deformation field saved as an image
Nii = nifti(field);
Def = single(Nii.dat(:,:,:,1,:));
d   = size(Def);
if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
Def = reshape(Def,[d(1:3) d(5)]);
mat = Nii.mat;

vx  = sqrt(sum(Nii.mat(1:3,1:3).^2));
if det(Nii.mat(1:3,1:3))<0, vx(1) = -vx(1); end

o   = Nii.mat\[0 0 0 1]';
o   = o(1:3)';
dm  = size(Nii.dat);
bb  = [-vx.*(o-1) ; vx.*(dm(1:3)-o)];

%_______________________________________________________________________
function Def = identity(d,M)
[y1,y2]   = ndgrid(single(1:d(1)),single(1:d(2)));
Def       = zeros([d 3],'single');
for y3=1:d(3)
    Def(:,:,y3,1) = y1*M(1,1) + y2*M(1,2) + (y3*M(1,3) + M(1,4));
    Def(:,:,y3,2) = y1*M(2,1) + y2*M(2,2) + (y3*M(2,3) + M(2,4));
    Def(:,:,y3,3) = y1*M(3,1) + y2*M(3,2) + (y3*M(3,3) + M(3,4));
end

%_______________________________________________________________________
function [Def,mat] = get_comp(field,vox)
% Return the composition of two deformation fields.

[Def,mat,vx,bb] = get_def(field);

% only estimate composite if job field is given
if nargin > 1
  % only move on if any vox or bb field is not NaN
  if any(isfinite(vox))
    Def1         = Def;
    mat1         = mat;
    vox(~isfinite(vox)) = vx(~isfinite(vox));
    
    [mat, dim]   = spm_get_matdim('', vox, bb);
    Def          = identity(dim, mat);
    M            = inv(mat1);
    tmp          = zeros(size(Def),'single');
    tmp(:,:,:,1) = M(1,1)*Def(:,:,:,1)+M(1,2)*Def(:,:,:,2)+M(1,3)*Def(:,:,:,3)+M(1,4);
    tmp(:,:,:,2) = M(2,1)*Def(:,:,:,1)+M(2,2)*Def(:,:,:,2)+M(2,3)*Def(:,:,:,3)+M(2,4);
    tmp(:,:,:,3) = M(3,1)*Def(:,:,:,1)+M(3,2)*Def(:,:,:,2)+M(3,3)*Def(:,:,:,3)+M(3,4);
    Def(:,:,:,1) = single(spm_diffeo('bsplins',Def1(:,:,:,1),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,2) = single(spm_diffeo('bsplins',Def1(:,:,:,2),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,3) = single(spm_diffeo('bsplins',Def1(:,:,:,3),tmp,[1,1,1,0,0,0]));
    clear tmp
  end
end

%_______________________________________________________________________
function out = apply_def(Def,vol,intrp,mat)
% Warp an image or series of images according to a deformation field
intrp = [intrp*[1 1 1], 0 0 0];
M = inv(mat);

C  = spm_bsplinc(vol,intrp);
dat = [];
for j=1:size(Def,3)
    d0    = {double(Def(:,:,j,1)), double(Def(:,:,j,2)), double(Def(:,:,j,3))};
    d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
    d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
    d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
    dat   = [dat spm_bsplins(C,d{:},intrp)];
end

sz = size(Def);
out = reshape(dat,sz(1:3));
