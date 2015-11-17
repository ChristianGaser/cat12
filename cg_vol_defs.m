function cg_vol_defs(job)
% Apply deformations to images. In contrast to spm_defs images are saved
% in the original directory.
%_______________________________________________________________________
% Christian Gaser
% $Id$
%
% based on spm_deformations.m

many_images = 0;

try
  PU = job.field;
catch
  PU = job.field1;
  many_images = 1;
end

PI = job.images;
intp = job.interp;

for i=1:numel(PU),

  [pth,nam,ext,num] = spm_fileparts(PU{i});
  PU{i} = fullfile(pth,[nam ext]);

  [Def,mat] = get_def(PU{i});
 
  for m=1:numel(PI)
    
    if many_images % many images
      apply_def(Def,mat,char(PI{m}),intp,job.modulate);
    else % many subjects
      apply_def(Def,mat,char(PI{m}{i}),intp,job.modulate);
    end
  end
end
return

%_______________________________________________________________________
function [Def,mat] = get_def(job)
% Load a deformation field saved as an image
Nii = nifti(job);
Def = single(Nii.dat(:,:,:,1,:));
d   = size(Def);
if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
Def = reshape(Def,[d(1:3) d(5)]);
mat = Nii.mat;

%_______________________________________________________________________
function apply_def(Def,mat,fnames,intrp,modulate)
% Warp an image or series of images according to a deformation field

intrp = [intrp*[1 1 1], 0 0 0];
ofnames = cell(size(fnames,1),1);

for i=1:size(fnames,1),
    V = spm_vol(fnames(i,:));
    M = inv(V.mat);
    [pth,nam,ext,num] = spm_fileparts(deblank(fnames(i,:))); ext = '.nii'; 
    switch modulate
    case 0
        ofnames{i} = fullfile(pth,['w',nam,ext]);
    case 1
        ofnames{i} = fullfile(pth,['mw',nam,ext]);
    case 2
        ofnames{i} = fullfile(pth,['m0w',nam,ext]);
    end
    Vo = struct('fname',ofnames{i},...
                'dim',[size(Def(:,:,:,1),1) size(Def(:,:,:,1),2) size(Def(:,:,:,1),3)],...
                'dt',[spm_type('uint16') spm_platform('bigend')],...
                'pinfo',[1/255 0 0]',...
                'mat',mat,...
                'n',V.n,...
                'descrip',V.descrip);
    ofnames{i} = [ofnames{i} num];
    C  = spm_bsplinc(V,intrp);
    Vo = spm_create_vol(Vo);
    
    if modulate
      dt = spm_diffeo('def2det',Def)/det(mat(1:3,1:3));
      dt(:,:,[1 end]) = NaN;
      dt(:,[1 end],:) = NaN;
      dt([1 end],:,:) = NaN;

      % for modulation of non-linear parts only we have to remove the affine part
      % of the jacobian determinant
      if modulate == 2
        M1  = Vo.mat;
        [x1,x2,x3] = ndgrid(single(1:size(Def,1)),single(1:size(Def,2)),single(1:size(Def,3)));
        X = cat(4,x1,x2,x3);
        Ma = spm_get_closest_affine(X,Def)
        M3 = Ma\M1;
        dt = dt*abs(det(M3));
      end
    
    end
    
    for j=1:size(Def(:,:,:,1),3)
      d0    = {double(Def(:,:,j,1)), double(Def(:,:,j,2)),double(Def(:,:,j,3))};
      d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
      d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
      d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
      dat   = spm_bsplins(C,d{:},intrp);
      if modulate
        dat = dat.*dt(:,:,j);
      end
      Vo    = spm_write_plane(Vo,dat,j);
    end;
end;
return;
