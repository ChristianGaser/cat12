function cg_vbm8_defs(job)
% Apply deformations to images. In contrast to spm_defs images are saved
% in the original directory.
%_______________________________________________________________________
% Christian Gaser
% $Id$

% remove potential file number at the end
[pth,nam,ext,num] = spm_fileparts(job.field{1});
job.field{1} = fullfile(pth,[nam ext]);

[Def,mat] = get_def(job.field);
apply_def(Def,mat,strvcat(job.fnames),job.interp,job.modulate);
return

%_______________________________________________________________________
function [Def,mat] = get_def(job)
% Load a deformation field saved as an image

P      = [repmat(job{:},3,1), [',1,1';',1,2';',1,3']];
V      = spm_vol(P);
Def    = cell(3,1);
Def{1} = spm_load_float(V(1));
Def{2} = spm_load_float(V(2));
Def{3} = spm_load_float(V(3));
mat    = V(1).mat;
%_______________________________________________________________________
function apply_def(Def,mat,fnames,intrp,modulate)
% Warp an image or series of images according to a deformation field

intrp = [intrp*[1 1 1], 0 0 0];

for i=1:size(fnames,1),
    V = spm_vol(fnames(i,:));
    M = inv(V.mat);
    [pth,nam,ext] = spm_fileparts(fnames(i,:));
    if modulate
        ofname = fullfile(pth,['mw',nam,ext]);
    else
        ofname = fullfile(pth,['w',nam,ext]);
    end
    Vo = struct('fname',ofname,...
                'dim',[size(Def{1},1) size(Def{1},2) size(Def{1},3)],...
                'dt',V.dt,...
                'pinfo',V.pinfo,...
                'mat',mat,...
                'n',V.n,...
                'descrip',V.descrip);
    C  = spm_bsplinc(V,intrp);
    Vo = spm_create_vol(Vo);
    if modulate
      dt = spm_def2det(Def{1},Def{2},Def{3},V.mat);
      abs(det(V.mat(1:3,1:3))/det(Vo.mat(1:3,1:3)))
      dt = dt*abs(det(V.mat(1:3,1:3))/det(Vo.mat(1:3,1:3)));

if 0
        x      = affind(rgrid(V.dim(1:3)),V.mat);
        y = zeros([Vo.dim(1:3) 3]);
        for i=1:3
          y(:,:,:,i) = Def{i};
        end
        size(y)
        size(x)
        whos
        y1     = affind(y,Vo.mat);
        
        [M3,R]  = spm_get_closest_affine(x,y1)

      M2 = Vo.mat\M3*V.mat;
      abs(det(M2(1:3,1:3)))
end
    end
    for j=1:size(Def{1},3)
        d0    = {double(Def{1}(:,:,j)), double(Def{2}(:,:,j)),double(Def{3}(:,:,j))};
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

%=======================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
    y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%=======================================================================

