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
apply_def(Def,mat,strvcat(job.fnames),job.interp);
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
function apply_def(Def,mat,fnames,intrp)
% Warp an image or series of images according to a deformation field

intrp = [intrp*[1 1 1], 0 0 0];

for i=1:size(fnames,1),
    V = spm_vol(fnames(i,:));
    M = inv(V.mat);
    [pth,nam,ext] = spm_fileparts(fnames(i,:));
    ofname = fullfile(pth,['w',nam,ext]);
    Vo = struct('fname',ofname,...
                'dim',[size(Def{1},1) size(Def{1},2) size(Def{1},3)],...
                'dt',V.dt,...
                'pinfo',V.pinfo,...
                'mat',mat,...
                'n',V.n,...
                'descrip',V.descrip);
    C  = spm_bsplinc(V,intrp);
    Vo = spm_create_vol(Vo);
    for j=1:size(Def{1},3)
        d0    = {double(Def{1}(:,:,j)), double(Def{2}(:,:,j)),double(Def{3}(:,:,j))};
        d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
        d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
        d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
        dat   = spm_bsplins(C,d{:},intrp);
        Vo    = spm_write_plane(Vo,dat,j);
    end;
end;
return;

