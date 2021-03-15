function vol = cat_vol_defs(job)
% Apply deformations to images. In contrast to spm_deformations images are saved
% in the original directory.
% FORMAT vol = cat_vol_defs(job)
% job - a job created via cat_conf_tools.m
% vol - cell of deformed output volumes (deformed images are only returned, but not saved as file)
%_______________________________________________________________________
% Christian Gaser
% $Id$
%
% based on spm_deformations.m

many_images = 0;

if isfield(job,'field')
  PU = job.field;
else
  PU = job.field1;
  many_images = 1;
end

def.verb = cat_get_defaults('extopts.verb'); 
job      = cat_io_checkinopt(job,def);

PI     = job.images;
interp = job.interp;

if interp < 0 && job.modulate
  warning('Modulation in combination with categorical interpolation is not meaningful and possible. Disable modulation.');
  job.modulate = 0;
end

if nargout == 1
  vol = cell(numel(PU),numel(PI));
end

for i=1:numel(PU)

  % external call with PU as deformation field
  if isnumeric(PU{i}) & isfield(job,'mat')
    Def = PU{i};
    mat = job.mat;
  else
    [pth,nam,ext] = spm_fileparts(PU{i});
    PU{i} = fullfile(pth,[nam ext]);
  
    if isfield(job,'vox') & isfield(job,'bb')
      [Def,mat] = get_comp(PU{i},job);
    else
      [Def,mat] = get_comp(PU{i});
    end
  end
  
  for m=1:numel(PI)
    
    if many_images % many images
      PIi = char(PI{m}); 
    else % many subjects
      PIi = char(PI{m}{i}); 
    end
    
    if nargout == 1
      [PIri, vol{i,m}] = apply_def(Def,mat,PIi,interp,job.modulate,job.verb);
    else
      PIri = apply_def(Def,mat,PIi,interp,job.modulate,job.verb);
    end
    
    if job.verb & ~isempty(PIri)
      fprintf('Display resampled %s\n',spm_file(PIri,'link','spm_image(''Display'',''%s'')'));
    end
  end
end
return

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
function [Def,mat] = get_comp(field,job)
% Return the composition of two deformation fields.

[Def,mat,vx,bb] = get_def(field);

% only estimate composite if job field is given
if nargin > 1
  % only move on if any vox or bb field is not NaN
  if any(isfinite(job.vox)) | any(isfinite(job.bb))
    Def1         = Def;
    mat1         = mat;
    job.vox(~isfinite(job.vox)) = vx(~isfinite(job.vox));
    job.bb(~isfinite(job.bb))   = bb(~isfinite(job.bb));
    
    [mat, dim]   = spm_get_matdim('', job.vox, job.bb);
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
function [out, wvol] = apply_def(Def,mat,filenames,interp0,modulate,verb)
% Warp an image or series of images according to a deformation field

interp = [interp0*[1 1 1], 0 0 0];
dim    = size(Def);
dim    = dim(1:3);
if nargout == 2
    wvol = cell(size(filenames,1),1); 
    out = '';
end

for i=1:size(filenames,1)

    % Generate headers etc for output images
    %----------------------------------------------------------------------
    [pth,nam,ext] = spm_fileparts(deblank(filenames(i,:)));  %#ok<ASGLU>
    NI = nifti(fullfile(pth,[nam ext]));
    j_range = 1:size(NI.dat,4);
    k_range = 1:size(NI.dat,5);
    l_range = 1:size(NI.dat,6);
    
    if nargout < 2
      NO = NI;
      ext = '.nii'; 
      
      % use float for modulated images
      if modulate
          NO.dat.scl_slope = 1.0;
          NO.dat.scl_inter = 0.0;
          NO.dat.dtype     = 'float32-le';
      end
      
      % set slope to 1 for categorical interpolation
      if interp0 < 0
          NO.dat.scl_slope = 1.0;
          NO.dat.scl_inter = 0.0;
          
          % select data type w.r.t. maximum value
          f0  = single(NI.dat(:,:,:,:,:,:));
          max_val = max(f0(:)); clear f0
          if max_val < 2^8
            NO.dat.dtype   = 'uint8-le';
            if verb, fprintf('Set data type to uint8\n'); end
          elseif max_val < 2^16
            NO.dat.dtype   = 'uint16-le';
            if verb, fprintf('Set data type to uint16\n'); end
          else 
            NO.dat.dtype   = 'float32-le';
            if verb, fprintf('Set data type to float32\n'); end
          end
      end
  
      NO.dat.dim     = [dim NI.dat.dim(4:end)];
      NO.dat.offset  = 0; % For situations where input .nii images have an extension.
      NO.mat         = mat;
      NO.mat0        = mat;
      NO.mat_intent  = 'Aligned';
      NO.mat0_intent = 'Aligned';
  
      switch modulate
      case 0
          NO.dat.fname = fullfile(pth,['w',nam,ext]);
          NO.descrip   = sprintf('Warped');
      case 1
          NO.dat.fname = fullfile(pth,['mw',nam,ext]);
          NO.descrip   = sprintf('Warped & Jac scaled');
      case 2
          NO.dat.fname = fullfile(pth,['m0w',nam,ext]);
          NO.descrip   = sprintf('Warped & Jac scaled (nonlinear only)');
      end
      out = NO.dat.fname; 
      
      NO.extras      = [];
      create(NO);
    end
    
    if modulate
      dt = spm_diffeo('def2det',Def)/det(mat(1:3,1:3));
      dt(:,:,[1 end]) = NaN;
      dt(:,[1 end],:) = NaN;
      dt([1 end],:,:) = NaN;

      % for modulation of non-linear parts only we have to remove the affine part
      % of the jacobian determinant
      if modulate == 2
        [x1,x2,x3] = ndgrid(single(1:size(Def,1)),single(1:size(Def,2)),single(1:size(Def,3)));
        X = cat(4,x1,x2,x3);
        Ma = spm_get_closest_affine(X,Def);
        M3 = Ma\mat;
        dt = dt*abs(det(M3));
      end
    
    end
    
    if nargout == 2, wvol{i} = zeros([dim NI.dat.dim(4:end)]); end
    for j=j_range

        M0 = NI.mat;
        if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
            M1 = NI.extras.mat;
            if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                M0 = M1(:,:,j);
            end
        end
        M  = inv(M0);
        % Generate new deformation (if needed)
        Y     = affine(Def,M);
        % Write the warped data for this time point
        %------------------------------------------------------------------
        for k=k_range
            for l=l_range
                f0  = single(NI.dat(:,:,:,j,k,l));
                if interp0>=0
                    c  = spm_diffeo('bsplinc',f0,interp);
                    f1 = spm_diffeo('bsplins',c,Y,interp);
                else
                    % Warp labels
                    U  = unique(f0(:));
                    if numel(U)>1000
                        error('Too many label values.');
                    end
                    f1   = zeros(dim(1:3),class(f0));
                    p1   = zeros(size(f1),'single');
                    for u=U'
                        g0       = single(f0==u);
                        tmp      = spm_diffeo('bsplins',g0,Y,[abs(interp(1:3)) interp(4:end)]);
                        msk1     = (tmp>p1);
                        p1(msk1) = tmp(msk1);
                        f1(msk1) = u;
                    end
                end
                if modulate
                  f1 = f1.*double(dt);
                end
                
                if nargout < 2
                  NO.dat(:,:,:,j,k,l) = f1;
                else
                  wvol{i}(:,:,:,j,k,l) = f1; 
                end
            end
        end
    end

end;
return;

%==========================================================================
% function Def = affine(y,M)
%==========================================================================
function Def = affine(y,M)
Def          = zeros(size(y),'single');
Def(:,:,:,1) = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def(:,:,:,2) = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def(:,:,:,3) = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);
