function vbm_tst_reduceRes(P,resi)
% Reduce and Reinterpolate images for resolution test.
% ______________________________________________________________________
% Create further resolution levels of the input images P to test the 
% influence of the resolution for preprocessing and the res_RMS quality 
% measure. Only 0.5 steps for reduced leves are allowed, because finer 
% substeps will lead to interpolation problems. Default for the BWP is 
% therefore 1.5:0.5:3.0 mm.  
% 
% For each resolution 2 images were created, one in reduced and one in
% reinterpolated (sinc interpolation) resolution.
% 
% As far as also the resampling depend strongly on the used grid this 
% script use first an nearst neighbor interpolation to double resolution
% and an than a PVE resampling that estimate the mean of the matrix
% block that form the new voxel. For example the reduction to 1.5 times
% of the original resolutions means an interpolation to 0.5 mm and a
% groupping of 3x3x3 voxel as new voxel. 
% 
% vbm_tst_reduceRes(P,resi)
%  P    = filenames of input images
%  resi = resolution of the image n x 3 matrix
% ______________________________________________________________________
% Robert Dahnke 2014_10
% Structural Brain Mapping Group
% University Jena
%  
% $Id$
% ______________________________________________________________________


if ~exist('P','var') || isempty(P)
  P = cellstr(spm_select(Inf,'image','Select images to reduce')); 
else
  P = cellstr(P);
end
if isempty(P), return; end


% resolutions
if ~exist('resi','var')
  res.res  = (1.5:0.5:3)';
  res.res1 = [res.res res.res res.res];                                   % isotropic
  res.res2 = [res.res res.res ones(numel(res.res),1)];                    % aniso 1 - lower inplane resolution
  res.res3 = [ones(numel(res.res),1) ones(numel(res.res),1) res.res];     % aniso 2 - lower slice resolution
  res.all  = [res.res1;res.res2;res.res3];
else
  if size(resi,2)~=3 
    if size(resi,1)==1
      res.res  = resi';
      res.res1 = [res.res res.res res.res];                                   % isotropic
      res.res2 = [res.res res.res ones(numel(res.res),1)];                    % aniso 1 - lower inplane resolution
      res.res3 = [ones(numel(res.res),1) ones(numel(res.res),1) res.res];     % aniso 2 - lower slice resolution
      res.all  = [res.res1;res.res2;res.res3]; 
    else
      error('Bad resolution matrix resi - has to be a nx3 matrix or a 1xn matrix.'); 
    end
  else
    res.all  = resi;
  end
end
res.all = max(1,round(res.all*2)/2); 
res.all = unique(res.all,'rows');



for pi=1:numel(P)
  for ri=1:size(res.all,1)
    % filename
    [pp,ff,ee] = spm_fileparts(P{pi});
    V          = spm_vol(P{pi});
    vx_vol     = sqrt(sum(V.mat(1:3,1:3).^2));
    vx         = round(vx_vol.*res.all(ri,:)*100);
     
    fplace  = strfind(ff,'_vx'); 
    if isempty(fplace)
      fnamer  = fullfile(pp,sprintf('%s_vx%dx%dx%dr%s',ff,vx,ee));
      fnamei  = fullfile(pp,sprintf('%s_vx%dx%dx%di%s',ff,vx,ee));
    else
      fplace2 = strfind(ff(fplace+1:end),'_');
      if isempty(fplace2) 
        fnamer  = fullfile(pp,sprintf('%s_vx%dx%dx%dr%s',ff(1:fplace-1),vx,ee));
        fnamei  = fullfile(pp,sprintf('%s_vx%dx%dx%di%s',ff(1:fplace-1),vx,ee));
      else
        fplace2 = fplace + fplace2; 
        fnamer  = fullfile(pp,sprintf('%s_vx%dx%dx%dr%s%s',ff(1:fplace-1),vx,ff(fplace2:end),ee));
        fnamei  = fullfile(pp,sprintf('%s_vx%dx%dx%di%s%s',ff(1:fplace-1),vx,ff(fplace2:end),ee));
      end
    end
    fprintf('%s',fnamer)

    if 1
      % load original image
      Y = single(spm_read_vols(V)); 


      % interpolate
      method = 'nearest'; % nearest, linear, or cubic.
      [Rx,Ry,Rz] = meshgrid(single(1-0.25:0.5:size(Y,2)),...
                            single(1-0.25:0.5:size(Y,1)),...
                            single(1-0.25:0.5:size(Y,3)));
      Yi = vbm_vol_interp3f(Y,Rx,Ry,Rz,method); 
      clear Rx Ry Rz Yo Yg;

      % reduce
      si = size(Yi); 
      ss = ceil(si./(res.all(ri,:)*2));

      Yi(end+20,end+20,end+20) = 0;  %#ok<AGROW>

      Yr = zeros(ss,'single'); 
      c  = 0;
      for x=1:res.all(ri,1)*2
        for y=1:res.all(ri,2)*2
          for z=1:res.all(ri,3)*2
            Yr = Yr + Yi(x:res.all(ri,1)*2:si(1)+x,y:res.all(ri,2)*2:si(2)+y,z:res.all(ri,3)*2:si(3)+z);
            c  = c + 1;
          end
        end
      end
      Yr = Yr / c; 
      clear Yi c si ss;
      
          
      % write low resolution file
      Vr        = V;
      Vr        = rmfield(Vr,'private');
      Vr.dim    = size(Yr);
      imat      = spm_imatrix(V.mat); 
      imat(1:3) = imat(1:3) -  (res.all(ri,:)-1)/2;
      imat(7:9) = imat(7:9) .* res.all(ri,:);
      Vr.mat    = spm_matrix(imat);
      Vr.fname  = fnamer;
      Vr        = spm_write_vol(Vr,Yr);
      clear Y Yr imat;
    else
      Vr        = V;
      Vr.dim    = ceil(Vr.dim ./ res.all(ri,:));
      imat      = spm_imatrix(V.mat); 
      imat(7:9) = imat(7:9) .* res.all(ri,:);
      Vr.mat    = spm_matrix(imat);
      Vr.fname  = fnamer;
      Vr        = vbm_vol_imcalc(V,Vr,'i1',struct('interp',6,'verb',0));
    end

    
    % reinterpolate
    Vo = V; Vo.fname = fnamei;
    vbm_vol_imcalc(Vr,Vo,'i1',struct('interp',6,'verb',0));
    fprintf(' done.\n'); 
    clear Vo Vr V;
    
  end
end


