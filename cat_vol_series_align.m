function out = cat_vol_series_align(job)
% Longitudinal rigid registration of image series
% FORMAT out = cat_vol_series_align(job)
%_______________________________________________________________________
%
% modified version of
% John Ashburner
% spm_series_align.m 5044 2012-11-09 13:40:35Z john
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

N = numel(job.data);

if numel(job.noise)==1
  noise = repmat(job.noise,[N,1]);
elseif numel(job.noise) ~= N
  error('Incompatible numbers of noise estimates and scans.');
else
  noise = job.noise(:);
end

prec   = noise.^(-2);

if isfield(job,'reg') && isfield(job.reg,'nonlin')
  cat_io_cprintf('blue','Non-linear Registration!\n');
  tim = job.reg.nonlin.times(:);
  if all(isfinite(tim))
    if numel(tim) == 1
        tim = (1:N)';
    elseif numel(tim) ~= N
        error('Incompatible numbers of times and scans.');
    end
    if any(abs(diff(tim)) > 50)
        error('Time differences should be in years.');
    end
    wparam0   = job.reg.nonlin.wparam;
    
    midtim = median(tim);
    tim    = tim - midtim;
    w_settings = kron(wparam0,1./(abs(tim)+1/365));
    s_settings = round(3*abs(tim)+2);
  else % use default regularization if tim is set to NAN
    w_settings = job.reg.nonlin.wparam;
    s_settings = 6; %round( job.reg.nonlin.wparam(5) / 25);
  end
else
  w_settings = [Inf Inf Inf Inf Inf];
  s_settings = Inf;
end

b_settings = [0 0 job.bparam];
Nii = nifti(strvcat(job.data));
ord = [3 3 3 0 0 0];

output = {};
if job.write_avg,  output = [output, {'wavg'}]; end
if job.write_rimg, output = [output, {'wimg'}]; end

if isfield(job.reg,'nonlin')
  if job.reg.nonlin.write_jac, output = [output, {'wjac'} ]; end
  if job.reg.nonlin.write_def, output = [output, {'wdef'} ]; end
end

if ~isfield(job,'use_brainmask')
  use_brainmask = 1;
else
  use_brainmask = job.use_brainmask;
end

if ~isfield(job,'reduce')
  reduce = 1;
else
  reduce = job.reduce;
end

if ~isfield(job,'sharpen')
  sharpen = 2;
else
  sharpen = job.sharpen;
end

if ~isfield(job,'setCOM')
  setCOM = 1;
else
  setCOM = job.setCOM;
end

% force isotropic average resolution (0-default,1-best,2-worst,3-optimal)
if ~isfield(job,'isores')
  isores = 0;
else
  isores = -job.isores;
end

% RD202510: Test voxel resolution and force minimal resolution 
%           this is relevant in case with high slice thickness 
%           for instance in the ultra-low-field data (Frantisek et al. 2025)
%             https://openneuro.org/datasets/ds006557/versions/1.0.0
% RD202511: Moreover, this would also allow superinterpolation in case of 
%           many rescans (similar timepoint)
minres   = min([ 1.2 abs( isores(isores~=0) ) ]); % 
vx_vol   = nan(numel(Nii),3); 
tempimgs = cell(numel(Nii),1); 
for vi = 1:numel(Nii)
  vx_vol(vi,:) = sqrt(sum(Nii(vi).mat(1:3,1:3).^2));

  %%
  if sharpen || any( minres  < vx_vol(vi,:) * .8 )
  % only for average estimation and not long
    Vx = spm_vol(job.data{vi}); 
    Yx = spm_read_vols(Vx);
  end

  if sharpen
    con   = abs(diff([ prctile(Yx(:),10) prctile(Yx(:),90) ]));  
    sharp = min(.5,sharpen * min(con, numel(Nii).^1.5 / 20));
    Yx = min(max(Yx(:)),max(min(Yx(:)),Yx + ... 
      (Yx~=0) .* max(-con/6.*Yx,min(con/6.*Yx, ...
      sharp/2 .* (Yx - smooth3(Yx)) + ...
      sharp/4 .* (Yx - cat_vol_smooth3X(Yx,2)))))); 
  end

  if any( minres  < vx_vol(vi,:) * .8 )
    % interpolate data 
    [Yx,Vx] = cat_vol_resize( Yx , 'interphdr', spm_vol(job.data{vi}), min(minres,vx_vol(vi,:)), 5);
    vx_vol(vi,:) = Vx.resN; 
    Vx = Vx.hdrN; 

    if 0
      % post correction is not working and cause issues in the longitudinal bias modeling
      Yx = min(max(Yx(:)),max(min(Yx(:)),Yx + ... 
        (Yx~=0) .* max(-con/6.*Yx,min(con/6.*Yx, ...
        sharp/2 .* (Yx - smooth3(Yx))))));
    end
  end

  if sharpen || any( minres  < vx_vol(vi,:) * .8 )
    % Write image to reload by the nifti routine
    % We use the tempdir to keep the file name but this means that we have
    % to move and cleanup data later 
    Vo = Vx; Vo.fname = spm_file(Vo.fname,'path',tempdir); 
    spm_write_vol(Vo,Yx);
    Nii(vi) = nifti( Vo.fname );
    tempimgs{vi} = Vo.fname;
  end
end


% sometimes for quite anisotropic data long. registration fails and will be
% called again with more isotropic spatial resolution using isores = 3
% (optimal)
try
  out = cat_vol_groupwise_ls(Nii, output, prec, w_settings, b_settings, s_settings, ord, use_brainmask, reduce, setCOM, isores);
catch
  fprintf('Recall cat_vol_groupwise_ls again with more isotropic spatial resolution.\n')
  out = cat_vol_groupwise_ls(Nii, output, prec, w_settings, b_settings, s_settings, ord, use_brainmask, reduce, setCOM, 3);
end

% RD202510: move and cleanup interpolated data
if any( cellfun(@isempty,tempimgs) == 0 )
  % if the first file had to be interpolated we have to move it from the tmp dir
  if ~strcmp( spm_file( out.avg{1} ,'path') , fileparts( job.data{1} ) )
    movefile( out.avg{1} , spm_file( out.avg{1} ,'path', fileparts( job.data{1} ) ) ); 
  end
  if isfield( out , 'rimg' )
    for vi=1:numel(tempimgs)
      % move realigned images if they where template
      movefile( out.rimg{vi} , spm_file( out.rimg{vi} ,'path', fileparts( job.data{vi} ) ) ); 
      delete(tempimgs{vi}); 
    end
  end
end


return

