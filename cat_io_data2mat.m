function out = cat_io_data2mat(opt,par,scaling)
% Save spatially registered volume or resampled surface data as Matlab data matrix for further 
% use with machine learning tools such as relevance/support vector approaches or Gaussian Process
% models. Spatial structure of the data is not considered. 
% Volume data will be optionally masked to remove non-brain areas.
% 
% FORMAT cat_io_data2mat(opt,scaling)
%
% data format:
% opt.use_double - optional parameter to save Y as double instead of single
%
% volume or surface data:
% opt.data       - cell of char array of filenames
% opt.c          - confounds data to be removed
% opt.fname      - filename for saving mat-file
% opt.outdir     - output directory of saved mat-file
%
% additional parameters for volume data only:
% opt.resolution - resampling spatial resolution for volume data
% opt.mask       - optional brainmask for volume data
% opt.mask_th    - optional threshold for brainmask for volume data
% opt.fwhm       - optional Gaussian smoothing kernel size in FWHM
%
% additional parameters to be saved with mat file:
% par            - structure with parameter as name and the values
%
% scaling        - optionally define either a vector for user-specified scaling data or a 
%                  constant (2 for global scaling)
%
% saved parameters:
% Y          - data matrix with size number of subjects x number of voxels/vertices
% label      - label of samples 
% ind        - index for volume or surface data inside mask
% dim        - dimension of original data
% V          - structure array containing data information of surface or resampled volume
% sample     - structure array containing sample information
%
% Examples:
% Select recursively all gray matter segments from folder1
% for the 1st sample and folder 2 from the 2nd sample and save resampled data 
% with 4mm resampling spatial resolution after filtering with fwhm of 8mm.
% Additionally save age and male as parameter in the mat-file.
% the parameter "files" should be defined as "{files}")
% files{1} = spm_select('FPListRec',folder1,'^mwp1.*\.nii$');
% files{2} = spm_select('FPListRec',folder2,'^mwp1.*\.nii$');
% cat_io_data2mat(struct('data',{files},'resolution',4,'fwhm',8,'mask',cat_get_defaults('extopts.brainmask'),...
%    'fname','Data.mat'),struct('age',age,'male',male));
%
% Select recursively all 12mm smoothed and resampled thickness data from current folder
% and save Data.mat in subfolder test
% files = spm_select('FPListRec',pwd,'^s12.mesh.thickness.resampled_32k.*\.gii$');
% cat_io_data2mat(struct('data',files,'fname','Data.mat','outdir','test'));
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*NASGU,*AGROW>

if ~exist('opt','var') || (~isfield(opt,'data_type') && ~isfield(opt,'data') )
  error('No data defined.');
end

if ~exist('opt','var'), opt = struct(); end
def.c      = [];
def.fname  = 'Data.mat';
def.outdir = {'.'};
def.use_double = 0;
opt = cat_io_checkinopt(opt,def);

if isfield(opt,'data') % use simpler field structure for call as script
  if ischar(opt.data)
    sample.data{1} = opt.data;
  else
    sample.data = opt.data;
  end
  if isfield(opt,'mask')
    brainmask = opt.mask;
  end
  if isfield(opt,'mask_th')
    mask_th = opt.mask_th;
  else
    mask_th = 0.5;
  end
  if isfield(opt,'resolution')
    resolution = opt.resolution;
  end
  if isfield(opt,'fwhm')
    fwhm = opt.fwhm;
  end
elseif isfield(opt.data_type,'vol_data') % volume data
  brainmask = char(opt.data_type.vol_data.mask);
  sample = opt.data_type.vol_data;
  resolution = sample.resolution;
elseif isfield(opt.data_type,'surf_data') % surface data
  sample = opt.data_type.surf_data;
else
  error('No data defined.');
end
  
n_samples = numel(sample.data);
n_subjects = zeros(n_samples,1);
label = []; 

[pth, name, ext] = spm_fileparts(char(sample.data{1}(1,:)));
if strcmp(ext,'.gz')
  is_gz = 1;
else
  is_gz = 0;
end

V = cell(n_samples,1);
label = [];
for i = 1:n_samples
  if is_gz
    V{i} = spm_vol(char(sample.data{i}));
  else
    V{i} = spm_data_hdr_read(char(sample.data{i}));
  end
  n_subjects(i) = size(V{i},1);
  label = [label; i*ones(n_subjects(i),1)]; 
end
n_all_subjects = numel(label);

% cell structure is expected
if ~isempty(opt.c) && isnumeric(opt.c), opt.c = {opt.c}; end

n_confounds = numel(opt.c);
confounds = [];
for i = 1:n_confounds
  [m,n] = size(opt.c{i});
  
  % transpose if necessary
  if m ~= n_all_subjects && n == n_all_subjects
    opt.c{i} = opt.c{i}';
    m = size(opt.c{i},1);
  end 
  
  % check size of confounds      
  if m ~= n_all_subjects
    error('Length of nuisance parameters (m=%d) differs from number of subjects (n=%d)',m,n_all_subjects);
  end

  confounds = [confounds opt.c{i}];
end

outname = opt.fname;
out.fname{1} = outname; 
outdir = opt.outdir{1};
if ~isempty(outdir)
  if ~exist(outdir,'dir')
    mkdir(outdir);
  end
  outname = fullfile(outdir,outname);
end

% 3D data
if ~spm_mesh_detect(V{1}(1))
  % this is probably not the most elegant way but used here for compatibility with my BrainAGE tools...
  % 1mm reference fields
  Vres.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
  Vres.dim = [181 217 181];
  
  if ~isempty(brainmask)
    Vm  = spm_vol(brainmask);
  end
  
  Vres.dim = round(Vres.dim/resolution);
  Vres.mat(1:3,1:3) = resolution*Vres.mat(1:3,1:3);
  Vres.mat(1:3,4) = Vres.mat(1:3,4) - [resolution resolution resolution]';
  dim = Vres.dim(1:3);

else
  % check that mesh contains data
  if isfield(V{1}(1),'private') && ~isfield(V{1}(1).private,'cdata')
    if ~isfield(V{1}(1),'cdata')
      error('No data found in mesh')
    end
  end
  
  % use mask from first mesh and assume this holds for all data
  if isfield(V{1}(1),'private')
    ind = find(isfinite(V{1}(1).private.cdata(:)));
  else
    ind = find(isfinite(V{1}(1).cdata(:)));
  end
  dim = V{1}(1).dim(1);
end

% global scaling
if nargin > 2
  % user specified scaling
  if isvector(scaling)
    % check size of scaling vector
    if numel(scaling) ~= n_all_subjects
      error('Size of scaling parameter (n=%d) does not fit to number of subjects (n=%d)',numel(scaling), n_all_subjects)
    end
    count = 1;
    for j=1:n_samples     
      for i = 1:n_subjects(j)
        % compute mean voxel value (within per image fullmean/8 mask)
        if ~spm_mesh_detect(V{1}(1))
          g = spm_global(V{j}(i));
        else
          if isfield(V{1}(1),'private')
            y = V{j}(i).private.cdata(:);
          else
            y = V{j}(i).cdata(:);
          end
          g = mean(y(isfinite(y)));
        end
        V{j}(i).pinfo(1:2,:) = V{j}(i).pinfo(1:2,:)/g(count); % <<< the g is not existing otherwise!
        count = count + 1;
      end
    end
  % mean scaling
  elseif scaling == 2
    fprintf('Calculating globals\n')
    for j=1:n_samples     
      for i = 1:n_subjects(j)
        % compute mean voxel value (within per image fullmean/8 mask)
        if ~spm_mesh_detect(V{1}(1))
          g = spm_global(V{j}(i));
        else
          if isfield(V{1}(1),'private')
            y = V{j}(i).private.cdata(:);
          else
            y = V{j}(i).cdata(:);
          end
          g = mean(y(isfinite(y)));
        end
        V{j}(i).pinfo(1:2,:) = V{j}(i).pinfo(1:2,:)/g;
      end
    end
  end
end

Y = [];
Ymean = [];
C = zeros(n_all_subjects);

% 3D data
if ~spm_mesh_detect(V{1}(1))

  if ~exist('mask_th') || isempty(mask_th)
    mask_th = 0.5;
  end
  M1     = cell(Vres.dim(3),1);
  
  if isempty(brainmask)
    mask_ind = ones(Vres.dim,'logical');
  else
    mask_ind = zeros(Vres.dim,'logical');
  end
  ind = find(mask_ind);
  
  for sl=1:Vres.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    
    if ~isempty(brainmask)
      Mm  = Vres.mat\Vm.mat\M;
      mask_slice = spm_slice_vol(Vm,Mm,Vres.dim(1:2),1);
      mask_ind(:,:,sl) = mask_slice > mask_th;
    end
    M1{sl} = Vres.mat\V{1}(1).mat\M;
  end
  
    
  for j=1:n_samples
    if n_subjects(j) > 500, cat_progress_bar('Init',n_subjects(j),'reading...','subjects completed','cmd%'); end
    yi = zeros(n_subjects(j), sum(mask_ind(:)), 'single');
    for i = 1:n_subjects(j)

      vol = spm_read_vols(V{j}(i));

      % optional smoothing if fwhm is defined
      if exist('fwhm','var')
        if isscalar(fwhm) 
          fwhm = repmat(fwhm,1,3);
        end
        if sum(fwhm) > 0
          spm_smooth(vol,vol,fwhm,0);
        end
      end
   
      ysl = []; 
      for sl=1:Vres.dim(3)
    
        % read data inside mask
        if ~isempty(mask_ind(:,:,sl))
          try
            d = spm_slice_vol(vol,M1{sl},Vres.dim(1:2),1);
          catch
            fprintf('File %s could not be read\n',V{j}(i).fname);
            return
          end
        end
        ysl = [ysl; d(mask_ind(:,:,sl))];
      end
      yi(i,:) = single(ysl);
      
      if n_subjects(j) > 500, cat_progress_bar('Set',i); end
    end
    if n_subjects(j) > 500, cat_progress_bar('Clear'); end
    
    Y = [Y; yi];
  end

else % meshes
  Y = zeros(n_all_subjects,numel(ind));
  count = 1;
  for j=1:n_samples     
    for i = 1:n_subjects(j)
      if dim ~= V{j}(i).dim(1)
        error('Mesh size of %s differs (%d vs. %d)',V{j}(i).fname,V{j}(i).dim(1),dim)
      end
      try
        if isfield(V{1}(1),'private')
          Y(count,:) = V{j}(i).private.cdata(ind);
        else
          Y(count,:) = V{j}(i).cdata(ind);
        end
      catch
        fprintf('File %s could not be read\n',V{j}(i).fname);
        return
      end
      count = count + 1;
    end
  end
end

% remove confounds
if ~isempty(confounds)
  beta = pinv(confounds)*Y;
  Y = Y - confounds*beta;
end

% convert to single if not otherwise defined
if ~opt.use_double
  Y = single(Y);
end

% save surface or resampled volume structure
if spm_mesh_detect(V{1}(1))
  if isfield(V{1}(1),'private')
    V = V{1}(1).private;
  else
    V = V{1}(1);
  end
  try V = rmfield(V,'cdata'); end %#ok<TRYNC>
else
  V = Vres;
end

save(outname,'Y','label','dim','V','ind','sample','-v7.3');
fprintf('Save data (Y,label,V,dim,ind) in %s.\n',outname);

% add additional parameters if defined
if nargin > 1
  fnames = fieldnames(par);
  str = [];
  fprintf('Append:');
  for i = 1:numel(fnames)
    name = fnames{i};
    eval([name '=par.(name);']);
    save(outname,'-append',name);
    fprintf(' %s',name);
  end
  fprintf('\n');
end
