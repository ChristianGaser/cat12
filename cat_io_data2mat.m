function varargout = cat_io_data2mat(opt)
% Save spatially registered volume or resampled surface data as Matlab data matrix for further 
% use with machine learning tools.
% Volume data will be optionally masked to remove non-brain areas
% 
% FORMAT cat_io_data2mat(opt)
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
%
% saved parameters:
% Y          - data matrix with size number of subjects x number of voxels/vertices
% label      - label of samples 
% ind        - index for volume or surface data inside mask
% dim        - dimension of original data
%
% Examples:
% Select recursively all 8mm smoothed gray matter segments from folder1
% for the 1st sample and folder 2 from the 2nd sample and save resampled data 
% with 4mm resampling spatial resolution (please note that for several samples
% the parameter "files" should be defined as "{files}")
% files{1} = spm_select('FPListRec',folder1,'^s8mwp1.*\.nii$');
% files{2} = spm_select('FPListRec',folder2,'^s8mwp1.*\.nii$');
% cat_io_data2mat(struct('data',{files},'resolution',4,'mask',cat_get_defaults('extopts.brainmask'),'fname','Data.mat'));
%
% Select recursively all 12mm smoothed and resampled thickness data from current folder
% and save Data.mat in subfolder test
% files = spm_select('FPListRec',pwd,'^s12.mesh.thickness.resampled_32k.*\.gii$');
% cat_io_data2mat(struct('data',files,'fname','Data.mat','outdir','test'));
%_______________________________________________________________________
% Christian Gaser
% $Id$

if ~exist('opt','var') || (~isfield(opt,'data_type') && ~isfield(opt,'data') )
  error('No data defined.');
end

if ~exist('opt','var'), opt = struct(); end
def.c      = [];
def.fname  = 'Data.mat';
def.outdir = {'.'};
opt        = cat_io_checkinopt(opt,def);

if isfield(opt,'data') % use simpler field structure for call as script
  if ischar(opt.data)
    sample.data{1} = opt.data;
  else
    sample.data = opt.data;
  end
  if isfield(opt,'mask')
    brainmask = opt.mask;
  end
  if isfield(opt,'resolution')
    resolution = opt.resolution;
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

V = cell(n_samples,1);
label = [];
for i = 1:n_samples;
  V{i} = spm_data_hdr_read(char(sample.data{i}));
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
    [m,n] = size(opt.c{i});
  end
  
  % check size of confounds      
  if m ~= n_all_subjects,
    error('Length of nuisance parameters (m=%d) differs from number of subjects (n=%d)',m,n_all_subjects);
  end

  confounds = [confounds opt.c{i}];
end

fname = opt.fname;
outdir = opt.outdir{1};
if ~isempty(outdir)
  fname = fullfile(outdir,fname);
end

% 3D data
if ~spm_mesh_detect(V{1}(1))
  % 1mm reference fields
  V0.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
  V0.dim = [181 217 181];
  
  if ~isempty(brainmask)
    Vm  = spm_vol(brainmask);
  end
  
  V0.dim = round(V0.dim/resolution);
  V0.mat(1:3,1:3) = resolution*V0.mat(1:3,1:3);
  V0.mat(1:3,4) = V0.mat(1:3,4) - [resolution resolution resolution]';
  dim = V0.dim(1:3);

else
  % check that mesh contains data
  if isfield(V{1}(1),'private') & ~isfield(V{1}(1).private,'cdata')
    if ~isfield(V{1}(1),'cdata')
      error('No data found in mesh')
    end
  end
  
  % use mask from first mesh and assume this holds for all data
  if isfield(V{1}(1),'private')
    ind = find(isfinite(V{1}(1).private.cdata));
  else
    ind = find(isfinite(V{1}(1).cdata));
  end
  dim = V{1}(1).dim(1);
end

% global scaling
if nargin > 3
  % user specified scaling
  if isvector(scaling)
    % check size of scaling vector
    if numel(scaling) ~= n_all_subjects
      error('Size of scaling parameter (n=%d) does not fit to number of subjects (n=%d)',numel(scaling), n_all_subjects)
    end
    count = 1;
    for j=1:n_samples     
      for i = 1:n_subjects(j)
        V{j}(i).pinfo(1:2,:) = V{j}(i).pinfo(1:2,:)/g(count);
        count = count + 1;
      end
    end
  % mean scaling
  elseif scaling == 2
    fprintf('Calculating globals\n')
    for j=1:n_samples     
      for i = 1:n_subjects(j)
        % compute as mean voxel value (within per image fullmean/8 mask)
        if ~spm_mesh_detect(V{1}(1))
          g = spm_global(V{j}(i));
        else
          if isfield(V{1}(1),'private')
            y = V{j}(i).private.cdata;
          else
            y = V{j}(i).cdata;
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
  spm_progress_bar('Init',V0.dim(3),'reading...','planes completed');
  ind = [];
  for sl=1:V0.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    if ~isempty(brainmask)
    
      Mm  = V0.mat\Vm.mat\M;
      mask_slice = spm_slice_vol(Vm,Mm,V0.dim(1:2),1);
      ind0 = find(mask_slice > 0.5);
    else
      ind0 = find(ones(V0.dim(1:2)));
    end
    M1  = V0.mat\V{1}(1).mat\M;
    
    ind = [ind; ind0];
    clear mask_slice
  
    % read data inside mask
    if ~isempty(ind0)
      yslice = [];
      for j=1:n_samples
        y = zeros(n_subjects(j), length(ind0));
          for i = 1:n_subjects(j)
              try
                d = spm_slice_vol(V{j}(i),M1,V0.dim(1:2),1);
              catch
                fprintf('File %s could not be read\n',V{j}(i).fname);
                return
              end
              y(i,:) = d(ind0);
          end
        yslice = [yslice; y];
      end
      Y = [Y yslice];
    end
    spm_progress_bar('Set',sl)
  end
  spm_progress_bar('Clear')
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

save(fname,'Y','label','dim','ind');
fprintf('Save data (Y,label,dim,ind) in %s.\n',fname);

if nargout == 1
  varargout{1}.fname = fname;
end
