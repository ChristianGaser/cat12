function out = cat_vol_avg(job)
% Calculate average of images of same dimension and voxel size
% This function also works for 5D images (e.g. deformations).
%
% FORMAT out = cat_vol_avg(job)
% job.data      - data array of input files
% job.weighting - weighting vector for input files
%                 If you want to use weighting maps checkto mimcalc.
% job.write_var - Write also the variance map with suffix '_var'.
%                 We use the suffix as it is connected to the avg map
%                 and both represent single maps in a dataset.
% job.output    - output name
% job.outdir    - output directory
%
% out           - output name
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
%

% interactive call function
if ~nargin
  job.outdir{1} = '';
  job.data  = cellstr(spm_select(Inf,'image','Select (spatially registered) images for average'));
  [~, name] = spm_str_manip(spm_str_manip(job.data,'t'),'C');
  pos = strfind(name.e,',1');
  if ~isempty(pos)
    name.e = name.e(1:pos-1);
  end
  Q = ['avg_' name.s name.e];
  job.output = spm_input('Output filename',1,'s',Q);
end

[p,nam,ext] = spm_fileparts(job.output);
if isempty(p)
    if isempty(job.outdir{1})
        p = spm_fileparts(job.data{1});
    else
        p = job.outdir{1};
    end
end
if isempty(nam)
    nam = ['avg_' spm_file(job.data{1},'basename')];
    ext = ['.' spm_file(job.data{1},'ext')];
end
if isempty(ext)
    ext = spm_file_ext;
end
out.files = { fullfile(p,[nam ext]) };

N = nifti(char(job.data));

if length(N)>1 && any(any(diff(cat(1,N.dat.dim),1,1),1))
  error('images don''t all have same dimensions')
end
if max(max(max(abs(diff(cat(3,N.mat),1,3))))) > 1e-8
  error('images don''t all have same orientation & voxel size')
end

d = N(1).dat.dim;
% extend dimensions if necessary
d = [d ones(1, 5-length(d))];


if isfield(job,'weighting') && ~isempty(job.weighting)
  if numel(job.weighting)~=length(N)
    error('number of weights has to be identical to the number of images')
  end
  weighting = job.weighting ./ sum(job.weighting); 
else
  weighting = repmat( 1/length(N) , 1, length(N));
end 


% write average map
avg = zeros(d);

for i = 1:length(N)
  for j = 1:d(4)
    for k = 1:d(5)
      avg(:,:,:,j,k) = avg(:,:,:,j,k) + N(i).dat(:,:,:,j,k) * weighting(i) ;
    end
  end
end

Nout = N(1);
Nout.dat.fname = out.files{1};
create(Nout);
Nout.dat(:,:,:,:,:) = avg;

% cmd line output
cmd = 'spm_image(''display'',''%s'')';
fprintf('Average Image:  %s\n', spm_file(out.files{1}, 'link', cmd) );


% write variance map
if isfield(job,'write_var') && job.write_var
  % estimate variance
  var = zeros(d);
  for i = 1:length(N)
    for j = 1:d(4)
      for k = 1:d(5)
        var(:,:,:,j,k) = var(:,:,:,j,k)  +  (avg(:,:,:,j,k) - N(i).dat(:,:,:,j,k)).^2 * weighting(i) ;
      end
    end
  end
  
  for i=1:min(2,job.write_var)
    % write output
    Nout = N(1);
    if i==1
      out.files{i+1} = spm_file( out.files{1}, 'suffix', '_var');
      var2 = var;
    else
      out.files{i+1} = spm_file( out.files{1}, 'suffix', '_vardivavg');
      var2 = var ./ max(eps,abs(avg));
    end
    Nout.dat.fname = out.files{i+1};
    create(Nout);
    Nout.dat(:,:,:,:,:) = var2;
    
    % cmd line output
    cmd = 'spm_image(''display'',''%s'')';
    fprintf('Variance Image: %s\n', spm_file( out.files{i+1} , 'link', cmd) );
  end
end

if ~nargout
  clear out
end