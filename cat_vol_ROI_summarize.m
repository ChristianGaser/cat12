function Y = cat_vol_ROI_summarize(job)
% Summarise data within a Region of Interest
% FORMAT Y = cat_vol_ROI_summarize(job)
%
% job fields:
% atlases  - cell structure of atlas files
% field1   - cell structure of one deformation file
% field    - cell structure of multiple deformation files
% images   - cell structure of value files for ROI estimation
% fhandle  - function handle to be applied on image data within ROI
%            Default is @mean. Use keyword 'volume' to estimate regional 
%            volume in ml.
%
%            Examples
%            calculate absolute amplitude between 10-90% percentile:
%            @(x) abs(diff(spm_percentile(x,[10 90]))) 
%
%            get mean inbetween 10-90% percentile:
%            @(x) mean(x>spm_percentile(x,10) & x<spm_percentile(x,90)
%
% Y        - output cell structure of ROI values with size n_atlases,n_values
%__________________________________________________________________________
%
% Example:
% opt = struct('atlases',{fullfile(spm('dir'),'toolbox','cat12','templates_volumes','cobra.nii')},...
%              'field1',{spm_select('FPList','mri','^y_.*nii')},...
%              'images',{spm_select('FPList','mri','^p1_.*nii')},...
%              'fhandle','@median');
% Y = cat_vol_ROI_summarize(opt)
%__________________________________________________________________________
%  $Id$

if nargin == 0
  atlas_file = cellstr(spm_select([1 Inf],'image','Select atlas file(s)',{},...
    fullfile(spm('dir'),'toolbox','cat12','templates_volumes')));
  def_file = cellstr(spm_select(1,'image','Select forward deformation field',{},'','^y_.*\.nii$'));
  value_file = cellstr(spm_select([1 Inf],'image','Select co-registered files for ROI estimation'));
  
  fhandle_sel = spm_input('Function to summarize?', 1, 'm','volume|avg|std|customized');
  switch fhandle_sel
    case 1, fhandle ='volume';
    case 2, fhandle = @(x) mean(x);
    case 3, fhandle = @(x) std(x);
    case 4, fhandle = str2func(spm_input('Function handle to summarize', 1, 's','@median'));
  end  
end

many_images = 1;

% consider different methods from conf-file
if isfield(job,'Method')
  if isfield(job.Method,'ManySubj')
    job = job.Method.ManySubj;
    def_file = job.field;
    many_images = 0;
  elseif isfield(job.Method,'ManyImages')
    job = job.Method.ManyImages;
    def_file = job.field1;
    many_images = 1;
  end

  value_file = job.images;

  if isfield(job.fhandle,'fun')
    fhandle = job.fhandle.fun;
  else
    fhandle = job.fhandle.cfun;
  end

  AN  = fieldnames(job.atlases);
  fai = 1;
  
  for ai=1:numel(AN)
    if ~iscell(job.atlases.(AN{ai}))
      if job.atlases.(AN{ai})
        atlas_file{fai} = fullfile(spm('dir'),'toolbox','cat12','templates_volumes',[AN{ai} '.nii']);
        fai = fai+1;
      end
    elseif ~isempty(char(job.atlases.(AN{ai})))
      atlas_file{fai} = char(job.atlases.(AN{ai}));
      fai = fai+1;
    end
  end
  
elseif nargin == 1
  atlas_file = job.atlases;
  
  if isfield(job,'field1')
    def_file = job.field1;
    many_images = 1;
  else
    def_file = job.field;
    many_images = 0;
  end
  
  value_file = job.images;
  fhandle = job.fhandle;
end

% recognize 'volume' keyword to estimate volume in ml
get_volume = 0;
if ischar(fhandle)
  if strcmpi(fhandle,'volume')
    fhandle = @(x) sum(x);
    get_volume = 1;
  else
    fhandle = str2func(fhandle);
  end
end

Vvalue = spm_vol(value_file);

% ignore NaNs
dropNaNs = @(x) x(~isnan(x));

n_atlas = numel(atlas_file);
n_value = numel(Vvalue);
n_def   = numel(def_file);

% check whether deformation file contains inverse deformations
if ~isempty(strfind(def_file{1},'iy_'))
  error('You have to select only forward deformations!');
end

Y = cell(n_atlas,n_value,n_def,1);
csv_arr = cell(n_atlas,n_def);
 
for ai = 1:n_atlas

  atlas = round(spm_read_vols(spm_vol(atlas_file{ai})));
  [pp,atlas_name] = spm_fileparts(atlas_file{ai});
  csv = get_atlas_csv(atlas_file{ai},atlas);
  
  % remove unnecessary columns and transpose to have ROI estimates 
  % for each value file in lines
  if csv{2,1} == 0, csv(2,:) = []; end % remove BG if found

  [bb,vox] = spm_get_bbox(atlas_file{ai});
  vox = abs(vox);
    
  % go through all deformation fields
  for di = 1:n_def
    def = struct('field1',{{def_file{di}}},'images',{value_file},'interp',1,'modulate',0,'vox',vox,'bb',bb,'verb',0);
  
    % use modulated data to get correct volumes
    if get_volume
      def.modulate = 1;
    end
    wvol = cat_vol_defs(def);

    csv_arr{ai,di} = transpose(csv(:,1:3));
    n_structures = size(csv_arr{ai,di},2)-1; % 1st column contains name
    
    Ytmp = zeros(n_structures,1);
    
    % go through all value files
    for vi = 1:n_value
      value = wvol{vi}{di};

      if many_images % many images
        [pth,nam] = spm_fileparts(value_file{vi});
      else % many subjects
        [pth,nam] = spm_fileparts(value_file{vi}{di});
      end
  
      % also consider 4D data
      for vii = 1:size(value,4)
  
        % add numbers for 4D data
        if size(value,4) > 1
          nam1 = sprintf('%s_%04d',nam,vii);
        else
          nam1 = nam;
        end
        
        sel_value = value(:,:,:,vii);
        csv_arr{ai,di}{end+1,1} = nam1;
        
        % go through all ROI structures
        for ri = 1:n_structures
          % 1st column of csv_array contains name, thus add 1 to index
          val = fhandle(dropNaNs(sel_value(atlas==csv_arr{ai,di}{1,ri+1})));
          if isempty(val)
            val = NaN;
          elseif get_volume
            val = val * prod(vox)/1000;
          end
          
          % save ROI values
          csv_arr{ai,di}{end,ri+1} = val;
          Ytmp(ri) = val;
        end
        
        Y{ai,vi,di,vii} = Ytmp(1:end);
      end
    end
  end
end

% write csv-file in folder of deformation file and use its name
for ai = 1:n_atlas
  [pp,atlas_name] = spm_fileparts(atlas_file{ai});
  
  for di=1:n_def
    [pp,def_name] = spm_fileparts(def_file{di});
    % get output name from deformation field and add atlas name
    def_name = strrep(def_name,'y_','');
    out_file = fullfile(pp,[atlas_name '_' def_name '.csv']);
    fprintf('ROI estimation %s saved.\n',out_file);
    cat_io_csv(out_file,csv_arr{ai,di},'','',struct('format','%g'));
  end
  
end

%=======================================================================
function csv = get_atlas_csv(atlas_file,atlas)
% load csv atlas file

[pp,ff] = fileparts(atlas_file);
csvf = fullfile(pp,[ff '.csv']);
if exist(csvf,'file')
  csv = cat_io_csv(csvf,'','',struct('delimiter',';')); 
else
  IDs = unique(atlas(isfinite(atlas)));
  csv = [{'ROIid','ROIabbr','ROIname'}; ...
    num2cell(IDs) ...
    cellstr([repmat('ROI',numel(IDs),1) num2str(IDs,'%03d')]) ...
    cellstr([repmat('ROI',numel(IDs),1) num2str(IDs,'%03d')])];
end

