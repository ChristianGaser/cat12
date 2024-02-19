function OV = cat_vol_slice_overlay(OV)
% Extension/wrapper to slice_overlay
% Call help for slice_overlay for any additional help
% 
% Additional fields to slice_overlay:
% OV.name       - char array of filenames for overlay that can be interactively
%                 selected
% OV.slices_str - char array of slice values (e.g. '-32:2:20')
%                 use empty string for automatically estimating slices with
%                 local maxima
% OV.xy         - define number of columns and rows
%                 comment this out for interactive selection or set the values
%                 to [Inf 1] for using one row and automatically estimate number
%                 of columns or use [1 Inf] for using one column
% OV.atlas      - define atlas for labeling (e.g. 'cat12_cobra')
%                 comment this out for interactive selection
%                 or use 'none' for no atlas information
% OV.save       - save result as png/jpg/pdf/tif
%                 comment this out for interactive selection or use '' for not 
%                 saving any file or use just file extension (png/jpg/pdf/tif) to 
%                 automatically estimate filename to save
% OV.FS         - normalized font size (default 0.08)
% OV.name_subfolder
%               - if result is saved as image use up to 2 subfolders to add their 
%                 names to the filename (default 1)
% OV.overview   - use empty brackets to not suppress slice overview (.e.g []);
% OV.pos        - define first two numbers of image position
% OV.bkg_col    - color of background ([0 0 0] as default)
% OV.fig        - figure number (default 22)
% OV.cbar       - show colorbar (leave empty for no colorbar)
% OV.labels     - show slice label text (leave empty for no label text)
%
% see cat_vol_slice_overlay_ui.m for an example
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

clear global SO
global SO

if ~nargin
    
  imgs = spm_select(2, 'image', 'Select additional overlay image', cat_get_defaults('extopts.shootingT1'));
  if isempty(imgs)
    return;
  end
  OV = pr_basic_ui(imgs, 0);
  
  % set options
  OV.opacity = Inf;
  OV.reference_image = deblank(imgs(1, :));
  OV.reference_range = OV.img(1).range;
  OV.name = imgs(2:end, :);
  OV.cmap = OV.img(2).cmap;
  OV.range = OV.img(2).range;
  OV.slices_str = '';
end

% get fontsize
if isfield(OV,'FS')
  FS = OV.FS;
else
  FS = 0.08;
end

if ~isfield(OV,'fig')
  OV.fig = 22;
end

% hot colormap by default
if ~isfield(OV,'cmap')
  OV.cmap = hot(256);
end

% clip colorbar
if isfield(OV,'clip') && all(isfinite(OV.clip)) && OV.clip(2) ~= OV.clip(1)
  ncol = length(OV.cmap);
  col_step = (OV.range(2) - OV.range(1)) / ncol;
  cmin = max([1, ceil((OV.clip(1) - OV.range(1)) / col_step)]);
  cmax = min([ncol, floor((OV.clip(2) - OV.range(1)) / col_step)]);
  OV.cmap(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
  OV.func = sprintf('i1(i1>%f & i1<%f)=NaN;',OV.clip(1),OV.clip(2));
else
  OV.clip = [Inf -Inf];
end

% black background by default
if ~isfield(OV,'bkg_col')
  SO.bkg_col = [0 0 0];
else
  SO.bkg_col = OV.bkg_col;
end

% check filename whether log. scaling was used
OV.logP = zeros(size(OV.name, 1));
for i = 1:size(OV.name, 1)
  if ~isempty(strfind(OV.name(i, :), 'logP')) || ~isempty(strfind(OV.name(i, :), 'log_'))
    OV.logP(i) = 1;
  end
end

% check fields of OV structure
fieldnames = char('reference_image', 'reference_range', ...
  'opacity', 'cmap', 'name', 'range', 'logP', 'slices_str', 'transform');
for i = 1:size(fieldnames, 1)
  str = deblank(fieldnames(i, :));
  if ~isfield(OV, str)
    error([str ' not defined']);
  end
end

cmap_bivariate = [1 - (hot); hot]; % colormap if range(1) < 0

if isfield(OV, 'labels')
  SO.labels = OV.labels;
  if ~isempty(SO.labels)
    SO.labels.size = FS*0.65;
  end
end

if isfield(OV, 'overview')
  SO.overview = OV.overview;
end

if isfield(OV, 'cbar')
  SO.cbar = OV.cbar;
else
  SO.cbar = 2; % colorbar
end

n = size(OV.name, 1);

str = deblank(OV.name(1, :));
for i = 2:n, str = [str '|' deblank(OV.name(i, :))]; end

if n > 1
  sel = spm_input('Select image', 1, 'm', str);
else
  sel = 1;
end

OV.name = deblank(OV.name(sel, :));

% if only one argument is given assume that parameters are the same for all files
if size(OV.range, 1) > 1
  range = OV.range(sel, :);
else
  range = OV.range;
end

if size(OV.logP, 1) > 1
  logP = OV.logP(sel);
else
  logP = OV.logP;
end

% for log-scaled p-values we should rather use gt than ge for comparison with threshold
if logP
  compare_to_threshold = @(a,b) gt(a,b);
else
  compare_to_threshold = @(a,b) ge(a,b);
end

img = OV.name;

n_slice = size(OV.slices_str, 1);
if n_slice > 0
  for i = 1:n_slice
    try
      slices{i} = eval(OV.slices_str(i, :));
    catch
      slices{i} = [];
    end
  end
else
  if isfield(OV,'slices')
    slices{1} = OV.slices;
  else
    SO.img(2).vol = spm_vol(OV.name);
    [mx, mn, XYZ, img2] = volmaxmin(SO.img(2).vol);

    % apply function to img2
    i1 = img2; eval(OV.func); img2 = i1;
    
    % threshold map and restrict coordinates
    Q = find(compare_to_threshold(img2,range(1)) & le(img2,range(2)));
    XYZ = XYZ(:, Q);
    img2 = img2(Q);
    
    M = SO.img(2).vol.mat;
    XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
    
    XYZ_unique = get_xyz_unique(XYZ, XYZmm, img2);
    orientn = strcmpi(OV.transform, {'sagittal', 'coronal', 'axial'});
    slices{1} = XYZ_unique{orientn};
  end
end

sl_name = [];
for i = 1:size(OV.transform, 1)
  if n_slice > 0
    sl_name = strvcat(sl_name, [OV.transform(i, :) ': ' OV.slices_str(i, :)]);
  else
    sl_name = strvcat(sl_name, OV.transform(i, :));
  end
end

str_select = deblank(sl_name(1, :));
for i = 2:n_slice, str_select = [str_select '|' deblank(sl_name(i, :))]; end
if n_slice
  ind = spm_input('Select slices', '+1', 'm', str_select);
else
  ind = 1;
end
OV.transform = deblank(OV.transform(ind, :));
slices = slices{ind};

SO.img(1).vol = spm_vol(OV.reference_image);
SO.img(1).prop = 1;
SO.img(1).cmap = gray;
SO.img(1).range = OV.reference_range;
SO.img(1).background = 0;

SO.img(2).vol = spm_vol(img);
SO.img(2).hold = 0; % use NN interpolation
SO.img(2).prop = OV.opacity; % transparent overlay
SO.img(2).cmap = OV.cmap; % colormap

if ~isfield(OV, 'func')
  SO.img(2).func = 'i1(i1==0)=NaN;';
else
  SO.img(2).func = OV.func;
end

if ~isfield(OV, 'range')
  [mx, mn] = volmaxmin(OV.img(2).vol);
  SO.img(2).range = spm_input('Intensity range for colormap', '+1', 'e', [mn mx], 2)';
else
  SO.img(2).range = range;
end

if range(1) >= 0
  SO.img(2).outofrange = {0, size(SO.img(2).cmap, 1)};
else
  SO.img(2).outofrange = {1, 1};
  % use bivariate colormap if OV was not defined
  if ~nargin, SO.img(2).cmap = cmap_bivariate; end
end

SO.transform = OV.transform;
SO.slices = slices;

if isempty(SO.slices)
  [mx, mn, XYZ, vol] = volmaxmin(SO.img(2).vol);
  
  % threshold map and restrict coordinates
  Q = find(compare_to_threshold(vol,SO.img(2).range(1)));
  XYZ = XYZ(:, Q);
  vol = vol(Q);
  
  M = SO.img(2).vol.mat;
  XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
  orientn = strcmpi(SO.transform, {'sagittal', 'coronal', 'axial'});
  
  XYZ_unique = get_xyz_unique(XYZ, XYZmm, vol);
  SO.slices = XYZ_unique{orientn};
  
  % update OV.slices_str for cat_surf_results to estimate available
  % rows/columns
  OV.slices_str = num2str(SO.slices);
end

n_images = length(SO.slices) + length(SO.cbar);
xy = get_xy(n_images);

n = size(xy, 1);
xy_name = num2str(xy);
str = deblank(xy_name(1, :));
for i = 2:n, str = [str '|' deblank(xy_name(i, :))]; end

% either interactively select columns/rows or use the defined values
if ~isfield(OV, 'xy')
  indxy = spm_input('Select number of columns/rows', '+1', 'm', str);
  xy = xy(indxy, :);
else
  if ~isfinite(OV.xy(1))
    ind = find(xy(:,2) == OV.xy(2));
    if isempty(ind)
      ind = n;
    end
    xy = xy(ind,:);
  elseif ~isfinite(OV.xy(2))
    ind = find(xy(:,1) == OV.xy(1));
    if isempty(ind)
      ind = 1;
    end
    xy = xy(ind,:);
  else
    xy = OV.xy;
  end
end

% get position of graphic figure
pos1 = spm('Winsize', 'Graphics');

screensize = get(0, 'screensize');

if ~isfield(SO,'overview')
  % prepare overview of slices
  V = SO.img(1).vol;
  ref_vol = spm_read_vols(V);
  ref_vol = 64 * (ref_vol - SO.img(1).range(1)) / (SO.img(1).range(2) - SO.img(1).range(1));
  vx = sqrt(sum(V.mat(1:3, 1:3).^2));
  Orig = round(V.mat \ [0 0 0 1]');
  
  h0 = figure(11); clf
  axes('Position', [0 0 1 1]);
  
  hold on
  dim = SO.img(1).vol.dim(1:3);
  switch lower(OV.transform)
    case 'sagittal'
      ref_img = ref_vol(:, :, Orig(3))';
      slices_vx = SO.slices / vx(1) + Orig(1);
      image(fliplr(ref_img))
      for i = slices_vx
        h = line([i i], [1 dim(2)]);
        set(h, 'Color', 'r')
      end
    case 'coronal'
      ref_img = squeeze(ref_vol(Orig(1), :, :))';
      slices_vx = SO.slices / vx(2) + Orig(2);
      image(ref_img)
      for i = slices_vx
        h = line([i i], [1 dim(3)]);
        set(h, 'Color', 'r')
      end
    case 'axial'
      ref_img = squeeze(ref_vol(Orig(1), :, :))';
      slices_vx = SO.slices / vx(3) + Orig(3);
      image(ref_img)
      for i = slices_vx
        h = line([1 dim(2)], [i i]);
        set(h, 'Color', 'r')
      end
  end

  set(h0, 'Position', [10 10 2 * size(ref_img, 2), 2 * size(ref_img, 1)], ...
      'MenuBar', 'none', ...
      'Resize', 'off', ...
      'PaperType', 'A4', ...
      'PaperUnits', 'normalized', ...
      'NumberTitle', 'off', ...
      'Name', 'Slices', ...
      'PaperPositionMode', 'auto');
  
  hold off
  axis off xy image
  colormap(gray)

end

SO.xslices = xy(:, 1);
switch lower(OV.transform)
  case 'sagittal'
    dim = xy .* SO.img(1).vol.dim(2:3);
  case 'coronal'
    dim = xy .* SO.img(1).vol.dim([1 3]);
  case 'axial'
    dim = xy .* SO.img(1).vol.dim(1:2);
end

% use double size
dim = 2 * dim;

scale = screensize(3:4) ./ dim;
% scale image only if its larger than screensize
if min(scale) < 1
  fig_size = min(scale) * dim * 0.975;
else
  fig_size = dim;
end

[pt, nm] = spm_fileparts(img);

if isfield(OV,'pos') && ishandle(OV.fig)
  pos0 = OV.pos(1:2);
else
  pos0 = [10 1200];
end

if ishandle(OV.fig)
  h = figure(OV.fig);
  pos = get(h, 'Position');
  set(h, 'Position', [pos(1:2) fig_size]);
else
  h = figure(OV.fig);
  set(h, 'Position', [pos0 fig_size]);
end

set(h, ...
    'MenuBar', 'none', ...
    'Resize', 'off', ...
    'PaperType', 'A4', ...
    'PaperUnits', 'normalized', ...
    'PaperPositionMode', 'auto', ...
    'NumberTitle', 'off', ...
    'Name', nm, ...
    'Visible', 'on');

SO.figure = h;
SO.area.units = 'pixels';

slice_overlay;

% remove remaining gray colored border
ax = get(SO.figure,'Children');
for i = 1:numel(ax)
  set(ax(i),'YColor',SO.bkg_col,'XColor',SO.bkg_col);
end

% change labels of colorbar for log-scale
H = gca;

if ~isempty(SO.cbar) && SO.cbar == 2 && logP

  YTick = get(H, 'YTick');

  % check whether lower threshold is P=0.05 and change values for YTick at
  % threshold
    if ~isempty(OV.clip) && abs(OV.clip(2)) >= 1.3 && abs(OV.clip(2)) <= 1.4 && OV.range(2) > OV.range(1)
      YTick_step = ceil((OV.range(2) - OV.range(1)) / numel(YTick));
      if OV.clip(1) <= - 1.3 && OV.clip(1) >= - 1.4 && OV.range(1) < 0
        values = [round(OV.range(1)):YTick_step:round(OV.range(2))];
        mid = find(YTick==0);
        if ~isempty(mid)
          values(mid-1:mid+1) = [log10(0.05) 0 -log10(0.05)];
        else
          values(values == -1) =  log10(0.05);
          values(values == 1)  = -log10(0.05);
        end
      else
        values = [0:YTick_step:round(OV.range(2))];
        values(2) = -log10(0.05);
      end
      
    else
      mn = floor(min(YTick));
      mx = ceil(max(YTick));

      % only allow integer values
      values = floor(mn:mx);
    end

     
  pos = get(get(gca, 'YLabel'), 'position');
  pos(1) = 2.5;
  
  set(H, 'YTick', values);
  YTick = get(H, 'YTick');
  
  YTickLabel = cell(length(YTick),1);
  for i = 1:length(YTick)
    if YTick(i) > 0 && YTick(i) >= OV.clip(2)
      if YTick(i) > 7
        % use 1E-x notation
        YTickLabel{i} = sprintf('%g', 10^(-YTick(i)));
      else
        % use 0.000x notation
        YTickLabel{i} = remove_zeros(sprintf('%3.7f', 10^(-YTick(i))));
      end
    elseif YTick(i) < 0 && YTick(i) <= OV.clip(1)
      if YTick(i) < -7
        % use 1E-x notation
        YTickLabel{i} = sprintf('-%g', 10^(YTick(i)));
      else
        % use 0.000x notation
        YTickLabel{i} = remove_zeros(sprintf('-%3.7f', 10^(YTick(i))));
      end
    else
      YTickLabel{i} = '';
    end
  end
  
  % update YTickLabel
  set(H, 'YTickLabel', YTickLabel,'YAxisLocation','left')
  set(get(gca, 'YLabel'), 'string', 'p-value', 'position', pos)
  
end

set(H, 'FontSize', FS, 'YColor', 1 - SO.bkg_col)
set(get(H, 'YLabel'), 'FontUnits', 'normalized', 'FontSize', 1.5*FS, 'Color', 1 - SO.bkg_col)

% we have to get rid off that annoying axis and simply draw a black box
% with 1 pixel width
posc = get(H,'Position');
posc(3) = 1;
a=axes(...
      'Parent',SO.figure,...
      'XTick',[],...
      'XTickLabel',[],...
      'YTick',[],...
      'YTickLabel',[],...
      'Units', 'pixels',...
      'YColor',[0 0 0],...
      'Color',[0 0 0],...
      'Box', 'off',...
      'Position',posc);

% select atlas for labeling
if isfield(OV, 'atlas')
  atlas_name = OV.atlas;
  if strcmpi(atlas_name,'none') || isempty(atlas_name)
    xA = [];
  else
    xA = spm_atlas('load',atlas_name);
  end
else
  list = spm_atlas('List','installed');
  atlas_labels{1} = 'None';
  j = 1;
  for i=1:numel(list)
    if ~strcmp(list(i).name,'Neuromorphometrics')
      atlas_labels{j+1} = list(i).name;
      j = j + 1;
    end
  end
    atlas = spm_input('Select atlas?', '1', 'm', atlas_labels);
    atlas_name = atlas_labels{atlas};
    if atlas > 1
      xA = spm_atlas('load',atlas_name);
    else
      xA = [];
    end
end

% atlas labeling
if ~isempty(xA)
  [mx, mn, XYZ, vol] = volmaxmin(SO.img(2).vol);
  
  % threshold map and restrict coordinates
  if SO.img(2).range(1) >= 0
    Q = find(compare_to_threshold(vol,SO.img(2).range(1)));
    XYZ = XYZ(:, Q);
    vol = vol(Q);
  end
  M = SO.img(2).vol.mat;
  XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
  
  % apply func that is defined for "i1"
  i1 = vol;
  eval(SO.img(2).func)
  
  % remove NaN values
  Q = find(isfinite(i1));
  XYZ = XYZ(:, Q);
  i1 = i1(Q);
  
  % find clusters
  A = spm_clusters(XYZ);
    
  labk   = cell(max(A)+2,1);
  Pl     = cell(max(A)+2,1);
  Zj     = cell(max(A)+2,1);
  maxZ   = zeros(max(A)+2,1);
  XYZmmj = cell(max(A)+2,1);
  Q      = [];
  
  for k = 1:min(max(A))
    j = find(A == k);
    Q = [Q j];
    
    [labk{k}, Pl{k}]  = spm_atlas('query',xA,XYZmm(:,j));
    Zj{k} = i1(j);
    XYZmmj{k} = XYZmm(:,j);
    maxZ(k) = sign(Zj{k}(1))*max(abs(Zj{k}));
  end

  % sort T/F values and print from max to min values
  [tmp, maxsort] = sort(maxZ,'descend');

  % use ascending order for neg. values
  indneg = find(tmp<0);
  maxsort(indneg) = flipud(maxsort(indneg));

  if ~isempty(maxsort)
    found_neg = 0;
    found_pos = 0;
    for l=1:length(maxsort)
      j = maxsort(l); 
      [tmp, indZ] = max(abs(Zj{j}));
    
      if ~isempty(indZ)
        if maxZ(j) < 0,  found_neg = found_neg + 1; end
        if maxZ(j) >= 0, found_pos = found_pos + 1; end
        
        if logP, valname = 'p-value'; else valname = 'Value'; end
        
        % print header if the first pos./neg. result was found
        if found_pos == 1
          fprintf('\n______________________________________________________');
          fprintf('\n%s: Positive effects\n%s',SO.img(2).vol.fname,atlas_name);
          fprintf('\n______________________________________________________\n\n');
          fprintf('%7s\t%12s\t%15s\t%s\n\n',valname,'Cluster-Size','    xyz [mm]   ','Overlap of atlas region');
        end
        if found_neg == 1
          fprintf('\n______________________________________________________');
          fprintf('\n%s: Negative effects\n%s',SO.img(2).vol.fname,atlas_name);
          fprintf('\n______________________________________________________\n\n');
          fprintf('%7s\t%12s\t%15s\t%s\n\n',valname,'Cluster-Size','    xyz [mm]   ','Overlap of atlas region');
        end
        
        if logP, val = 10^(-maxZ(j)); else val = maxZ(j); end
        fprintf('%7.2g\t%12d\t%4.0f %4.0f %4.0f',val,length(Zj{j}),XYZmmj{j}(:,indZ));
        for m=1:numel(labk{j})
          if Pl{j}(m) >= 1,
            if m==1, fprintf('\t%3.0f%%\t%s\n',Pl{j}(m),labk{j}{m});
            else     fprintf('%7s\t%12s\t%15s\t%3.0f%%\t%s\n','       ','       ','               ',...
              Pl{j}(m),labk{j}{m});
            end
          end
        end
      end
    end
  end
  fprintf('\n');
end

auto_savename = 0;
% save image
if ~isfield(OV, 'save')
  image_ext = spm_input('Save image file?', '+1', 'none|png|jpg|pdf|tif', char('none', 'png', 'jpg', 'pdf', 'tiff'), 2);
else
  if isempty(OV.save)
    image_ext = spm_input('Save image file?', '+1', 'none|png|jpg|pdf|tif', char('none', 'png', 'jpg', 'pdf', 'tiff'), 2);
  else
    [pp, nn, ee] = spm_fileparts(OV.save);
    if ~isempty(ee)
      image_ext = ee(2:end);
    else
      % if only the extension is given then automatically estimate filename for saving
      image_ext = OV.save;
      auto_savename = 1;
    end
  end
end

if ~strcmp(image_ext, 'none')
    
  [pt, nm] = spm_fileparts(img);
  if isempty(pt)
    pt1 = '';
  else
    [tmp,nm2] = spm_fileparts(pt);
    if isempty(nm2)
      pt1 = [pt '_']; 
    else
      pt1 = [nm2 '_']; 
    end
    
    [tmp,nm3] = spm_fileparts(spm_fileparts(pt));
    if isempty(nm3)
      pt2 = [pt1 '_']; 
    else
      pt2 = [nm3 '_']; 
    end
  end
  
  if ~isfield(OV, 'save')
    imaname = spm_input('Filename', '+1', 's', [pt1 nm '_' lower(OV.transform) '.' image_ext]);
  else
    if auto_savename
    
      % use up to 2 subfolders for getting filename
      if isfield(OV,'name_subfolder')
        % subfolder should be 0..2
        subfolder = max(min(OV.name_subfolder,3),0);
      else
        subfolder = 1;
      end
      
      % use up to 2 subfolders for getting filename
      switch subfolder
        case 0, pt1 = '';
        case 2, try, pt1 = [pt2 pt1]; end 
      end
      if numel(slices) == 1
        imaname = [pt1 nm '_' lower(OV.transform) num2str(slices) '.' image_ext];
      else
        imaname = [pt1 nm '_' lower(OV.transform) '.' image_ext];
      end
    else
      imaname = OV.save;
    end
  end
  
  % jpg needs full name to be accepted
  if strcmp(image_ext, 'jpg')
    image_ext = 'jpeg';
  end

  % and print
  H = findobj(get(SO.figure, 'Children'), 'flat', 'Type', 'axes');
  set(H, 'Units', 'normalized')
  
  saveas(SO.figure, imaname, image_ext);
  
  % read image, remove white border and save it again
  tmp = imread(imaname);
  sz = size(tmp);
  wborderx = 1;
  wbordery = sz(2);
  % we search for a white border inside a width of 4 pixels using
  % effect size. A line would be indicated by large mean and very low std
  for k=1:4
    if mean(double(tmp(k,:,1)))./(eps+std(double(tmp(k,:,1)))) > 25 
      wborderx = wborderx + 1;
    end
    if mean(double(tmp(:,sz(2)-k+1,1)))./(eps+std(double(tmp(:,sz(2)-k+1,1)))) > 25 
      wbordery = wbordery - 1;
    end
  end
  
  % save the image without borders
  imwrite(tmp(wborderx:sz(1),1:wbordery,:),imaname);
  fprintf('Image %s saved.\n', imaname);
  
  if ~isfield(SO,'overview')
    if n_slice > 0
      imaname = [lower(OV.transform) '_' replace_strings(OV.slices_str(ind, :)) '.' image_ext];
    else
      imaname = [lower(OV.transform) '.' image_ext];
    end
    
    saveas(h0, imaname, image_ext);
    fprintf('Image %s saved.\n', imaname);
  end
end

% check whether bounding box is from previous version that is not compatible
% with new template
BB = spm_get_bbox(SO.img(2).vol);
if sum(sum(BB-[-90 -126 -72;90 90 108])) == 0
  fprintf('WARNING: Check that %s is really compatible to new template in MNI152NLin2009cAsym template space. If not, you should use the old T1-template from CAT12.7 or older for overlay.\n',SO.img(2).vol.fname);
end

if ~nargout
  clear OV
end

% --------------------------------------------------------------------------
function xy = get_xy(n)

nn = round(n^0.4);
if n > 8, x = nn:round(n / nn); else x = 1:n; end
xy = [];
for i = 1:length(x)
    y = round(n / x(i));
  % check whether y is to small
  while y * x(i) < n, y = y + 1; end
  if i > 2
    if y * x(i - 1) < n, xy = [xy; [x(i) y]]; end
  else xy = [xy; [x(i) y]]; end
end

% change order of x and y
for i = 1:size(xy, 2)
  yx = [xy(i, 2) xy(i, 1)];
  xy = [xy; yx];
end

% remove duplicates
xy = [[n 1];xy];
xy = unique(xy, 'rows');
return

% --------------------------------------------------------------------------
function s = replace_strings(s)

s = deblank(s);
% replace spaces with "_" and characters like "<" or ">"
s(strfind(s, ' ')) = '_';
s(strfind(s, ':')) = '_';
s = spm_str_manip(s, 'v');

return

% --------------------------------------------------------------------------
function [mx, mn, XYZ, img] = volmaxmin(vol)

if nargout > 2
  XYZ = [];
end
if nargout > 3
  img = [];
end

mx = -Inf; mn = Inf;
for i = 1:vol.dim(3)
  tmp = spm_slice_vol(vol, spm_matrix([0 0 i]), vol.dim(1:2), [0 NaN]);
  tmp1 = tmp(isfinite(tmp(:)) & (tmp(:) ~= 0));
  if ~isempty(tmp1)
    if nargout > 2
      [Qc Qr] = find(isfinite(tmp) & (tmp ~= 0));
      if size(Qc, 1)
        XYZ = [XYZ; [Qc Qr i * ones(size(Qc))]];
        if nargout > 3
          img = [img; tmp1];
        end
      end
    end
    mx = max([mx; tmp1]);
    mn = min([mn; tmp1]);
  end
end

if nargout > 2
    XYZ = XYZ';
end

return

% --------------------------------------------------------------------------
function SO = pr_basic_ui(imgs, dispf)
% GUI to request parameters for slice_overlay routine
% FORMAT SO = pr_basic_ui(imgs, dispf)
%
% GUI requests choices while accepting many defaults
%
% imgs  - string or cell array of image names to display
%         (defaults to GUI select if no arguments passed)
% dispf - optional flag: if set, displays overlay (default = 1)
%
% $Id$

if nargin < 1
  imgs = '';
end
if isempty(imgs)
  imgs = spm_select(Inf, 'image', 'Image(s) to display');
end
if ischar(imgs)
  imgs = cellstr(imgs);
end
if nargin < 2
  dispf = 1;
end

clear global SO
global SO %#ok<REDEF> this is print as error

spm_clf('Interactive'); 
spm_input('!SetNextPos', 1);

% load images
nimgs = length(imgs);

% process names
nchars = 20;
imgns = spm_str_manip(imgs, ['rck' num2str(nchars)]);

SO.transform = deblank(spm_input('Image orientation', '+1', ['Axial|' ...
    ' Coronal|Sagittal'], strvcat('axial', 'coronal', 'sagittal'), 1));
orientn = find(strcmpi(SO.transform, {'sagittal', 'coronal', 'axial'}));

% identify image types
SO.cbar = [];
XYZ_unique = cell(3, 1);
for i = 1:nimgs
  SO.img(i).vol = spm_vol(imgs{i});
  if i == 1
    SO.img(i).cmap = gray;
    %[mx, mn] = volmaxmin(SO.img(i).vol);
    [tmp, th]=cat_stat_histth(spm_read_vols(SO.img(i).vol),0.95,0);

    SO.img(i).range = th;
  else
    [mx, mn, XYZ, img] = volmaxmin(SO.img(i).vol);
    if ~isempty(strfind(SO.img(i).vol.fname, 'logP')) || ~isempty(strfind(SO.img(i).vol.fname, 'log_'))
      logP = 1;
    else
      logP = 0;
    end
    
    SO.img(i).func = 'i1(i1==0)=NaN;';
    SO.img(i).prop = Inf;
    SO.cbar = [SO.cbar i];
    SO.img(i).cmap = return_cmap('Colormap:', 'jet');
    if logP
      % only ask for threshold if images is probably not thresholded 
      if mn < -log10(0.05)
        thresh = spm_input('Threshold P','+1','b','0.05|0.01|0.001',[0.05 0.01 0.001],1);
        mn = -log10(thresh);
        % use slightly larger maximum value to ensure that YTickLabel fits
        if thresh == 0.05
          mx = floor(mx) + 0.3011;
        end
      end
    end
        
    SO.img(i).range = spm_input('Image range for colormap', '+1', 'e', [mn mx], 2)';
    define_slices = spm_input('Slices', '+1', 'm', 'Estimate slices with local maxima|Define slices', [0 1], 1);
    
    if ~define_slices

      % for log-scaled p-values we should rather use gt than ge for comparison with threshold
      if logP
        compare_to_threshold = @(a,b) gt(a,b);
      else
        compare_to_threshold = @(a,b) ge(a,b);
      end
    
      % threshold map and restrict coordinates
      Q = find(compare_to_threshold(img,SO.img(i).range(1)) & le(img,SO.img(i).range(2)));
      XYZ = XYZ(:, Q);
      img = img(Q);
      
      M = SO.img(i).vol.mat;
      XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
      
      XYZ_unique = get_xyz_unique(XYZ, XYZmm, img);
    end
    
  end
end

% slices for display
ts = [0 0 0 pi / 2 0 -pi / 2 -1 1 1; ...
    0 0 0 pi / 2 0 0 1 -1 1; ...
    0 0 0 0 0 0 1 1 1];

V = SO.img(2).vol;
D = V.dim(1:3);
T = spm_matrix(ts(orientn, :)) * V.mat;
vcorners = [1 1 1; D(1) 1 1; 1 D(2) 1; D(1:2) 1; ...
    1 1 D(3); D(1) 1 D(3); 1 D(2:3); D(1:3)]';
corners = T * [vcorners; ones(1, 8)];

SO.slices = spm_input('Slices to display (mm)', '+1', 'e', XYZ_unique{orientn});
SO.figure = figure(22);

% and do the display
if dispf, slice_overlay; end

return

% --------------------------------------------------------------------------
function cmap = return_cmap(prompt, defmapn)
cmap = [];
while isempty(cmap)
  cmap = slice_overlay('getcmap', spm_input(prompt, '+1', 's', defmapn));
end
return

% --------------------------------------------------------------------------
function XYZ_unique = get_xyz_unique(XYZ, XYZmm, img)

xyz_array = [];

% cluster map
A = spm_clusters(XYZ);
for j = 1:max(A)
  ind = find(A == j);
  xyz = XYZmm(:, ind);
  xyz_array = [xyz_array xyz(:, img(ind) == max(img(ind)))];
end

% only keep unique coordinates
XYZ_unique = cell(3, 1);
for j = 1:3
  XYZ_unique{j} = unique(xyz_array(j, :));
end

return

%==========================================================================
function s = remove_zeros(s)

s = deblank(s);
pos = length(s);
while pos > 1
  if strcmp(s(pos), '0')
    s(pos) = '';
    pos = pos - 1;
  else break
  end
end

return