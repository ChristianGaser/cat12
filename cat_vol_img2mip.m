function cat_vol_img2mip(OV)
% Visualise up to 3 images as RGB colored MIP (glass brain)
% ______________________________________________________________________
%
% OV - either char array of 1..3 nifti filenames or structure with
%              the following fields:
%  name       - char array of 1..3 nifti filenames
%  cmap       - colormap for single MIP wth just one result (default
%               jet(64))
%  func       - function to apply to image before scaling to cmap
%               (and therefore before min/max thresholding. E.g. a func of
%               'i1(i1==0)=NaN' would convert zeros to NaNs
%               (default 'i1(i1==0)=NaN')
%  range      - 2x1 vector of values for image to distribute colormap across
%               the first row of the colormap applies to the first
%               value in 'range', and the last value to the second
%               value in 'range' (default [-Inf Inf])
%  gamma_scl  - gamma value to provide non-linear intensity
%               scaling (default 0.7)
%  save_image - save mip as png image (default '')
%  RGB_order  - array of RGB order (default [1 2 3])
%  bkg_col    - color of background ([0 0 0] as default)
%  cbar       - if empty skip showing colorbar
%  fig_mip    - figure number (default 12)
%
%  style      - MIP style:
%               0 - use old glassbrain
%               1 - use cat_vol_glassbrain
%               options 1 is only available for a single filename
%
% If < 3 arguments are given, you can save the png-file by using the
% context menu (right mouse click in the image)
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin < 1
  OV = spm_select([1 3],'image','Select images');
end

def = struct('name',OV,'func','i1(i1==0)=NaN;','cmap',jet(64),'range',...
  [-Inf Inf],'gamma_scl',0.7,'save_image','','RGB_order',1:3,'Pos',...
  [10 10],'bkg_col',[0 0 0],'fig_mip',12,'cbar',2,'roi',[]);

style = 1;
    
if ischar(OV)
  OV = def;
else
  OV = cat_io_checkinopt(OV,def); 
end

if ischar(OV.name)
  V = spm_vol(OV.name);
  dim = V(1).dim;
  fname = V(1).fname;
elseif isa(OV.name,'nifti')
  V = OV.name;
  dim = V(1).dat.dim;
  fname = V(1).dat.fname;
else
  error('You have to either define a nifti structure or a filename as input.');
end

n = numel(V);
if ~isempty(OV.roi)
  if ~isa(OV.name,'nifti')
    error('The ROI option cannot only be used in conjunction with nifti structure as OV.');
  end
  if n > 1
    error('The ROI option cannot be used for multiple images.');
  end
  if any(size(OV.roi) ~= dim)
    error('Dimension of ROI and image data differs.');
  end
end

% new glassbrain does not yet support RGB MIP
if n > 1
  style = 0;
end

if n > 3
  error('At maximum 3 images are allowed.');
end

% we use different affine corrections to better fit MN2009cAsym into the
% existing MIPs
if style
  Affine = spm_matrix([0.5 0 -1.5 0 0 0 1 1 0.98 0 0 0]);
else
  Affine = spm_matrix([0.5 2.5 3.5 0 0 0 0.95 0.95 0.95 0.1 0.1 0.1]);     
end
for i=1:n
  V(i).mat = Affine * V(i).mat;
end

M   = V(1).mat;

col = {'R','G','B'};
for i=1:n
  if n > 1
    fprintf('%s color: %s\n',col{OV.RGB_order(i)},V(i).fname);
  end
  XYZ{i} = [];
  Y{i}   = [];
  if ~isempty(OV.roi)
    ROI{i} = logical([]);
  end
end

Origin = V(1).mat\[0 0 0 1]';
vox    = sqrt(sum(V(1).mat(1:3,1:3).^2));

mnI = 1e15; mxI = -1e15;
cat_progress_bar('Init',dim(3),'Mip',' ');
for j = 1:dim(3)
  B  = spm_matrix([0 0 -j 0 0 0 1 1 1]);
  M1 = inv(B);
  
  for i=1:n
    % read slice and flip for MIP
    if isa(V(i),'nifti')
      i1  = V(i).dat(:,:,j);
    else
      i1  = spm_slice_vol(V(i),M1,V(i).dim(1:2),1);
    end
    i1 = flipud(i1);

    % apply defined function
    eval(OV.func)
        
    % find indices
    if ~isempty(OV.roi)
      roi = flipud(OV.roi(:,:,j));
      [Qc Qr] = find((isfinite(i1) & i1 ~=0) | roi ~= 0);
    else
      [Qc Qr] = find(isfinite(i1) & i1 ~=0);
    end
    
    Q = sub2ind(size(i1),Qc,Qr);
        
    if ~isempty(Q)
      Qc = (Qc - Origin(1))*vox(1);
      Qr = (Qr - Origin(2))*vox(2);
      XYZ{i} = [XYZ{i}; [Qc Qr ones(size(Qc,1),1)*(j - Origin(3))*vox(3)]];
      
      mnI = min(mnI, min(i1(Q)));
      mxI = max(mxI, max(i1(Q)));
      
      Y{i} = [Y{i}; i1(Q)];
      if ~isempty(OV.roi)
        ROI{i} = [ROI{i}; roi(Q)];
      end
    end
  end
  
  cat_progress_bar('Set',j);
end
cat_progress_bar('Clear');

[pt,nm,xt] = fileparts(fname);

if style

  % show new glassbrain
  if ishandle(OV.fig_mip)
    fig = figure(OV.fig_mip);
    set(fig, 'Name',nm);
  else
    fig = figure(OV.fig_mip);
    png_name = ['mip' num2str(n) '_' nm '.png'];
    set(fig, 'MenuBar', 'none','Position',[10 10 2*182 2*200],'Name',nm,'NumberTitle','off');
  end
  
  if isempty(OV.cbar), OV.cbar = 0; end
  S = struct('dark',all(OV.bkg_col==0),'cmap',OV.cmap,'grid',false,'colourbar',OV.cbar,...
    'order',style);
  if ~isempty(OV.roi)
    S.roi = ROI{1};
  end
  cat_vol_glassbrain(Y{1},XYZ{1},S);

  set(gca,'units','pixels'); x = get(gca,'position');
  set(gcf,'units','pixels'); y = get(gcf,'position');
  set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
  set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
  
  if isempty(OV.save_image)
    cmenu = uicontextmenu(fig);
    m2 = uimenu(cmenu, 'Label','Save png image','Callback',@(src,event)save_png(OV.save_image));
    try, p.ContextMenu = cmenu; end
  else
    save_png(OV.save_image);
  end
  
  return
else
  load MIP
  mip96 = double(mip96);
  mip  = repmat(full(rot90(mip96/max(mip96(:)))),[1 1 3]);
  c0   = [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1
          1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1] -0.5;
  c   = (M(1:3,1:3)*c0')';
  dim = [(max(c)-min(c)) size(mip96)];
  
  % compute colored spheres for colorbar
  if n > 1
    s = 35;
    coord = [270 310];
    [x,y,z] = sphere(4*s);
    yy = y'*y;
    yy = yy(1:2*s,1:2*s);
    yy = yy/max(yy(:));
    zx = coord(1)-s:coord(1)+s-1;
    zy = coord(2)-s:coord(2)+s-1;

    col = {'R','G','B'};
    col = col(OV.RGB_order);
  end
end

for i=1:n
   % apply range
  if isfinite(OV.range(1))
    Y{i}(Y{i}<OV.range(1)) = OV.range(1);
  end

  if isfinite(OV.range(2))
    Y{i}(Y{i}>OV.range(2)) = OV.range(2);
  end

  % for some range min combinations we have to subtract the minimum
  if isfinite(OV.range(1)) && ( (mnI < 0 && OV.range(1) <= 0) || mnI > 0)
    Y{i} = Y{i} - OV.range(1);
%    Y{i} = Y{i} - mnI;
  end

end

sz = size(mip);
for i=1:3
  rgb{i} = zeros(sz(1:2));
end

mx = 0;
for i=1:n
  XYZ{i} = XYZ{i}';
  if ~isempty(Y{i})
    Y{i} = Y{i}';
    rgb{i}   = rot90(spm_project(Y{i},round(XYZ{i}),dim));
  end
  if OV.gamma_scl ~= 1
  rgb{i} = rgb{i}.^(1/OV.gamma_scl);
  end
  mx = max([mx; rgb{i}(:)]);
end

for i=1:n
  rgb{i} = rgb{i}/mx;
end

% just draw colorbar for a single image
if n == 1 && ~isempty(OV.cbar)
  rgb{1}(230:329,305:315) = repmat((100:-1:1)'/100,1,11);
end

% use RGB color-spheres
if n > 1 && ~isempty(OV.cbar)
  rgb{1}(zx,zy-20) = rgb{1}(zx,zy-20) + yy;
  rgb{2}(zx,zy) = rgb{2}(zx,zy) + yy;
end
if n > 2 && ~isempty(OV.cbar)
  rgb{3}(zx-14,zy-10) = rgb{3}(zx-14,zy-10) + yy;
end

col = reshape([rgb{OV.RGB_order(1)},rgb{OV.RGB_order(2)},rgb{OV.RGB_order(3)}],size(mip));
old_mip = mip;
mip = max(cat(3,col),mip);

sz = size(mip(:,:,1));

% show mip image
if ishandle(OV.fig_mip)
  fig = figure(OV.fig_mip);
else
  fig = figure(OV.fig_mip);
  set(fig, 'MenuBar', 'none','Name',nm,'Position',[10 10 sz([2 1])],'NumberTitle','off');
end

% single MIP of one result supports colored MIP
if n == 1
  % only use 1st channel
  mip = mip(:,:,1);
  
  % add MIP grid with white color
  mip(old_mip(:,:,1)>0) = 1.1; 
  
  p = imagesc(mip);
  
  % force defined background and inverted MIP grid
  OV.cmap(1,:)     = OV.bkg_col;
  OV.cmap(end+1,:) = 1 - OV.bkg_col;
  
  if mxI < 0
    colormap(OV.cmap(1:round(size(OV.cmap,1)/2)-1,:));
  else
    colormap(OV.cmap);
  end
  
  if ~isempty(OV.cbar)
    
    if mxI < 0
      t1 = text(320,230,'0');
      t2 = text(320,329,'-max');
    elseif mnI > 0
      t1 = text(320,230,'max');
      t2 = text(320,329,'0');
    else
      t1 = text(320,230,'max');
      t2 = text(320,329,'-max');
    end
  
    set(t1,'Color',1 - OV.bkg_col);
    set(t2,'Color',1 - OV.bkg_col);
  end
  
else
  p = image(mip);
end

axis image; axis off;

% remove border
set(gca,'units','pixels'); x = get(gca,'position');
set(gcf,'units','pixels'); y = get(gcf,'position');
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

if isempty(OV.save_image)
  cmenu = uicontextmenu(fig);
  m2 = uimenu(cmenu, 'Label','Save png image','Callback',@(src,event)save_png(OV.save_image));
  try, p.ContextMenu = cmenu; end
else
  save_png(mip,OV.save_image);
end
           
function save_png(png_name)

hh = getframe(gcf);
img = frame2im(hh);
imwrite(img,png_name,'png','BitDepth',8);

fprintf('Image %s saved.\n',png_name);
