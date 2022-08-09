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

load MIP
mip96 = double(mip96);

if nargin < 1
  OV = spm_select([1 3],'image','Select images');
end

def = struct('name',OV,'func','i1(i1==0)=NaN;','cmap',jet(64),'range',...
  [-Inf Inf],'gamma_scl',0.7,'save_image','','RGB_order',1:3,'Pos',...
  [10 10], 'bkg_col',[0 0 0]);
    
if ischar(OV)
  OV = def;
else
  OV = cat_io_checkinopt(OV,def); 
  if ~isfield(OV,'cbar')
    OV.cbar = 1;
  end
end
P = OV.name;

n = size(P,1);
if n > 3
  error('At maximum 3 images are allowed.');
end

V = spm_vol(P);
M   = V(1).mat;

mip  = repmat(full(rot90(mip96/max(mip96(:)))),[1 1 3]);
c0   = [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1
        1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1] -0.5;
c   = (M(1:3,1:3)*c0')';
dim = [(max(c)-min(c)) size(mip96)];

% compute colored spheres for colorbar
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
for i=1:n
  if n > 1
    fprintf('%s color: %s\n',col{i},V(i).fname);
  end
  XYZ{i} = [];
  Y{i} = [];
end

Origin = V(1).mat\[0 0 0 1]';
vox = sqrt(sum(V(1).mat(1:3,1:3).^2));

spm_progress_bar('Init',V(1).dim(3),'Mip',' ');
for j = 1:V(1).dim(3)
  B  = spm_matrix([0 0 -j 0 0 0 1 1 1]);
  M1 = inv(B);
  
  for i=1:n
    % read slice and flip for MIP
    i1  = spm_slice_vol(V(i),M1,V(i).dim(1:2),1);
    i1 = flipud(i1);
    
    % apply defined function
    eval(OV.func)

    % find indices in defined range
    [Qc Qr] = find(i1 >= OV.range(1) & i1 <= OV.range(2));
    Q = sub2ind(size(i1),Qc,Qr);
        
    if ~isempty(Q)
      Qc = (Qc - Origin(1))*vox(1);
      Qr = (Qr - Origin(2))*vox(2);
      XYZ{i} = [XYZ{i}; [Qc Qr ones(size(Qc,1),1)*(j - Origin(3))*vox(3)]];
      
      % if finite lower range is defined this should be subtracted from
      % image
      if isfinite(OV.range(1))
        i1(Q) = i1(Q) - OV.range(1);
      end
      
      Y{i} = [Y{i}; i1(Q)];
    end
  end
  
  spm_progress_bar('Set',j);
end

spm_progress_bar('Clear');

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

% show mip image
fig = figure(12);
[pt,nm,xt] = fileparts(P(1,:));
png_name = ['mip' num2str(n) '_' nm '.png'];
sz = size(mip(:,:,1));
set(fig, 'MenuBar', 'none','Name',nm,'Position',[10 10 sz([2 1])]);

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
  colormap(OV.cmap)
  
  if ~isempty(OV.cbar)
    t = text(320,230,'max');
    set(t,'Color',1 - OV.bkg_col);
    if OV.range(1) > 1
      t = text(320,329,'min');
    else
      t = text(320,329,'-max');
    end
    set(t,'Color',1 - OV.bkg_col);
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
  m2 = uimenu(cmenu, 'Label','Save png image','Callback',@(src,event)save_png(mip,OV.save_image));
  try, p.ContextMenu = cmenu; end
else
  save_png(mip,OV.save_image);
end
           
function save_png(mip,png_name)

hh = getframe(gcf);
img = frame2im(hh);
imwrite(img,png_name,'png','BitDepth',8);

fprintf('Image %s saved.\n',png_name);
