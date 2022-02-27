function cat_vol_img2mip(P, gamma_scl, save_image, RGB_order)
% Visualise up to 3 images as RGB colored MIP (glass brain)
% ______________________________________________________________________
%
% P          - char array of 1..3 nifti filenames
% gamma_scl  - gamma value to provide non-linear intensity scaling
% save_image - save mip as png image
% RGB_order  - array of RGB order (default [1 2 3])
%
% If < 3 arguments are given, you can save the png-file by using the
% context menu (right mouse click in the image)
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

load MIP
mip96 = double(mip96);

if nargin < 1
  P = spm_select([1 3],'image','Select images');
end

n = size(P,1);
if n > 3
  error('At maximum 3 images are allowd.');
end

% define gamma to provide non-linear scaling
if nargin < 2
  gamma_scl = 1;
end

if nargin < 3
  save_image = 0;
end

if nargin < 4
  RGB_order = 1:3;
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
col = col(RGB_order);
for i=1:n
  fprintf('%s color: %s\n',col{i},V(i).fname);
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
	  d{i}  = spm_slice_vol(V(i),M1,V(i).dim(1:2),1);
    d{i} = flipud(d{i});
	  Q = find(d{i} > 0);
    
	  if ~isempty(Q)
		  [Qc Qr] = find(d{i} > 0);
		  Qc = (Qc - Origin(1))*vox(1);
		  Qr = (Qr - Origin(2))*vox(2);
		  XYZ{i} = [XYZ{i}; [Qc Qr ...
			ones(size(Qc,1),1)*(j - Origin(3))*vox(3)]];
		  Y{i} = [Y{i}; d{i}(Q)];
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
  Y{i} = Y{i}';
  rgb{i}   = rot90(spm_project(Y{i},round(XYZ{i}),dim));
  if gamma_scl ~= 1
	rgb{i} = rgb{i}.^(1/gamma_scl);
  end
  mx = max([mx; rgb{i}(:)]);
end

for i=1:n
  rgb{i} = rgb{i}/mx;
end

rgb{1}(zx,zy-20) = rgb{1}(zx,zy-20) + yy;
if n > 1
  rgb{2}(zx,zy) = rgb{2}(zx,zy) + yy;
end
if n > 2
  rgb{3}(zx-14,zy-10) = rgb{3}(zx-14,zy-10) + yy;
end

col = reshape([rgb{RGB_order(1)},rgb{RGB_order(2)},rgb{RGB_order(3)}],size(mip));
mip = max(cat(3,col),mip);

% show mip image
fig = figure(11);
[pt,nm,xt] = fileparts(P(1,:));
png_name = ['mip' num2str(n) '_' nm '.png'];
set(fig, 'MenuBar', 'none','Name',nm);

p = image(mip);
axis image; axis off;

if nargin < 3

  cmenu = uicontextmenu(fig);        
  m2 = uimenu(cmenu, 'Label','Save png image','Callback',@(src,event)save_png(mip,png_name));
  p.ContextMenu = cmenu;
else
  if save_image
    save_png(mip,png_name);
  end
end
              
function save_png(mip,png_name)
imwrite(mip,png_name,'png','BitDepth',8)
fprintf('Image %s saved.\n',png_name);
