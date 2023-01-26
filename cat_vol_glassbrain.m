function fig = cat_vol_glassbrain(X,pos,varargin)
% Glass brain plot
% FORMAT fig = cat_glass(X,pos,S)
%   X               - (REQUIRED) values to be painted
%   pos             - (REQUIRED) coordinates in MNI head (not voxel) space
%   S               - (optional) config structure
% Fields of S:
%   S.cmap          - colormap of plot             - Default: 'gray'
%   S.dark          - dark mode                    - Default: false
%   S.detail        - glass brain detail level:
%                     0=LOW, 1=NORMAL, 2=HIGH      - Default: 1
%   S.grid          - overlay grid                 - Default: false
%   S.colourbar     - show colorbar                - Default: false
%                     0 - no, 1 - with string min/max, 2 - with min/max value 
% Output:
%   fig             - Handle for generated figure
%__________________________________________________________________________

% modified version of spm_glass
% George O'Neill & Guillaume Flandin
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging
%---------------------------------------------------------------------
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

switch nargin
  case 1
    error('need at least two arguments, values, and positions!')
  case 2
    S = [];
  case 3
    S = varargin{1};
end

assert(length(X) == length(pos), ['number of values do not match '...
  'number of positions!']);

if ~isfield(S, 'dark'),     S.dark = false; end
if ~isfield(S, 'cmap'),     S.cmap = 'gray'; end
if ~isfield(S, 'detail'),   S.detail = 1; end
if ~isfield(S, 'grid'),     S.grid = false; end
if ~isfield(S, 'colourbar'),S.colourbar = 0; end

M = [-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1];
dim = [91 109 91];
pos = ceil(M \ [pos';ones(1,size(pos,1))])';

% exclude positions that exceed image dimensions
ind = find(pos(:,1) < 1 | pos(:,1) > dim(1) | pos(:,2) < 1 | pos(:,2) > dim(2) | pos(:,3) < 1 | pos(:,3) > dim(3));
if ~isempty(ind)
  pos(ind,:) = [];
  X(ind) = [];
end

if any(X<0) && any(X>0)
  div = 1;
else
  div = 0;
end

[~,id] = sort(abs(X),'ascend');
if div
  [~,bin] = histc(X,linspace(-max(abs(X)),max(abs(X)),65));
else
  [~,bin] = histc(X,linspace(min(abs(X)),max(abs(X)),65));
end

% saggital plane
%----------------------------------------------------------------------
p_sag = NaN(dim(2),dim(3));

for ii = 1:length(id)
    
  p1 = pos(id(ii),2);
  p2 = pos(id(ii),3);
  
  if p1 > 0 && p1 <= dim(2) && p2 > 0 && p2 <= dim(3)
    p_sag(p1,p2) = bin(id(ii));
  end
  
end

% coronal plane
%----------------------------------------------------------------------
p_cor = NaN(dim(1),dim(3));

for ii = 1:length(id)
    
  p1 = pos(id(ii),1);
  p2 = pos(id(ii),3);
  
  if p1 > 0 && p1 <= dim(1) && p2 > 0 && p2 <= dim(3)
    p_cor(p1,p2) = bin(id(ii));
  end  
end


% axial plane
%----------------------------------------------------------------------
p_axi = NaN(dim(2),dim(1));

for ii = 1:length(id)
    
  p1 = pos(id(ii),2);
  p2 = pos(id(ii),1);
  
  if p1 > 0 && p1 <= dim(2) && p2 > 0 && p2 <= dim(1)
    p_axi(p1,p2) = bin(id(ii));
  end
end

% optional colorbar
%---------------------------------------------------------------------
p_col = NaN(dim(1),dim(3));
if S.colourbar
  for ii = 42:52
    p_col(ii,14:78) = linspace(2,65,numel(14:78));
  end
  if div
    rmin = -max(abs(X));
  else
    rmin = min(abs(X));
  end
  rmax = max(abs(X));
end


% combine and plot
%---------------------------------------------------------------------
p_all = [rot90(p_sag,1) fliplr(rot90(p_cor,1));...
rot90(p_axi,1) rot90(p_col,1)];
p_all(isnan(p_all)) = 0;
imagesc(p_all)
set(gca,'XTickLabel',{},'YTickLabel',{});
axis image

caxis([0 64])
load(fullfile(fileparts(mfilename('fullpath')),'glass_brain.mat'));
overlay_glass_brain(glass,'side',S.dark,S.detail);
overlay_glass_brain(glass,'back',S.dark,S.detail);
overlay_glass_brain(glass,'top', S.dark,S.detail);

if ischar(S.cmap)
  c = feval(S.cmap,64);
else
  c = S.cmap;
end

if S.dark
  c(1,:) = [0 0 0];
else
  c(1,:) = [1 1 1];
end
colormap(c);

if S.dark
  set(gcf,'color','k');
else
  set(gcf,'color','w');
end

if S.colourbar > 1
  text(175,170,sprintf('%0.2f',rmin),'color',~c(1,:),'fontsize',12,'horizontalalignment','center');
  text(175,105,sprintf('%0.2f',rmax),'color',~c(1,:),'fontsize',12,'horizontalalignment','center');
elseif S.colourbar == 1
  if rmin < 0
    text(175,170,'-max','color',~c(1,:),'fontsize',12,'horizontalalignment','center');
  else
    text(175,170,'min','color',~c(1,:),'fontsize',12,'horizontalalignment','center');
  end
    
  text(175,105,'max','color',~c(1,:),'fontsize',12,'horizontalalignment','center');
end

if S.grid
  grid on
else
  axis off
end

fig = gcf;

end

% supporting functions
%---------------------------------------------------------------------

function overlay_glass_brain(glass,orient,dark,detail)

dat = glass.(orient);

switch orient
  case 'top'
    xform = [0 -1 0; 1 0 0; 0 0 1]*[0.185 0 0; 0 0.185 0; 10.5 173 1];
  case 'back'
%    xform = [0.185 0 0; 0 -0.185 0; 120 89 1];
    % xform was distorted for coronal slice
    xform = [1.05 0 0; 0 1.05 0; 0 0 1]*[0.185 0 0; 0 -0.185 0; 118 92 1];
  case 'side'
    xform = [0.185 0 0; 0 -0.185 0; 10.5 89 1];
end

for ii = 1:length(dat.paths)
  pth = dat.paths(ii);
  % see if we need to draw based on the complexity option
  switch detail
    case 0
      draw = pth.linewidth > 1 & sum(hex2rgb(pth.edgecolor))==0;
    case 1
      draw = sum(hex2rgb(pth.edgecolor))==0;
    otherwise
      draw = 1;
  end
  
  if draw
    for jj = 1:length(pth.items)
      pts = pth.items(jj).pts;
      v = [generate_bezier(pts) ones(10,1)];
      v2 = v*xform;
      if dark
        c = 1 - hex2rgb(pth.edgecolor);
      else
        c = hex2rgb(pth.edgecolor);
      end
      line(v2(:,1),v2(:,2),'LineWidth',pth.linewidth,'Color',c);
    end
  end
end

end

function [points, t] = generate_bezier(controlPts, varargin)

% bezier generation from control points based on code by
% Adrian V. Dalca, https://www.mit.edu/~adalca/
% https://github.com/adalca/bezier

% estimate nDrawPoints
if nargin == 1
  nCurvePoints = 10;
else
  nCurvePoints = varargin{1};
end

% curve parametrization variable
t = linspace(0, 1, nCurvePoints)';

% detect the type of curve (linear, quadratic, cubic) based on the
% number of points given in controlPts.
switch size(controlPts, 1)
  case 1
    error('Number of Control Points should be at least 2');
    
  case 2
    % linear formula
    points = (1 - t) * controlPts(1, :) + ...
      t * controlPts(2, :);
    
  case 3
    % quadratic formula
    points = ((1 - t) .^ 2) * controlPts(1, :) + ...
      (2 * (1 - t) .* t) * controlPts(2, :) + ...
      (t .^ 2) * controlPts(3, :);
    
  case 4
    % cubic formula
    points =  ((1 - t) .^ 3) * controlPts(1, :) + ...
      (3 * (1 - t) .^ 2 .* t) * controlPts(2, :) + ...
      (3 * (1 - t) .* t .^ 2) * controlPts(3, :) + ...
      (t .^ 3) * controlPts(4, :);
end

% verify dimensions
assert(size(points, 2) == size(controlPts, 2));
end

function rgb = hex2rgb(hex)
% converts hex string to matlab rgb triplet
if strcmpi(hex(1,1),'#')
  hex(:,1) = [];
end
rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
end
