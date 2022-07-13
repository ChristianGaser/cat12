function cat_plot_cov(data,opt)
% Plot rotated covariance matrix
% FORMAT cat_plot_cov(data,opt)
% data        - covariance/correlation matrix
% opt
%   .ax       - currect axes
%   .name     - tick labels
%   .group    - group coding (1..k) for different samples (indicated with
%               different colors)
%   .color    - background color
%   .pos_cbar - position of colorbar
%   .ctick    - number of tick labels for colorbar
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% sample size
[m,n] = size(data);

% data should be quadratic matrix
if m ~= n
  error('Data should be a covariance or correlation matrix');
end

% get axes
if isfield(opt,'ax')
  ax = opt.ax;
else
  ax = gca;
end

% get group for coloring names for different samples
if isfield(opt,'group')
  group = opt.group;
else
  group = ones(n,1);
end

% get position of axes
pos0 = get(ax,'Position');

% tick labels
if isfield(opt,'name')
  name = opt.name;
else
  name = num2str(1:n);
end

% background color
if isfield(opt,'color')
  col = opt.color;
else
  col = [0.8 0.8 0.8];
end

% position of colorbar
if isfield(opt,'pos_cbar')
  pos = opt.pos_cbar;
else
  pos = [pos0(1) 0.950 pos0(3) 0.025*pos0(4)];
end

% number of tick labels for colorbar
if isfield(opt,'ctick')
  ctick = opt.ctick;
else
  ctick = 5;
end

shift = min(max(0.08,1/(13 - 45/n)),0.4);
shift = pos0(3)/(0.0692*n+8.871);

% scale data to min..max
mn = min(data(:));
mx = max(data(data~=1));

% show only lower left triangle
data_scaled = tril((data - mn)/(mx - mn));

image(64*data_scaled);
set(gca,'XTickLabel',[],'YTickLabel',[],'Color',col);
view([-45 90])
set(get(ax,'children'),'AlphaData',tril(isfinite(data)));

axis image
axis off

% colorbar
H.cbar = axes('Position',pos,'Parent',gcf);
image(1:64);

set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,ctick), 'XTickLabel',...
  round(100*linspace(min(data(:)),max(data(data~=1)),ctick))/100,'TickLength',[0 0]);

% axes of rotated tick labels
pos = [pos0(1)+shift pos0(2)+pos0(4)/1.8 pos0(3)-2*shift pos0(4)/2];
ax = axes('Position',pos,'Parent',gcf,'Color',col);

% position of tick labels
pos_text = linspace(0,1,n);

% limit ticks to maximum of 50
gap = round(n/min(n,50));

% colormap for groups
if exist('lines')
  cm = lines(max(group));
else
  cm = jet(max(group));
end

for i=1:gap:n
  t = text(pos_text(i),0.15,name(i,:),'parent',ax,'interpreter','none','Color',0.8*cm(group(i),:));
  set(t,'Units','normalized','VerticalAlignment','middle','HorizontalAlignment','left','Rotation',90);
  if isfield(opt,'fontsize')
    set(t,'FontSize',opt.fontsize);
  end
end
axis off
