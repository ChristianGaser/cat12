function h = cat_plot_circular(data, opt)

% cat_plot_curcular Plots correlation matrix as a cat_plot_circular
%
% usage: vargout = cat_plot_boxplot(data,opt);
%
%   cat_plot_curcular(data) data is a square numeric matrix with values in [0,1].
%
%                 NOTE: only the off-diagonal lower triangular section of data is
%                       considered, i.e. tril(data,-1).
%                                               
%     opt.lbls      Plot a cat_plot_circular with custom labels at each node.
%                   LBLS is either a cellstring of length M, where
%                   M = size(r,1), or a M by N char array, where each
%                   row is a label.
%
%     opt.ccolor    Supply an RGB triplet or a colormap that specifies the color of
%                   the curves. CURVECOLOR can also be a 2 by 3 matrix
%                   with the color in the first row for negative
%                   correlations and the color in the second row for
%                   positive correlations.
%
%     opt.ncolor    Change color of the nodes with an RGB triplet.
%
%
%   H = cat_plot_curcular(...) Returns a structure with handles to the graphic objects
%
%       h.l     handles to the curves (line objects), one per color shade. 
%               If no curves fall into a color shade that handle will be NaN.
%       h.s     handle  to the nodes (scattergroup object)
%       h.t     handles to the node text labels (text objects)
%
%
% Examples
%
%   % Base demo
%   cat_plot_circular
%
%   % Supply your own correlation matrix (only lower off-diagonal triangular part is considered)
%   x = rand(10).^3;
%   x(:,3) = 1.3*mean(x,2);
%   cat_plot_circular(x)
%
%   % Supply custom labels as ['aa'; 'bb'; 'cc'; ...] or {'Hi','how','are',...}
%   cat_plot_circular(x, repmat(('a':'j')',1,2))
%   cat_plot_circular(x, {'Hi','how','is','your','day?', 'Do','you','like','cat_plot_circulars?','NO!!'})
%
%   % Customize curve colors
%   cat_plot_circular([],[],[1,0,1;1 1 0])
%
%   % Customize node color
%   cat_plot_circular([],[],[],[0,1,0])
%
%   % Customize manually other aspects
%   h   = cat_plot_circular;
%   set(h.l(~isnan(h.l)), 'LineWidth',1.2)
%   set(h.s, 'MarkerEdgeColor','red','LineWidth',2,'SizeData',100)
%
%
% Additional features:
% - <a href="matlab: web('http://www.mathworks.com/matlabcentral/fileexchange/42279-schemaball','-browser')">FEX schemaball page</a>
% - <a href="matlab: web('http://www.stackoverflow.com/questions/17038377/how-to-visualize-correlation-matrix-as-a-schemaball-in-matlab/17111675','-browser')">Origin: question on Stackoverflow.com</a>
% - <a href="matlab: web('https://github.com/GuntherStruyf/matlab-tools/blob/master/schemaball.m','-browser')">Schemaball by Gunther Struyf</a>
%
% See also: CORR, CORRPLOT
%
% modified version of schemaball.m
% Author: Oleg Komarov (oleg.komarov@hotmail.it)
% Tested on R2013a Win7 64 and Vista 32
% 15 jun 2013 - Created
% _________________________________________________________________________
% $Id: cat_plot_boxplot.m 831 2016-01-14 14:00:03Z gaser $

% Number of color shades/buckets (large N simply creates many perceptually indifferent color shades)
N      = 20;

% Points in [0, 1] for bezier curves: leave space at the extremes to detach a bit the nodes. 
% Smaller step will use more points to plot the curves.
t      = (0.025: 0.05 :1)';

% Nodes edge color
ecolor = [.25 .103922 .012745];

% Text color
tcolor = [.7 .7 .7];

% Some defaults
if nargin < 1 || isempty(data);        data      = (rand(50)*2-1).^29;                                  end
sz = size(data);

% default parameter
if ~exist('opt','var'), opt = struct(''); end
def.ncolor         = [0 0 1];
def.maxlinewidth   = 3;
def.lbls           = cellstr(reshape(sprintf('%-4d',1:sz(1)),4,sz(1))');
def.ccolor         = hsv2rgb([[linspace(.8333, .95, N); ones(1, N); linspace(1,0,N)],...
                              [linspace(.03, .1666, N); ones(1, N); linspace(0,1,N)]]');
  
opt = cat_io_checkinopt(opt,def);

% data
if ~isnumeric(data) || any(abs(data(:)) > 1) || sz(1) ~= sz(2) || numel(sz) > 2 || sz(1) == 1
    error('cat_plot_circular:validR','data should be a square numeric matrix with values in [0, 1].')
end

% lbls
if (~ischar(opt.lbls) || size(opt.lbls,1) ~= sz(1)) && (~iscellstr(opt.lbls) || ~isvector(opt.lbls) || length(opt.lbls) ~= sz(1))
    error('cat_plot_circular:validLbls','LBLS should either be an M by N char array or a cellstring of length M, where M is size(R,1).')
end
if ischar(opt.lbls)
    opt.lbls = cellstr(opt.lbls);
end

% ccolor
if isempty(opt.ccolor);
    opt.ccolor = hsv2rgb([[linspace(.8333, .95, N); ones(1, N); linspace(1,0,N)],...
                      [linspace(.03, .1666, N); ones(1, N); linspace(0,1,N)]]');
else
    szC = size(opt.ccolor);
    if ~isnumeric(opt.ccolor) || szC(2) ~= 3  
        error('cat_plot_circular:validCcolor','CCOLOR should be a 1 by 3, 2 by 3 or N by 3 numeric matrix with RGB colors.')
    elseif szC(1) == 1
        opt.ccolor = [opt.ccolor; opt.ccolor];
    end
    if szC(1) < 3
      opt.ccolor = rgb2hsv(opt.ccolor);
      opt.ccolor = hsv2rgb([repmat(opt.ccolor(1,1:2),N,1), linspace(opt.ccolor(1,end),0,N)';
                            repmat(opt.ccolor(2,1:2),N,1), linspace(0,opt.ccolor(2,end),N)']);
    end
end

% ncolor
szN = size(opt.ncolor);
if ~isnumeric(opt.ncolor) || szN(2) ~= 3  || szN(1) > 1
    error('cat_plot_circular:validNcolor','NCOLOR should be a single RGB color, i.e. a numeric row triplet.')
end
opt.ncolor = rgb2hsv(opt.ncolor);
%% Engine

% Create figure
figure('renderer','zbuffer','visible','off')
axes('NextPlot','add')

% Index only low triangular matrix without main diag
tf        = tril(true(sz),-1);

% Index correlations into bucketed colormap to determine plotting order (darkest to brightest)
N2        = 2*N;
data(data == 0) = NaN;
[n, isrt] = histc(data(tf), linspace(-1,1 + eps(100),N2 + 1));
plotorder = reshape([N:-1:1; N+1:N2],N2,1);

% create varying linewidth
linewidth = [linspace(opt.maxlinewidth,1,N) linspace(1,opt.maxlinewidth,N)];

% Retrieve pairings of nodes
[row,col] = find(tf);

% Use tau http://tauday.com/tau-manifesto
tau   = 2*pi;

% Positions of nodes on the circle starting from (0,-1), useful later for label orientation
step  = tau/sz(1);
theta = -.25*tau : step : .75*tau - step;

% Get cartesian x-y coordinates of the nodes
x     = cos(theta);
y     = sin(theta);

% PLOT BEZIER CURVES 
% Calculate Bx and By positions of quadratic Bezier curves with P1 at (0,0)
% B(t) = (1-t)^2*P0 + t^2*P2 where t is a vector of points in [0, 1] and determines, i.e.
% how many points are used for each curve, and P0-P2 is the node pair with (x,y) coordinates.
t2  = [1-t, t].^2;
s.l = NaN(N2,1);

% LOOP per color bucket
for c = 1:N2
    pos = plotorder(c);
    idx = isrt == pos;
    if nnz(idx)
        Bx     = [t2*[x(col(idx)); x(row(idx))]; NaN(1,n(pos))];
        By     = [t2*[y(col(idx)); y(row(idx))]; NaN(1,n(pos))];
        s.l(c) = plot(Bx(:),By(:),'Color',opt.ccolor(pos,:),'LineWidth',linewidth(pos));
    end
end

% PLOT NODES
% Do not rely that r is symmetric and base the mean on lower triangular part only
[row,col]  = find(tf(end:-1:1,end:-1:1) | tf);
subs       = col;
iswap      = row < col;
tmp        = row(iswap);
row(iswap) = col(iswap);
col(iswap) = tmp;

% Plot in brighter color those nodes which on average are more absolutely correlated
[Z,isrt]   = sort(accumarray(subs,abs(data( row + (col-1)*sz(1) )),[],@mean));
Z          = (Z-min(Z)+0.01)/(max(Z)-min(Z)+0.01);
opt.ncolor = hsv2rgb([repmat(opt.ncolor(1:2), sz(1),1) Z*opt.ncolor(3)]);
s.s        = scatter(x(isrt),y(isrt),[], opt.ncolor,'fill','MarkerEdgeColor',ecolor,'LineWidth',1);

% PLACE TEXT LABELS such that you always read 'left to right'
ipos       = x > 0;
s.t        = zeros(sz(1),1);
s.t( ipos) = text(x( ipos)*1.1, y( ipos)*1.1, opt.lbls( ipos),'Color',tcolor);
set(s.t( ipos),{'Rotation'}, num2cell(theta(ipos)'/tau*360))
s.t(~ipos) = text(x(~ipos)*1.1, y(~ipos)*1.1, opt.lbls(~ipos),'Color',tcolor);
set(s.t(~ipos),{'Rotation'}, num2cell(theta(~ipos)'/tau*360 - 180),'Horiz','right')

% ADJUST FIGURE height width to fit text labels
xtn        = cell2mat(get(s.t,'extent'));
post       = cell2mat(get(s.t,'pos'));
sg         = sign(post(:,2));
posfa      = cell2mat(get([gcf gca],'pos'));

% Calculate xlim and ylim in data units as x (y) position + extension along x (y)
ylims      = post(:,2) + xtn(:,4).*sg;
ylims      = [min(ylims), max(ylims)];
xlims      = post(:,1) + xtn(:,3).*sg;
xlims      = [min(xlims), max(xlims)];

% Stretch figure
posfa(1,3) = (( diff(xlims)/2 - 1)*posfa(2,3) + 1) * posfa(1,3);
posfa(1,4) = (( diff(ylims)/2 - 1)*posfa(2,4) + 1) * posfa(1,4);

% Position it a bit lower (movegui slow)
posfa(1,2) = 100;

% Axis settings
set(gca, 'Xlim',xlims,'Ylim',ylims, 'color', 'k', 'layer','bottom', 'Xtick',[],'Ytick',[])
set(gcf, 'pos' ,posfa(1,:),'Visible','on')
axis equal

if nargout == 1
    h = s;
end
end