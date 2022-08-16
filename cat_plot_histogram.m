function out = cat_plot_histogram(data,opt)
% Format out = cat_plot_histogram(data,opt);
% Show histogram of one or more image or surface data 
% If data are spmT-files also print mean, std, effect size,
% and upper 5%-tail cutoff
%
% data           - char array of input files or data matrix
% opt fields:
% color          = cat_io_colormaps('nejm',size(data,1)) as default, see cat_io_colormaps for more categorical colormaps
% norm_frequency = true;
% winsize        = [750 500];
% xrange         = [];
% xlim           = [];
% ylim           = [];
% dist           = 'kernel' (normal, gamma, rician, rayleigh, poisson, weibull)
%                  see fitdist for all distributions
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% default parameter
if ~exist('opt','var'), opt = struct(''); end

if nargin == 0 || (nargin > 0 && isempty(data))
  data = spm_select([1 Inf],{'txt','image','gii','(lh|rh).*'},'Select images to work on');
end

if isempty(data), error('no input images specified'), end

if ischar(data)
  n = size(data,1);
else
  if ~iscell(data)
    [n, ind] = min(size(data));
    if n==1
      data0{1} = data(:);
    else
      if ind == 2
        data = data';
      end
      for i=1:n
        data0{i} = data(i,:);
      end
    end
    data = data0; clear data0;
  end
  n = numel(data);
end

def.color          = cat_io_colormaps('nejm',n);
def.norm_frequency = true;
def.winsize        = [750 500];
def.xrange         = [];
def.xlim           = [];
def.ylim           = [];
def.dist           = 'kernel';

opt = cat_io_checkinopt(opt,def);

if ~isempty(opt.dist) && exist('fitdist') ~= 2
  fprintf('Function fitdist not found: Disable curve fitting.');
  opt.dist = '';
end

% ignore NaNs
dropNaNs = @(x) double(x(~isnan(x)));

mn = zeros(n,1);
mx = zeros(n,1);
cdata = cell(n,1);

for i = 1:n

  if ~ischar(data)
    cdata{i} = single(data{i}(:));
    mn(i) = min(data{i}(:));
    mx(i) = max(data{i}(:));
  else
  
    [pth,nam,ext] = spm_fileparts(deblank(data(i,:)));

    % 1 - txt; 2 - volume; 3 - mesh; 4 - Freesurfer
    if strcmp(ext,'.txt')
      filetype = 1;
    elseif strcmp(ext,'.nii') || strcmp(ext,'.img')
      filetype = 2;
    elseif strcmp(ext,'.gii')
      filetype = 3;
    else
      filetype = 4;
    end

    switch filetype
    case 1
      [cdata{i}, mn(i), mx(i)] = loadsingle_txt(deblank(data(i,:)));
    case 2
      [cdata{i}, mn(i), mx(i)] = loadsingle(nifti(data(i,:)));
    case 3
      [cdata{i}, mn(i), mx(i)] = loadsingle(spm_data_hdr_read(data(i,:)));
    case 4
      try
        [cdata{i}, mn(i), mx(i)] = loadsingleFS(deblank(data(i,:)));
      catch
        error('Unknown data format');
      end
    end
  end
end

if n == 2
  if length(cdata{1}(:)) == length(cdata{2}(:))
    
    fig = figure(11);
    set(fig,'MenuBar', 'none', 'Position',[100,0,500,500]);

    cat_plot_scatter(cdata{1}(:), cdata{2}(:), 'fig', fig);

    ax  = gca;
    set(ax,'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal');

    title('Histogram','Parent',ax);
    if ischar(data(1,:))
      xlabel(spm_str_manip(data(1,:),'a90'),'Parent',ax,'Interpreter','none');
      ylabel(spm_str_manip(data(2,:),'a90'),'Parent',ax,'Interpreter','none');
    else
      xlabel('Data 1','Parent',ax,'Interpreter','none');
      ylabel('Data 2','Parent',ax,'Interpreter','none');
    end
    if ischar(data(1,:)) && (filetype == 2) && (size(cdata{1},2) ~= 1) && (size(cdata{2},2) ~= 1)

      d = double(cdata{2}) - double(cdata{1});
      d(isnan(d)) = 0;
  
      d1 = squeeze(sum(d,1));
      d2 = squeeze(sum(d,2));
      d3 = squeeze(sum(d,3));
      
      mx1 = max(abs(d(:)));
      
      if mx1 > 0.0

        fprintf('Blue: i1>i2; Red: i2>i1\ni1 = %s\ni2 = %s\n',data(1,:),data(2,:));
        fig = figure(12);
        cm = hot(64);
        set(fig,'Menubar','none');
        colormap([1-(cm); cm])

        mx2 = max(abs([d1(:); d2(:); d3(:)]));
        if mx2 == 0, mx2 = eps; end
  
        subplot(2,2,1)
        imagesc(rot90(d1),[-mx2 mx2])
        axis off image
    
        subplot(2,2,2)
        imagesc(rot90(d2),[-mx2 mx2])
        axis off image

        subplot(2,2,3)
        imagesc(d3,[-mx2 mx2])
        axis off image

        subplot(2,2,4)
        colorbar
        set(gca,'CLim',[-mx1 mx1]);
        axis off image
      else
        disp('Images are identical!');
      end
    % display surface rendering of difference
    elseif ischar(data(1,:)) && filetype == 3
      d = double(cdata{2}) - double(cdata{1});
      d(isnan(d)) = 0;
      sinfo = cat_surf_info(data(1,:));
      Pmesh = '';
      if ~isempty(strfind(sinfo.Pmesh,'templates_surfaces'))
        Pmesh = sinfo.Pmesh;
      elseif sinfo.resampled
        if sinfo.resampled_32k
          templates_surfaces = 'templates_surfaces_32k';
        else
          templates_surfaces = 'templates_surfaces';
        end
        Pmesh = fullfile(spm('dir'),'toolbox','cat12',templates_surfaces,[sinfo.side '.central.freesurfer.gii']);
      end
      if exist(Pmesh,'file')
        % scale in 2..98% range
        range = cat_vol_iscaling(d,[0.02 0.98]);
        absmx = max(abs(range));
        S = gifti(Pmesh);
        cat_surf_render2(struct('vertices',S.vertices,'faces',S.faces,'cdata',d));
        cat_surf_render2('colorbar'); 
        cat_surf_render2('view','top');
        cat_surf_render2('clim',[-absmx absmx]); 
        fprintf('Blue: i1>i2; Red: i2>i1\ni1 = %s\ni2 = %s\n',data(1,:),data(2,:));
      end
    end
  else
    disp('No 2D histogram plotted because size differs between data.');
  end   
end

if isempty(opt.xrange)
  X0 = linspace(min(mn), max(mx), max(min(round(numel(dropNaNs(cdata{1}))/100),500),10));
elseif numel(opt.xrange) == 2
  X0 = linspace(opt.xrange(1), opt.xrange(2), 500);
else
  error('Parameter xrange does not consist of two entries');
end

fig = figure(13);
set(fig,'MenuBar', 'none', 'Position',[100, 0, opt.winsize]);

for j = 1:n
  y = dropNaNs(cdata{j});
  if ~isempty(opt.dist)
    n = numel(y);
    H0 = hist(y,X0);
    pd = fitdist(y,opt.dist);

    [bincounts,binedges] = histcounts(y,X0);

    % Normalize the density to match the total area of the histogram
    Hfit(j,:) = n * (binedges(2)-binedges(1)) * pdf(pd,X0);
  else
    H0 = hist(y,X0);
    Hfit(j,:) = H0;
  end

  if opt.norm_frequency
    if ~isempty(opt.dist)
      Hfit(j,:) = Hfit(j,:)/sum(H0);
    end
    H0 = H0/sum(H0);
  end
  H(j,:) = H0;
  X(j,:) = X0;
  if ischar(data)
    legend_str{j} = char(spm_str_manip(data(j,:),'a90'));
  
    % give some specific output for (normally distributed) T-values or
    % effect size (D)
    [pth,nam] = spm_fileparts(deblank(data(j,:)));
    spmT_found = ~isempty(strfind(nam,'spmT')) || strcmp(nam(1),'D');
    if spmT_found
      mn = mean(y);
      sd = std(y);
      ES = mn/sd;
      TH5 = X0(min(find(cumsum(H0)/sum(H0) > 0.95)));
      fprintf('%s\tmean=%g\tSD=%g\tES=%g\tTH5=%g\n',data(j,:),mn,sd,ES,TH5);
      legend_str{j} = sprintf('TH5=%.4f %s',TH5,legend_str{j}); 
    else
      fprintf('%s\tSD=%g\n',data(j,:),std(y));
    end
  else
    legend_str{j} = num2str(j);
  end
end

% check whether there are (almost) identical data
Hcorr = corrcoef(H');
Hcorr = Hcorr - triu(Hcorr);
[xc,yc] = find(Hcorr > 0.99999);
for i = 1:numel(xc)
%  fprintf('File %s and %s are (almost) identical\n',deblank(data(xc(i),:)),deblank(data(yc(i),:)));
end

if ~isempty(opt.dist)
  HP  = plot(X(:,2:end-1)', Hfit(:,2:end-1)');
  hold on
  HP0 = plot(X(:,2:end-1)', H(:,2:end-1)');
  hold off
else
  HP = plot(X(:,2:end-1)', H(:,2:end-1)');
end

for i = 1:length(HP)
  set(HP(i),'LineWidth',1);
  if ~isempty(opt.dist)
    set(HP0(i),'LineWidth',1,'Linestyle',':');
  end
  if ~isempty(opt.color)
    set(HP(i),'Color',opt.color(i,:));
    if ~isempty(opt.dist)
      set(HP0(i),'Color',opt.color(i,:));
    end
  end
end

h = legend(legend_str);
set(h,'Interpreter','none');
grid on
if opt.norm_frequency
  ylabel('Normalized Frequency');
else
  ylabel('Frequency');
end

if ~isempty(opt.xlim) && numel(opt.xlim) == 2
  xlim(opt.xlim)
end

if ~isempty(opt.ylim) && numel(opt.ylim) == 2
  xlim(opt.ylim)
end

if nargout
  out = HP;
end

%_______________________________________________________________________
function [udat, mn, mx] = loadsingle(V)
% Load surface or volume data from file indicated by V into an array of floats.

% use fast method for file reading for nifti files
if isa(V,'nifti')
  udat(:,:,:) = V.dat(:,:,:);
else
  udat = spm_data_read(V);
end

% remove zero background
ind0 = find(udat == 0);
if ~isempty(ind0)
  if length(ind0) > 0.01*numel(udat)
    udat(ind0) = NaN;
  end
end

mx = max(udat(:));
mn = min(udat(:));

udat = single(udat);

%_______________________________________________________________________
function [udat, mn, mx] = loadsingle_txt(P)
% Load txt data from file indicated by V into an array of floats.

udat = spm_load(P);

% remove zero background
ind0 = find(udat == 0);
if ~isempty(ind0)
  if length(ind0) > 0.01*numel(udat)
    udat(ind0) = NaN;
  end
end

mx = max(udat(:));
mn = min(udat(:));

udat = single(udat);

%_______________________________________________________________________
function [udat, mn, mx] = loadsingleFS(P)
%
% [udat, mn, mx] = loadsingleFS(P)
% reads a binary curvature file into a vector with single data type
%

udat = cat_io_FreeSurfer('read_surf',filename);

mx = max(udat(:));
mn = min(udat(:));

udat = single(udat);

