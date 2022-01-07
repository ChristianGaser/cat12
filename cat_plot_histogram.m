function out = cat_plot_histogram(data,opt)
%CG_HIST2D
% Format out = cat_plot_histogram(data,opt);
% Show histogram of one or more image or surface data 
% If data are spmT-files also print mean, std, effect size,
% and upper 5%-tail cutoff
%
% data           - input files
% opt fields:
% color          = nejm(size(data,1)) as default, see cat_io_colormaps for more categorical colormaps
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
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
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
    mx(i) = min(data{i}(:));
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
    
    % Create the joint histogram
    [edgesX,edgesY, H2] =   ndhist([cdata{2}(:), cdata{1}(:)]);

    H2 = log(H2+1);

    figure
    h = pcolor(edgesX,edgesY,H2);
    set(h,'EdgeColor','none')

    set(gcf,'MenuBar', 'none', 'Position',[100,0,500,500]);
    ax  = gca;
    set(ax,'PlotBoxAspectRatioMode','auto','XDir','normal','YDir','normal');

    title('Histogram','Parent',ax);
    xlabel(spm_str_manip(data(2,:),'a90'),'Parent',ax,'Interpreter','none');
    ylabel(spm_str_manip(data(1,:),'a90'),'Parent',ax,'Interpreter','none');
        
    if (filetype == 2) && (size(cdata{1},2) ~= 1) && (size(cdata{2},2) ~= 1)

      d = double(cdata{2}) - double(cdata{1});
      d(isnan(d)) = 0;
  
      d1 = squeeze(sum(d,1));
      d2 = squeeze(sum(d,2));
      d3 = squeeze(sum(d,3));
      
      mx1 = max(abs(d(:)));
      
      if mx1 > 0.0

        fprintf('Blue: i1>i2; Red: i2>i1\ni1 = %s\ni2 = %s\n',data(1,:),data(2,:));
        h = figure;
        cm = hot(64);
        set(h,'Menubar','none');
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
    elseif filetype == 3
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

fig = figure('MenuBar', 'none', 'Position',[100, 0, opt.winsize]);

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
  
    % give some specific output for (normally distributed) T-values
    [pth,nam] = spm_fileparts(deblank(data(j,:)));
    spmT_found = ~isempty(strfind(nam,'spmT'));
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
    legend_str{j} = num2str(i);
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

%_______________________________________________________________________
% Add logx,logy,loglog
% 

% Displays a 2d histogram for your data, it will set the bins appropriately
% 
% NDHIST(x,y); where x and y are equal length vectors. It will choose
%              reasonable axis bounds and bins to present your data. The
%              default parameters may leave some data off the chart.
% 
% NDHIST(XY);  where XY = [x y] and x,y are vectors
% 
% NDHIST(z);   where z is a vector of complex numbers, (x+1i*y) or amplitude*exp(1i*theta)
% 
% NDHIST(y);   where y is a vector of real numbers will plot a 2d histogram
%              of the points, as if it were a line-chart. It is equivalent
%              to calling ndhist(1:length(y),y);
%  
% N = NDHIST(x,y); returns a matrix N containing the counts for each bin
%              determined by the histogram.
% 
% [edgesX2,edgesY2,N,h] = NDHIST(x,y); which returns a matrix N containing
%              the counts for each bin determined by the histogram. You can
%              plot it with sanePColor(edgesX2,edgesY2,N); (from Matlabcentral)
%              h is the plot handle.
% 
% NDHIST(...,'param','value','param','value', ... ); Run ndhist with specific
%              parameters
%           
% List of special parameters: 
% 
%  'filter' : This will apply a gaussian filter to the final histogram data.
%             The default filter width is 5 bins wide. If you pass a number
%             then that will be used. Even numbered filter parameters will be
%             changed to odd numbers to keep the filter perfectly symetrical.
%             'filt','filtering','smooth'
% 
%    'log' : Change the colormap to be on a log scale to represent data
%            over a large dynamic range. 
%            'logplot'
% 
%   'bins' : Change the size of the bins. For example '2' will create a
%            plot with twice the default number of bins; 0.5 will have half
%            the default number of bins. The default uses Scott's normal
%            reference rule. Unclear if it is ideal for 2d histograms...
%            If you are looking for a histogram with specific bins, use the
%            subfunction hist3. Feel free to implement it as an additional
%            parameter 'edgdes','edgesx' or 'edgesy'
%            'f','numbins'
% 
%  'binsx' : Change the size of only the x bins. 'fx'
%  'binsy' : Change the size of only the y bins. 'fy'
%            
%     axis : This is to set the range of the plot, [xmin xmax ymin ymax]
%            The default range is set to 3*std(x) and 3*std(y) where the
%            parameter stdTimes = 3 is hard-coded in this version and
%            potentially added as a parameter in a later version.
% 
%      max : This is to set the range of the plot to be such that every
%            point will be contained within the plot.
%            'themax'
%  
%  intbins : Set the bins to be intiger widths. For both x and y
%            'int'
% 
% intbinsx : Set the x bins to be intiger widths. 'intx'
% intbinsy : Set the y bins to be intiger widths. 'inty'
% 
%normalizex: Normalize the plot so that the sum of all the y values in each
%            x bin sum to one.
%            'normx','nx' 
% 
%normalizey: Normalize the plot so that the sum of all the x values in each
%            y bin sum to one.
%            'normy','ny'
% 
%normalizeR: Normalize the plot so that the you can clearly see how the
%            distribution vary's over angle. It weights points in the outer
%            radius by the diameter at that radius.
%            'nr'
%    points: Plot the points on top of the colored histogram.
% 
%     
% 
% NOT IMPLEMENTED YET
%'samebins': NOT IMPLEMENTED YET. Would set the width of the x and y bins
%            to be equal to each other and the axis equal too.

% 
% user parameters:
% filter: This will filter the data, you may choose to follow it with a
%         number. This number will represent the radius of the circular
%         gaussian filter. Other ways to call it: 'filt','filtering','f'
%  
% 
% examples
% 
% To test the function you may use this example:
% z = 2*randn(1,100000)+1i*(randn(1,100000));
% 
% If you have amplitude and angle measures then pass this:
% z = amp*exp(1i*ang);
% 
% NDHIST(z)
% NDHIST(z,'lansey')
% NDHIST(z,'filter')
% 
% % Note
% The name of this function comes because really its a 2d hist, but since I
% already have an 'nhist' I thought I might name it this.
% 
% SEE ALSO: HIST, HIST3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Lansey May 2009,     questions to Lansey at gmail.com          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan C. Lansey (2021). Efficient 2D histogram, no toolboxes needed 
% https://www.mathworks.com/matlabcentral/fileexchange/45325-efficient-2d-histogram-no-toolboxes-needed
% MATLAB Central File Exchange. Retrieved December 13, 2021.
%%
%_______________________________________________________________________
function [edgesX2,edgesY2,N,h] = ndhist(z,varargin)

%% Errors and warning for bad data
if ~isnumeric(z)
  error('you must take the histogram of numeric data silly');
end
if numel(z)<10
  warning(['you have only ' num2str(numel(z)) ' being plotted, are you sure a histogram is such a good idea?']);
end
if isempty(z)
  warning('no data, returning');
  return;
end

% Also re-call the function in case of things
if isreal(z)
  if size(z,2) == 2 % we have a [x y] case
    [edgesX2,edgesY2,N] = ndhist(z(:,1)+1i*z(:,2),varargin{:});
  else % we've got either 
    if isempty(varargin) % just one value passed
      disp('your data are all ''real'' so we intepreted it as a timeseries');
      idx = (1:length(z))';
      [edgesX2,edgesY2,N] = ndhist(idx+1i*z(:),'fx',4,'fy',2);
      colormap(linspecer('blue'));
    else% at least two values passed
      y = varargin{1};
      if isnumeric(y) % we've got the (x,y) case
        %         if do something about nargout here
        if length(y)~= length(z)
          error('x and y must be the same length');
        end
        [edgesX2,edgesY2,N] = ndhist(z(:)+1i*y(:),varargin{2:end});
      else % we've got just one value, but with special arguments passed
        idx = (1:length(z))';
        [edgesX2,edgesY2,N] = ndhist(idx+1i*z(:),'fx',4,'fy',1,varargin{:});
      end
    end
  end
  if nargout == 1
    edgesX2 = N;
  end
  if nargout == 2
    warning('you are being passed out ''edgesX2,edgesY2'' is that really what you want from this function?');
  end
  
  return;
end
% great we can continue


%% Standardize the data
% Pcolor does not work with 'single' type so this makes it is a double.
% Making it into a vertical vector is also important (see below).
z = double(z(:));

% remove nan data
I = isnan(z);
if sum(I)>0
  z = z(~I);
end

% separate x and y
x = real(z); y = imag(z);

%% Set program defaults
filtering = 0;
userParam = 0;
stdTimes = 3; % the number of times the standard deviation to set the upper end of the axis to be.
minBins = 10;
maxBins = 1000;
binFactorX = 1;
binFactorY = 1;
intbinsFlagX = isdiscrete(x);
intbinsFlagY = isdiscrete(y);

normalizeFlagX = 0;
normalizeFlagY = 0;

pointsFlag = 0; % to plot the points also

maxisFlag = 0; % If you choose this then no data will fall off the screen
logColorFlag = 0; % plot the logarithm of the counts instead to get more contrast

S = length(x);

stdX = std(x); meanX = mean(x);
stdY = std(y); meanY = mean(y);

if stdX>(10*eps) % just checking there are more than two different points to the data, checking for rounding errors.
% we include some padding in here
  leftEdge = max(meanX-stdX*stdTimes,min(x));
  riteEdge = min(meanX+stdX*stdTimes,max(x));
else % stdV == 0, wow, all your data points are equal
  leftEdge = min(x)-1000*eps; % padd it by 100, seems reasonable
  riteEdge = max(x)+1000*eps;
end

if stdY>(10*eps) % just checking there are more than two different points to the data, checking for rounding errors.
  botEdge = max(meanY-stdY*stdTimes,min(y));
  topEdge = min(meanY+stdY*stdTimes,max(y));
else % stdV == 0, wow, all your data points are equal
  botEdge = min(y)-1000*eps; % padd it by 100, seems reasonable
  topEdge = max(y)+1000*eps;
end

padX = (riteEdge - leftEdge)*.01;
padY = (topEdge - botEdge)*.01;

axisXY = [leftEdge riteEdge botEdge topEdge]+[-padX padX -padY padY]; % do we need this much padding really?
rangeX = riteEdge-leftEdge;
rangeY =  topEdge -botEdge;


%% interperet user parameters
k = 1;
while k <= nargin-1
  if ischar(varargin{k})
  switch lower(varargin{k})
    case {'filter','filt','filtering','smooth'}
      filtering = 5; % default filter radius
      if k<= nargin-2% -1+1, lol
        if ~ischar(varargin{k+1})
          filtering = varargin{k+1};
          k = k + 1;
        end
      end
    case {'axis'}
      axisXY = varargin{k+1};
      k = k + 1;
    case {'themax','max'} % make it so all of your data are plotted
      maxisFlag = 1;

    case {'range'} % Note: this feature does nothing yet.
      rangeR = varargin{k+1};
      k = k + 1;
    case {'bins','numbins','f'} % 'f' comes from the binfactor of nhist
      binFactorY = varargin{k+1};
      binFactorX = varargin{k+1};
      k = k + 1;       
    case {'binsx','fx'} % 'f' comes from the binfactor of nhist
      binFactorX = varargin{k+1};
      k = k + 1;       
    case {'binsy','fy'} % 'f' comes from the binfactor of nhist
      binFactorY = varargin{k+1};
      k = k + 1;       
    case {'intbins','int'} % 'f' comes from the binfactor of nhist
      intbinsFlagX = 1;
      intbinsFlagY = 1;
      if k+1<= length(varargin)
        temp = varargin{k + 1};
        if ~ischar(temp) % if its a number then we want to use it.
          intbinsFlagX = temp;
          intbinsFlagY = temp;
          k = k+1;
        end
      end
    case {'intbinsx','intx'} % 'f' comes from the binfactor of nhist
      intbinsFlagX = 1;
      if k+1<= length(varargin)
        temp = varargin{k + 1};
        if ~ischar(temp) % if its a number then we want to use it.
          intbinsFlagX = temp;
          k = k+1;
        end
      end
    case {'intbinsy','inty'} % 'f' comes from the binfactor of nhist
      intbinsFlagY = 1;
      if k+1<= length(varargin)
        temp = varargin{k + 1};
        if ~ischar(temp) % if its a number then we want to use it.
          intbinsFlagY = temp;
          k = k+1;
        end
      end
    case {'lansey','normalizer','nr'} % this guy weights the values based on radius
      userParam = 1;     
    case {'log','logplot'}
      logColorFlag = 1;
    case {'normalizex','normx','nx'}
      normalizeFlagX = 1;
    case {'normalizey','normy','ny'}
      normalizeFlagY = 1;
    case {'points','.'}
      pointsFlag = 1;
    case 'stdtimes' % the number of times the standard deviation to set the upper end of the axis to be.
      stdTimes = varargin{k + 1};
      k = k + 1;
      if ischar(stdTimes)
        fprintf(['\nstdTimes set to: ' stdTimes]);
        error('stdTimes must be a number')
      end

    otherwise 
      warning(['you have passed a strange parameter: ''' varargin{k} ''' please roll again']);
  end
  else
    warning(['input parameter ''' num2str(varargin{k}) ''' not understood, please use ''param'',''value'' pairs']);
  end
  k = k + 1; % increment it so it goes to the next one
end
%%

if normalizeFlagX && normalizeFlagY
  warning('Only normalize X was used');
  normalizeFlagY = 0;
end
  

%% set the bin widths
% Using Scott's normal reference rule, unclear if it is ideal for 2D histograms ...
binWidthX = 3.5*stdX/(binFactorX*S^(1/3));
binWidthY = 3.5*stdY/(binFactorY*S^(1/3));

% Instate a mininum and maximum number of bins
numBinsX = rangeX/binWidthX; % Approx number of bins
numBinsY = rangeY/binWidthY; % Approx number of bins

if numBinsX<minBins % if this will imply less than 10 bins
  binWidthX = rangeX/(minBins); % set so there are ten bins
end
if numBinsX>maxBins % if there would be more than the max bins
  binWidthX = rangeX/maxBins;
end

if numBinsY<minBins % if this will imply less than 10 bins
  binWidthY = rangeY/(minBins); % set so there are ten bins
end
if numBinsY>maxBins % if there would be more than the max bins
  binWidthY = rangeY/maxBins;
end

% check for maxis
if maxisFlag
  temp = [max(x)-min(x), max(y)-min(y)]*[-1 1 0 0 ; 0 0 -1 1]*.05;
  axisXY = [min(x) max(x) min(y) max(y)]+temp; % do some padding with one matrix equation woo!
end

% round the edges if intbins are a thing
if intbinsFlagX
  axisXY(1) = round(axisXY(1))-.5; % subtract 1/2 to make the bin peaks appear on the numbers.
  axisXY(2) = round(axisXY(2))+.5;
  binWidthX = max(round(binWidthX),1);
    
end
if intbinsFlagY
  axisXY(3) = round(axisXY(3))-.5; % subtract 1/2 to make the bin peaks appear on the numbers.
  axisXY(4) = round(axisXY(4))+.5;
  binWidthY = max(round(binWidthY),1);

end

% finally set the bins
edgesX = axisXY(1):binWidthX:axisXY(2);
edgesY = axisXY(3):binWidthY:axisXY(4);

%% start the real histogram computation

N =  hist3(x,y,edgesX,edgesY);

[X,Y] = meshgrid(edgesX,edgesY);

%% user parameters to adjust the results
if filtering
  N2 = smooth2(N,filtering);
  N = N2;
  
end

if 0 % IMPLEMENTATION NOT TESTED if you wanted real probabilities instead of counts then you would use this one.
%   what exactly is 'x' here?
  N = N/(binWidthX*binWidthY*length(x)); % normalize by area to get probability distribution
end

if userParam % to do the radial normalizing
  R2 = X.^2+Y.^2; % R.^2
  R = sqrt(R2); % R.^2
  N = N.*R;
end

if logColorFlag
% Note: we are safe from negative numbers because nothing can be less than zero here
  N = log10(N+1);
end

if normalizeFlagX
%   N = N./repmat(sum(N,1)+eps,size(N,1),1);
  N = N./repmat(max(N,[],1)+eps,size(N,1),1);
end
if normalizeFlagY
%   N = N./repmat(sum(N,2)+eps,1,size(N,2));
  N = N./repmat(max(N,[],2)+eps,1,size(N,2));
end

%%
edgesX2 = edgesX+binWidthX/2;
edgesY2 = edgesY+binWidthY/2;

if nargout == 1
  edgesX2 = N;
end
if nargout == 2
  warning('you are being passed out ''edgesX2,edgesY2'' is that really what you want from this function?');
end


%% tweak the colormap to be lighter
% colormap((2*jet(64)/3+1/3));
% colormap(linspecer(128))
% C = colormap;

%% plot the points if we need to
if pointsFlag
  hold on;
  plot(real(z),imag(z),'k.');
end


%%
% 2D histogram which is actually kind of fast making use of matlab's histc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Lansey May 2009,     questions to Lansey at gmail.com          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%_______________________________________________________________________
function allN = hist3(x,y,edgesX,edgesY)

allN = zeros(length(edgesY),length(edgesX));
[~,binX] = histc(x,edgesX);
for ii = 1:length(edgesX)
  I = (binX == ii);
  N = histc(y(I),edgesY);
  allN(:,ii) = N';
end


%%
% This will tell if the data is an integer or not.
% first it will check if matlab says they are integers, but even so, they
% might be integers stored as doubles!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Lansey May 2009,     questions to Lansey at gmail.com          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%_______________________________________________________________________
function L = isdiscrete(x,varargin) % L stands for logical

minError = eps*100; % the minimum average difference from integers they can be.
L = 0; % double until proven integer
if ~isempty(varargin)
  minError = varargin{1};
end 
if isinteger(x)||islogical(x) % your done, duh, its an int.
  L = 1;
  return; 
else
  if sum(abs(x-round(x)))/length(x)<minError
    L = 1;
  end
end


% This is a fast smooth function that will return a smoothed
% version of the original that you pass it.
% 
% Hey now the default is actually going to be a gaussian filter not a
% moving average filter
% 
% It will be the same length as the original, and it will be centered on
% the original function. To padd the edges it extends the average of the
% last 'n' values on the end out further.
% 
% Note that it uses a convolution which uses the
% cool fft trick to do it effeciently.
%% set things up if you want to test it as a script 
% N = 301;
% % y = rand(N,1);
% y = -1./[1:N];
% y = sin(linspace(1,5*pi,N));
% n = 10;

% This guys is going to smooth it for 2 dimensions I hope;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Lansey May 2009,     questions to Lansey at gmail.com          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%_______________________________________________________________________
function yout = smooth2(A,n)

if length(A)<3
  warning('Sorry bud, you can''t smooth less than 3 points, thats silly');
  yout = A;
  return;
end

if length(n) == 1
  if n<2
    yout = A;
    return;
  end
  
%   forcing you to have an odd 'n'
  if double(~logical(round(n/2)-n/2))
    n = n+1;
  end
  
  bee = linspace(-1.96,1.96,n); % normal distribution with 95% confidence bounds
  [BX, BY] = meshgrid(bee);
  R2 = BX.^2+BY.^2;
  toConvolve = exp(-R2)/sum(exp(-R2(:)));
  
else
  toConvolve = n;
  if round(sum(toConvolve)*100000)/100000 ~= 1
    warning('the sum here does not equal to one.');
  end
  n = length(toConvolve);
end

if min(size(A))<= n
  warning('Sorry bud, you can''t smooth that, pick a smaller n or use more points');
  yout = A;
  return;
end


% padding on the left
padLeft = repmat(mean(A(:,1:n),2),1,n);
padRite = repmat(mean(A(:,end-n:end),2),1,n);
A = [padLeft A padRite];

padTop = repmat(mean(A(1:n,:),1),n,1);
padBot = repmat(mean(A(end-n:end,:),1),n,1);

A = [padTop; A; padBot];

% the main event
As = conv2(A,toConvolve);

% outputting a centered subset
isEven = double(~logical(round(n/2)-n/2));
if isEven
  yout = As(n+(n/2):end-n-(n/2),n+(n/2):end-n-(n/2));
else % it is odd then
  yout = As(n+(n/2)+.5:end-n-(n/2)+.5,n+(n/2)+.5:end-n-(n/2)+.5);
end

