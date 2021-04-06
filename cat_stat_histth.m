function varargout = cat_stat_histth(src,percent,opt)
% ______________________________________________________________________
% Remove outliers based on the histogram and replace them by the new
% limits. E.g. some MRI images have some extremely high or low values 
% that can trouble other functions that try to work on the full given 
% input range of the data. Removing only 0.2% of the data often often
% helps to avoid problems without removing important information. 
% 
% The function can also print a histogram, box- or violin plot of 
% the given data (using cat_plot_boxplot) and give some basic values.
%
% [res,ths] = cat_stat_histth(src,percent,verb)
%  
%   src     .. input data
%   res     .. limited input data 
%   ths     .. estimated thresholds
%   percent .. included values (default) = 0.998; 
%              can be defined with upper and lower limit, e.g., to be more 
%              aggressive in the background that has more values in general 
%   verb    .. 0 - none 
%              1 - histogram
%              2 - histogram without boundary values
%              3 - box-plot
%              4 - violin-plot
%              5 - violin- and box-plot 
%
% Expert options: 
%   
%  [res,ths] = cat_stat_histth(src,percent, opt )
%
%   opt      .. structur with further fields
%    .verb   .. see above 
%    .fs     .. font size 
%    .hbins  .. bins for histogram estimation 
%    .vacc   .. limit number of elements used in the violin plot
%               e.g. vacc=100 means src(1:100:end)
%    .scale  .. [low high], default [] - no scaling
%               update ths to opt.scale !
%
% Examples: 
%   s=10; b = randn(s,s,s); cat_stat_histth(b,0.9,4);
%   s=10; b = rand(s,s,s);  cat_stat_histth(b,0.9,5);
%   s=10; b = rand(s,s,s);  cat_stat_histth(b,[0.8,0.999],5);
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


  %% check input
  if nargin==0, help cat_stat_histth; return; end
  if ~exist('src','var') || isempty(src) 
    varargout{1} = src;
    varargout{2} = nan(1,2); 
    return; 
  end
  
  if ~exist('opt','var'), opt = struct(); end
  if ~isstruct(opt)
    verb = opt; clear opt; opt.verb = verb; 
  end
  def.verb        = 0;
  def.fs          = 16; 
  def.hbins       = 10000; 
  def.vacc        = max(1,min(100000,round(numel(src)/1000))); % reduce elements in violin plot
  def.scale       = [];
  opt = cat_io_checkinopt(opt,def); 

  tol = [ 0.002 0.002 ]; 
  if nargin==0, help cat_stat_histth; return; end
  if ~exist('percent','var') || isempty(percent)
    tol = [ 0.002 0.002 ]; 
  else
    if numel(percent)==1, percent(2) = percent(1); end
    for pi = 1:2
      if percent(pi)<=1
        tol(pi) = 1 - percent(pi);
      elseif percent<=100
        tol(pi) = 1 - percent(pi)/100; 
      else 
        error('cat_stat_histth:percent','Percent has to be in the range of 1 to 100');  
      end
    end
  end
  
  if ~isreal(src)
    % RD202004: This should not happen but it does in surprise some chimp 
    %           images with strong negative intensities.  At some point in 
    %           cat_vol_sanlm or here the image becomes complex and create 
    %           an error in the histogram fucntion - so I convert it back. 
    src = real(src);
  end

  % histogram
  % use adaptive number of bins to 
  hsrc = zeros(1,opt.hbins); hbins = opt.hbins;
  while (hbins == opt.hbins) || (sum(hsrc>0)/numel(hsrc)<0.3 && hbins<intmax/10 && numel(src)/hbins>2^4)
    [hsrc,hval] = hist(src(~isinf(src(:)) & ~isnan(src(:)) & src(:)<3.4027e+38 & src(:)>-3.4027e+38),hbins);
    hbins = hbins*2; 
  end
  hp = cumsum(hsrc)./sum(hsrc); 

  
  % lower limit
  if opt.verb, srco=src; end
  % Use the last value below tol rather than the first value above tol to 
  % avoid problems 
  ind    = max( [1,find(hp<tol(1),1,'last')]); 
  ths(1) = mean( hval( ind:min(ind+1,numel(hval) ) ) ); 
  src(src<ths(1)) = ths(1); 
  % upper limit
  ind    = min( [numel(hval) , find(hp>(1-tol(2)),1,'first') ] ); 
  ths(2) = mean( hval( max(1,ind-1):ind ) ); 
  src(src>ths(2)) = ths(2); 
  
  if ~isempty(opt.scale) && opt.scale(1)~=opt.scale(2) && diff(ths)~=0
    src = ( src - ths(1) ) ./ diff(ths); 
    src = src * diff(opt.scale) + opt.scale(1); 
    ths = opt.scale; 
  end
  
  %% display
  if opt.verb
    % get figure
    fh = findobj('name','cat_stat_histth');
    if isempty(fh), figure('name','cat_stat_histth'); else, figure(fh); clf; end
    
    % create main plot
    subplot('Position',[0.10 0.06 0.64 0.86]); 
    if opt.verb == 1
      hist(src(:),100); 
    elseif opt.verb == 2
      hist(src(src(:)>ths(1) & src(:)<ths(2)),100); 
    else
      src2 = src(src(:)>ths(1) & src(:)<ths(2));
      if opt.verb == 3
        src2     = src2(1:opt.vacc:end); 
        boxwidth = 0.8;
      else
        src2     = src2(:); 
        boxwidth = 1.0;
      end
      cat_plot_boxplot( { src2 } , struct('violin',opt.verb - 3,'boxwidth',boxwidth)); 
    end
    if opt.verb <= 2
      title('Histogram','fontsize',opt.fs); 
    elseif opt.verb == 3
      title('Boxplot','fontsize',opt.fs); 
    elseif opt.verb == 4
      title('Violinplot without boundary values','fontsize',opt.fs); 
    elseif opt.verb == 5
      title('Violin/Boxplot without boundary values','fontsize',opt.fs); 
    end
    set(gca,'fontsize',opt.fs);
    
    
    % print some values
    annotation('textbox',[.75 0.06 0.23 0.86],'String',sprintf([...
      'max:  %10.4f \nmin:  %10.4f \nmedian:%9.4f \nmean: %10.4f \nstd:  %10.4f \n\n' ...
      'lth:  %10.4f \nhth:  %10.4f \n\n' ...
      'lth_{99}:%10.4f \nhth_{99}:%10.4f \n\n' ...
      'lth_{95}:%10.4f \nhth_{95}:%10.4f \n\n' ...
      ],[max(srco(:)),min(srco(:)),cat_stat_nanmedian(srco(:)),...
      cat_stat_nanmean(srco(:)),cat_stat_nanstd(srco(:))],...
      ths, ...
      [hval(find(hp>0.005,1,'first')),hval(find(hp<0.995,1,'last'))], ...
      [hval(find(hp>0.025,1,'first')),hval(find(hp<0.975,1,'last'))]), ...
      'FitBoxToText','On','fontsize',opt.fs*0.85,'linestyle','none','fontname','fixedwidth'); 
  end
  
  if nargout>=1, varargout{1} = src; end
  if nargout>=2, varargout{2} = ths; end
end