function varargout = cat_stat_histth(src,percent,verb)
% ______________________________________________________________________
% Remove outliers based on the histogram. 
% Print a histogram of your data and give some threshholds
%
% [res,ths] = cat_stat_histth(src,percent,verb)
%  
%   src     .. input data
%   res     .. limited input data 
%   ths     .. estimated thresholds
%   percent .. included values (default) = 0.998; 
%   verb    .. verbose (default = 0); 
% ______________________________________________________________________
% $Id$

  %% check input
  if nargin==0, help cat_stat_histth; return; end
  if ~exist('verb','var'), verb = 0; end
  if ~exist('percent','var') || isempty(percent)
    tol = 0.002; 
  else
    if percent<=1
      tol = 1 - percent; 
    elseif percent<=100
      tol = 1 - percent/100; 
    else 
      error('cat_stat_histth:percent','Percent has to be in the range of 1 to 100');  
    end
  end
  
  
  % histogram
  [hsrc,hval] = hist(src(:),10000);
  hp          = cumsum(hsrc)./sum(hsrc); 

  
  % lower limit
  if verb, srco=src; end
  if min(src(:))~=0
    ths(1) = hval(find(hp>tol,1,'first')); 
    src(src<ths(1)) = ths(1); 
  else
    ths(1) = 0; 
  end
  % upper limit
  ths(2) = hval(find(hp<(1-tol),1,'last')); 
  src(src>ths(2)) = ths(2); 
  
  % display
  if verb
    %%
    fs = 16;
    
    fh = findobj('name','cat_stat_histth');
    if isempty(fh), figure('name','cat_stat_histth'); else figure(fh); clf; end
    
    subplot('Position',[0.06 0.06 0.68 0.86]); 
    hist(src(src(:)>ths(1)),100); title('Histogram','fontsize',fs); 
    set(gca,'fontsize',fs); 
    
    annotation('textbox',[.75 0.06 0.23 0.86],'String',sprintf([...
      'max:  %10.4f \nmin:  %10.4f \nmedian:%9.4f \nmean: %10.4f \nstd:  %10.4f \n\n' ...
      'lth:  %10.4f \nhth:  %10.4f \n\n' ...
      'lth_{99}:%10.4f \nhth_{99}:%10.4f \n\n' ...
      'lth_{95}:%10.4f \nhth_{95}:%10.4f \n\n' ...
      ],[max(srco(:)),min(srco(:)),cat_stat_nanmedian(srco(:)),cat_stat_nanmean(srco(:)),cat_stat_nanstd(srco(:))],...
      ths, ...
      [hval(find(hp>0.005,1,'first')),hval(find(hp<0.995,1,'last'))], ...
      [hval(find(hp>0.025,1,'first')),hval(find(hp<0.975,1,'last'))]), ...
      'FitBoxToText','On','fontsize',fs*0.85,'linestyle','none','fontname','fixedwidth'); 
  end
  
  if nargout>=1, varargout{1} = src; end
  if nargout>=2, varargout{2} = ths; end
end