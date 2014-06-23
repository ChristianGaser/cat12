function [out,s] = vbm_plot_boxplot(data,opt)
% _________________________________________________________________________
%
% usage: vargout = vbm_plot_boxplot(data,opt);
%
%  opt.notched     = 0;             % thinner at median [0 1] with 1=0.5
%  opt.symbol      = '+o';          % outlier symbols
%  opt.vertical    = 1;             % boxplot orientation 
%  opt.maxwhisker  = 1.5;           % 
%  opt.sort        = 0;             % no sorting
%                  = 1;             % sort groups (ascending)
%                  = 2;             % sort groups (descending)[inactive]
%                  = [index];       % or by a index matrix
%  opt.names       = {};            & cell of group names
%  opt.fill        = 1;             % filling of boxes
%  opt.groupnum    = 1;             % add number of elements
% [opt.groupmin    = 5;]            % minimum number of non-nan-elements
%                                     in a group [inactive]
%  opt.ylim        = [-inf inf];    % y-axis scaling
%  opt.ygrid       = 1;             % activate y-grid-lines
%  opt.groupcolor  = [R G B];       % matrix with (group)-bar-color(s) 
%                                     use jet(numel(data)) 
%                                     or other color functions
%  opt.fontsize    = [];            % achsis fontsize 
%                                     important for ygird size!
%
% The box plot is a graphical display that simultaneously describes several 
% important features of a data set, such as center, spread, departure from 
% symmetry, and identification of observations that lie unusually far from
% the bulk of the data.
%
% data is a matrix with one column for each dataset, or data is a cell
% vector with one cell for each dataset.
% opt.notched = 1 produces a notched-box plot. Notches represent a robust 
% estimate of the uncertainty about the median.
% opt.notched = 0 (default) produces a rectangular box plot. 
% opt.notched in (0,1) produces a notch of the specified depth.
% opt.notched values outside [0,1] are amusing if not exactly practical.
% opt.notched sets the notched for the outlier values, default notched for
% points that lie outside 3 times the interquartile range is 'o',
% default opt.notched for points between 1.5 and 3 times the interquartile
% range is '+'. 
%
% Examples
% opt.notched = '.' points between 1.5 and 3 times the IQR is marked with
% '.' and points outside 3 times IQR with 'o'.
% opt.notched = ['x','*'] points between 1.5 and 3 times the IQR is marked with
% 'x' and points outside 3 times IQR with '*'.
% opt.vertical = 0 makes the boxes horizontal, by default opt.vertical = 1.
% maxwhisker defines the length of the whiskers as a function of the IQR
% (default = 1.5). If maxwhisker = 0 then boxplot displays all data  
% values outside the box using the plotting opt.notched for points that lie
% outside 3 times the IQR.   
%
% The returned matrix s has one column for each dataset as follows:
%
%    1  minimum
%    2  1st quartile
%    3  2nd quartile (median)
%    4  3rd quartile
%    5  maximum
%    6  lower confidence limit for median
%    7  upper confidence limit for median
%
% Example
%
%   title("Grade 3 heights");
%   tics("x",1:2,["girls";"boys"]);
%   axis([0,3]);
%   boxplot({randn(10,1)*5+140, randn(13,1)*8+135});
%
% _________________________________________________________________________
%
% Author: Alberto Terruzzi <t-albert@libero.it>
% Version: 1.4
% Created: 6 January 2002
% Copyright (C) 2002 Alberto Terruzzi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% modified by Christian Gaser (christian.gaser@uni-jena.de) and
% Robert Dahnke (robert.dahnke@uni-jena.de)
% original version was written for octave by Alberto Terruzzi
% _________________________________________________________________________
% $Id$


  % default parameter
  if ~exist('opt','var'), opt = struct(''); end
  def.notched     = 0;
  def.symbol      = '+o';
  def.vertical    = 1;
  def.maxwhisker  = 1.5;
  def.sort        = 0; 
  def.names       = num2str( (1:numel(data))' );
  def.fill        = 1;
  def.groupcolor  = repmat([0.6 0.6 0.8],numel(data),1);
  def.groupnum    = 1;
  def.groupmin    = 5;
  def.ylim        = [-inf inf];
  def.ygrid       = 1;  
  def.fontsize    = []; % empty = default font size
  

  opt = checkinopt(opt,def);
  opt.notched = max(0,min(1,opt.notched));
  
  % figure out how many data sets we have
  if iscell(data), 
    nc = length(data);
  else
    if isvector(data), data = data(:); end
    nc = columns(data);
  end
  opt.names = cellstr(opt.names);
  if numel(opt.names) < nc
    error('ERROR:vbm_plot_boxplot:names','ERROR: To short name list.'); 
  end
  
  % update colortable
  if size(opt.groupcolor,1)==1
    if size(opt.groupcolor,1)<nc 
      warning('WARNING:vbm_plot_boxplot:groupcolor','WARNING: To short colortable.'); 
    end
    opt.groupcolor = repmat(opt.groupcolor(1,:),numel(data),1);
  end
  if numel(opt.sort)>1 && numel(opt.sort) ~= nc
    error('ERROR:vbm_plot_boxplot:sort','ERROR: To sort list.'); 
  end
  
  groupnr = cellfun(@(x) sum(~isnan(x)),data);
  out.sortj = 1:length(data);
  rmdata = zeros(1,nc);
  % remove groups with to few elemnts
  % ... require addaption for many other fields like names, color, ...
  if 0 && opt.groupmin>0
    rmdata=cellfun('isempty',data) | groupnr<opt.groupmin;
    if numel(opt.sort)==numel(data)
      opt.sort(rmdata) = [];
    end
    opt.names(rmdata) = [];
    opt.groupcolor(rmdata) = [];
    data(rmdata) = [];
  end
  if isempty(data), 
    error('ERROR:vbm_plot_boxplot:data','ERROR: Not enought (non-NaN) data (may change opt.groupmin).'); 
  end
    
  
  % add number of group elements
  if opt.groupnum
    for ni=1:numel(opt.names)
      opt.names{ni}=sprintf('%s[%d]',opt.names{ni},groupnr(ni));
    end
  end
  
  
  % sort groups by their median value or a specific order
  if isfield(opt,'sort') 
    if numel(opt.sort)==1 && opt.sort
    % sort by median
      mdata = zeros(1,numel(data));
      for i=1:numel(data), mdata(i) = vbm_stat_nanmedian(data{i}(:)); end
      [mdata,sorti] = sort(mdata);
      clear mdata;
      data          = data(sorti);
      opt.names     = opt.names(sorti); 
      opt.groupcolor = opt.groupcolor(sorti,:);
      tmp = out.sortj(~rmdata);
      [tmp,out.sortj(~rmdata)] = sort(tmp(sorti));
    elseif ~opt.sort
      % noting to do, just to avoid an error
    elseif numel(opt.sort)==numel(data)
    % sort by given order
      sorti = opt.sort;
      data          = data(sorti);
      opt.names     = opt.names(sorti); 
      opt.groupcolor = opt.groupcolor(sorti,:);
      tmp = out.sortj(~rmdata);
      [tmp,out.sortj(~rmdata)] = sort(tmp(sorti));
    end
  end  
  out.sorti = 1:nc;
  out.sortj = 1:nc;

  if length(opt.symbol)==1, opt.symbol(2)=opt.symbol(1); end

  if opt.notched==1, opt.notched=0.5; end
  a = 1-opt.notched;

  
  
  %% anova, if statistical toolbox exist
  %{
  if exist('kruskalwallis','file') 
    try
      gr = []; for gi=1:nc, gr=[gr gi*ones(1,numel(data{gi}))]; end %#ok<AGROW>
      [out.kw.p,out.kw.tn,out.kw.st] = kruskalwallis(cell2mat(data),gr,'off');
      out.kw.cn05   = multcompare(out.kw.st,'display','off','alpha',0.05);
      out.kw.cn01   = multcompare(out.kw.st,'display','off','alpha',0.01);
      out.kw.cn001  = multcompare(out.kw.st,'display','off','alpha',0.001);
      out.kw.cn0001 = multcompare(out.kw.st,'display','off','alpha',0.0001);
    %{
      out.anova.F=ones(nc,nc); out.tt2p=ones(nc,nc);
      for i=1:nc
        for j=1:nc
          out.rs(i,j) = ranksum(data{i},data{j});
          out.kw.p(i,j) = kruskalwallis(cell2mat(data([i j])),[zeros(1,numel(data{i})),ones(1,numel(data{j}))],'off');
        end
      end
      %}

      out.kw.pg = ones(nc);
      FN={'cn05','cn01','cn001','cn0001'};val=[0.05,0.01,0.001,0.0001];
      for fi=1:numel(FN)
        for i=1:size(out.kw.(FN{fi}),1)
          if out.kw.(FN{fi})(i,3)>0
            out.kw.pg(out.kw.(FN{fi})(i,1),out.kw.(FN{fi})(i,2)) = val(fi);
            out.kw.pg(out.kw.(FN{fi})(i,2),out.kw.(FN{fi})(i,1)) = val(fi);
          end
        end
      end
    end
  end
  %}
  
  

  %% compute statistics
  % s will contain
  %    1,5    min and max
  %    2,3,4  1st, 2nd and 3rd quartile
  %    6,7    lower and upper confidence intervals for median
  s = zeros(7,nc);
  box = zeros(1,nc);
  whisker_x   = ones(2,1)*[1:nc,1:nc];
  whisker_y   = zeros(2,2*nc);
  outliers_x  = [];
  outliers_y  = [];
  outliers2_x = [];
  outliers2_y = [];

  for i=1:nc
    % Get the next data set from the array or cell array
    if iscell(data)
      col = data{i}(:);
    else
      col = data(:,i);
    end
    % Skip missing data
    col(isnan(col)) = [];
    % Remember the data length
    nd = length(col);
    box(i) = nd;
    if (nd > 1)
      % min,max and quartiles
      s(1:5,i) = [min(col) prctile(col,[25 50 75]) max(col)]';
      % confidence interval for the median
      est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
      s(6,i) = max([s(3,i)-est, s(2,i)]);
      s(7,i) = min([s(3,i)+est, s(4,i)]);
      % whiskers out to the last point within the desired inter-quartile range
      IQR = opt.maxwhisker*(s(4,i)-s(2,i));
      whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
      whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
      % outliers beyond 1 and 2 inter-quartile ranges
      ol = (col < s(2,i)-IQR & col >= s(2,i)-2*IQR) | (col > s(4,i)+IQR & col <= s(4,i)+2*IQR);
      ol2 = col < s(2,i)-2*IQR | col > s(4,i)+2*IQR;
      outliers = col(ol);
      outliers2 = col(ol2);
      
      oll1 = (col < s(2,i)-IQR & col >= s(2,i)-2*IQR);
      olh1 = (col > s(4,i)+IQR & col <= s(4,i)+2*IQR);
      oll2 = col < s(2,i)-2*IQR;
      olh2 = col > s(4,i)+2*IQR;
      
      ind = 1:numel(col);
      out.names   = opt.names;
      out.indn.l1{i}  = ind(oll1);
      out.indn.l2{i}  = ind(oll2);
      out.indn.h1{i}  = ind(olh1);
      out.indn.h2{i}  = ind(olh2);
      out.matn.l1{i}  = oll1;
      out.matn.l2{i}  = oll2;
      out.matn.h1{i}  = olh1;
      out.matn.h2{i}  = olh2;

      if exist('sorti','var'), out.sorti   = sorti; end
      out.indo.l1{out.sortj(i)} = ind(oll1);
      out.indo.l2{out.sortj(i)} = ind(oll2);
      out.indo.h1{out.sortj(i)} = ind(olh1);
      out.indo.h2{out.sortj(i)} = ind(olh2);
      out.mato.l1{out.sorti(i)} = oll1;
      out.mato.l2{out.sorti(i)} = oll2;
      out.mato.h1{out.sorti(i)} = olh1;
      out.mato.h2{out.sorti(i)} = olh2;
      
      outliers_x = [outliers_x; i*ones(size(outliers))];
      outliers_y = [outliers_y; outliers];
      outliers2_x = [outliers2_x; i*ones(size(outliers2))];
      outliers2_y = [outliers2_y; outliers2];
    elseif (nd == 1)
      % all statistics collapse to the value of the point
      s(:,i) = col;
      % single point data sets are plotted as outliers.
      outliers_x = [outliers_x; i];
      outliers_y = [outliers_y; col];
    else
      % no statistics if no points
      s(:,i) = NaN;
    end
  end

  % Note which boxes don't have enough stats
  chop = find(box <= 1);

  % Draw a box around the quartiles, with width proportional to the number of
  % items in the box. Draw notches if desired.
  box = 0.15 + box*(0.30/max(box));
  quartile_x = ones(11,1)*[1:nc] + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
  quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

  % Draw a line through the median
  median_x = ones(2,1)*[1:nc] + [-a;+a]*box;
  median_y = s([3,3],:);

  % Chop all boxes which don't have enough stats
  quartile_x(:,chop) = [];
  quartile_y(:,chop) = [];
  whisker_x(:,[chop,chop+nc]) = [];
  whisker_y(:,[chop,chop+nc]) = [];
  median_x(:,chop) = [];
  median_y(:,chop) = [];
  
  % Add caps to the remaining whiskers
  cap_x = whisker_x;
  cap_x(1,:) = cap_x(1,:) - 0.1;
  cap_x(2,:) = cap_x(2,:) + 0.1;
  cap_y = whisker_y([1,1],:);

  
  
  
  
  %% Do the plot
  children0=get(gca,'Children');
  
  if opt.vertical
    if opt.fill
      qn=size(quartile_x,2);
      for i=1:qn
        fill(quartile_x(:,i), quartile_y(:,i),'b-','FaceColor',1-max(0,(1-opt.groupcolor(i,:))/3),'EdgeColor',[0 0 0.5]); 
        if i==1, hold on; end
        plot(cap_x(:,[i i+qn]), cap_y(:,[i i+qn]),'Color',[0 0 0.5]); %max(0,opt.groupcolor(i,:)-0.5))
        plot(whisker_x(:,[i i+qn]), whisker_y(:,[i i+qn]),'Color',[0 0 0.5]); %max(0,opt.groupcolor(i,:)-0.5))
        if opt.groupcolor(i,1)>0.2 && opt.groupcolor(i,2)<0.5 && opt.groupcolor(i,3)<0.5
          plot(median_x(:,i), median_y(:,i),'Color',[0.5 0 0])
        else
          plot(median_x(:,i), median_y(:,i),'Color',[1 0 0])
        end
        plot(quartile_x(:,i), quartile_y(:,i),'Color',[0 0 0.5]); 
      end
      plot(outliers_x,  min(opt.ylim(2),max(opt.ylim(1),outliers_y)) ,'MarkerSize',...
        max(4,min(8,80/nc)),'Marker',opt.symbol(1),'MarkerEdgeColor','r','LineStyle','none')
      plot(outliers2_x, min(opt.ylim(2),max(opt.ylim(1),outliers2_y)),'MarkerSize',...
        max(4,min(8,80/nc)),'Marker',opt.symbol(2),'MarkerEdgeColor','r','LineStyle','none');
    else
      plot(quartile_x, quartile_y, 'b-')
      hold on
      plot(whisker_x, whisker_y, 'b-')
      plot(cap_x, cap_y, 'b-')
      plot(median_x, median_y, 'r-')
      plot(outliers_x,  min(opt.ylim(2),max(opt.ylim(1),outliers_y)) , [opt.symbol(1),'r'])
      plot(outliers2_x, min(opt.ylim(2),max(opt.ylim(1),outliers2_y)), [opt.symbol(2),'r']);
      plot(median_x, median_y, 'r-') % yes, two times... otherwise the last median is below the box ...
    end
    
    
    % add labels
    linecolor = [0.8 0.8 0.8];
    set(gca,'XTick',1:numel(opt.names),'XTickLabel',opt.names,'TickLength',[0 0],'xlim',[0.4 numel(opt.names)+0.6]);
    if ~isempty(opt.ylim)
      ylim(gca,opt.ylim);
    end

    if ~isempty(opt.fontsize)
      set(gca,'FontSize',opt.fontsize);
    end

    % plot yticks
    if opt.ygrid
      ytick=get(gca,'YTick');
      if numel(ytick)<5, ytick=interp1(ytick,1:0.5:numel(ytick)); elseif numel(ytick)>10, ytick=ytick(1:2:end); end
      if ytick(1)==opt.ylim(1),   ytick(1)=[];   end
      if ytick(end)==opt.ylim(2), ytick(end)=[]; end
      h1=plot(repmat([0;numel(opt.names)+1],1,numel(ytick)),[ytick;ytick],'Color',linecolor);
      uistack(h1,'bottom')
    end


  else
    if opt.fill
      qn=size(quartile_x,2);
      for i=1:qn
        fill(quartile_y(:,i), quartile_x(:,i),'b-','FaceColor',1-max(0,(1-opt.groupcolor(i,:))/3),'EdgeColor',[0 0 0.5]); 
        if i==1, hold on; end
        plot(cap_y(:,[i i+qn]), cap_x(:,[i i+qn]),'Color',[0 0 0.5]);
        plot(whisker_y(:,[i i+qn]), whisker_x(:,[i i+qn]),'Color',[0 0 0.5]); 
        if opt.groupcolor(i,1)>0.2 && opt.groupcolor(i,2)<0.5 && opt.groupcolor(i,3)<0.5
          plot(median_y(:,i), median_x(:,i),'Color',[0.5 0 0])
        else
          plot(median_y(:,i), median_x(:,i),'Color',[1 0 0])
        end
        plot(quartile_y(:,i), quartile_x(:,i),'Color',[0 0 0.5]); 
      end
      plot(min(opt.ylim(2),max(opt.ylim(1),outliers_y)) , outliers_x ,'MarkerSize',...
        max(4,min(8,80/nc)),'Marker',opt.symbol(1),'MarkerEdgeColor','r','LineStyle','none')
      plot(min(opt.ylim(2),max(opt.ylim(1),outliers2_y)), outliers2_x ,'MarkerSize',...
        max(4,min(8,80/nc)),'Marker',opt.symbol(2),'MarkerEdgeColor','r','LineStyle','none');
      
    else
      plot(quartile_y, quartile_x, 'b-')
      hold on
      plot(whisker_y, whisker_x, 'b-')
      plot(cap_y, cap_x, 'b-')
      plot(median_x, median_y, 'r-')
      plot(min(opt.ylim(2),max(opt.ylim(1),outliers_y)),  outliers_x,  [opt.symbol(1),'r'])
      plot(min(opt.ylim(2),max(opt.ylim(1),outliers2_y)), outliers2_x, [opt.symbol(2),'r']);
      plot(median_y, median_x, 'r-')     
    end
    
    
    % add labels
    linecolor = [0.8 0.8 0.8];
    set(gca,'YTick',1:numel(opt.names),'YTickLabel',opt.names,'TickLength',[0 0],'ylim',[0.4 numel(opt.names)+0.6]);
    if ~isempty(opt.ylim)
      xlim(gca,opt.ylim);
    end

    if ~isempty(opt.fontsize)
      set(gca,'FontSize',opt.fontsize);
    end

    % plot yticks
    if opt.ygrid
      ytick=get(gca,'XTick');
      if numel(ytick)<5, ytick=interp1(ytick,1:0.5:numel(ytick)); elseif numel(ytick)>10, ytick=ytick(1:2:end); end
      if ytick(1)==opt.ylim(1),   ytick(1)=[];   end
      if ytick(end)==opt.ylim(2), ytick(end)=[]; end
      h1=plot([ytick;ytick],repmat([0;numel(opt.names)+1],1,numel(ytick)),'Color',linecolor);
      uistack(h1,'bottom')
    end

  end

  %%
  hold off

end