function cg_check_mesh_cov(vargin)
%cg_check_mesh_cov to check covriance across sample
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_check_mesh_cov.m 475 2013-04-19 08:32:39Z gaser $

global fname jY h1 h2 YpY
rev = '$Rev: 475 $';

if nargin == 1
  P = char(vargin.data_surf);
  measure = char(vargin.measure);
else
  P = spm_select(Inf,'gifti','Select surfaces');
  measure = spm_input('Which data to check?',1,'m',{'Surface values','Surface coordinates (shape)'},[1 2],1);
end

V = gifti(P);
n = size(P,1);

%-Start progress plot
%-----------------------------------------------------------------------
if measure == 1
  data_array = zeros(length(V(1).cdata),n);
else
  data_array = zeros(prod(size(V(1).vertices)),n);
end
if measure == 1
  for i = 1:n
    data_array(:,i) = V(i).cdata;
  end
else
  for i = 1:n
    data_array(:,i) = V(i).vertices(:);
  end
end

YpY = (data_array'*data_array)/n;

% normalize YpY
d = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd = d*d';
YpY = YpY./(dd+eps);
t = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n+1:end) = sign(diag(YpY));

YpYsum = sum(YpY,1);
[iY, jY] = sort(YpYsum, 2, 'descend');
YpYsorted = YpY(jY,jY);
Nsorted = P(jY,:);

% extract mean covariance
mean_cov = zeros(n,1);
for i=1:n
  % extract row for each subject
  cov0 = YpY(i,:);
  % remove cov with its own
  cov0(i) = [];
  mean_cov(i) = mean(cov0);
end

threshold_cov = mean(mean_cov) - 2*std(mean_cov);

fprintf('Mean covariance: %3.2f\n',mean(mean_cov));

[tmp fname] = spm_str_manip(char(P),'C');
fprintf('Compressed filenames: %s  \n',tmp);

% print suspecious files with cov>0.9
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.9);
if ~isempty(indx) & (sqrt(length(indx)) < 0.5*n)
  fprintf('\nUnusual large covariances (check that subjects are not identical):\n');
  for i=1:length(indx)
    % exclude diagonal
    if indx(i) ~= indy(i)
      % report file with lower mean covariance first
      if mean_cov(indx(i)) < mean_cov(indy(i))
        fprintf('%s and %s: %3.3f\n',fname.m{indx(i)},fname.m{indy(i)},YpY(indx(i),indy(i)));
      else
        fprintf('%s and %s: %3.3f\n',fname.m{indy(i)},fname.m{indx(i)},YpY(indy(i),indx(i)));
      end
    end
  end
end

% sort files
fprintf('\nMean covariance for data below 2 standard deviations:\n');
[mean_cov_sorted, ind] = sort(mean_cov,'descend');
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

for i=n_thresholded:n
  fprintf('%s: %3.3f\n',P(ind(i),:),mean_cov_sorted(i));
end

Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
FS    = spm('FontSizes');

xpos = 2*(0:n-1)/(n-1);
for i=1:n
  text(xpos(i),mean_cov(i),fname.m{i},'FontSize',FS(8),'HorizontalAlignment','center')
end

hold on
cg_boxplot({mean_cov});
set(gca,'XTick',[],'XLim',[-.25 2.25]);

title(sprintf('Boxplot: mean covariance  \nCommon filename: %s*%s',fname.s,fname.e),'FontSize',FS(12),'FontWeight','Bold');
ylabel('<----- low (poor quality) --- mean covariance --- large (good quality)------>  ','FontSize',FS(10),'FontWeight','Bold');
xlabel('<----- first --- file order --- last ------>  ','FontSize',FS(10),'FontWeight','Bold');
hold off

% covariance
f = figure(4);
ws = spm('Winsize','Graphics');

set(f,'Name','Click in image to get file names','NumberTitle','off');
h = datacursormode(f);
set(h,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
set(f,'MenuBar','none','Position',[10 10 ws(3) ws(3)]);

cmap = [gray(64); hot(64)];

% scale YpY to 0..1
mn = min(YpY(:));
mx = max(YpY(:));
YpY_scaled = (YpY - mn)/(mx - mn);
YpYsorted_scaled = (YpYsorted - mn)/(mx - mn);

% show upper right triangle in gray
ind_tril = find(tril(ones(size(YpY))));
ima = YpY_scaled;
ima(ind_tril) = 0;
ima(ind_tril) = 1 + 1/64 + YpY_scaled(ind_tril);
image(64*ima)
a = gca;
set(a,'XTickLabel','','YTickLabel','');
axis image
xlabel('<----- first --- file order --- last ------>  ','FontSize',10,'FontWeight','Bold');
ylabel('<----- last --- file order --- first ------>  ','FontSize',10,'FontWeight','Bold');
title('Covariance','FontSize',12,'FontWeight','Bold');
colormap(cmap)

% Close button
hCloseButton = uicontrol(f,...
        'position',[round(ws(3)/2)-75 10 150 20],...
        'style','Pushbutton',...
        'string','Close windows',...
        'callback','for i=2:20, try close(i); end; end;',...
        'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on');

% ordered covariance
f = figure(5);
set(f,'Name','Click in image to get file names','NumberTitle','off');
h = datacursormode(f);
set(h,'UpdateFcn',@myupdatefcn_ordered,'SnapToDataVertex','on','Enable','on');
set(f,'MenuBar','none','Position',[11+ws(3) 10 ws(3) ws(3)]);

% show upper right triangle in gray
ind_tril = find(tril(ones(size(YpY))));
ima = YpYsorted_scaled;
ima(ind_tril) = 0;
ima(ind_tril) = 1 + 1/64 + YpYsorted_scaled(ind_tril);
image(64*ima)
if n_thresholded <= n
  hold on
  line([n_thresholded-0.5, n_thresholded-0.5], [0.5,n_thresholded-0.5])
  line([0.5,n_thresholded-0.5],[n_thresholded-0.5, n_thresholded-0.5])
  hold off
end
a = gca;
set(a,'XTickLabel','','YTickLabel','');
axis image
xlabel('<----- high --- mean covariance --- low ------>  ','FontSize',10,'FontWeight','Bold');
ylabel('<----- low --- mean covariance --- high ------>  ','FontSize',10,'FontWeight','Bold');
title({'Sorted Covariance','Blue line indicates 2-SD threshold'},'FontSize',12,'FontWeight','Bold');
colormap(cmap)

% check for replicates
for i=1:n
  for j=1:n
  if (i>j) & (mean_cov(i) == mean_cov(j))
    [s,differ] = unix(['diff ' P(i) ' ' P(j)]);
    if (s==0), fprintf(['\nWarning: ' P(i) ' and ' P(j) ' are same files?\n']); end
  end
  end
end

%-End
%-----------------------------------------------------------------------

show = spm_input('Show files with poorest cov?',1,'yes|no',[1 0],2);
if show

  if measure == 1
    data = data_array(:);
    data(isnan(data) | isinf(data)) = [];
    mn_data = min(data);
    mx_data = max(data);
  end
    
  number = min([n 16]);
  number = spm_input('How many files ?','+1','e',number);
  
  list = str2mat(P(ind(n:-1:1),:));
  list2 = list(1:number,:);

  for i=1:number
    h = spm_mesh_render(deblank(list2(i,:)));
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_str_manip(list2(i,:),'k40d'),'NumberTitle','off');
    spm_mesh_render('ColourMap',h.axis,jet);
    if measure == 1
      h2 = spm_mesh_render('ColourBar',h.axis,'on');
      set(h2.colourbar,'Ylim',[1.25*mn_data 0.8*mx_data]);
    end
  end
end
return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global fname YpY
pos = get(event_obj, 'Position');
h = gca;

x = pos(1);
y = pos(2);

txt = {sprintf('Covariance: %3.3f',YpY(x,y)),fname.m{x},fname.m{y}};

return

%-----------------------------------------------------------------------
function txt = myupdatefcn_ordered(obj, event_obj)
%-----------------------------------------------------------------------
global fname jY YpY
pos = get(event_obj, 'Position');
h = gca;

x = jY(pos(1));
y = jY(pos(2));

txt = {sprintf('Covariance: %3.3f',YpY(x,y)),fname.m{x},fname.m{y}};

return

