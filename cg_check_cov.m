function cg_check_cov(vargin)
%cg_check_cov to check covriance across sample
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%_______________________________________________________________________
% Christian Gaser
% $Id$

global fname H YpY YpYsorted data_array pos ind_sorted mean_cov FS P issurf mn_data mx_data V sample xml_files sorted
rev = '$Rev$';

% show data by fileorder
sorted = 0;

if nargin == 1
  P      = [];
  sample = [];
  G      = [];

  % read filenames for each sample and indicate sample parameter
  if isfield(vargin,'data_vbm')
    for i=1:numel(vargin.data_vbm)
      P = char([P; vargin.data_vbm{i}]);
      sample = [sample, i*ones(1,size(vargin.data_vbm{i},1))];
    end
    sep = vargin.gap;
  else
    for i=1:numel(vargin.data_surf)
      P = char([P; vargin.data_surf{i}]);
      sample = [sample, i*ones(1,size(vargin.data_surf{i},1))];
    end
  end
  
  n_subjects = size(P,1);
  
  if ~isempty(vargin.c)
    for i=1:numel(vargin.c)
      G = [G vargin.c{i}];
    end
  end

  if isempty(vargin.data_xml)
    xml_files = [];
  else
    xml_files = char(vargin.data_xml{1});
  end
  
else
  error('No argument give.');
end

if ~isempty(xml_files)

  if size(xml_files,1) ~= n_subjects
    error('XML-files must have the same number as sample size');
  end
  
  QM = ones(n_subjects,3);
  QM_names = str2mat('Noise','Bias','PQ processibility');

  spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
    % get basename for xml- and data files
    [pth, xml_name] = fileparts(deblank(xml_files(i,:)));
    [pth, data_name] = fileparts(deblank(P(i,:)));
    
    % remove leading 'vbm_'
    xml_name = xml_name(5:end);
    
    % check for filenames
    if isempty(strfind(data_name,xml_name))
      warning(sprintf('Please check file names because of deviating subject names\n: %s vs. %s\n',P(i,:),xml_files(i,:)));
    end
    
    xml = convert(xmltree(deblank(xml_files(i,:))));
    QM(i,:) = [str2double(xml.qam.QM.NCR) str2double(xml.qam.QM.ICR) str2double(xml.qam.QM.rms)];
    
    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');

else
  QM_names = '';
end

[pth,nam,ext] = spm_fileparts(deblank(P(1,:)));

if strcmp(ext,'.gii')
  issurf = 1;
else
  issurf= 0;
end

if issurf
  % load surface texture data
  spm_progress_bar('Init',n_subjects,'Load surfaces','subjects completed')

  V = gifti(deblank(P(1,:)));
  sz = length(V.cdata);

  Y = zeros(n_subjects,sz);
  tmp = V.cdata;
  tmp(isnan(tmp)) = 0;
  Y(1,:) = tmp';

  for i = 2:n_subjects
    V = gifti(deblank(P(i,:)));
    if length(V.cdata) ~= sz
      error(sprintf('File %s has different surface size than %',P(i,:),P(1,:)));
    end
    tmp = V.cdata;
    tmp(isnan(tmp)) = 0;
    Y(i,:) = tmp';
    spm_progress_bar('Set',i);  
  end

  spm_progress_bar('Clear');

else
  % load volume data
  V = spm_vol(deblank(P));
  
  % voxelsize and origin
  vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
  Orig = V(1).mat\[0 0 0 1]';

  if length(V)>1 & any(any(diff(cat(1,V.dim),1,1),1))
    error('images don''t all have same dimensions')
  end
  if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
    error('images don''t all have same orientation & voxel size')
  end
end

% add constant to nuisance parameter
if ~isempty(G)
  if size(G,1) ~= n_subjects
    G = G';
  end
  G = [ones(n_subjects,1) G];
end

% positions & font size
ws = spm('Winsize','Graphics');
FS = spm('FontSizes');

pos = struct(...
    'fig',   [10  10  1.2*ws(3) ws(3)],... % figure
    'cbar',  [0.240 0.950 0.300 0.020],... % colorbar for correlation matrix
    'corr',  [-0.02 0.050 0.825 0.825],... % correlation matrix
    'close', [0.775 0.900 0.200 0.050],... % close button
    'show',  [0.775 0.850 0.200 0.050],... % button to show worst cases
    'boxp',  [0.775 0.795 0.200 0.050],... % button to display boxplot
    'sort',  [0.775 0.750 0.200 0.050],... % button to enable ordered matrix
    'text',  [0.775 0.550 0.200 0.200],... % textbox
    'slice', [0.775 0.050 0.200 0.400],... % two single images according to position of mouse pointer
    'slider',[0.775 0.000 0.200 0.030]);   % slider for z-slice   

if issurf
  % rescue unscaled data min/max
  mn_data = min(Y(:));
  mx_data = max(Y(:));
  Y = Y - repmat(mean(Y,2), [1 length(V.cdata)]);

  % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
  if ~isempty(G) 
    Ymean = repmat(mean(Y), [n_subjects 1]);
    Y = Y - G*(pinv(G)*Y) + Ymean;
  end

  data_array = Y';
  YpY = (Y*Y')/n_subjects;

  % calculate residual mean square of mean adjusted Y
  Y = Y - repmat(mean(Y,1), [n_subjects 1]);
  MSE = sum(Y.*Y,2);

  clear Y
else
  % consider image aspect ratio
  pos.slice = [0.775 0.050 0.20 0.40*V(1).dim(2)/V(1).dim(1)];

  slices = 1:sep:V(1).dim(3);

  dimx = length(1:sep:V(1).dim(1));
  dimy = length(1:sep:V(1).dim(2));
  Y = zeros(n_subjects, prod(dimx*dimy));
  YpY = zeros(n_subjects);
  MSE = zeros(n_subjects,1);
  data_array = zeros([V(1).dim(1:2) n_subjects]);

  %-Start progress plot
  %-----------------------------------------------------------------------
  spm_progress_bar('Init',V(1).dim(3),'Check correlation','planes completed')

  for j=slices

    M  = spm_matrix([0 0 j 0 0 0 sep sep sep]);

    for i = 1:n_subjects
      img = spm_slice_vol(V(i),M,[dimx dimy],[1 0]);
      img(isnan(img)) = 0;
      Y(i,:) = img(:);
    end

    % make sure data is zero mean
    Y = Y - repmat(mean(Y,2), [1 prod(dimx*dimy)]);

    % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
    if ~isempty(G) 
      Ymean = repmat(mean(Y), [n_subjects 1]);
      Y = Y - G*(pinv(G)*Y) + Ymean;
    end
    YpY = YpY + (Y*Y')/n_subjects;

    % calculate residual mean square of mean adjusted Y
    Y = Y - repmat(mean(Y,1), [n_subjects 1]);
    MSE = MSE + sum(Y.*Y,2);

    spm_progress_bar('Set',j);  

  end

  spm_progress_bar('Clear');
end


% normalize YpY
d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd     = d*d';
YpY    = YpY./(dd+eps);
t      = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n_subjects+1:end) = sign(diag(YpY));

% extract mean correlation for each data set
mean_cov = zeros(n_subjects,1);
for i=1:n_subjects
  % extract row for each subject
  cov0 = YpY(i,:);

  % remove cov with its own
  cov0(i) = [];
  mean_cov(i) = mean(cov0);
end

[tmp fname] = spm_str_manip(char(P),'C');
fprintf('\nCompressed filenames: %s  \n',tmp);

% print suspecious files with cov>0.9
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.9);
if ~isempty(indx) & (sqrt(length(indx)) < 0.5*n_subjects)
  fprintf('\nUnusual large correlation (check that subjects are not identical):\n');
  for i=1:length(indx)
    % exclude diagonal
    if indx(i) ~= indy(i)
      % report file with lower mean correlation first
      if mean_cov(indx(i)) < mean_cov(indy(i))
        fprintf('%s and %s: %3.3f\n',fname.m{indx(i)},fname.m{indy(i)},YpY(indx(i),indy(i)));
      else
        fprintf('%s and %s: %3.3f\n',fname.m{indy(i)},fname.m{indx(i)},YpY(indy(i),indx(i)));
      end
    end
  end
end

% sort data
[mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
YpYsorted = YpY(ind_sorted,ind_sorted);

threshold_cov = mean(mean_cov) - 2*std(mean_cov);
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

if ~isempty(n_thresholded)
  fprintf('\nMean correlation for data below 2 standard deviations.\n');
  fprintf('This does not necessarily mean that you have to exclude these data. However, these data have to be carefully checked:\n');
  for i=n_thresholded:n_subjects
    fprintf('%s: %3.3f\n',P(ind_sorted(i),:),mean_cov_sorted(i));
  end
end

% create figure
H.figure = figure(2);
clf(H.figure);

set(H.figure,'MenuBar','none','Position',pos.fig,...
    'Name','Click in correlation matrix to display slices','NumberTitle','off');
    
cm = datacursormode(H.figure);
set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
try, set(cm,'NewDataCursorOnClick',false); end

show_matrix(YpY, sorted);

% create two colormaps
cmap = [hot(64); gray(64)];
colormap(cmap)

% add colorbar
axes('Position',pos.cbar,'Parent',H.figure);
image((1:64));

% display YTick with 5 values (limit accuracy for floating numbers)
set(gca,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
  round(100*linspace(min(YpY(:)),max(YpY(:)),5))/100,'TickLength',[0 0]);

% add button for closing all windows
H.close = uicontrol(H.figure,...
        'string','Close','Units','normalized',...
        'position',pos.close,...
        'style','Pushbutton','HorizontalAlignment','center',...
        'callback','for i=2:26, try close(i); end; end;',...
        'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on');

% check button
H.show = uicontrol(H.figure,...
        'string','Check worst data','Units','normalized',...
        'position',pos.show,...
        'style','Pushbutton','HorizontalAlignment','center',...
        'callback',@check_worst_data,...
        'ToolTipString','Display worst files',...
        'Interruptible','on','Enable','on');

% create popoup menu 
if isempty(xml_files)
  str  = { 'Boxplot...','Mean correlation','Mean squared error'};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation', 1},...
           {@show_mean_boxplot, MSE, 'Mean squared error', -1} };
else
  str  = { 'Boxplot...','Mean correlation','Mean squared error',QM_names};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation', 1},...
           {@show_mean_boxplot, MSE, 'Mean squared error', -1},...;
           {@show_mean_boxplot, QM(:,1), QM_names(1,:), -1},...
           {@show_mean_boxplot, QM(:,2), QM_names(2,:), -1},...
           {@show_mean_boxplot, QM(:,3), QM_names(3,:), -1} };
end

H.boxp = uicontrol(H.figure,...
        'string',str,'Units','normalized',...
        'position',pos.boxp,'UserData',tmp,...
        'style','PopUp','HorizontalAlignment','center',...
        'callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Display boxplot',...
        'Interruptible','on','Visible','on');

str  = { 'Correlation matrix...','Order by selected filename','Sorted by mean correlation'};
tmp  = { {@show_matrix, YpY, 0},...
         {@show_matrix, YpYsorted, 1} };

H.sort = uicontrol(H.figure,...
        'string',str,'Units','normalized',...
        'position',pos.sort,'UserData',tmp,...
        'style','PopUp','HorizontalAlignment','center',...
        'callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Sort matrix',...
        'Interruptible','on','Visible','on');

H.text = uicontrol(H.figure,...
        'Units','normalized','position',pos.text,...
        'String','Click in correlation matrix to display slices',...
        'style','text','HorizontalAlignment','center',...
        'ToolTipString','Select slice for display');

% add slider only for volume data
if ~issurf
  H.mm = uicontrol(H.figure,...
        'Units','normalized','position',pos.slider,...
        'Min',(1 - Orig(3))*vx(3),'Max',(V(1).dim(3) - Orig(3))*vx(3),...
        'style','slider','HorizontalAlignment','center',...
        'callback',@update_slices_array,...
        'ToolTipString','Select slice for display',...
        'SliderStep',[0.005 0.05],'Visible','off');

  update_slices_array;
end
   
show_mean_boxplot(mean_cov,'Mean correlation',1);

% check for replicates
for i=1:n_subjects
  for j=1:n_subjects
    if (i>j) & (mean_cov(i) == mean_cov(j))
      try
        [s,differ] = unix(['diff ' P(i,:) ' ' P(j,:)]);
        if (s==0), fprintf(['\nWarning: ' P(i,:) ' and ' P(j,:) ' are same files?\n']); end
      end
    end
  end
end

%-End
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj)
%-----------------------------------------------------------------------
global P ind_sorted issurf mn_data mx_data H

n = size(P,1);
number = min([n 24]);
number = spm_input('How many files ?',1,'e',number);
number = min([number 24]);
number = min([number size(P,1)]);
  
list = str2mat(P(ind_sorted(n:-1:1),:));
list2 = list(1:number,:);

if issurf
  % display single meshes and correct colorscale of colorbar
  for i=1:number
    h = spm_mesh_render(deblank(list2(i,:)));
    
    % shift each figure slightly
    if i==1
        pos = get(h.figure,'Position');
    else
        pos = pos - [20 20 0 0];
    end
    
    % remove menubar and toolbar, use filename as title
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(list2(i,:),'short50'),...
         'NumberTitle','off','Position',pos);
    spm_mesh_render('ColourMap',h.axis,jet);
    spm_mesh_render('ColourBar',h.axis,'on');
    spm_mesh_render('CLim',h,[mn_data mx_data]);
  end
else
  spm_check_registration(list2)
  set(H.boxp,'Visible','on');
end
return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
global P H FS pos sorted

% get sorting order
sorted = order;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 0.75 1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',pos.corr,'Parent',H.figure);

% scale data to 0..1
mn = min(data(:));
mx = max(data(:));
data_scaled = (data - mn)/(mx - mn);

% show only lower left triangle
ind_tril = find(tril(ones(size(data))));
ima = zeros(size(data));
ima(ind_tril) = data_scaled(ind_tril);
image(64*ima)
set(gca,'XTickLabel','','YTickLabel','');
axis image

if sorted
  xlabel('<----- Best ---      File Order      --- Worst ------>  ','FontSize',FS(8),'FontWeight','Bold');
  ylabel('<----- Worst ---      File Order      --- Best ------>  ','FontSize',FS(8),'FontWeight','Bold');
  title('Sorted Sample Correlation Matrix','FontSize',FS(10),'FontWeight','Bold');
else
  xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(8),'FontWeight','Bold');
  ylabel('<----- Last ---      File Order      --- First ------>  ','FontSize',FS(8),'FontWeight','Bold');
  title('Sample Correlation Matrix','FontSize',FS(10),'FontWeight','Bold');
end

return

%-----------------------------------------------------------------------
function show_mean_boxplot(data_boxp, name, quality_order)
%-----------------------------------------------------------------------
global fname H FS sample ind_sorted
 
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);

n_samples = max(sample);

xpos = cell(1,n_samples);
data = cell(1,n_samples);

for i=1:n_samples
  ind = find(sample == i);
  data{i} = data_boxp(ind);
  
  if n_samples == 1
    xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
  else
    xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
  end

  for j=1:length(ind)
    text(xpos{i}(j),data{i}(j),fname.m{ind(j)},'FontSize',FS(7),'HorizontalAlignment','center')
  end
end

hold on

opt = struct('groupnum',0,'ygrid',0,'groupcolor',jet(n_samples));

vbm_plot_boxplot(data,opt);

set(gca,'XTick',[],'XLim',[-.25 n_samples+1.25]);
if max(data_boxp) > min(data_boxp)
  yamp = max(data_boxp) - min(data_boxp);
  ylim_min = min(data_boxp) - 0.05*yamp;
  ylim_max = max(data_boxp) + 0.05*yamp;
  set(gca,'YLim',[ylim_min ylim_max]);
end

% add colored labels and title
title(sprintf('Boxplot: %s  \nCommon filename: %s*%s',name,spm_file(fname.s,'short25'),fname.e),'FontSize',FS(10),'FontWeight','Bold');
xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(10),'FontWeight','Bold');

xpos = -0.35 - n_samples*0.1;

if quality_order > 0
  text(xpos, ylim_min,'<----- Low rating (poor quality)','Color','red','Rotation',90,'HorizontalAlignment','left','FontSize',FS(10),'FontWeight','Bold')
  text(xpos, ylim_max,'High rating (good quality) ------>','Color','green','Rotation',90,'HorizontalAlignment','right','FontSize',FS(10),'FontWeight','Bold')
else
  text(xpos, ylim_max,'Low rating (poor quality) ------>','Color','red','Rotation',90,'HorizontalAlignment','right','FontSize',FS(10),'FontWeight','Bold')
  text(xpos, ylim_min,'<----- High rating (good quality)','Color','green','Rotation',90,'HorizontalAlignment','left','FontSize',FS(10),'FontWeight','Bold')
end
text(xpos, (ylim_max+ylim_min)/2,sprintf('%s',name),'Color','black','Rotation',90,'HorizontalAlignment','center','FontSize',FS(10),'FontWeight','Bold')

hold off

% estimate sorted index new for displaying worst files
if quality_order > 0
  [tmp, ind_sorted] = sort(data_boxp,'descend');
else
  [tmp, ind_sorted] = sort(data_boxp,'ascend');
  
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global V fname data_array H YpY pos sorted ind_sorted

if isfield(H,'mm')
  slice_mm = get(H.mm,'Value');
else
  slice_mm = 0;
end

vx   =  sqrt(sum(V(1).mat(1:3,1:3).^2));
Orig = V(1).mat\[0 0 0 1]';
sl   = round(slice_mm/vx(3)+Orig(3));

% if slice is outside of image use middle slice
if (sl>V(1).dim(3)) | (sl<1)
  sl = round(V(1).dim(3)/2);
end

M  = spm_matrix([0 0 sl]);

for i = 1:length(V)
  img = spm_slice_vol(V(i),M,V(1).dim(1:2),[1 0]);
  img(isnan(img)) = 0;
  
  % scale image according to mean
  data_array(:,:,i) = img/mean(img(find(img ~= 0)));
end

% enhance contrast and scale image to 0..64
mn = min(data_array(:));
mx = max(data_array(:));
data_array = 64*((data_array - mn)/(mx-mn));

if sorted
  if isfield(pos,'x')
    x = ind_sorted(pos.x);
    y = ind_sorted(pos.y);
  end
else
  if isfield(pos,'x')
    x = pos.x;
    y = pos.y;
  end
end

% check whether mouse position is defined
if isfield(pos,'x')
  img = [data_array(:,:,y) data_array(:,:,x)]';

  % use gray scale colormap for values > 64
  axes('Position',pos.slice);
  image(65 + flipud(img))
  set(gca,'XTickLabel','','YTickLabel','');

  txt = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Top: ',spm_file(fname.m{x},'short25')],...
      ['Bottom: ',spm_file(fname.m{y},'short25')],[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  set(H.text,'String',txt);
end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global fname H YpY YpYsorted data_array pos issurf mn_data mx_data ind_sorted sorted

pos_mouse = get(event_obj, 'Position');

% check for valid mouse position
if pos_mouse(1) > pos_mouse(2) || pos_mouse(1)>length(fname.m) || pos_mouse(2)>length(fname.m)
  txt = {''};
  return
end

% save position of mouse
pos.x = pos_mouse(1);
pos.y = pos_mouse(2);

if sorted
  if isfield(pos,'x')
    x = ind_sorted(pos.x);
    y = ind_sorted(pos.y);
  end
else
  if isfield(pos,'x')
    x = pos.x;
    y = pos.y;
  end
end

% text info for data cursor window
if issurf
  txt = {sprintf('Correlation: %3.3f',YpY(x,y)),['Left: ',fname.m{x}],...
    ['Right: ',fname.m{y}]};
else
  txt = {sprintf('Correlation: %3.3f',YpY(x,y)),['Top: ',fname.m{x}],...
    ['Bottom: ',fname.m{y}]};
end

% text info for textbox
if issurf
  txt2 = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Left: ',...
    spm_file(fname.m{x},'short25')],['Right: ',spm_file(fname.m{y},'short25')]};
else
  txt2 = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Top: ',...
    spm_file(fname.m{x},'short25')],['Bottom: ',spm_file(fname.m{y},'short25')],...
    [],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
end      

set(H.text,'String',txt2);
axes('Position',pos.slice);

if issurf 
  % use indexed 2D-sheet to display surface data as image
  % check surface size to use indexed 2D map
  if (length(data_array(:,x)) == 163842)
    ind = spm_load(fullfile(spm('dir'),'toolbox','vbm12','templates_surfaces','fsavg.index2D_256x128.txt'));
    img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
    img = circshift(img,128);
  else
    img = [data_array(:,y) data_array(:,x)]';
  end
  
  % scale img to 0..64
  mn = min(data_array(:));
  mx = max(data_array(:));
  img = 64*((img - mn)/(mx-mn));
else
  % add slider for colume data
  set(H.mm,'Visible','on');
  img = [data_array(:,:,y) data_array(:,:,x)]';
end

% display image with 2nd colorbar (gray)
image(65 + flipud(img));
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

if issurf
  xlabel('2D surface maps');
end

return
