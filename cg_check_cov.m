function cg_check_cov(vargin)
%cg_check_cov to check covriance across sample
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. normalized images)
%
% Surfaces have to be same size
%_______________________________________________________________________
% Christian Gaser
% $Id$

global fname H YpY data_array pos ind_sorted mean_cov FS P issurf mn_data mx_data V sample xml_files
rev = '$Rev$';

if nargin == 1
  if isfield(vargin,'data_vbm')
    P = char(vargin.data_vbm);
    sep = vargin.gap;
  else
    P = char(vargin.data_surf);
  end
  if isempty(vargin.nuisance)
    G = [];
  else
    G = vargin.nuisance.c;
  end
  if isempty(vargin.coding)
    sample = ones(1,size(P,1));
  else
    sample = vargin.coding.c;
  end
  if isempty(vargin.xml)
    xml_files = [];
  else
    xml_files = char(vargin.xml.data_xml);
  end
end


if nargin < 1
  P = spm_select(Inf,'image','Select images');
end

n_subjects = size(P,1);

if size(sample,2) == 1
  sample = sample';
end

if size(sample,2) ~= n_subjects
  error('Vector for sample coding must have same length as sample');
end

if ~isempty(xml_files)
  if size(xml_files,1) ~= n_subjects
    error('Number of xml-files must have same as sample size');
  end
  
  QM = ones(n_subjects,3);
  QM_names = str2mat('Noise','Bias','PQ processibility');
  spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
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

if nargin < 1
  def_nuis = spm_input('Variable to covariate out (nuisance parameter)?','+1','yes|no',[1 0],2);
  if def_nuis
    G = spm_input('Nuisance parameter:','+1','r',[],n_subjects);
  else
    G = [];
  end
  sep = spm_input('Separation between points to speed up','+1','e',0,1);
  sample = [];
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
FS    = spm('FontSizes');

pos = struct(...
    'fig',   [10 10 1.2*ws(3) ws(3)],... % figure
    'cbar',  [0.715 0.050 0.02 0.50],... % colorbar for correlation matrix
    'corr',  [-0.05 0.050 0.80 0.80],... % correlation matrix
    'quit',  [0.775 0.900 0.20 0.05],... % quit button
    'show',  [0.775 0.850 0.20 0.05],... % button to show worst cases
    'boxp',  [0.775 0.800 0.20 0.05],... % button to display boxplot
    'text',  [0.775 0.625 0.20 0.15],... % textbox
    'slice', [0.775 0.050 0.20 0.40],... % two single images according to position of mouse pointer
    'slider',[0.775 0.000 0.20 0.03]);   % slider for z-slice

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
  clear Y
else
  % consider image aspect ratio
  pos.slice = [0.775 0.050 0.20 0.40*V(1).dim(2)/V(1).dim(1)];

  slices = 1:sep:V(1).dim(3);

  vol = zeros(n_subjects, prod(length(1:sep:V(1).dim(1))*length(1:sep:V(1).dim(2))));
  YpY = zeros(n_subjects);
  data_array = zeros([V(1).dim(1:2) n_subjects]);

  %-Start progress plot
  %-----------------------------------------------------------------------
  spm_progress_bar('Init',V(1).dim(3),'Check correlation','planes completed')

  for j=slices

    M  = spm_matrix([0 0 j]);

    for i = 1:n_subjects
      img = spm_slice_vol(V(i),M,V(1).dim(1:2),[1 0]);
      img = img(1:sep:V(1).dim(1),1:sep:V(1).dim(2));
      img(isnan(img)) = 0;
      vol(i,:) = img(:);
    end

    % find mask with non-zeros voxels in all images
    mask = all(vol ~= 0);

    if sum(mask)>0
      % make sure data is zero mean
      Y = vol(:,mask);
      Y = Y - repmat(mean(Y,2), [1 sum(mask)]);

      % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
      if ~isempty(G) 
        Ymean = repmat(mean(Y), [n_subjects 1]);
        Y = Y - G*(pinv(G)*Y) + Ymean;
      end
      YpY = YpY + (Y*Y')/n_subjects;
    end 
    spm_progress_bar('Set',j);  
  end

  spm_progress_bar('Clear');

end


% normalize YpY
d = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd = d*d';
YpY = YpY./(dd+eps);
t = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n_subjects+1:end) = sign(diag(YpY));

YpYsum = sum(YpY,1);
[iY, jY] = sort(YpYsum, 2, 'descend');

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

% sort files
[mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
threshold_cov = mean(mean_cov) - 2*std(mean_cov);
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

if ~isempty(n_thresholded)
  fprintf('\nMean correlation for data below 2 standard deviations:\n');
  for i=n_thresholded:n_subjects
    fprintf('%s: %3.3f\n',P(ind_sorted(i),:),mean_cov_sorted(i));
  end
end

% create figure
H.figure = figure(2);
clf(H.figure);

set(H.figure,'MenuBar','none','Position',pos.fig,...
    'Name','Click in correlation matrix to display slices','NumberTitle','off');
    
H.ax = axes('Position',pos.corr,'Parent',H.figure);
cm = datacursormode(H.figure);
set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
try, set(cm,'NewDataCursorOnClick',false); end

% create two colormaps
cmap = [hot(64); gray(64)];

% scale YpY to 0..1
mn = min(YpY(:));
mx = max(YpY(:));
YpY_scaled = (YpY - mn)/(mx - mn);

% show only lower left triangle
ind_tril = find(tril(ones(size(YpY))));
ima = zeros(size(YpY_scaled));
ima(ind_tril) = YpY_scaled(ind_tril);
image(64*ima)
set(gca,'XTickLabel','','YTickLabel','');
axis image
xlabel('<----- First ---      File order      --- Last ------>  ','FontSize',FS(8),'FontWeight','Bold');
ylabel('<----- Last ---      File order      --- First ------>  ','FontSize',FS(8),'FontWeight','Bold');
title('Sample Correlation Matrix','FontSize',FS(10),'FontWeight','Bold');
colormap(cmap)

% add colorbar
axes('Position',pos.cbar,'Parent',H.figure);
image((64:-1:1)');

% display YTick with 5 values (limit accuracy for floating numbers)
set(gca,'XTickLabel','','YTickLabel','','YTick',linspace(1,64,5), 'YTickLabel',...
  round(100*linspace(mx,mn,5))/100,'TickLength',[0 0]);
xlabel(gca,'Correlation')


% add button for closing all windows
H.quit = uicontrol(H.figure,...
        'string','Quit','Units','normalized',...
        'position',pos.quit,...
        'style','Pushbutton','HorizontalAlignment','center',...
        'callback','for i=2:26, try close(i); end; end;',...
        'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on');

H.show = uicontrol(H.figure,...
        'string','Check the worst data','Units','normalized',...
        'position',pos.show,...
        'style','Pushbutton','HorizontalAlignment','center',...
        'callback',@check_worst_data,...
        'ToolTipString','Display worst files',...
        'Interruptible','on','Enable','on');

if isempty(xml_files)
  H.boxp = uicontrol(H.figure,...
        'string','Boxplot of mean correlation','Units','normalized',...
        'position',pos.boxp,...
        'style','Pushbutton','HorizontalAlignment','center',...
        'callback',{@show_mean_boxplot, mean_cov, 'Mean correlation', 1},...
        'ToolTipString','Display boxplot of mean correlation',...
        'Interruptible','on','Visible','off');
else
  data_boxp = mean_cov;
  str  = { 'Boxplot...','Mean correlation',QM_names};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation', 1},...
           {@show_mean_boxplot, QM(:,1), QM_names(1,:), -1},...
           {@show_mean_boxplot, QM(:,2), QM_names(2,:), -1},...
           {@show_mean_boxplot, QM(:,3), QM_names(3,:), -1}};

  H.boxp = uicontrol(H.figure,...
        'string',str,'Units','normalized',...
        'position',pos.boxp,'UserData',tmp,...
        'style','PopUp','HorizontalAlignment','center',...
        'callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Display boxplot',...
        'Interruptible','on','Visible','on');
end

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
function show_mean_boxplot(data_boxp, name, quality_order)
%-----------------------------------------------------------------------
global fname H YpY pos FS sample xml_files ind_sorted

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
  ylim_min = 0.99*min(data_boxp);
  ylim_max = 1.01*max(data_boxp);
  if ylim_min > min(data_boxp)
    ylim_min = 1.01*min(data_boxp);
  end
  if ylim_max < max(data_boxp)
    ylim_max = 0.99*max(data_boxp);
  end
  set(gca,'YLim',[ylim_min ylim_max]);
end

title(sprintf('Boxplot: %s  \nCommon filename: %s*%s',name,spm_file(fname.s,'short25'),fname.e),'FontSize',FS(10),'FontWeight','Bold');
if quality_order > 0
  ylabel(sprintf('<----- Low (poor quality) --- %s --- High (good quality)------>  ',name),'FontSize',FS(10),'FontWeight','Bold');
else
  ylabel(sprintf('<----- High rating (good quality) --- %s --- Low rating (poor quality)------>  ',name),'FontSize',FS(10),'FontWeight','Bold');
end
xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(10),'FontWeight','Bold');
hold off

% hide button again
if isempty(xml_files)
  set(H.boxp,'Visible','off');
end

% estimate sorted index new fo displaying worst files
if quality_order > 0
  [tmp, ind_sorted] = sort(data_boxp,'descend');
else
  [tmp, ind_sorted] = sort(data_boxp,'ascend');
  
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global V fname data_array H YpY pos

if isfield(H,'mm')
  slice_mm = get(H.mm,'Value');
else
  slice_mm = 0;
end

vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
Orig = V(1).mat\[0 0 0 1]';
sl = round(slice_mm/vx(3)+Orig(3));

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

% check whether mouse position is defined
if isfield(pos,'x')
  img = [data_array(:,:,pos.y) data_array(:,:,pos.x)]';

  % use gray scale colormap for values > 64
  axes('Position',pos.slice);
  image(65 + flipud(img))
  set(gca,'XTickLabel','','YTickLabel','');

  txt = {sprintf('Correlation: %3.3f',YpY(pos.x,pos.y)),[],['Top: ',spm_file(fname.m{pos.x},'short25')],...
      ['Bottom: ',spm_file(fname.m{pos.y},'short25')],[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  set(H.text,'String',txt);
end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global fname H YpY data_array pos issurf mn_data mx_data

pos_mouse = get(event_obj, 'Position');

% check for valid mouse position
if pos_mouse(1) > pos_mouse(2) || pos_mouse(1)>length(fname.m) || pos_mouse(2)>length(fname.m)
  txt = {''};
  return
end

% save position of mouse
pos.x = pos_mouse(1);
pos.y = pos_mouse(2);

% text info for data cursor window
if issurf
  txt = {sprintf('Correlation: %3.3f',YpY(pos.x,pos.y)),['Left: ',fname.m{pos.x}],...
    ['Right: ',fname.m{pos.y}]};
else
  txt = {sprintf('Correlation: %3.3f',YpY(pos.x,pos.y)),['Top: ',fname.m{pos.x}],...
    ['Bottom: ',fname.m{pos.y}]};
end

% text info for textbox
if issurf
  txt2 = {sprintf('Correlation: %3.3f',YpY(pos.x,pos.y)),[],['Left: ',...
    spm_file(fname.m{pos.x},'short25')],['Right: ',spm_file(fname.m{pos.y},'short25')]};
else
  txt2 = {sprintf('Correlation: %3.3f',YpY(pos.x,pos.y)),[],['Top: ',...
    spm_file(fname.m{pos.x},'short25')],['Bottom: ',spm_file(fname.m{pos.y},'short25')],...
    [],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
end      

set(H.text,'String',txt2)
axes('Position',pos.slice);

if issurf 
  % use indexed 2D-sheet to display surface data as image
  % check surface size to use indexed 2D map
  if (length(data_array(:,pos.x)) == 163842)
    ind = spm_load(fullfile(spm('dir'),'toolbox','vbm12','fsaverage','fsavg.index2D_256x128.txt'));
    img = [reshape(data_array(ind,pos.x),[256,128]) reshape(data_array(ind,pos.y),[256,128])];
    img = circshift(img,128);
  else
    img = [data_array(:,pos.y) data_array(:,pos.x)]';
  end
  
  % scale img to 0..64
  mn = min(data_array(:));
  mx = max(data_array(:));
  img = 64*((img - mn)/(mx-mn));
else
  % add slider for colume data
  set(H.mm,'Visible','on');
  img = [data_array(:,:,pos.y) data_array(:,:,pos.x)]';
end

% display image with 2nd colorbar (gray)
image(65 + flipud(img));
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

if issurf
  xlabel('2D surface maps');
end

return
