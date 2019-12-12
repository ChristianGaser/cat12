% [signal_change0, xyz] = cat_stat_confplot_spm(SPM,xSPM,hReg,scale,names,Ic)
%
% SPM, xSPM, hReg - parameters saved in workspace
% scale    - scale factor for percent signal change of data
% name     - optional names of columns given as {'name1','name2'...}
% Ic       - number of contrast (usually 1 for effects of interest)   
%
% signal_change0 - scaled beta to obtain percent signal change
% xyz            - coordinates of local cluster maximum     
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_stat_confplot_spm.m 38 2014-04-09 09:01:05Z gaser $

global xY H

try
    [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg);
catch
    [hReg xSPM SPM] = spm_results_ui('Setup');
    [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg);
end

CI    = 1.6449;         % = spm_invNcdf(1 - 0.05);

%-Colour specifications
%-----------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 0 0];

Cplot = 'Parameter estimates';

%-Specify VOI
%-----------------------------------------------------------------------
xY.def = spm_input('VOI definition...',1,'b',...
      {'sphere','box','cluster','voxel'},[],4);
Q      = ones(1,size(xSPM.XYZmm,2));

% default font size and boxplot
if ~exist('H','var') | (exist('H','var') & isempty(H))
  H.fs = 18;
  H.boxplot = 1;
  H.medianplot = 0;
  H.rawdata = 0;
end

if ~exist('scale','var')
    scale = spm_input('Scaling factor','+1', 'm',['1 (no scaling)|100 (for percent signal change)'],[1 100],2);
end

if ~exist('colored','var')
    colored = spm_input('Boxplot','+1', 'm',['Colored|Define colors|Grey'],[2 1 0],1);
  switch colored
  case 0
    groupcolor = [0.7 0.7 0.7];
  case 1
    groupcolor = spm_input('Colors','!+0','r',[],[n_effects 3]);
  case 2
    groupcolor = [];
  end
else
  groupcolor = [];
end

switch xY.def

  case 'sphere'
  %---------------------------------------------------------------
  xY.spec = spm_input('VOI radius (mm)','!+0','r',0,1,[0,Inf]);
  d     = [xSPM.XYZmm(1,:) - xyz(1);
  xSPM.XYZmm(2,:) - xyz(2);
  xSPM.XYZmm(3,:) - xyz(3)];
  Q     = find(sum(d.^2) <= xY.spec^2);
  XYZstr = sprintf(' averaged in sphere (radius %d mm)', xY.spec);
  xY.string = sprintf('sphere_%dmm_at_%g_%g_%gmm',xY.spec,xyz);

  case 'box'
  %---------------------------------------------------------------
  xY.spec = spm_input('box dimensions [x y z] {mm}',...
      '!+0','r','0 0 0',3);
  Q     = find(all(abs(xSPM.XYZmm - xyz*Q) <= xY.spec(:)*Q/2));
  XYZstr = sprintf(' averaged in box dimensions (%3.2f %3.2f %3.2f)', xY.spec);
  xY.string = sprintf('box_%g_%g_%gmm_at_%g_%g_%gmm',xY.spec,xyz);

  case 'cluster'
  %---------------------------------------------------------------
  [x i] = spm_XYZreg('NearestXYZ',xyz,xSPM.XYZmm);
  A     = spm_clusters(xSPM.XYZ);
  Q     = find(A == A(i));
  XYZstr = sprintf(' averaged in cluster');
  xY.string = sprintf('cluster_at_%g_%g_%gmm',x);

  case 'voxel'
  %---------------------------------------------------------------
  d     = [xSPM.XYZmm(1,:) - xyz(1);
  xSPM.XYZmm(2,:) - xyz(2);
  xSPM.XYZmm(3,:) - xyz(3)];
  d2 = sum(d.^2);
  Q = find(d2==min(d2));
  XYZstr = sprintf(' in voxel');
  xY.string = sprintf('voxel_at_%g_%g_%gmm',xyz);
end

XYZ     = xSPM.XYZ(:,Q);    % coordinates

%-Parameter estimates:   beta = xX.pKX*xX.K*y;
%-Residual mean square: ResMS = sum(R.^2)/xX.trRV
%---------------------------------------------------------------

beta0  = spm_get_data(SPM.Vbeta, XYZ);
beta   = mean(beta0,2);

try
  fprintf('Read raw data...');
  y = spm_get_data(SPM.xY.VY, XYZ);
  fprintf(sprintf('%s',repmat('\b',1,150)));
  fprintf(sprintf('%s',repmat(' ',1,150)));
  H.y_found = 1;
catch
  warning('No raw data found! Please check that you have not moved your data.\n');
  H.y_found = 0;
end

ResMS  = spm_get_data(SPM.VResMS,XYZ);
ResMS  = mean(ResMS,2);
Bcov   = ResMS*SPM.xX.Bcov;
Bcov   = Bcov;

% determine which contrast
%---------------------------------------------------------------
if ~exist('Ic','var')
    Ic    = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});
end

TITLE = {Cplot XYZstr};

% find contrast and related columns in design matrix
%-------------------------------------------------------------- 
c = SPM.xCon(Ic).c';
[ind_x, ind_y] = find(c~=0);
ind_y = unique(ind_y);
X = SPM.xX.X;
X = X(:,ind_y);
n_effects = size(X,2);

if ~exist('names','var')
    define_names = spm_input('Define names?',1,'yes|use numbers',[1 0],1);
    if define_names
        names = [];
        for i=1:n_effects
            new_name = spm_input(['Name for parameter ' num2str(i)],1,'s');
            names = strvcat(names,new_name);
        end
    else
        names = num2str((1:n_effects)');
    end
end

% compute contrast of parameter estimates and 90% C.I.
%-------------------------------------------------------------- 
signal_change0 = SPM.xCon(Ic).c'*beta0;
signal_change  = SPM.xCon(Ic).c'*beta;
CI    = CI*sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));

if ~exist('repeated_anova','var')
  repeated_anova = ~isempty(SPM.xX.iB);
  [rw,cl] = find(SPM.xX.I == length(SPM.xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
  % expect that subject factors are 2nd colum, group 3rd column, time 4th column
  if cl(1) == 2
    n_groups = max(SPM.xX.I(:,3));
    count = 0;
    for i=1:n_groups
      ind_times{i} = count + (1:max(SPM.xX.I(find(SPM.xX.I(:,3)==i),4)));
      count = max(SPM.xX.I(find(SPM.xX.I(:,3)==i),4));
    end
    n_time = max(SPM.xX.I(:,4));
    n_groupsxtime = n_groups*n_time;
    if n_groupsxtime ~= n_effects
      repeated_anova = [];
    end
  else
    repeated_anova = [];
  end
end

% GUI figure
%--------------------------------------------------------------
H.h10 = figure(10);
%clf
set(H.h10,'Position',[0 800 150 550],'MenuBar','none','NumberTitle','off');
hNewButton = uicontrol(H.h10,...
    'Position',[20 500 110 20],...
    'Callback','cat_stat_confplot_spm',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Plot');
hClearButton = uicontrol(H.h10,...
    'position',[20 460 110 20],...
    'Callback','clear names Ic scale colored groupcolor repeated_anova H',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Reset variables');
hSaveButton = uicontrol(H.h10,...
    'position',[20 420 110 20],...
    'Callback',{@save_image},...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Save images');
hCloseButton = uicontrol(H.h10,...
    'position',[20 380 110 20],...
    'Callback','close(10,11,12)',...
    'Interruptible','on',...
    'Style','Pushbutton',...
    'String','Close windows');
hShowBoxplot = uicontrol(H.h10,...
    'position',[20 340 110 20],...
    'Callback',(@show_boxplot),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',H.boxplot,...
    'String','Show Boxplot');
hShowRawdata = uicontrol(H.h10,...
    'position',[20 300 110 20],...
    'Callback',(@show_rawdata),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',H.rawdata,...
    'String','Show Raw Data');
hShowMedianplot = uicontrol(H.h10,...
    'position',[20 260 110 20],...
    'Callback',(@show_medianplot),...
    'Interruptible','on',...
    'Style','CheckBox',...
    'Visible','off',...
    'Value',H.medianplot,...
    'String','Show Medianplot');
htext = uicontrol(H.h10,...
    'position',[20 200 60 20],...
    'Style','Text',...
    'String','Font Size');
hedit = uicontrol(H.h10,...
    'position',[80 200 50 20],...
    'Callback',{@set_font_size},...
    'Interruptible','on',...
    'Style','Edit',...
    'String',num2str(H.fs));

if H.y_found
  set(hShowBoxplot,'Visible','on');
  set(hShowRawdata,'Visible','on');
  if ~isempty(repeated_anova)
    set(hShowMedianplot,'Visible','on');
  end
end

% % signal change plot
%--------------------------------------------------------------

if ~exist('H','var') | (exist('H','var') & ~isfield(H,'h11'))
  H.h11 = figure(11);
  set(H.h11,'Position',[150 800 800 550],'NumberTitle','off','MenuBar','none');
else
  H.h11 = figure(11);
end

cla
hold on

% estimates
%--------------------------------------------------------------
h = bar(signal_change');
set(h,'FaceColor',Col(2,:));

% standard error
%--------------------------------------------------------------
for j = 1:length(signal_change)
  line([j j],([CI(j) 0 - CI(j)] + signal_change(j)),...
        'LineWidth',2,'Color',Col(3,:))
end

title(TITLE,'FontSize',14,'FontWeight','bold')
ylabel('parameter estimate','FontSize',12)
set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));

if exist('names','var')
  if size(names,1) == length(signal_change)
    set(gca,'XTickLabel',names,'TickLabelInterpreter','none');
  end
end

hold off

% prepare raw values for boxplot
%--------------------------------------------------------------
if H.y_found
	if scale == 1
			y_label = 'raw signal';
	else
			y_label = 'percent signal change';
	end
	
  y2 = cell(1,n_effects);
  for i=1:n_effects
    y2{i} = scale*y(find(X(:,i)==1),:);
  end
  
  if ~exist('H','var') | (exist('H','var') & ~isfield(H,'h12'))
     H.h12 = figure(12);
    set(H.h12,'Position',[950 800 800 550],'NumberTitle','off','MenuBar','none');
  else
     H.h12 = figure(12);
  end
  cla
  
  if H.rawdata
    vshowdata = 1;
  else
    vshowdata = 0;
  end
  
  if H.boxplot
    vbox = 1;
    voutliers = 1;
  else
    vbox = 0;
    voutliers = 0;
  end

  vstruct = struct('showdata',vshowdata,'box',vbox,'outliers',voutliers);
  if ~isempty(groupcolor)
    vstruct = setfield('groupcolor',groupcolor);
  end

  cat_plot_boxplot(y2,vstruct);
    
  TITLE = {'Boxplot of raw data ' XYZstr};
  title(TITLE,'FontSize',14,'FontWeight','bold')
  ylabel(y_label,'FontSize',12)
  set(gca,'XLim',[0.4 (length(signal_change) + 0.6)],'XTick',1:length(signal_change));
  
  if exist('names','var')
    if size(names,1) == length(signal_change)
      set(gca,'XTickLabel',names,'TickLabelInterpreter','none');
    end
  end  

  if ~isempty(repeated_anova) & H.medianplot
    hold on

    plot_data = zeros(n_effects,1);
    count = 1;
    for i=1:n_groups
      for j=1:length(ind_times{i})
        plot_data(count) = median(y2{ind_times{i}(j)});
        count = count + 1;
      end
      plot(ind_times{i},plot_data(ind_times{i}),'r','Linewidth',2);
    end
    hold off
  end
    
  set(gca(H.h12),'FontSize',H.fs);
end

set(gca(H.h11),'FontSize',H.fs);

%==========================================================================
function set_font_size(obj, event_obj)

global H

H.fs = str2num(get(obj,'String'));

if isempty(H.fs) | numel(H.fs)>1
  fprintf('Error: Please enter a single number for defining font size\n');
else
  set(gca(H.h11),'FontSize',H.fs);
  if H.y_found
    set(gca(H.h12),'FontSize',H.fs);
  end
end

end

%==========================================================================
function show_rawdata(obj, event_obj, filename)

global H

if H.boxplot | H.medianplot
  H.rawdata = get(obj, 'Value');
end

end

%==========================================================================
function show_medianplot(obj, event_obj, filename)

global H

if H.rawdata | H.boxplot
  H.medianplot = get(obj, 'Value');
end

end

%==========================================================================
function show_boxplot(obj, event_obj, filename)

global H

if H.rawdata | H.medianplot
  H.boxplot = get(obj, 'Value');
end

end

%==========================================================================
function save_image(obj, event_obj, filename)

global xY H

if ~exist('filename', 'var')
    
  filename = xY.string;
    
  [filename, newpth] = uiputfile({ ...
      '*.png' 'PNG files (*.png)'}, 'Save as', filename);
else
    [pth, nam, ext] = fileparts(filename);
    if isempty(pth), pth = cd; end
    if isempty(nam)
        [filename, newpth] = uiputfile({ ...
            '*.png' 'PNG files (*.png)'}, 'Save as', nam);
    else
        filename = fullfile(pth, nam);
        newpth = pth;
    end
end

% remove potential .png
filename = regexprep(filename,'.png','');

try
  % keep background color
  set(H.h10, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
  hh = getframe(H.h11);
  img = hh.cdata;
  col = colormap;
  saved_file = fullfile(newpth,['estimates_' filename '.png']);
  imwrite(img,col,saved_file);
  fprintf('File %s saved.\n',saved_file);
catch
  fprintf('File %s could not be saved.\n',saved_file);
end

if H.y_found
  try
    % keep background color
    set(H.h12, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
    hh = getframe(H.h12);
    img = hh.cdata;
    col = colormap;
    saved_file = fullfile(newpth,['boxplot_' filename '.png']);
    imwrite(img,col,saved_file);
    fprintf('File %s saved.\n',saved_file);
  catch
    fprintf('File %s could not be saved.\n',saved_file);
  end
end

end