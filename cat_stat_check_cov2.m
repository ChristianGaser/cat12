function varargout = cat_stat_check_cov2(job)
% cat_stat_check_cov to check covariance and image quality across sample
%
% Images have to be in the same orientation with same voxel size and
% dimension (e.g. spatially registered images), whereas surfaces have 
% to be same size (number of vertices).
%
% varargout = cat_stat_check_cov(job)
%
% job           .. spm job structure
%   .data_vol   .. volume input files
%   .data_surf  .. surface input files 
%   .c          ..
%  [.gap]       .. gap between slices (in case of volume input)
%
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: cat_stat_check_cov.m 1286 2018-03-05 11:44:37Z gaser $

%#ok<*AGROW,*ASGLU,*TRYNC,*MINV,*INUSD,*INUSL>

cat_io_cprintf('err','\nWARNING: cat_stat_check_cov2 is in an ealy development stage!\n\n')

global filename ... structure from spm_str_manip with grouped filenames
       H        ... main structure with object/button handles
       pos      ... structure with position values for GUI objects/buttons
       YpY YpYsorted mean_cov         ... (sorted) covariance matrix
       mask1d mask2d                  ... masks for files on the trash list 
       trashlist                      ... index list of subjects to remove
       ind_sorted ind_sorted_display  ... index lists 
       QM                             ... Quality measures from xml_files
       data_files org_files surf_files xml_files log_files pdf_files ...  different lists of filenames 
       data_array data_array_diff     ... slices of all subjects 
       img                            ... slice image(s) for GUI display
       FS FSi                         ... SPM fontsize 
       mesh_detected isscatter isxml sorted show_name useicons ... binary variables for GUI control 
       mn_data mx_data                ... minimum/maximum value of the mesh
       V Vo                           ... volume header structures
       Vchanged names_changed         ... modified volume header for 4D-structures
       sample                         ... sample groups of each scan
       dataprefix                     ... prefix of the input files
       inorm                          ... normalize slice intensity even in normalized data
       X X2 MD MD2; 
        % cbar img_alpha

% create default SPM windows if required        
if isempty(spm_figure('FindWin','Interactive')), spm('createintwin'); end
if isempty(spm_figure('FindWin','Graphics')),    spm_figure('Create','Graphics'); end        
        

% positions & global font size option
%  ------------------------------------------------------------------------
trashlist     = [];                             % start with empty list 
sorted        = 0;                              % show data by file order
isscatter     = 0;
show_name     = 0;                              % show filenames in boxplot rather small dots
inorm         = 1;                              % normalize slice intensity even in normalized data
ws            = spm('Winsize','Graphics');
FS            = spm('FontSizes');
FSi           = 8; 
useicons      = 1; 


popb = [0.038 0.035]; 
popm = 0.780;
posp = struct('naviui',0.745,'trashui',0.690,'checkui',0.600);
pos = struct(...
  'fig',            [10  10  1.4*ws(3) 1.2*ws(3)],... % figure
  'popup',          [10  10  200       100  ],... % popup in case of closing with non-empty trash list
  'corr',           [-0.015 0.050 0.820 0.820],... % correlation matrix
  'scat',           [0.045 0.050 0.700 0.820],... % scatter plot
  'slice',          [0.780 0.060 0.190 0.450],... % image plot
  'cbar',           [0.050 0.950 0.580 0.020],... % colorbar for correlation matrix
  ... 
  'boxplot',        [0.100 0.055 0.880 0.915],... % boxplot axis
  'fnamesbox',      [0.830 0.003 0.160 0.038],... % show filenames in boxplot 
  ...
  'cbarfix',        [0.660 0.935 0.100 0.045],... % colorbar fix/auto option
  'close',          [0.775 0.925 0.100 0.045],... % close button
  'help',           [0.875 0.925 0.100 0.045],... % help button
  'show',           [0.775 0.880 0.200 0.045],... % button to show worst cases
  'sort',           [0.772 0.825 0.210 0.045],... % list to use ordered matrix or Maha-distance 
  'boxp',           [0.772 0.785 0.210 0.045],... % list to display different variables as boxplot
  'alphabox',       [0.775 0.000 0.200 0.035],... % show filenames in boxplot 
  'sslider',        [0.780 0.030 0.193 0.040],... % slider for z-slice  
  ...
  ... == navigation unit ==
  'scSelect',       [popm+popb(1)*0 posp.naviui popb],... % select (default) 
  'scZoomReset',    [popm+popb(1)*1 posp.naviui popb],... % standard zoom
  'scZoomIn',       [popm+popb(1)*2 posp.naviui popb],... % zoom in 
  'scZoomOut',      [popm+popb(1)*3 posp.naviui popb],... % zoom out
  'scPan',          [popm+popb(1)*4 posp.naviui popb],... % pan (moving hand)
  ...
  ... == tashlist unit ==
  'newtrash',       [popm+popb(1)*0 posp.trashui popb],... % new trash list
  'trash',          [popm+popb(1)*1 posp.trashui popb],... % add data to trash list
  'detrash',        [popm+popb(1)*2 posp.trashui popb],... % remove data from trash list
  'disptrash',      [popm+popb(1)*3 posp.trashui popb],... % print trash list
  'ziptrash',       [popm+popb(1)*4 posp.trashui popb],... % pack data from trash list
  ... second row?
  'autotrash',      [popm+popb(1)*0 posp.trashui-popb(2) popb],... % button to mark data with low IQR
  'undo',           [popm+popb(1)*1 posp.trashui-popb(2) popb],... % undo last trash list operation
  'redo',           [popm+popb(1)*2 posp.trashui-popb(2) popb],... % redo last trash list operation
  'trashrow',       [popm+popb(1)*3 posp.trashui-popb(2) popb],... % add data to trash list
  'unziptrash',     [popm+popb(1)*4 posp.trashui-popb(2) popb],... % remove data from trash list
  ...
  ... == checklist unit ==
  'checkvol',       [popm+popb(1)*0 posp.checkui popb],... % open checkvol 
  'checksurf',      [popm+popb(1)*1 posp.checkui popb],... % open checksurf
  'checklog',       [popm+popb(1)*2 posp.checkui popb],... % open log-txt
  'checkxml',       [popm+popb(1)*3 posp.checkui popb],... % open xml-txt
  'checkpdf',       [popm+popb(1)*4 posp.checkui popb]);   % open pdf in external viewer
  ... 'checklow',     
     

        
if nargin == 0, error('No argument given.'); end



%% get all filenames from the data_vol/surf input
%  ------------------------------------------------------------------------
if isfield(job,'data_vol')
  datafield     = 'data_vol'; 
  datadir       = 'mri'; 
  mesh_detected = 0;
elseif isfield(job,'data_surf')
  datafield     = 'data_surf'; 
  datadir       = 'surf'; 
  mesh_detected = 1;
end
% get all input scans/surfaces
data_files = {}; 
for i = 1:numel(job.(datafield))
  data_files = [data_files;job.(datafield){i}];
end
% number of samples and scans, trash mask arrays, sample array
n_subjects = numel(data_files);
n_samples  = numel(job.(datafield));
mask1d     = true(n_subjects,1);           % trash list mask 1D matrix
mask2d     = true(n_subjects,n_subjects);  % trash list mask 2D matrix
sample     = [];
for i=1:n_samples
  sample = [sample, i*ones(1,size(job.(datafield){i},1))];
end



%% get the different files
%  ------------------------------------------------------------------------
spm_progress_bar('Init',n_samples,'Search files','subects completed')
[filenames,fparts] = spm_str_manip(data_files,'trC'); 
out_files  = data_files;
org_files  = data_files;
pdf_files  = data_files;
xml_files  = data_files; 
log_files  = data_files;
surf_files = data_files;
dataprefix = fparts.s;
for i = 1:numel(data_files)
  [pp,ff,ee] = spm_fileparts(data_files{i});
  [pp1,pp2]  = spm_fileparts(pp); 
  
  % output files
  out_files{i} = fullfile(pp,ff,ee); 
  
  % set subdirectories
  if strcmp(pp2,datadir)
    reportdir = 'report';
    surfdir   = 'surf';
  else
    reportdir = ''; 
    surfdir   = '';
  end

  % set original input files of the CAT preprocessing
  org_files{i} = fullfile(pp1,[cat_io_strrep(ff,fparts.s,'') '.nii']); 
  if ~exist(org_files{i},'file')
    org_files{i} = fullfile(pp1,[cat_io_strrep(ff,fparts.s,'') '.img']);
    if ~exist(org_files{i},'file')
      org_files{i} = ''; 
    end
  end
  
  % try to find the XML file if not given
  if isempty( char(job.data_xml) ) 
    xml_files{i} = fullfile(pp1,reportdir,...
      ['cat_' cat_io_strrep(ff,fparts.s,'') ...
        cat_io_strrep(ee,{'.nii','.img','.gii'},'.xml')]);
    if ~exist(xml_files{i},'file')
      xml_files{i} = ''; 
    end
  else
    xml_files = cellstr(job.data_xml);
  end
  
  % set report pdf
  pdf_files{i} = fullfile(pp1,reportdir,...
    ['catreport_' cat_io_strrep(ff,fparts.s,'') '.pdf']); 
  if ~exist(pdf_files{i},'file')
    pdf_files{i} = ''; 
  end
  
  % log files
  log_files{i} = fullfile(pp1,reportdir,...
    ['catlog_' cat_io_strrep(ff,fparts.s,'') '.txt']);
  if ~exist(log_files{i},'file')
    log_files{i} = ''; 
  end

  % surface files
  surf_files{i} = fullfile(pp1,surfdir,...
    ['lh.thickness.' cat_io_strrep(ff,fparts.s,'') ]);
  if ~exist(surf_files{i},'file')
    surf_files{i} = ''; 
  end

  spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');



%% load header 
%  ------------------------------------------------------------------------
V  = spm_data_hdr_read(char(data_files));
Vo = spm_data_hdr_read(char(org_files));



%% load XML data
%  ------------------------------------------------------------------------
if isempty( xml_files )
  isxml     = 0;
  QM_names  = '';
else
  isxml     = 1;
  
  if size(xml_files,1) ~= n_subjects
    error('XML-files must have the same number as sample size');
  end
  
  QM = ones(n_subjects,3);
  QM_names = char(...
    'Noise rating (NCR)',...
    'Bias Rating (ICR)',...
    'Weighted overall image quality rating (IQR)');

  spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
    % get basename for xml- and data files
    [pth, xml_name]  = fileparts(deblank(xml_files{i}));
    [pth, data_name] = fileparts(V(i).fname);
    
    % remove leading 'cat_'
    xml_name = xml_name(5:end);
    
    % check for filenames
    if isempty(strfind(data_name,xml_name))
      warning('Please check file names because of deviating subject names\n: %s vs. %s\n',...
        V(i).fname,xml_files{i});
    end
    
    xml = cat_io_xml(deblank(xml_files{i}));
    try
      QM(i,:)   = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR];
      site(i,1) = xml.qualityratings.res_RMS; 
    catch % also try to use old version
      QM(i,:)   = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
      site(i,1) = xml.QAM.res_RMS;
    end
    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');
 
  
  % added protocol depending QA parameter
  if cat_get_defaults('extopts.expertgui')>1
    [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_cleaner_intern(QM(:,3),struct('site',site,'figure',0));
    QM_names = char([cellstr(QM_names);{'Protocol IQR difference (IQRD)'}]);
    QM(:,4) = rth(:,3)   - QM(:,3);
    %QM(:,5) = rthsc(:,3) - QM(:,3);
  end
  
  
  % convert marks into rps rating
  mark2rps   = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  markd2rpsd = @(mark) ( mark*10) + isnan(mark).*mark;
  QM(:,1:3)  = mark2rps(QM(:,1:3));
  if cat_get_defaults('exptops.expertgui')>1
    QM(:,4)    = markd2rpsd(QM(:,4));
    %QM(:,5)    = markd2rpsd(QM(:,5));
  end
end



%% add constant to nuisance parameter
%  ------------------------------------------------------------------------
G = [];
if ~isempty(job.c)
  for i=1:numel(job.c)
    G = [G job.c{i}];
  end
  if size(G,1) ~= n_subjects
    G = G';
  end
  G = [ones(n_subjects,1) G];
end



%% load surface data, prepare volume data loading
%  ------------------------------------------------------------------------
if mesh_detected
  % load surface texture data
  Y = spm_data_read(V)';
  
  % optional global scaling
  if isfield(job,'gSF')
    for i=1:numel(V)
      Y(:,2) = Y(:,2)*job.gSF(i);
    end
  end
  
  Y(isnan(Y)) = 0;
  
  % rescue unscaled data min/max
  mn_data = min(Y(:));
  mx_data = max(Y(:));
  Y = Y - repmat(mean(Y,2), [1 size(Y,2)]);

  % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
  if ~isempty(G) 
    Ymean = repmat(mean(Y), [n_subjects 1]);
    Y = Y - G*(pinv(G)*Y) + Ymean;
  end

  data_array = Y';
  YpY = (Y*Y')/n_subjects;

  % calculate residual mean square of mean adjusted Y
  Y = Y - repmat(mean(Y,1), [n_subjects 1]);
  data_array_diff = Y';
  
  %MSE = sum(Y.*Y,2);
else
  % voxelsize and origin
  vx   = sqrt(sum(V(1).mat(1:3,1:3).^2));
  Orig = V(1).mat\[0 0 0 1]';

  if length(V)>1 && any(any(diff(cat(1,V.dim),1,1),1))
    error('images don''t all have same dimensions')
  end
  if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
    error('images don''t all have same orientation & voxel size')
  end
  
  % consider image aspect ratio
  pos.slice(4) = pos.slice(4) * V(1).dim(2)/V(1).dim(1);
  
  slices = 1:job.gap:V(1).dim(3);

  dimx = length(1:job.gap:V(1).dim(1));
  dimy = length(1:job.gap:V(1).dim(2));
  Y    = zeros(n_subjects, prod(dimx*dimy));
  YpY  = zeros(n_subjects);
  %MSE  = zeros(n_subjects,1);
  data_array = zeros([V(1).dim(1:2) n_subjects]);

  
  
  %-Start progress plot
  %-----------------------------------------------------------------------
  spm_progress_bar('Init',V(1).dim(3),'Check correlation','planes completed')

  for j=slices

    M  = spm_matrix([0 0 j 0 0 0 job.gap job.gap job.gap]);

    for i = 1:n_subjects
      img = spm_slice_vol(V(i),M,[dimx dimy],[1 0]);
      img(isnan(img)) = 0;
      Y(i,:) = img(:);
      if isfield(job,'gSF')
        Y(i,:) = Y(i,:)*job.gSF(i);
      end
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
    
    %MSE = MSE + sum(Y.*Y,2);

    spm_progress_bar('Set',j);  

  end

  % correct filenames for 4D data
  if strcmp(V(1).fname, V(2).fname)
    names_changed = 1;
    Vchanged      = V;
    for i=1:n_subjects
      [pth,nam,ext] = spm_fileparts(V(i).fname);
      V(i).fname    = fullfile(pth, [nam sprintf('%04d',i) ext]);
    end
  else
    names_changed = 0;
  end
  
  spm_progress_bar('Clear');
end
clear Y



%% normalize YpY and estimate mean_cov
%  ------------------------------------------------------------------------
d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd     = d*d';
YpY    = YpY./(dd+eps);
t      = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n_subjects+1:end) = sign(diag(YpY));
clear t d dd;

% extract mean correlation for each data set
mean_cov = zeros(n_subjects,1);
for i=1:n_subjects
  cov0        = YpY(i,:);     % extract row for each subject
  cov0(i)     = [];           % remove cov with its own
  mean_cov(i) = mean(cov0);
end
clear cov0;



%% output compressed filenames structure
%  ------------------------------------------------------------------------
fprintf('\n');
fname_m = [];
fname_tmp = cell(n_samples,1);
fname_s   = cell(n_samples,1);
fname_e   = cell(n_samples,1);
for i=1:n_samples
  [tmp, fname_tmp{i}] = spm_str_manip(char(V(sample == i).fname),'C');
  fname_m = [fname_m; fname_tmp{i}.m];
  fname_s{i} = fname_tmp{i}.s;
  cat_io_cprintf('n','Compressed filenames sample %d: ',i);
  cat_io_cprintf('b',sprintf('%s %s \n',spm_str_manip(tmp,'f120'),...
    repmat('.',1,3*(numel(tmp)>120))));
end
filename = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});
clear fname_e fname_m fname_s fname_tmp tmp



%% print suspecious files with cov>0.925
%  ------------------------------------------------------------------------
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.925);
[siv,si] = sort(YpY(sub2ind(size(YpY),indx,indy)),'descend');
% if more than 25% of the data this points to longitudinal data of one subject and no warning will appear
if ~isempty(indx) && (sqrt(length(indx)) < 0.25*n_subjects)
  fprintf('\nUnusual large correlation (check that subjects are not identical):\n');
  for i=si'
    % exclude diagonal
    if indx(i) ~= indy(i)
      % report file with lower mean correlation first
      if mean_cov(indx(i)) < mean_cov(indy(i))
        cat_io_cprintf('w',sprintf('  %0.4f',YpY(indx(i),indy(i)))); 
        cat_io_cprintf('n',' between ');
        cat_io_cprintf('b',filename.m{indx(i)}); cat_io_cprintf('n',' and ');
        cat_io_cprintf('b',filename.m{indy(i)}); fprintf('\n');
      else
        cat_io_cprintf('w',sprintf('  %0.4f',YpY(indy(i),indx(i)))); 
        cat_io_cprintf('n',' between ');
        cat_io_cprintf('b',filename.m{indy(i)}); cat_io_cprintf('n',' and ');
        cat_io_cprintf('b',filename.m{indx(i)}); fprintf('\n');
      end
    end
  end
end



%% sort data and estimate critical files
%  ------------------------------------------------------------------------
[mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
YpYsorted          = YpY(ind_sorted,ind_sorted);
ind_sorted_display = ind_sorted;
threshold_cov      = mean(mean_cov) - 2*std(mean_cov);
n_thresholded      = find(mean_cov_sorted < threshold_cov,1,'first');
if ~isempty(n_thresholded)
  fprintf('\nThese data have a mean correlation below 2 standard deviations. \n');
  fprintf('This does not necessarily mean that you have to exclude these data. \n');
  fprintf('However, these data have to be carefully checked:\n');
  for i=n_thresholded:n_subjects
    cat_io_cprintf('r',sprintf('  %0.4f ',mean_cov_sorted(i)));
    cat_io_cprintf('b',V(ind_sorted(i)).fname); fprintf('\n');
  end
end



%% output structure
%  ------------------------------------------------------------------------
if nargout>0
  varargout{1} = struct('table',[out_files,num2cell(mean_cov)],...
                        'covmat',YpY,...
                        'sorttable',[cellstr(V(ind_sorted).fname),num2cell(mean_cov_sorted)],...
                        'sortcovmat',YpYsorted, ...
                        'cov',mean_cov,...
                        'threshold_cov',threshold_cov);
end
clear mean_cov_sorted threshold_cov;



%% create figure
%  ------------------------------------------------------------------------

H.graphics = spm_figure('GetWin','Graphics');
H.figure   = figure(2);
clf(H.figure);
set(H.figure,'MenuBar','none','Position',pos.fig,'NumberTitle','off','Resize','off');    

if mesh_detected
  set(H.figure,'Name','CAT Check Covarance: Click in image to display surfaces');
else
  set(H.figure,'Name','CAT Check Covarance: Click in image to display slices');
end

cm = datacursormode(H.figure);
set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
%try set(cm,'NewDataCursorOnClick',false); end 

% add colorbar
H.jet = axes('Position',pos.cbar,'Parent',H.figure);
cbarimg = image((1:64)); set(cbarimg,'HitTest','off','Interruptible','off');


show_matrix(YpY, sorted);

% create two colormaps
cmap = [jet(64); gray(64)];
colormap(cmap)

%%
%set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,7),...
%  'XTickLabel',linspace(0.7,1.0,7)); 

% add button for closing all windows
H.close = uicontrol(H.figure,...
  'Units','normalized','position',pos.close,'Style','Pushbutton','callback',@closeWindows,...
  'string','Close','ToolTipString','Close windows','FontSize',FS(FSi),'ForegroundColor',[0.8 0 0]);

% add button to open the help HTML
H.help = uicontrol(H.figure,...
  'Units','normalized','position',pos.help,'Style','Pushbutton',...
  'string','Help','ToolTipString','Open help window','ForegroundColor',[0 0 0.8],'FontSize',FS(FSi),...
  'callback',['spm_help(''!Disp'',''/Users/dahnke/Documents/MATLAB/' ...
    'spm12/toolbox/cat12/html/cat_methods_QA.html'','''',H.graphics);']);

% check button ... maybe use this one as small button?
H.show = uicontrol(H.figure,...
  'Units','normalized','position',pos.show,'Style','Pushbutton','callback',@check_worst_data,...
  'string','Check most deviating data','ToolTipString','Display most deviating files','FontSize',FS(FSi));
      


% create popoup menu for SPM grafix window
if isxml

  % estimate Mahalanobis distance between mean corr. and weighted overall quality
  X  = [mean_cov, QM(:,3)]; % mean correlation and IQR
  S  = cov(X);
  mu = mean(X);
  MD = (X-repmat(mu,[length(X),1]))*inv(S)*(X-repmat(mu,[length(X),1]))'; 
  MD = diag(MD);
  
  if cat_get_defaults('extopts.expertgui')>1
    X2  = [mean_cov, QM(:,4)]; 
    S2  = cov(X2);
    mu2 = mean(X2);
    MD2 = (X2-repmat(mu2,[length(X2),1]))*inv(S2)*(X2-repmat(mu2,[length(X2),1]))';
    MD2 = diag(MD2);
  end  
  
  str  = { 'Boxplot...','Mean correlation',QM_names,'Mahalanobis distance'};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1},...
           {@show_mean_boxplot, QM(:,1), QM_names(1,:), 1},...
           {@show_mean_boxplot, QM(:,2), QM_names(2,:), 1},...
           {@show_mean_boxplot, QM(:,3), QM_names(3,:), 1},...
           {@show_mean_boxplot, MD, 'Mahalanobis distance  ', -2} };
  if cat_get_defaults('extopts.expertgui')>1
    tmp = [ tmp(1:4),...
        {{@show_mean_boxplot, QM(:,4), QM_names(4,:), 1}},...
        tmp(5:end),...
        {{@show_mean_boxplot, MD2, 'Mahalanobis distance 2  ', -2}}];
  end
else
  str  = { 'Boxplot...','Mean correlation'};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1} };
end
H.boxp = uicontrol(H.figure,...
  'Units','normalized','position',pos.boxp,'Style','PopUp','callback','spm(''PopUpCB'',gcbo)',...
  'string',str,'ToolTipString','Display boxplot','FontSize',FS(FSi),'UserData',tmp);


% create popoup menu for main check_cov window
if isxml
  str  = { 'Image...','Mean Correlation: Order by selected filenames', ...
           'Mean Correlation: Sorted by mean correlation','Mahalanobis distance'};
  tmp  = { {@show_matrix, YpY, 0},...
           {@show_matrix, YpYsorted, 1},...
           {@show_mahalanobis, X}};
  if cat_get_defaults('extopts.expertgui')>1
    str = char([cellstr(str),{'Mahalanobis distance IQRp'}]);
    tmp = [ tmp , ...
           {{@show_mahalanobis, X2}}]; 
  end
else
  str  = { 'Correlation matrix...','Order by selected filename',...
           'Sorted by mean correlation'};
  tmp  = { {@show_matrix, YpY, 0},...
           {@show_matrix, YpYsorted, 1} };
end

H.sort = uicontrol(H.figure,...
  'Units','normalized','position',pos.sort,'Style','PopUp','UserData',tmp,...
  'callback','spm(''PopUpCB'',gcbo)','string',str,'ToolTipString','Sort matrix','FontSize',FS(FSi));

H.alphabox = uicontrol(H.figure,...
  'Units','normalized','position',pos.alphabox,'Style','CheckBox','callback',@update_alpha,...
  'string','Colorize diff. to sample mean','Value',1,...
  'ToolTipString','Colorize difference to sample mean (pos=green;neg=red)',...
  'Visible','off','BackgroundColor',[0.8 0.8 0.8],'FontSize',FS(FSi));
        
H.cbarfix = uicontrol(H.figure,...
  'Units','normalized','Style','CheckBox','position',pos.cbarfix,'callback',{@checkbox_cbarfix},...
  'string','Fixed range','ToolTipString','Switch between fixed and auto-scaled colorbar',...
  'Value',1,'BackgroundColor',[0.8 0.8 0.8],'FontSize',FS(FSi));

      
% add slider only for volume data
if ~mesh_detected
  H.mm = uicontrol(H.figure,...
    'Units','normalized','position',pos.sslider,...
    ...'Min',(1 - Orig(3))*vx(3) ,'Max',(V(1).dim(3) - Orig(3))*vx(3),...
    'Min', -sum(slices<Orig(3)) * job.gap * vx(3),...
    'Max',  sum(slices>Orig(3)) * job.gap * vx(3),...
    'Style','slider','HorizontalAlignment','center',...
    'callback',@update_slices_array,...
    'ToolTipString','Select slice for display',...
    'SliderStep',[1 job.gap] / (V(1).dim(3)-1),'Visible','off');

  H.mm_txt = uicontrol(H.figure,...
    'Units','normalized','HorizontalAlignment','center',...
    'Style','text','BackgroundColor',[0.8 0.8 0.8],...
    'Position',[pos.sslider(1) pos.sslider(2)-0.005 0.2 0.02],...
    'String','Slice [mm]','Visible','off','FontSize',FS(FSi));
        
  update_slices_array;
end
   

if isxml
  
  % == navigation unit ==
  H.naviui.text = uicontrol(H.figure,...
    'Units','normalized',...
    'Style','text','BackgroundColor',[0.8 0.8 0.8],...
    'Position',[pos.scSelect(1) pos.scSelect(2)+0.035 0.2 0.02],...
    'String','Navigation options','FontSize',FS(FSi));
      
  H.naviui.select = uicontrol(H.figure,...
    'Units','normalized','position',pos.scSelect,'callback','datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Data selection');
    buttonicon(H.naviui.select,'DC',fullfile(matlabroot,'toolbox','matlab','icons','tool_data_cursor.png'));
    
  H.naviui.zoomReset = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomReset,'callback','zoom out; datacursormode(''on'')',...
    'Style','Pushbutton','enable','on','ToolTipString','Reset zoom'); 
    buttonicon(H.naviui.zoomReset,'Zo',fullfile(matlabroot,'toolbox','shared','dastudio','resources','glue','zoom_fit_view_mo.png'));

  H.naviui.zoomIn = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomIn,'callback',@scZoomIn,...
    'Style','Pushbutton','enable','on','ToolTipString','Zoom in');
    buttonicon(H.naviui.zoomIn,'Z+',fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_in.png'));
      
  H.naviui.zoomOut = uicontrol(H.figure,...
    'Units','normalized','position',pos.scZoomOut,'callback',@scZoomOut,...
    'Style','Pushbutton','enable','on','ToolTipString','Zoom out');
    buttonicon(H.naviui.zoomOut,'Z-',fullfile(matlabroot,'toolbox','matlab','icons','tool_zoom_out.png'));
     
  H.naviui.pan = uicontrol(H.figure,...
    'Units','normalized','position',pos.scPan,'Enable','off','callback','pan on',...
    'Style','Pushbutton','enable','on','ToolTipString','Hand');
    buttonicon(H.naviui.pan,'H',fullfile(matlabroot,'toolbox','matlab','icons','tool_hand.png'));
     
  
  
  % == check unit ==
  H.checkui.text = uicontrol(H.figure,...
    'Units','normalized','HorizontalAlignment','center','Style','text',...
    'BackgroundColor',[0.8 0.8 0.8],'ForegroundColor',[0.4 0.4 0.4],...
    'Position',[pos.checkvol(1) pos.checkvol(2)+0.035 0.2 0.02],...
    'String','View selected data','FontSize',FS(FSi));
  
  % add button to open one image with SPM check_reg
  H.checkui.vol = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkvol,'callback',@checkvol,...
    'string','VOL','ToolTipString','Display original volume in SPM Graphics',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  
    %buttonicon(H.trashui.disptrash,'',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
      
  % add button to open one image with SPM check_reg
  H.checkui.surf = uicontrol(H.figure,...
    'Units','normalized','position',pos.checksurf,'callback',@checksurf,...
    'string','SURF','ToolTipString','Display processed surfaces in own figure',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  
    %buttonicon(H.trashui.disptrash,'',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
      
  % add button to open one image with SPM check_reg
  H.checkui.log = uicontrol(H.figure,...
    'Units','normalized','position',pos.checklog,'callback',@checklog,...
    'string','LOG','ToolTipString','Display log-file in SPM Graphics',...
    'Style','Pushbutton','FontSize',FS(FSi),'Enable','off');  
    %buttonicon(H.trashui.disptrash,'',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  
  % add button to open one image with SPM check_reg
  H.checkui.xml = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkxml,'callback',@checkxml,...
    'string','XML','FontSize',FS(FSi),'ToolTipString','Display xml-file in SPM Graphics',...
    'Style','Pushbutton','Enable','off'); 
    %buttonicon(H.trashui.disptrash,'',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
  
  % add button to open the pdf in an external viewer 
  H.checkui.pdf = uicontrol(H.figure,...
    'Units','normalized','position',pos.checkpdf,'callback',@checkpdf,...
    'ToolTipString','Display PDF report in external viewer',...
    'Style','Pushbutton','Enable','off');
    buttonicon(H.checkui.pdf,'PDF',fullfile(spm('dir'),'toolbox','cat12','html','icons','ico_pdf.png'));
      
        
        
  % == trashlist unit ==
  H.trashui.text = uicontrol(H.figure,...
    'Units','normalized','Style','text','Position',[pos.newtrash(1) pos.newtrash(2)+0.035 0.2 0.02],...
    'BackgroundColor',[0.8 0.8 0.8],'ForegroundColor',[0.4 0.4 0.4],...
    'String','Trashlist operations','FontSize',FS(FSi));

  % add button for new garbage mask
  H.trashui.new = uicontrol(H.figure,...
    'Units','normalized','position',pos.newtrash,'callback',@newtrash,...
    'string','NT','ForegroundColor',[ 0 0 0.8],'FontSize',FS(FSi),...
    'ToolTipString','Reset trash list','Style','Pushbutton','Enable','off');  
        
  % add button to set the active image as garbage
  H.trashui.trash = uicontrol(H.figure,...
    'Units','normalized','position',pos.trash,'callback',@trash,...
    'string','T+','ForegroundColor',[0.8 0 0],'FontSize',FS(FSi),...
    'ToolTipString','Trash selected subject','Style','Pushbutton','Enable','off');
  %/Applications/MATLAB_R2016a.app/toolbox/shared/comparisons/+comparisons/+internal/+text/deleted.png
        
  % add button to remove the active image from garbage
  H.trashui.detrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.detrash,'callback',@detrash,...
    'string','T-','ForegroundColor',[0 0.8 0],'FontSize',FS(FSi),...
    'ToolTipString','Trash selected subject','Style','Pushbutton','Enable','off');
  %RestoreOrphanedSignals_16
  
  % add button to remove the active image from garbage
  H.trashui.ziptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.ziptrash,'callback',@ziptrash,...
    'string','ZT','FontSize',FS(FSi),'ToolTipString','ZIP selected subject',...
    'Style','Pushbutton','Enable','off');
        
  % add button for mask below threshold as garbage
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.disptrash,'callback',@disptrash,...
    'string','PT','FontSize',FS(FSi),'ToolTipString','Display trash',...
    'Style','Pushbutton','Enable','off'); 
  
  % == second row ==
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.autotrash,'callback',@disptrash,...
    'string','AT','FontSize',FS(FSi),'ToolTipString','Automatic IQR thresholding',...
    'Style','Pushbutton','Enable','off'); 
   
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.undo,'callback',@disptrash,...
    'Style','Pushbutton','Enable','off','ToolTipString','Undo last trashlist operation'); 
    buttonicon(H.trashui.disptrash,'UD',fullfile(spm('dir'),'toolbox','cat12','html','icons','restore_24.png'));
  
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.redo,'callback',@disptrash,...
    'Style','Pushbutton','Enable','off','ToolTipString','Redo last trashlist operation'); 
    buttonicon(H.trashui.disptrash,'RD',fullfile(spm('dir'),'toolbox','cat12','html','icons','Refresh_16.png'));
  
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.trashrow,'callback',@disptrash,...
    'string','TR','FontSize',FS(FSi),'ToolTipString','Tresh row',...
    'Style','Pushbutton','Enable','off'); 
  
  H.trashui.disptrash = uicontrol(H.figure,...
    'Units','normalized','position',pos.trashrow,'callback',@disptrash,...
    'string','TR','FontSize',FS(FSi),'ToolTipString','Redo last trashlist operation',...
    'Style','Pushbutton','Enable','off'); 
 
end




show_mean_boxplot(mean_cov,'Mean correlation  ',1);


% check for replicates
for i=1:n_subjects
  for j=1:n_subjects
    if (i>j) && (mean_cov(i) == mean_cov(j))
      try
        nami = deblank(V(i).fname);
        namj = deblank(V(j).fname);
        if strcmp(nami(end-1:end),',1')
          nami = nami(1:end-2);
        end 
        if strcmp(namj(end-1:end),',1')
          namj = namj(1:end-2);
        end 
        s = unix(['diff ' nami ' ' namj]);
        if (s==0), fprintf(['\nWarning: ' nami ' and ' namj ' are same files?\n']); end
      end
    end
  end
end

%-End
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function buttonicon(h,str,Picon) 
%-----------------------------------------------------------------------
% Function to print an image file Picon on the button with handle h. Use
% the string str if the global variable useicons<1 or other errors.
%-----------------------------------------------------------------------
  global FS FSi useicons

  usethisicon = useicons;
  if ~exist(Picon,'file')
    warning('Button Icon "%s" does not exist!',Picon); 
  end
  if ~exist('findjobj','file')
    warning('JAVA Function for  "%s" does not exist!',Picon);
    useicons = 0; 
  end
  if usethisicon
    try
      jButton = findjobj(h);
      jButton.setIcon(javax.swing.ImageIcon(Picon));
      jButton.setHorizontalTextPosition(javax.swing.SwingConstants.LEFT);
      jButton.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
    catch 
      usethisicon = 0;
    end
  end
  if usethisicon<1
    set(h,'string',str,'FontSize',FS(FSi));
  end
return

%-----------------------------------------------------------------------
function trash(obj, event_obj) 
%-----------------------------------------------------------------------
% Put a record on the trash list and mark them with a red cross in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist mask1d mask2d
  if isfield(pos,'x') && all( trashlist~=pos.x ) 
    trashlist = unique([trashlist pos.x]);
    set(H.trashui.trash  ,'Enable','off');
    set(H.trashui.detrash,'Enable','on' );
    set(pos.tar_mouse,'sizedatasource',get(pos.tar_mouse,'marker'),...
      'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
      'ZDataSource','trash','MarkerFaceAlpha',0);
   
    mask1d(pos.x)    = 0;
    mask2d(pos.x,:)  = 0;
    mask2d(:,pos.x)  = 0;
  end
return

%-----------------------------------------------------------------------
function detrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove a record from trash list and restore the old look like in the 
% Mahalanobis plot.
%-----------------------------------------------------------------------
  global H pos trashlist mask1d mask2d

  if isfield(pos,'x') && any( trashlist==pos.x ) 
    trashlist = setdiff(trashlist,pos.x);
    set(H.trashui.trash  ,'Enable','on' );
    set(H.trashui.detrash,'Enable','off');
    set(pos.tar_mouse,'marker',get(pos.tar_mouse,'sizedatasource'),... 
      'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
    
    mask1d(pos.x)    = 1;
    mask2d(pos.x,:)  = mask1d';
    mask2d(:,pos.x)  = mask1d;
  end
return

%-----------------------------------------------------------------------
function newtrash(obj, event_obj)
%-----------------------------------------------------------------------
% Create an empty trash list. 
%-----------------------------------------------------------------------
  global trashlist H 
  trashlist = [];
  sc = findobj('ZDataSource','trash');
  for sci=1:numel(sc)
   set(sc(sci),'marker',get(sc(sci),'sizedatasource'),... 
        'ZDataSource','','MarkerEdgeColor','flat','MarkerFaceAlpha',1/3);
  end
  set(H.trashui.disptrash,'Enable','off');
  set(H.trashui.ziptrash,'Enable','off');
  set(H.trashui.new,'Enable','off');
  set(H.trashui.detrash,'Enable','off');
  set(H.trashui.trash,'Enable','on');
return

%-----------------------------------------------------------------------
function disptrash(obj, event_obj)
%-----------------------------------------------------------------------
% List all records of the trash list in the command window.
%-----------------------------------------------------------------------
  global trashlist data_files 

  fprintf('Trashlist:\n')
  for fi=1:numel(trashlist)
    fprintf('  %s\n',data_files{trashlist(fi)});
  end
return

%-----------------------------------------------------------------------
function ziptrash(obj, event_obj)
%-----------------------------------------------------------------------
% Remove records and related files from the file system by zipping or
% storing in a separate directory (NOT READY). 
%-----------------------------------------------------------------------
  global X QM trashlist org_files 
  
  trashtime = datestr(clock,'yyyymmdd-HHMMSS');
  fprintf('Zip preprocessed data (trashtime = %s) of:\n',trashtime);
  for fi=numel(trashlist):-1:1
    fprintf('  Zip %s\n',org_files{trashlist(fi)});
    
    %% find preprocessed images 
    [pp,ff,ee] = spm_fileparts(org_files{trashlist(fi)});
    sim_files  = cat_vol_findfiles(pp,['*' ff '*.nii'],struct('maxdepth',1)); 
    sim_files  = [fullfile(pp,[ff ee]);setdiff(sim_files,fullfile(pp,[ff ee]))]; 
    for si=1:numel(sim_files); 
      [pps,ffs]    = spm_fileparts(sim_files{si});
      pp_files{si} = cat_vol_findfiles(pps,['*' ffs '*']);
      pp_files{si} = setdiff(pp_files{si},sim_files);
      for fsi = numel(pp_files{si}):-1:1
        ppfs = spm_str_manip(pp_files{si}{fsi},'hht');
        switch ppfs
          case 'err', pp_files{si}(fsi) = [];
        end
      end
      pp_files{si} = setdiff(pp_files{si},sim_files{si}); 
      if si==1
        pp_filescor = pp_files{si};
      else
        pp_filescor  = setdiff(pp_filescor,pp_files{si}); 
      end
    end
    
    
    %% zip the list and remove the files
    %  we need to go into the directory and use the short filenames 
    %  to opbtain save the relative path in the zip file!
    %odir = dir; cd(pp); 
    %pp_filescor1 = cat_io_strrep(pp_filescor,[pp filesep],''); 
    ffzip = fullfile(pp,sprintf('%s_cor%2.2f_IQR%2.2f_trashed%s',...
      ff,X(trashlist(fi),1),QM(trashlist(fi),3),trashtime));
    zip(ffzip,pp_filescor1,pp); 
    
    %for fii=1:numel(pp_filescor1), delete(pp_filescor1{fii}); end
    cd(odir); 
    
    
  end
  % update all variables :-/ 
  % or recreate job?

return

%-----------------------------------------------------------------------
function scZoomIn(obj, event_obj)
%-----------------------------------------------------------------------
  global H 
  hz = zoom(H.ax);
  set(hz,'enable','on','direction','in');
return

%-----------------------------------------------------------------------
function scZoomOut(obj, event_obj)
%-----------------------------------------------------------------------
  global H 
  hz = zoom(H.ax);
  set(hz,'enable','on','direction','out');
return 

%-----------------------------------------------------------------------
function closeWindows(obj, event_obj)
%-----------------------------------------------------------------------
% Close all windows. 
% Remove variables (NOT DONE).
%-----------------------------------------------------------------------
  global trashlist pos
  %%
  posx = get(get(event_obj.Source,'Parent'),'Position');
  pos.popup(1:2) = [posx(1) + posx(3)*0.8, posx(1) + posx(4)*0.9];  
  if ~isempty(trashlist)
    d = dialog('Position',pos.popup,'Name','Close Check Sample');
    uicontrol('Parent',d,'Style','text','Position',[20 60 160 20],...
       'String','Trashlist not empty!');
    uicontrol('Parent',d,'TooltipString','Sopt closing',...
       'Position',[25 20 70 25],'String','Cancel','Callback','delete(gcf)');    
    uicontrol('Parent',d,'TooltipString','Close windows without last changes',...
       'Position',[100 20 70 25],'String','Close','ForegroundColor',[0.8 0 0],...
       'Callback','for i=2:26, try close(i); end; end; delete(gcf)');   
  else
    for i=2:26, try close(i); end; end; 
    spm_clf
  end
return

%-----------------------------------------------------------------------
function checkpdf(obj, event_obj)
%-----------------------------------------------------------------------
% Open PDF report of selected subjects. 
% This is only possible for using an extern viewer. 
% Hence, it would be useful to save a JPG or HTML file in cat_main.
%-----------------------------------------------------------------------
  global pos pdf_files isscatter
  open(pdf_files{pos.x});
  if ~isscatter
    open(pdf_files{pos.y});
  end
return

%-----------------------------------------------------------------------
function checksurf(obj, event_obj)
%-----------------------------------------------------------------------
% Open surface files of selected subjects. 
% This is very slow and some feedback would be useful. 
%-----------------------------------------------------------------------
  global pos surf_files isscatter
  
  if isscatter
    cat_surf_display(struct('data',surf_files{pos.x},'multisurf',1));
  else
    cat_surf_display(struct('data',char(surf_files([pos.x,pos.y])),'multisurf',1));
  end
  
  % give some feedback
return

%-----------------------------------------------------------------------
function checkvol(obj, event_obj)
%-----------------------------------------------------------------------
% Load the original image of selected files in SPM graphics window.
% Some further information or legend would be helpful.
%-----------------------------------------------------------------------
  global H pos org_files isscatter
  
  spm_figure('Clear',H.graphics);
  
  if isscatter
    spm_check_registration(org_files{pos.x});
  else
    spm_check_registration(char(org_files([pos.x,pos.y])));
  end
  spm_orthviews('MaxBB')
return

%-----------------------------------------------------------------------
function checkxml(obj, event_obj)
%-----------------------------------------------------------------------
% Load XML report in SPM graphics window (see also checklog).
% This is just the first fast version of this function. 
% Finally, I want to use the xml structure from the file to print some
% specific informations similar to the CAT report in cat_main. 
%-----------------------------------------------------------------------
  global H pos xml_files isscatter
  
  spm_figure('Clear',H.graphics); axis off;
  
  if isscatter
    textbox = [0 0 1 1];
    files   = xml_files(pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = xml_files([pos.y,pos.x]); 
  end
  
  for fi=1:numel(files);
    fid = fopen(files{fi});
    ph  = uipanel(H.graphics,'Units','normalized','position',textbox(fi,:), ...
      'BorderWidth',0,'title',spm_str_manip(files{fi},'k100'),'ForegroundColor',[0 0 0.8]);
    lbh = uicontrol(ph,'style','listbox','Units','normalized',...
      'fontname','Fixedwidth','position',[ 0 0 1 1 ],'FontSize',9);
    indic = 1;
    while 1
     tline = fgetl(fid);
     if ~ischar(tline), 
       break
     end
     strings{indic}=tline; 
     indic = indic + 1;
    end
    fclose(fid);
    set(lbh,'string',strings);
    set(lbh,'Value',1);
    set(lbh,'Selected','on');
  end
return

%-----------------------------------------------------------------------
function checklog(obj, event_obj)
%-----------------------------------------------------------------------
% Load the log-file from cat_main of the selected subjects into the SPM
% graphics window.
%-----------------------------------------------------------------------
  global H pos log_files isscatter
  
  spm_figure('Clear',H.graphics); axis off;
  
  if isscatter
    textbox = [0 0 1 1];
    files   = log_files(pos.x); 
  else
    textbox = [0 0.5 1 0.5; 0 0 1 0.5];
    files   = log_files([pos.x,pos.y]); 
  end
  
  for fi=1:numel(files); 
    fid = fopen(files{fi});
    ph  = uipanel(H.graphics,'Units','normalized','position',textbox(fi,:), ...
      'BorderWidth',0,'title',spm_str_manip(files{fi},'k100'),'ForegroundColor',[0 0 0.8]);
    lbh = uicontrol(ph,'style','listbox','Units','normalized',...
      'fontname','Fixedwidth','position',[ 0 0 1 1 ],'FontSize',9);
    indic = 1;
    while 1
     tline = fgetl(fid);
     if ~ischar(tline), 
       break
     end
     strings{indic}=tline; 
     indic = indic + 1;
    end
    fclose(fid);
    set(lbh,'string',strings);
    set(lbh,'Value',1);
    set(lbh,'Selected','on');
  end
return

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj)
%-----------------------------------------------------------------------
% Old check worst function. The spm_input could be replaced by an popup 
% window. A specification of the data range would rather than the x worst 
% images would be useful. 
%-----------------------------------------------------------------------
global V ind_sorted_display mesh_detected mn_data mx_data H

if isempty(spm_figure('FindWin','Interactive')), spm('createintwin'); end

n = length(V);
number = min([n 24]);
number = spm_input('How many files ?',1,'e',number);
number = min([number 24]);
number = min([number length(V)]);
  
list = char(V(ind_sorted_display(n:-1:1)).fname);
list2 = list(1:number,:);

if mesh_detected
  % display single meshes and correct colorscale of colorbar
  for i=1:number
    h = cat_surf_render('Disp',deblank(list2(i,:)));
    
    % shift each figure slightly
    if i==1
        pos = get(h.figure,'Position');
    else
        pos = pos - [20 20 0 0];
    end
    
    % remove menubar and toolbar, use filename as title
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(list2(i,:),'short50'),...
         'NumberTitle','off','Position',pos);
    cat_surf_render('ColourMap',h.axis,jet);
    cat_surf_render('ColourBar',h.axis,'on');
    cat_surf_render('CLim',h,[mn_data mx_data]);
  end
else
  spm_check_registration(list2);
  spm_orthviews('Resolution',0.2);
  set(H.boxp,'Visible','on');
end
return

%-----------------------------------------------------------------------
function checkbox_names(obj, event_obj)
%-----------------------------------------------------------------------
  global H show_name
  show_name = get(H.chbox,'Value');
  show_mean_boxplot;
return
      
%-----------------------------------------------------------------------
function checkbox_cbarfix(obj, event_obj)
%-----------------------------------------------------------------------
global H show_name cbar

  show_name  = get(H.cbarfix,'Value');
  cbar.fixed = 1; 
  
return  
%-----------------------------------------------------------------------
function show_mahalanobis(X)
%-----------------------------------------------------------------------
global H FS pos isscatter ind_sorted_display MD sample trashlist YpY

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);
H.ax = axes('Position',pos.scat,'Parent',H.figure);
set(H.ax,'Color',[0.85 0.85 0.85]);
box on;

% get very similar scans
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.925);

groups        = unique(sample);
symbols       = repmat('.',1:numel(groups));  % default symbol
symbols(1:11) = 'o+^v<>ph*sd';                % need x for unset

% plot first object for the legend
hold on; 
for gi=1:numel(groups)
  txt{gi} = sprintf('sample %d \n',gi);
  Xt = X(sample==groups(gi),:); 
  scface{gi} = scatter(Xt(1,1),Xt(1,2),30,[0 0 0],symbols(gi),'Linewidth',2);
	if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
    set(scface{gi},'MarkerFaceColor','flat','markerFaceAlpha',1/3);
  end
end

txt{end+1} = 'trashlist';
scface{gi+1} = scatter(Xt(1,1),Xt(1,2),30,[1 0 0],'x','Linewidth',2,'Visible','off');
if numel(indx)/size(YpY,1)<0.5
  txt{end+1} = 'highly corr. scans'; 
  scface{gi+2} = plot([X(indx(1),1);X(indy(1),1)],[X(indx(1),2);X(indy(1),2)],'Color',[0 0 0],'Linewidth',2);
end
hl = legend(txt,'location','southwest');
set(get(hl,'title'),'string','Legend');
for gi=1:numel(groups)+2
  set(scface{gi},'visible','off');
end

if isfield(pos,'tar_mouse')
  delete(pos.tar_mouse); 
  pos = rmfield(pos,'tar_mouse');
  set(H.checkui.vol,'Enable','off'); 
  set(H.checkui.pdf,'Enable','off');
end
if isfield(pos,'x'), pos = rmfield(pos,'x'); end
set(H.alphabox,'Enable','off');
set([H.mm,H.mm_txt],'Visible','off');

%unit = struct2cell(H.trashui); set([unit{:}],'Visible','on');
%unit = struct2cell(H.checkui); set([unit{:}],'Visible','on');  


cmap = [jet(64); gray(64)];

S = cov(X);
mu = mean(X);
MD = (X-repmat(mu,[length(X),1]))*inv(S)*(X-repmat(mu,[length(X),1]))';
MD = diag(MD);

% because we use a splitted colormap we have to set the color
% values explicitely
MD2 = 63*MD/max(MD);
C = zeros(length(MD),3);
for i=1:length(MD)
  C(i,:) = cmap(round(MD2(i))+1,:);
end
hold on

% plot lines between similar objects
if numel(indx)/size(YpY,1)<0.5
  for i=1:numel(indx)
    plot([X(indx(i),1);X(indy(i),1)],[X(indx(i),2);X(indy(i),2)],...
      '-','Color',repmat(0.9 - 0.9*((YpY_tmp(indx(i),indy(i))-0.925)/0.0725),1,3),...
      'LineWidth',2,'HitTest','off','Interruptible','off');
  end
end

I   = 1:size(X,1);
for gi=1:numel(groups)
  It = I(sample==groups(gi)); 
  Xt = X(sample==groups(gi),:); 
  Ct = C(sample==groups(gi),:);
  sc = cell(size(Xt,1),1);
  for sci=1:size(Xt,1)
    sc{sci} = scatter( ...
      Xt(sci,1), ...
      Xt(sci,2), ...
      30,...
      Ct(sci,:),...
      symbols(gi), ...
      'Linewidth',2);
    if ~isempty( strfind('osd^v<>ph', symbols(gi) ) )
      set(sc{sci},'MarkerFaceColor','flat','markerFaceAlpha',1/3);
    end
    if any(It(sci)==trashlist)
      set(sc{sci},'sizedatasource',get(sc{sci},'marker'),...
      'marker','x','SizeData',40,'MarkerEdgeColor',[1 0 0.5],...
      'ZDataSource','trash','MarkerFaceAlpha',0);
    end
  end
end


xlabel('<----- Worst ---      Mean correlation       --- Best ------>  ','FontSize',FS(8),'FontWeight','Bold');
ylabel('<----- Worst ---      Weighted overall image quality rating     --- Best ------>  ','FontSize',FS(8),'FontWeight','Bold');
title('<--- Best -- Mahalanobis distance (Color) -- Worst ---->  ','FontSize',FS(10),'FontWeight','Bold');

% add colorbar
H.cbar = axes('Position',pos.cbar,'Parent',H.figure);
cbarimg = image((1:64)); set(cbarimg,'HitTest','off','Interruptible','off');

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,7),'TickLength',[0 0],...
  'HitTest','off','Interruptible','off','XTickLabel',linspace(0.7,1.0,7)); %round(100*linspace(min(MD),max(MD),5))/100);

% update index of worst files
[tmp, ind_sorted_display] = sort(MD,'ascend');

colormap(cmap)

isscatter = 1;
zoom reset
return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
global H FS FSi pos sorted isscatter mesh_detected mask1d mask2d ind_sorted

try
  if isfield(H,'alphabox'), set(H.alphabox,'Visible','off'); end
  unit = struct2cell(H.trashui); set([unit{:}],'Enable','off');
  unit = struct2cell(H.checkui); set([unit{:}],'Enable','off');  
end

if isfield(pos,'tar_mouse')
  delete(pos.tar_mouse); 
  pos = rmfield(pos,'tar_mouse'); 
end
if isfield(pos,'x'), pos = rmfield(pos,'x'); end

% get sorting order
sorted = order;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',pos.corr,'Parent',H.figure);

if order 
  mask1dt = mask1d(ind_sorted);
  mask2dt = mask2d(ind_sorted,:); 
  mask2dt = mask2dt(:,ind_sorted); 
else
  mask1dt = mask1d;
  mask2dt = mask2d;
end
data = reshape(data(mask2dt),sum(mask1dt),sum(mask1dt)); 

% scale data to 0..1
mn = 0.7; %min(data(:));
mx = 1.0; %max(data(:));
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
  title('Sorted Sample Correlation Matrix  ','FontSize',FS(10),'FontWeight','Bold');
else
  xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(8),'FontWeight','Bold');
  ylabel('<----- Last ---      File Order      --- First ------>  ','FontSize',FS(8),'FontWeight','Bold');
  title('Sample Correlation Matrix  ','FontSize',FS(FSi+2),'FontWeight','Bold');
end

H.cbar = axes('Position',pos.cbar,'Parent',H.figure);
cbarimg = image((1:64)); set(cbarimg,'HitTest','off','Interruptible','off');

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,7),'TickLength',[0 0],...
  'XTickLabel',linspace(0.7,1.0,7)); %round(100*linspace(min(YpYt(:)),max(YpYt(:)),5))/100);

cmap = [jet(64); gray(64)];
colormap(cmap)

isscatter = 0;

% remove sliders and text
try
  if ~mesh_detected
    set([H.mm,H.mm_txt],'Visible','off'); 
  end
  set([H.alpha,H.alpha_txt],'Visible','off');
end
zoom reset
return

%-----------------------------------------------------------------------
function show_mean_boxplot(data_boxp, name_boxp, quality_order)
%-----------------------------------------------------------------------
global H pos filename FS FSi sample ind_sorted_display bp mask1d show_name

if nargin == 0
  data_boxp     = bp.data;
  name_boxp     = bp.name;
  quality_order = bp.order;
end

spm_figure('Clear',H.graphics); figure(H.graphics)  
H.boxplot = axes('Position',pos.boxplot,'Parent',H.graphics);
set(H.graphics,'Renderer','OpenGL','color',[0.95 0.95 0.95]);

%if isfield(H,'chbox'), nval = ~get(H.chbox,'value'), else nval = 1; end
H.chbox = uicontrol(H.graphics,...
  'string','Show filenames','Units','normalized',...
  'position',pos.fnamesbox,'callback',@checkbox_names,...
  'Style','CheckBox','HorizontalAlignment','center',...
  'ToolTipString','Show filenames in boxplot','value',show_name,...
  'Interruptible','on','Visible','on','FontSize',FS(FSi));

n_samples = max(sample);

xpos = cell(1,n_samples);
data = cell(1,n_samples);

allow_violin = 2;

%% create filenames
hold on
for i=1:n_samples
  indtype   = { mask1d' 'k.' [0 0 0]; ~mask1d' 'rx' [1 0 0]}; 
  gnames{i} = sprintf('S%d',i);
  for ii=1:size(indtype,1)
    ind  = find(sample == i & indtype{ii,1});
    if numel(ind>0)
      data{i} = data_boxp(ind);

      if length(ind) < 8
        allow_violin = 0;
      end

      if n_samples == 1
        xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
      else
        xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
      end

      for j=1:length(ind)
        if get(H.chbox,'value')
          H.fnames{j,i} = text(xpos{i}(j),data{i}(j),filename.m{ind(j)},'Color',indtype{ii,3},...
            'FontSize',FS(FSi-1),'HorizontalAlignment','center');
        else
          H.fnames{j,i} = plot(xpos{i}(j),data{i}(j),indtype{ii,2});
        end
      end
    end
  end 
end

%% create boxplot
opt = struct('groupnum',0,'ygrid',1,'box',1,'violin',allow_violin,'median',2,...
             'groupcolor',jet(n_samples),'names',{gnames},'xlim',[-.25 n_samples+1.25]); 
if max(data_boxp) > min(data_boxp)
  ylim_add = 0.075;
  yamp = max(data_boxp) - min(data_boxp);
  ylim_min = min(data_boxp) - ylim_add*yamp;
  ylim_max = max(data_boxp) + ylim_add*yamp; 
  opt.ylim = [ylim_min ylim_max];
end    
cat_plot_boxplot(data,opt); box on;


%% add colored labels and title
if n_samples > 1
  [tmp,  tmp2] = spm_str_manip(char(filename.s),'C');
  title_str = sprintf('%s ',strrep(tmp,tmp2.s,''));
  fprintf('\nCommon filename: %s\n',tmp);
else
  title_str = sprintf('Common filename: %s*',spm_file(char(filename.s),'short25'));
end
title({['Boxplot: ' name_boxp],title_str},'FontSize',FS(FSi+1),'FontWeight','Bold');
xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(10),...
    'FontWeight','Bold');

xpos = -0.35 - n_samples*0.1;

if quality_order == -2 
  % reverse order to have the good things allways on the top
  set(gca, 'YDir','reverse');
  quality_order = 1; 
  t = ylim_min; ylim_min = ylim_max; ylim_max = t; 
end
if (length(data_boxp) > 2)
  if quality_order > 0 
    text(xpos, ylim_min,'<----- Low rating (poor quality)  ','Color','red','Rotation',...
        90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
    text(xpos, ylim_max,'High rating (good quality) ------>  ','Color',[0 0.8 0],'Rotation',...
        90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
  else
    text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
        90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
    text(xpos, ylim_min,'<----- High rating (good quality)  ','Color',[0 0.8 0],'Rotation',...
        90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
  end
  text(xpos, (ylim_max+ylim_min)/2,sprintf('%s',name_boxp),'Color','black','Rotation',...
        90,'HorizontalAlignment','center','FontSize',FS(9),'FontWeight','Bold')
end

hold off

% estimate sorted index new for displaying worst files
if quality_order > 0
  [tmp, ind_sorted_display] = sort(data_boxp,'descend');
else
  [tmp, ind_sorted_display] = sort(data_boxp,'ascend');
  
end

bp = struct('data',data_boxp,'name',name_boxp,'order',quality_order);

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
global H mesh_detected img img_alpha

alphaval = get(H.alphabox,'Value');

% display image with 2nd colorbar (gray)
image(65 + img);
if ~mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'HitTest','off','Interruptible','off');

% prepare alpha overlays for red and green colors
if alphaval > 0
  % get 2%/98% ranges of difference image
  range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

  hold on
  alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
  alpha_r = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
  hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hold off
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global V Vchanged data_array data_array_diff H pos dataprefix inorm ...
  sorted ind_sorted isscatter names_changed img_alpha img

alphaval = get(H.alphabox,'Value');

if isfield(H,'mm')
  slice_mm = get(H.mm,'Value');
else
  slice_mm = 0;
end

if names_changed
  P = Vchanged;
else
  P = V;
end

vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
Orig = P(1).mat\[0 0 0 1]';
sl   = round(slice_mm/vx(3)+Orig(3));

% if slice is outside of image use middle slice
if (sl>P(1).dim(3)) || (sl<1)
  sl = round(P(1).dim(3)/2);
end
set(H.mm,'Value',(sl-Orig(3))*vx(3));

M  = spm_matrix([0 0 sl]);
data_array_diff = data_array;

img = spm_slice_vol(P(round(length(V))),M,P(1).dim(1:2),[1 0]);
imgscale = mean(img(img ~= 0)); % scale image according to mean
  
for i = 1:length(V)
  img = spm_slice_vol(P(i),M,P(1).dim(1:2),[1 0]);
  img(isnan(img)) = 0;
  
  % rescue unscaled data
  data_array_diff(:,:,i) = img;

  % scale image according to mean
  if inorm==0 && ~isempty(strfind(dataprefix,'wp'))
    data_array(:,:,i) = img/0.4;
  else
    data_array(:,:,i) = img/imgscale; 
  end
end

% calculate individual difference to mean image
for i=1:size(data_array_diff,3)
  data_array_diff(:,:,i) = data_array_diff(:,:,i) - mean(data_array_diff,3);
end

% enhance contrast and scale image to 0..64
mn = min(data_array(:));
mx = max(data_array(:));
data_array = 64*((data_array - mn)/(mx-mn));

if sorted
  if isfield(pos,'x')
    x = ind_sorted(pos.x);
    if ~isscatter
      y = ind_sorted(pos.y);
    end
  end
else
  if isfield(pos,'x')
    x = pos.x;
    if ~isscatter
      y = pos.y;
    end
  end
end

% check whether mouse position is defined
if isfield(pos,'x')
  if isscatter
    img       = data_array(:,:,x)';
    img_alpha = data_array_diff(:,:,x)';
  else
    img       = [data_array(:,:,y) data_array(:,:,x)]';
    img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
  end
  
  % correct orientation
  img = rot90(img,2);
  img_alpha = rot90(img_alpha,2);
  
  % use gray scale colormap for values > 64
  axes('Position',pos.slice);
  image(65 + img);
  axis image
  set(gca,'XTickLabel','','YTickLabel','','HitTest','off','Interruptible','off');
  
  % prepare alpha overlays for red and green colors
  if alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

    hold on
    alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
    alpha_r = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
    hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    axis image
    hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),...
      'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
    axis image
    hold off
  end
  
  %{
  if isscatter
    txt = {sprintf('%s',spm_file(filename.m{x},'short25')),[],...
      ['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  else
    txt = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Top: ',...
      spm_file(filename.m{x},'short25')],...
      ['Bottom: ',spm_file(filename.m{y},'short25')],[],...
      ['Displayed slice: ', num2str(round(get(H.mm,'Value'))),' mm']};
  end
  set(H.text,'String',txt,'FontSize',FS(FSi));
  %}
  set(H.mm_txt,'String',sprintf('%0.1f mm',get(H.mm,'Value')));
end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj) 
%-----------------------------------------------------------------------
global trashlist filename sample H X YpY ...
  data_array data_array_diff pos mesh_detected ind_sorted sorted ...
  isscatter img img_alpha mask1d ...
  org_files pdf_files log_files surf_files xml_files

alphaval = get(H.alphabox,'Value');

pos_mouse     = get(event_obj, 'Position');
pos.tar_mouse = get(event_obj, 'Target');
      
if isscatter
  pos.x = find(X(:,1) == pos_mouse(1));
  if isempty(pos.x)
    pos.x = find(X(:,2) == pos_mouse(2));
  end

  % text info for data cursor window
  txt = {sprintf('S%d:%s',sample(pos.x),filename.m{pos.x})};

  
  %{
  % text info for textbox
  txt2 = {sprintf('%s',spm_file(filename.m{pos.x},'short25')),[],...
    'Difference to Sample Mean (red: - green: +)'};
  set(H.text,'String',txt2,'FontSize',FS(FSi));
  %}
  
  axes('Position',pos.slice .* [1 1 1 0.5]);

  x = pos.x;
   
else % covariance matrix
  
  % check for valid mouse position
  if pos_mouse(1) > pos_mouse(2) || pos_mouse(1)>length(sample) || pos_mouse(2)>length(sample)
    txt = {''};
    return
  end
  
  % save position of mouse
  pos.x = find(cumsum(mask1d)==pos_mouse(1),1,'first');
  pos.y = find(cumsum(mask1d)==pos_mouse(2),1,'first');

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
  if mesh_detected
    txt = {
      sprintf('Correlation: %3.3f',YpY(x,y)),...
      sprintf('Left:  S%d:%s',sample(x),filename.m{x}),...
      spirntf('Right: S%d:%s',sample(y),filename.m{y})};
  else
    txt = {
      sprintf('Correlation:  %3.3f',YpY(x,y)), ...
      sprintf('Column (Top): S%d:%s',sample(x),filename.m{x}), ...
      sprintf('Row (Bottom): S%d:%s',sample(y),filename.m{y})};
  end
  
  axes('Position',pos.slice);
end

% == check unit ==
onoff = {'on','off'};
if isscatter
  set(H.checkui.vol ,'Enable',onoff{ isempty(org_files{pos.x})+1  });
  set(H.checkui.surf,'Enable',onoff{ isempty(surf_files{pos.x})+1 });
  set(H.checkui.xml ,'Enable',onoff{ isempty(xml_files{pos.x})+1  });
  set(H.checkui.log ,'Enable',onoff{ isempty(log_files{pos.x})+1  });
  set(H.checkui.pdf ,'Enable',onoff{ isempty(pdf_files{pos.x})+1  });
else
  set(H.checkui.vol ,'Enable',onoff{ (isempty(org_files{pos.x}) | isempty(org_files{pos.y}))  + 1 });
  set(H.checkui.surf,'Enable',onoff{ (isempty(org_files{pos.x}) | isempty(surf_files{pos.y})) + 1 });
  set(H.checkui.xml ,'Enable',onoff{ (isempty(xml_files{pos.x}) | isempty(xml_files{pos.y}))  + 1 });
  set(H.checkui.log ,'Enable',onoff{ (isempty(log_files{pos.x}) | isempty(log_files{pos.y}))  + 1 });
  set(H.checkui.pdf ,'Enable',onoff{ (isempty(pdf_files{pos.x}) | isempty(pdf_files{pos.y}))  + 1 });
end
% == trash list unit ==
if ~isempty(pos.x) 
  if all( trashlist~=pos.x ) 
    set(H.trashui.trash  ,'Enable','on ');
    set(H.trashui.detrash,'Enable','off');
  else
    set(H.trashui.trash  ,'Enable','off');
    set(H.trashui.detrash,'Enable','on' );
  end
else
  set([H.trashui.trash,H.trashui.detrash],'Enable','off');
end
if ~isempty(trashlist)
  set([H.trashui.new,H.trashui.disptrash,H.trashui.ziptrash],'Enable','on');
else
  set([H.trashui.new,H.trashui.disptrash,H.trashui.ziptrash],'Enable','off');
end  
set(H.alphabox,'Visible','on');
  



if mesh_detected 
  % use indexed 2D-sheet to display surface data as image
  % check surface size to use indexed 2D map
  if (length(data_array(:,x)) == 163842)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
    if isscatter
      img = reshape(data_array(ind,x),[256,128]);
    else
      img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
    end
    img = circshift(img,128);
    % alpha overlay
    if isscatter
      img_alpha = reshape(data_array_diff(ind,x),[256,128]);
    else
      img_alpha = [reshape(data_array_diff(ind,x),[256,128]) reshape(data_array_diff(ind,y),[256,128])];
    end
    img_alpha = circshift(img_alpha,128);
  elseif (length(data_array(:,x)) == 327684)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
    data_array_x_lh = data_array(1:163842,x);
    data_array_x_rh = data_array(163843:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array(1:163842,y);
      data_array_y_rh = data_array(163843:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img = [circshift(img_lh,128); img_rh];
    % alpha overlay
    data_array_x_lh = data_array_diff(1:163842,x);
    data_array_x_rh = data_array_diff(163843:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array_diff(1:163842,y);
      data_array_y_rh = data_array_diff(163843:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img_alpha = [circshift(img_lh,128); img_rh];
  elseif (length(data_array(:,x)) == 32492)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
    if isscatter
      img = reshape(data_array(ind,x),[256,128]);
    else
      img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
    end
    img = circshift(img,128);
    % alpha overlay
    if isscatter
      img_alpha = reshape(data_array_diff(ind,x),[256,128]);
    else
      img_alpha = [reshape(data_array_diff(ind,x),[256,128]) reshape(data_array_diff(ind,y),[256,128])];
    end
    img_alpha = circshift(img_alpha,128);
  elseif (length(data_array(:,x)) == 64984)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
    data_array_x_lh = data_array(1:32492,x);
    data_array_x_rh = data_array(32493:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array(1:32492,y);
      data_array_y_rh = data_array(32493:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img = [circshift(img_lh,96); circshift(img_rh,96)];
    % alpha overlay
    data_array_x_lh = data_array_diff(1:32492,x);
    data_array_x_rh = data_array_diff(32493:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array_diff(1:32492,y);
      data_array_y_rh = data_array_diff(32493:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img_alpha = [circshift(img_lh,96); circshift(img_rh,96)];
  else
    if isscatter
      img = data_array(:,x)';
      % alpha overlay
      img_alpha = data_array_diff(:,x)';
    else
      img = [data_array(:,y) data_array(:,x)]';
      % alpha overlay
      img_alpha = [data_array_diff(:,y) data_array_diff(:,x)]';
    end
  end

  % scale img to 0..64
  mn = 0.7; %min(data_array(:));
  mx = 1.0; %max(data_array(:));
  img = 64*((img - mn)/(mx-mn));
else
  % add slider for colume data
  set(H.mm,'Visible','on');
  set(H.mm_txt,'Visible','on');
  if isscatter
    img = data_array(:,:,x)';
    % alpha overlay
    img_alpha = data_array_diff(:,:,x)';
  else
    img = [data_array(:,:,y) data_array(:,:,x)]';
    % alpha overlay
    img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
  end
end

set(H.cbar,'TickLength',[0 0],'XTickLabel',linspace(0.7,1.0,7)); %round(100*linspace(min(YpY(mask2d(:))),max(YpY(mask2d(:))),5))/100);


% correct orientation
img = rot90(img,2);
img_alpha = rot90(img_alpha,2);

% display image with 2nd colorbar (gray)
image(65 + img);
if ~mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0],'HitTest','off','Interruptible','off');

% prepare alpha overlays for red and green colors
if alphaval > 0
  % get 2%/98% ranges of difference image
  range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

  hold on
  alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
  alpha_r   = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
  hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),...
    'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),...
    'HitTest','off','Interruptible','off','AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hold off
end

if mesh_detected
  xlabel('2D surface maps');
end

return
function varargout = cat_tst_qa_cleaner_intern(data,opt)
%% THIS FUNCTION IS TO FAT - SEPARATE AND CLEAN IT! 
%  Do not forget to remove old external version from SVN if this is done.
%  _____________________________________________________________________
%  Estimate quality grades of given rating of one (or more) protocols
%  with 2 to 6 grads to separate passed, (unassignable) and failed 
%  images, by finding the first peak in the image quality histogram  
%  and using its width (standard deviation) in a limited range. 
%  If multiple protocols are used, than use the site variable opt.site 
%  and use the site depending output rths.
%
%  The passed range can be variated by opt.cf with lower values for harder 
%  and higher values for softer thresholds (more passed images), where 
%  opt.cf=1, describes a range that is similar to about 1% BWP noise that 
%  is equal to 5 rps.
%  ROC evaluation showed that opt.cf=0.72 allows the best separation of 
%  images without and with artifacts, but if the majority of your data 
%  include light artifacts (e.g. by movements in young children) that 
%  a softer weighing, e.g. opt.cf=2, is preferable (maximum is 4). 
%
%  Use the selftest with randomly generated data to get a first impression:
%    cat_tst_qa_cleaner('test')
%  _____________________________________________________________________
%
%  This tool is still in development / undert test:
%   * the combination of different sites is not finished
%   * multiside output required a 'stacked' output
%
%  [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_remover(data[,opt])
%
%    Pth      .. global threshold for passed images 
%                (for odd grades this is in the middle of the unassignable)
%    rth      .. all global threshold(s) between grads
%    sq       .. estimated first peak and its std, where the std depend on
%                the number of grades!
%    rths     .. site depending thresholds between grads of each input 
%    rthsc    .. site depending thresholds between grads of each input 
%                (global corrected, removed further low quality data)
%    sqs      .. site depending first peaks and stds of passed data 
%
%    data     .. array of quality ratings or xml-files
%    opt      .. option structure
%     .grads  .. number of grads (2:6, default=6, see below)
%     .cf     .. factor for harder/softer thresholds (defaults=0.72)
%     .figure .. display histogramm with colored ranges of grads
%                 1 - use current figure
%                 2 - create new figure (default)
%                 3 - use one test figure (default in the selftest)
%  _____________________________________________________________________
%
%  Grades:
%    2 grads:
%      P   passed
%      F   failed
%    3 grads:
%      P   passed
%      U   unassignable
%      F   failed
%    4 grads:
%      P+  clear passed 
%      P-  just passed
%      F+  just failed
%      F-  clear failed
%    5 grads:
%      P+  clear passed 
%      P-  just passed
%      U   unassignable
%      F+  just failed
%      F-  clear failed
%    6 grads (default):
%      P+  clear passed 
%      P   passed 
%      P-  just passed
%      F+  just failed
%      F   failed
%      F-  clear failed
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_tst_qa_cleaner.m 1338 2018-07-23 12:51:04Z dahnke $ 


  clear th; 
  if ~exist('opt','var'), opt = struct(); end
  def.cf        = 0.72;                 % normalization factor for rating 
  def.grads     = 6;                    % number of grads (default = 6)
  def.model     = 1;                    % model used for rating
  def.figure    = 2;                    % figure=2 for new/own figure
  def.smooth    = 0;                    % smoothing of output data
  def.siterf    = 1000000;              % round factor to identify similar resolution level 
  def.siteavgperc = [0.10 0.90];        % ?
  opt = cat_io_checkinopt(opt,def); 
  opt.cf = max( 0 , min( 4 , opt.cf )); % limit of cf
  
  % test options
  %opt.model = 2;
  %opt.grads = 6;
  
  % if no intput is given use SPM select to get some xml-files
  if ~exist('data','var') || isempty(data)
    data = cellstr(spm_select(inf,'XML','select qa XML-files',{},pwd,'^cat_.*')); 
  elseif ischar(data)
    data = cellstr(data);
  end
  if isempty(data) || (iscell(data) && all(cellfun('isempty',data)))
    if nargout>=1, varargout{1} = 3; end
    if nargout>=2, varargout{2} = 3; end
    if nargout>=3, varargout{3} = [2.5 0.5]; end
    if nargout>=4, varargout{4} = 3*ones(size(data)); end
    if nargout>=5, varargout{5} = 3*ones(size(data)); end
    if nargout>=6, varargout{6} = repmat([2.5 0.5],numel(data),1); end

    return;
  end
  if iscell(data) && ~strcmp(data{1},'test')
    fprintf('Load XML data');
    P = data; 
    xml = cat_io_xml(data,struct(),'read',1); clear data; 
    for di=1:numel(xml)
      opt.site(di,1) = xml(di).qualityratings.res_RMS; 
      data(di,1)     = xml(di).qualityratings.NCR; 
    end,
  end
  

  % --------------------------------------------------------------------
  % If a site variable is given (e.g. by the RMS resolution) then call
  % the cleanup for each subset. The threshold will be collected in a 
  % vector [markthss x opt.grads] with the same length as data. 
  % Nevertheless an average threshold will is estimated as average of 
  % the percentual range give by opt.siteavgperc with e.g. [0.1 0.9] to
  % concider 80% of the data.
  %  -------------------------------------------------------------------
  if isfield(opt,'site')
    if numel(opt.site)~=numel(data),
      error('cat_tst_qa_cleaner:numelsitedata','Numer of elements in data and opt.site have to be equal.\n');
    end
    opt.site = round(opt.site*opt.siterf)/opt.siterf; 
    sites    = unique(opt.site); 
    markth   = zeros(numel(sites),opt.grads-1); 
    markths  = zeros(numel(data),opt.grads-1); 
    siteth   = zeros(numel(data),2); 
    for si=1:numel(sites)
      sdatai = find(opt.site==sites(si));
      opts = opt; 
      opts = rmfield(opts,'site');
      opts.figure = 0; 
      [Sth,markth(si,:),out{1:4}] = cat_tst_qa_cleaner(data(sdatai),opts); %#ok<ASGLU>
      markths(sdatai,:) = repmat(markth(si,:),numel(sdatai),1); 
      siteth(sdatai,:)  = out{4}; 
    end
    % estimate global threshold
    markthss = sortrows(markth);
    th = cat_stat_nanmean(markthss(max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(1)))):...
                                   max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(2)))),:),1); 
    sd  = out{3}; 
    thx = out{4}; 
    % modify local rating based on the global one                 
    markths2 = markths;
    markths2 = min(markths2,1.2*repmat(th,size(markths2,1),1)); % higher thresholds even for sides with low rating 
    markths2 = max(markths2,0.8*repmat(th,size(markths2,1),1)); % lower  thresholds even for sides with high rating 
    d  = data; 

  else
    %  -----------------------------------------------------------------
    %  Simulate data, if no data is given by several normal distributed
    %  random numbers.
    %  -----------------------------------------------------------------
    if exist('data','var') && ~(iscell(data) && strcmp(data{1},'test'))
      d = data; 
      if numel(d)==0, 
        if nargout>=1, varargout{1} = nan; end
        if nargout>=2, varargout{2} = nan(1,opt.grads); end
        if nargout>=3, varargout{3} = nan(1,2); end
        if nargout>=4, varargout{4} = nan(size(data)); end
        if nargout>=5, varargout{5} = nan(size(data)); end
        if nargout>=6, varargout{6} = nan(size(data)); end
        return;
      end
    elseif iscell(data) && strcmp(data{1},'test')
      % Testcases with different quality ratings
      scans      = 100; % number of scans (per site) for simulation
      testcase   = round(rand(1)*10);
      randoffset = 0.5*randn(1,4);

      switch testcase
        case 0 % good quality, no outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.15)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.03)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.02))];
         case 1 % good quality, with average outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.40)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.15)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.05))];
        case 2 % good-average quality, with outlier group 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];
        case 3 % good-average quality, without outlier group 
          d = [2.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               3.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))]; 
        case 4 % average to low quality, with light falloff  
          d = [3.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               3.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 5 % high to good quality, with light falloff  
          d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               2.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 6 % high quality, no outlier
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.13)), ...
               3.0 + randoffset(3) + 0.3*randn(1,round(scans*0.05)), ...
               5.0 + randoffset(4) + 0.3*randn(1,round(scans*0.02))];
        case 7 % good quality with second average peak 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];   
        case 8 % good quality with second low quality peak 
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(2) + 0.2*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];    
        case 9 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
               3.0 + randoffset(2) + 0.3*randn(1,round(scans*0.20)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.0 + randoffset(4) + 0.8*randn(1,round(scans*0.50))];          
        case 10 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.10)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(4) + 1.0*randn(1,round(scans*0.60))];           
      end

      % remove high quality outlier and set them to normal
      cor = max(1,median(d)-std(d)/2);
      md= d<(cor); d(md) = cor + 0.05*randn(1,sum(md));
      
      % set selftest figure
      opt.figure = 3; 
    end

    
  %% Models
  %  -------------------------------------------------------------------
  %  I start with several ideas that all based on a similar idea: to 
  %  find the first peak that is given by the subset of images without
  %  inferences and to use the variance of this peak for further scaling
  %  of subsets for other grads. As far as IQR is already scaled, we 
  %  can limit the variance value ... e.g. the rating has an error of 
  %  0-2 rps (0.0-0.2 mark points) that is very low for high-quality data
  %  and higher for low-quality data. Due to our the general subdivion 
  %  of the rating scale in +,o, and - (e.g. B+,B,B-) we got a subrange 
  %  of 3.33 rps (1/3 mark points) that gives some kind of upper limit.
  %  -------------------------------------------------------------------
    thx = nan; sd = nan; th = zeros(1,opt.grads-1); 
    switch opt.model
      case 0
        % only global thresholding ... 
        % this is just to use the color bar output 
        thx = 3; 
        sd  = 1; 
        th  = 1.5:1:100;
        th(6:end) = []; 
      case 1 
        % kmeans model:
        % * estimate peaks based on the histogram
        % * mix the first and second peak until it fits to 30% of the data 
        %   or until the number of loops is similar the number of peaks 
        % * use the std give by one BWP noise level (0.5) to describe the 
        %   variance the passed interval.
        
        hx = hist(d,0.5:1:5.5);
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          if sum(d<thx(i))/numel(d) < 0.3
            thx(1) = cat_stat_nanmean(thx(1:2));
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sd    = 0.25 / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise
        th(1) = thx(1) - sdx(1) + 2*sd(1); %- mean(sdx(1:min(3,numel(sdx)))) 
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 
        end
      case 2
        % similar to case 1, but with std optimization based on the data 
        % ... surprisingly the simple model 1 works better
        
        hx = hist(d,0.5:1:5.5); 
        %for i=1:1, hx(2:end-1) = cat_stat_nanmean(cat(1,hx(1:end-2),hx(2:end-1),hx(3:end)),1); end
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          %if numel(thx)>i && sum(d<thx(i))/numel(d) < 0.05
          %  thx(1) = []; sdx(1) = [];
          if sum(d<thx(i))/numel(d) < 0.3 %numel(thx)>i && 
            thx(1) = cat_stat_nanmean(thx(1:2)); 
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sdx(1) = cat_stat_nanstd(d(d<thx(1)));
        [thx,sdx] = kmeans3D(d(d<=(max([min(d),thx(1)+sdx(1)]))),3); thx=thx(2); sdx=sdx(2);  %sdx = sdx./thx;
        sd    = min(1/3,max(1/6,sdx(1))) / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise*16
        th(1) = thx(1) - sdx(1) + 2*sd(1);
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 2*2/3*
        end
      
    end
    
    markths  = repmat(mean(th(floor(opt.grads/2):ceil(opt.grads/2))),size(data));
    markths2 = markths;
    siteth   = repmat([thx(1) sd],numel(data),1); 
  end
  
  
%% Print
%  ---------------------------------------------------------------------
%  This part is just for to plot a colorated histogram and the percents
%  of images in each group.
%  ---------------------------------------------------------------------
  if opt.figure
    if opt.figure==2
      f = figure;
      set(f,'color','w')
    elseif opt.figure==3
      f = findobj('type','figure','name','qa_cleaner_test');
      if isempty(f), figure('name','qa_cleaner_test'); else figure(f(1)); clf(f(1)); end
    end
    box on;
    
    %figure
    ss = 0.05; 
    [h,r]  = hist(d,0.5:ss:10.5); 
    for i=1:opt.smooth, h(2:end-1) = cat_stat_nanmean(cat(1,h(1:end-2),h(2:end-1),h(3:end)),1); end
    sh = 1; %sum(h);
    
    % background histogram (all data)
    %bar(r,h/sh,'facecolor',[0.8 0.8 0.8],'edgecolor','none');
    %fill(r,h/sh,[0.8 0.8 0.8],'edgecolor','none');
    hold on
    
    yl = [0 max(h)+1]; ylim(yl);
    % main grid
    for i=1.5:6,       plot([i i],ylim,'color',[0.8 0.8 0.8]); end
    switch numel(th)
      case 1
        hx = h; hx(r> th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(1))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 2
        hx = h; hx(r>=th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1) | r>th(2)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% unassignable' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 3
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(2))/numel(d)*100),'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% failed+',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed-',sum(d>=th(3))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 4
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(4)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(2))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% check' ,sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100)          ,'color',[0.7  0.0  0.0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.71,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.67,sprintf('%5.2f%% unassignable',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.63,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.59,sprintf('%5.2f%% failed-',sum(d>=th(4))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 5
        % plot
        testbar=0; % it would be cool to use bars but they failed at least in MATLAB R2013 and killed the axis positions...
        if testbar==1
          hx = h; hx(r>=th(1)+ss) = 0;           bar(r,hx/sh,'facecolor',[0.0  0.5  0.2],'edgecolor','none','barwidth',1);
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; bar(r,hx/sh,'facecolor',[0.4  0.7  0.1],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; bar(r,hx/sh,'facecolor',[0.7  0.8  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; bar(r,hx/sh,'facecolor',[0.9  0.6  0.4],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; bar(r,hx/sh,'facecolor',[0.75 0.3  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(5)-ss) = 0;           bar(r,hx/sh,'facecolor',[0.6  0.15 0.1],'edgecolor','none','barwidth',1);      
        else
          hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(5)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none'); 
        end
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(3))/numel(d)*100) ,'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% passed-',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.55,sprintf('%5.2f%% failed' ,sum(d>=th(4) & d<th(5))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.50,sprintf('%5.2f%% failed-',sum(d>=th(5))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
    end
    xlim([min(r),6.5]); 
    
    % subgrid
    for i=5/6:1/3:6.4, plot([i i],[0 0.03]*max(ylim),'color',[0.2 0.2 0.2]); end
        
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    
    
    % colored main grads
    FS = get(gca,'Fontsize')*1.3;
    set(gca,'XTick',0.5:1:6.5,'XTickLabel',{'100','90','80','70','60','50','40'},'TickLength',[0.02 0.02]);
    % further color axis objects...
    axA = copyobj(gca,gcf); axB = copyobj(axA,gcf); axC = copyobj(gca,gcf); 
    axD = copyobj(gca,gcf); axE = copyobj(gca,gcf); axF = copyobj(gca,gcf);
    % set colors...
    set(axA,'YTick',[],'XTickLabel',{},'XTick',1,'XColor',color(QMC,1),'Color','none','XTicklabel','A','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axB,'YTick',[],'XTickLabel',{},'XTick',2,'XColor',color(QMC,2),'Color','none','XTicklabel','B','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axC,'YTick',[],'XTickLabel',{},'XTick',3,'XColor',color(QMC,3),'Color','none','XTicklabel','C','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axD,'YTick',[],'XTickLabel',{},'XTick',4,'XColor',color(QMC,4),'Color','none','XTicklabel','D','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axE,'YTick',[],'XTickLabel',{},'XTick',5,'XColor',color(QMC,5),'Color','none','XTicklabel','E','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axF,'YTick',[],'XTickLabel',{},'XTick',6,'XColor',color(QMC,6),'Color','none','XTicklabel','F','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    hold off; 
    
    if isfield(opt,'site') && numel(sites>1);
      title(sprintf('Histogram (cf=%0.2f) - global treshold for multisite output (n=%d)',opt.cf,numel(sites)),'Fontsize',FS);
    else
      title(sprintf('Histogram (cf=%0.2f)',opt.cf),'Fontsize',FS);
    end
    xlabel('IQR (rps)','Fontsize',FS); 
    ylabel('number of scans','Fontsize',FS); 
  end
  %%
  MarkColor = cat_io_colormaps('marks+',40); 
  if isfield(opt,'site') && numel(sites)>1, globcorr = ' (global corrected)'; else globcorr = ''; end
  if exist('P','var')
    files = P(data<=markths2(:,3)); 
    fprintf('PASSED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 0
      iqrs  = [xml(data<=markths2(:,3)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    else
      
    end
    
    % bad files ...
    files = P(data>markths2(:,3) & data<=markths2(:,4)); 
    fprintf('FAILED+%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,3) & data<=markths2(:,4)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,4) & data<=markths2(:,5)); 
    iqrs  = [xml(data>markths2(:,4) & data<=markths2(:,5)).qualityratings];
    if 1
      fprintf('FAILED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,5)); 
    fprintf('FAILED-%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,5)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
  end
  
  
  %% create output
  if nargout>=1, varargout{1} = mean(th(floor(opt.grads/2):ceil(opt.grads/2))); end
  if nargout>=2, varargout{2} = th; end
  if nargout>=3, varargout{3} = [thx(1) sd(1)]; end
  if nargout>=4, varargout{4} = markths;  end
  if nargout>=5, varargout{5} = markths2; end
  if nargout>=6, varargout{6} = siteth; end
  
  if 0
    %%
    b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
                'value',zeta, 'min',0, 'max',1);
    bgcolor = f.Color;
    bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                    'String','0','BackgroundColor',bgcolor);
    bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                    'String','1','BackgroundColor',bgcolor);
    bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                    'String','Damping Ratio','BackgroundColor',bgcolor);
  end
return
