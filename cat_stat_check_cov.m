function varargout = cat_stat_check_cov(job)
%cat_stat_check_cov. To check covariance across sample.
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%
% varargout = cat_stat_check_cov(job)
%  
% job                .. SPM job structure
%  .data_vol         .. volume files
%  .data_surf        .. surface files
%  .gap              .. gap between slices (default=3)
%  .c                .. confounds
%  .data_xml         .. optional xml QC data
%  .verb             .. print figures
%
% varargout          .. output structure 
%  .table            .. covariance table
%  .covmat           .. covariance matrix
%  .sorttable        .. sorted covariance table
%  .sortcovmat       .. sorted covariance matrix
%  .threshold_cov    .. lower threshold for covariance (mean - 4*std)
%
% Example: 
%   cat_stat_check_cov(struct('data_vol',{{ files }} ,'gap',3,'c',[],'data_xml',{{}}));
%
% See also cat_stat_check_cov2
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if cat_get_defaults('extopts.expertgui')>2
  if nargout
    varargout = cat_stat_check_cov2(job);
  else
    cat_stat_check_cov2(job);
  end
  return
end

clearvars -GLOBAL H;       % clear old
global H

% show data by fileorder
H.sorted = 0;

if nargin == 0
  error('No argument given.');
end

H.sample   = [];
G          = [];
n_subjects = 0;
H.alphaval = 0.5;
H.names_changed = 0;

% check for global scaling
if isfield(job,'gSF')
  is_gSF = 1;
  gSF = job.gSF;
else
  is_gSF = 0;
end

if ~isfield(job,'verb')
  job.verb = true;
end

if isfield(job,'show_violin')
  H.show_violin = job.show_violin;
else
  H.show_violin = 0;
end

if isfield(job,'show_name')
  H.show_name = job.show_name;
else
  H.show_name = 0;
end

% read filenames for each sample and indicate sample parameter
if isfield(job,'data_vol')
  H.mesh_detected = 0;
  n_samples = numel(job.data_vol);
  for i=1:n_samples
    
    if size(job.data_vol{i},1) == 1 % 4D data
      [pth,nam,ext] = spm_fileparts(char(job.data_vol{i}));
      % remove ",1" at the end
      job.data_vol{i} = fullfile(pth,[nam ext]);
    end
    
    V0 = spm_data_hdr_read(char(job.data_vol{i}));
    n_subjects = n_subjects + length(V0);
      
    if i==1, H.V = V0;
    else,    H.V = [H.V; V0]; end

    H.sample = [H.sample, i*ones(1,length(V0))];
  end
  sep = job.gap;
else
  H.mesh_detected = 1;
  n_samples = numel(job.data_surf);
  sinfo = cat_surf_info(char(job.data_surf{1}(1,:)));
  H.Pmesh = gifti(sinfo.Pmesh);
  for i=1:n_samples
    V0 = spm_data_hdr_read(char(job.data_surf{i}));
    n_subjects = n_subjects + length(V0);
      
    if i==1, H.V = V0;
    else,    H.V = [H.V; V0]; end
    H.sample = [H.sample, i*ones(1,size(job.data_surf{i},1))];
  end
end
    
if ~isempty(job.c)
  for i=1:numel(job.c)
    G = [G job.c{i}];
  end
end

if isempty(char(job.data_xml))
  H.isxml = 0;
  QM_names = '';
  xml_files = [];
else
  xml_files = char(job.data_xml);
  H.isxml = 1;
end

pth = spm_fileparts(H.V(1).fname);
report_folder = fullfile(spm_fileparts(pth),'report');
subfolder = 1;
% check whether report subfolder exists
if ~exist(report_folder,'dir')
  report_folder = pth;
  subfolder = 0;
end

% search xml report files if not defined
H.found_xml = 0;
if ~H.isxml
  xml_files = spm_select('List',report_folder,'^cat_.*\.xml$');
  if ~isempty(xml_files)

    % find part of xml-filename in data files to get the prepending string
    % (e.g. mwp1)
    i = 1; j = 1;
    while i <= n_subjects
      while j <= size(xml_files,1)
        % remove "cat_" and ".xml" from name
        fname = deblank(xml_files(j,:));
        fname = fname(5:end-4);

        % and find that string in data filename
        ind = strfind(H.V(i).fname,fname);
        if ~isempty(ind)
          [pth, prep_str] = spm_fileparts(H.V(1).fname(1:ind-1));
          H.isxml = 1;
          H.found_xml = 1;
          fprintf('Corresponding xml-files were found.\n')
          i = n_subjects;
          j = size(xml_files,1);
          break
        else
          j = j + 1;
        end
      end
      i = i + 1;
    end
  end  
end

if H.isxml && (size(xml_files,1) < n_subjects)
    fprintf('WARNING: Less XML-files than data sets found. XML-files will be not used.\n');
    H.isxml = 0;
  end
  
if H.isxml
  if H.mesh_detected
    QM = ones(n_subjects,5);
    QM_names = char('Noise','Bias','Weighted overall image quality (IQR)','Euler number','Size of topology defects');
  else
    QM = ones(n_subjects,3);
    QM_names = char('Noise','Bias','Weighted overall image quality (IQR)');
  end
  
  spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
    % get basename for data files
    [pth, data_name] = fileparts(H.V(i).fname);
    
    % remove ending for rigid or affine transformed files
    data_name = strrep(data_name,'_affine','');
    data_name = strrep(data_name,'_rigid','');
    % use xml-file if found by name
    if H.found_xml
      % get report folder
      if subfolder
        report_folder = fullfile(spm_fileparts(pth),'report');
      else
        report_folder = pth;
      end

      % remove prep_str from name and use report folder and xml extension
      if H.mesh_detected
        % for meshes we aso have to remove the additional "." from name
        tmp_str = strrep(data_name,prep_str,'');
        xml_file = fullfile(report_folder,['cat_' tmp_str(2:end) '.xml']);
      else
        xml_file = fullfile(report_folder,['cat_' strrep(data_name,prep_str,'') '.xml']);
      end
    else
    
      [pth, xml_name] = fileparts(deblank(xml_files(i,:)));
      % remove leading 'cat_'
      xml_name = xml_name(5:end);

      % check for filenames
      if isempty(strfind(data_name,xml_name))
        fprintf('Please check file names because of deviating subject names:\n%s vs. %s\n',H.V(i).fname,xml_files(i,:));
      end
      xml_file = deblank(xml_files(i,:));
    end
    
    if exist(xml_file,'file')
      xml = cat_io_xml(xml_file);
    else
      fprintf('File %s not found. Skip use of xml-files for quality measures.\n',xml_file);
      H.isxml = 0;
      H.found_xml = 0;
      break
    end
    if ~isfield(xml,'qualityratings') && ~isfield(xml,'QAM')
      fprintf('Quality rating is not saved for %s. Report file %s is incomplete.\nPlease repeat preprocessing amd check for potential errors in the ''err'' folder.\n',H.V(i).fname,xml_files(i,:));    
      return
    end
    if H.mesh_detected
      if isfield(xml.qualityratings,'NCR')
      % check for newer available surface measures
        if isfield(xml.subjectmeasures,'EC_abs')
          QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR xml.subjectmeasures.EC_abs xml.subjectmeasures.defect_size];
        else
          QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR NaN NaN];
        end
      else % also try to use old version
        QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
      end
    else
      if isfield(xml.qualityratings,'NCR')
        QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR];
      else % also try to use old version
        QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
      end
    end
    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');
  
  % remove last two columns if EC_abs and defect_size are not defined
  if H.mesh_detected & all(isnan(QM(:,4))) & all(isnan(QM(:,5)))
    QM = QM(:,1:3);
  end
  
end

[pth,nam] = spm_fileparts(H.V(1).fname);

if H.mesh_detected
  % load surface texture data
  H.data = single(spm_data_read(H.V));
  Y = H.data';
  
  H.range_data98 =  cat_vol_iscaling(Y(:),[0.02 0.98]);
  
  % optional global scaling
  if is_gSF
    for i=1:numel(H.V)
      Y(:,2) = Y(:,2)*gSF(i);
    end
  end
  
  Y(isnan(Y)) = 0;
  
else
  % voxelsize and origin
  vx =  sqrt(sum(H.V(1).mat(1:3,1:3).^2));
  Orig = H.V(1).mat\[0 0 0 1]';

  if length(H.V)>1 && any(any(diff(cat(1,H.V.dim),1,1),1))
    error('images don''t all have same dimensions')
  end
  if max(max(max(abs(diff(cat(3,H.V.mat),1,3))))) > 1e-8
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
H.FS = cat_get_defaults('extopts.fontsize');

H.pos = struct(...
    'fig',    [10  10  1.3*ws(3) 1.1*ws(3)],... % figure
    'cbar',   [0.045 0.050 0.700 0.020],... % colorbar for correlation matrix
    'corr',   [-0.02 0.100 0.825 0.825],... % correlation matrix
    'scat',   [0.050 0.050 0.700 0.825],... % scatter plot
    'close',  [0.775 0.925 0.200 0.050],... % close button
    'show',   [0.775 0.875 0.200 0.050],... % button to show worst cases
    'boxp',   [0.775 0.820 0.200 0.050],... % button to display boxplot
    'sort',   [0.775 0.775 0.200 0.050],... % button to enable ordered matrix
    'fnambox',[0.775 0.735 0.200 0.050],... % show filenames?
    'text',   [0.775 0.600 0.200 0.150],... % textbox
    'aslider',[0.775 0.555 0.200 0.040],... % slider for alpha overlay
    'slice',  [0.775 0.050 0.200 0.400],... % two single images according to position of mouse pointer
    'sslider',[0.775 0.010 0.200 0.040],... % slider for z-slice   
    'plotbox',[0.875 0.735 0.200 0.050]);   % switch between boxplot and violin plot 

if H.mesh_detected

  Y = Y - repmat(mean(Y,2), [1 size(Y,2)]);

  % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
  if ~isempty(G) 
    Ymean = repmat(mean(Y), [n_subjects 1]);
    Y = Y - G*(pinv(G)*Y) + Ymean;
  end

  H.YpY = (Y*Y')/n_subjects;

  % calculate residual mean square of mean adjusted Y
  Y = Y - repmat(mean(Y,1), [n_subjects 1]);
  
%  MSE = sum(Y.*Y,2);

  clear Y
else
  % consider image aspect ratio
  H.pos.slice = [0.775 0.050 0.20 0.40*H.V(1).dim(2)/H.V(1).dim(1)];

  slices = 1:sep:H.V(1).dim(3);

  dimx = length(1:sep:H.V(1).dim(1));
  dimy = length(1:sep:H.V(1).dim(2));
  Y = zeros(n_subjects, prod(dimx*dimy));
  H.YpY = zeros(n_subjects);
%  MSE = zeros(n_subjects,1);
  H.data = zeros([H.V(1).dim(1:2) n_subjects],'single');

  %-Start progress plot
  %-----------------------------------------------------------------------
  spm_progress_bar('Init',H.V(1).dim(3),'Check correlation','planes completed')

  for j=slices

    M  = spm_matrix([0 0 j 0 0 0 sep sep sep]);

    for i = 1:n_subjects
      H.img = spm_slice_vol(H.V(i),M,[dimx dimy],[1 0]);
      H.img(isnan(H.img)) = 0;
      Y(i,:) = H.img(:);
      if is_gSF
        Y(i,:) = Y(i,:)*gSF(i);
      end
    end

    % make sure data is zero mean
    Y = Y - repmat(mean(Y,2), [1 prod(dimx*dimy)]);

    % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
    if ~isempty(G) 
      [ind_inf,tmp] = find(isinf(G) | isnan(G));
      if ~isempty(ind_inf)
        fprintf('Nuisance parameter for %s is Inf or NaN.\n',H.V(ind_inf).fname);
        return
      end
      Ymean = repmat(mean(Y), [n_subjects 1]);
      Y = Y - G*(pinv(G)*Y) + Ymean;
    end
    H.YpY = H.YpY + (Y*Y')/n_subjects;

    % calculate residual mean square of mean adjusted Y
    Y = Y - repmat(mean(Y,1), [n_subjects 1]);
    
    %MSE = MSE + sum(Y.*Y,2);

    spm_progress_bar('Set',j);  

  end

  % correct filenames for 4D data
  if strcmp(H.V(1).fname, H.V(2).fname)
    H.names_changed = 1;
    H.Vchanged = H.V;
    for i=1:n_subjects
      [pth,nam,ext] = spm_fileparts(H.V(i).fname);
      H.V(i).fname = fullfile(pth, [nam sprintf('%04d',i) ext]);
    end
  end
  
  clear Y
  spm_progress_bar('Clear');
end

% normalize YpY
d        = sqrt(diag(H.YpY)); % sqrt first to avoid under/overflow
dd       = d*d';
H.YpY    = H.YpY./(dd+eps);
t        = find(abs(H.YpY) > 1); 
H.YpY(t) = H.YpY(t)./abs(H.YpY(t));
H.YpY(1:n_subjects+1:end) = sign(diag(H.YpY));

% extract mean correlation for each data set
H.mean_cov = zeros(n_subjects,1);
for i=1:n_subjects
  % extract row for each subject
  cov0 = H.YpY(i,:);

  % remove cov with its own
  cov0(i) = [];
  H.mean_cov(i) = mean(cov0);
end

if job.verb, fprintf('\n'); end
fname_m = [];
fname_tmp = cell(n_samples,1);
fname_s   = cell(n_samples,1);
fname_e   = cell(n_samples,1);

for i=1:n_samples
  [tmp, fname_tmp{i}] = spm_str_manip(char(H.V(H.sample == i).fname),'C');
  fname_m = [fname_m; fname_tmp{i}.m];
  fname_s{i} = fname_tmp{i}.s;
  fname_e{i} = fname_tmp{i}.e;
  if job.verb
    fprintf('Compressed filenames sample %d: %s  \n',i,tmp);
  end
end

H.filename = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});

% print suspecious files with high cov
% use slightly higher threshold for (smoothed) mesh data
YpY_tmp = H.YpY - triu(H.YpY);
if H.mesh_detected
  [indx, indy] = find(YpY_tmp>0.950 & YpY_tmp < (1-eps));
else
  [indx, indy] = find(YpY_tmp>0.925);
end

% if more than 25% of the data were found, this points to longitudinal data of one subject
% and no warning will appear
if ~isempty(indx) && job.verb
 if (length(indx) < 0.25*n_subjects)
    fprintf('\nWARNING: Unusual large correlation (check that subjects are not identical):\n');
    for i = 1:length(indx)
      % exclude diagonal
      if indx(i) ~= indy(i)
        if n_samples > 1
          fprintf('%s (sample %d) and %s (sample %d): %g\n',H.filename.m{indx(i)},H.sample(indx(i)),H.filename.m{indy(i)},H.sample(indy(i)),H.YpY(indx(i),indy(i)));
        else
          fprintf('%s and %s: %g\n',H.filename.m{indx(i)},H.filename.m{indy(i)},H.YpY(indx(i),indy(i)));
        end
      end
    end
  else
    fprintf('\nMany unusual large correlations were found (e.g. common in longitudinal data).\n');
  end
end

[indx, indy] = find(YpY_tmp==1);

% give warning that data are identical
if ~isempty(indx) && job.verb
  fprintf('\nWARNING: Data of these subjects are identical!\n');
  for i = 1:length(indx)
    if n_samples > 1
      fprintf('%s (sample %d) and %s (sample %d)\n',H.filename.m{indx(i)},H.sample(indx(i)),H.filename.m{indy(i)},H.sample(indy(i)));
    else
      fprintf('%s and %s\n',H.filename.m{indx(i)},H.filename.m{indy(i)});
    end
  end
end

% sort data
[H.mean_cov_sorted, H.ind_sorted] = sort(H.mean_cov,'descend');
H.YpYsorted = H.YpY(H.ind_sorted,H.ind_sorted);

H.ind_sorted_display = H.ind_sorted;

threshold_cov = mean(H.mean_cov) - 2*std(H.mean_cov);
n_thresholded = min(find(H.mean_cov_sorted < threshold_cov));

if ~isempty(n_thresholded) && job.verb
  fprintf('\nThese data have a mean correlation below 2 standard deviations.\n');
  fprintf('This does not necessarily mean that you have to exclude these data. However, these data have to be carefully checked:\n');
  for i=n_thresholded:n_subjects
    fprintf('%s: %3.3f\n',H.V(H.ind_sorted(i)).fname,H.mean_cov_sorted(i));
  end
end

if nargout>0
  varargout{1} = struct('table',{[cellstr([{H.V.fname}]'),num2cell(H.mean_cov)]},...
                        'covmat',H.YpY,...
                        'sorttable',{[cellstr([{H.V(H.ind_sorted).fname}]'),num2cell(H.mean_cov_sorted)]},...
                        'sortcovmat',H.YpYsorted, ...
                        'cov',H.mean_cov,...
                        'threshold_cov',threshold_cov);
end

if job.verb
  % create figure
  H.figure = figure(2);
  clf(H.figure);

  set(H.figure,...
     'MenuBar','none',...
     'Position',H.pos.fig,...
  ...   'DefaultTextFontSize',H.FS,...
  ...   'DefaultUicontrolFontSize',H.FS,...
     'NumberTitle','off');

  if H.mesh_detected
    set(H.figure,'Name','Click in image to display surfaces');
  else
    set(H.figure,'Name','Click in image to display slices');
  end

  cm = datacursormode(H.figure);
  set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
  try set(cm,'NewDataCursorOnClick',false); end

  % add colorbar
  H.cbar = axes('Position',H.pos.cbar,'Parent',H.figure);
  try, image(H.cbar); end
  set(get(H.cbar,'children'),'HitTest','off','Interruptible','off');
  set(H.cbar,'Ytick',[],'YTickLabel',''); 

  H.isscatter = 0;
  show_matrix(H.YpY, H.sorted);

  % create two colormaps
  cmap = [jet(64); gray(64)];
  colormap(cmap)

  % display YTick with 5 values (limit accuracy for floating numbers)
  set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
    round(100*linspace(min(H.YpY(:)),max(H.YpY(H.YpY~=1)),5))/100,'TickLength',[0 0]);

  % add button for closing all windows
  H.close = uicontrol(H.figure,...
          'string','Close','Units','normalized',...
          'position',H.pos.close,...
          'Style','Pushbutton','HorizontalAlignment','center',...
          'callback','for i=2:26, try close(i); end; end;',...
          'ToolTipString','Close windows',...
          'Interruptible','on','Enable','on');

  % check button
  H.show = uicontrol(H.figure,...
          'string','Check most deviating data','Units','normalized',...
          'position',H.pos.show,...
          'Style','Pushbutton','HorizontalAlignment','center',...
          'callback',@check_worst_data,...
          'ToolTipString','Display most deviating files',...
          'Interruptible','on','Enable','on');

  % create popoup menu 
  if H.isxml

    % estimate ratio between weighted overall quality (IQR) and mean corr. 
    H.X = [H.mean_cov, QM(:,3)];
    H.IQRratio = (H.X(:,2)/std(H.X(:,2)))./(H.X(:,1)/std(H.X(:,1)));

    str  = { 'Boxplot...','Mean correlation',QM_names,'Norm. Ratio IQR/Mean Correlation'};

    if size(QM,2) == 5
      tmp  = { {@show_boxplot, H.mean_cov, 'Mean correlation  ', 1},...
               {@show_boxplot, QM(:,1), QM_names(1,:), -1},...
               {@show_boxplot, QM(:,2), QM_names(2,:), -1},...
               {@show_boxplot, QM(:,3), QM_names(3,:), -1},...
               {@show_boxplot, QM(:,4), QM_names(4,:), -1},...
               {@show_boxplot, QM(:,5), QM_names(5,:), -1},...
               {@show_boxplot, H.IQRratio, 'Norm. Ratio IQR/Mean Correlation  ', -1} };
    else
      tmp  = { {@show_boxplot, H.mean_cov, 'Mean correlation  ', 1},...
               {@show_boxplot, QM(:,1), QM_names(1,:), -1},...
               {@show_boxplot, QM(:,2), QM_names(2,:), -1},...
               {@show_boxplot, QM(:,3), QM_names(3,:), -1},...
               {@show_boxplot, H.IQRratio, 'Norm. Ratio IQR/Mean Correlation  ', -1} };
    end

  else
    str  = { 'Boxplot...','Mean correlation'};
    tmp  = { {@show_boxplot, H.mean_cov, 'Mean correlation  ', 1} };
  end

  H.boxp = uicontrol(H.figure,...
          'string',str,'Units','normalized',...
          'position',H.pos.boxp,'UserData',tmp,...
          'Style','PopUp','HorizontalAlignment','center',...
          'callback','spm(''PopUpCB'',gcbo)',...
          'ToolTipString','Display boxplot',...
          'Interruptible','on','Visible','on');

  if H.isxml
    str  = { 'Image...','Mean Correlation: Order by selected filenames','Mean Correlation: Sorted by mean correlation','Norm. Ratio IQR/Mean Correlation'};
    tmp  = { {@show_matrix, H.YpY, 0},...
             {@show_matrix, H.YpYsorted, 1},...
             {@show_IQRratio, H.X} };
  else
    str  = { 'Correlation matrix...','Order by selected filename','Sorted by mean correlation'};
    tmp  = { {@show_matrix, H.YpY, 0},...
             {@show_matrix, H.YpYsorted, 1} };
  end

  H.sort = uicontrol(H.figure,...
          'string',str,'Units','normalized',...
          'position',H.pos.sort,'UserData',tmp,...
          'Style','PopUp','HorizontalAlignment','center',...
          'callback','spm(''PopUpCB'',gcbo)',...
          'ToolTipString','Sort matrix',...
          'Interruptible','on','Visible','on');

  H.text = uicontrol(H.figure,...
          'Units','normalized','position',H.pos.text,...
          'String',{'','Click in image to display slices'},...
          'Style','text','HorizontalAlignment','center',...
          'ToolTipString','Select slice for display',...
          'FontSize',H.FS-2);

  % add slider and opacity control only for volume data
  if ~H.mesh_detected
    H.alpha = uicontrol(H.figure,...
          'Units','normalized','position',H.pos.aslider,...
          'Min',0,'Max',1,...
          'Style','slider','HorizontalAlignment','center',...
          'callback',@update_alpha,'Value',0.5,...
          'ToolTipString','Change Opacity of pos. (green colors) and neg. (red colors) differences to sample mean',...
          'SliderStep',[0.01 0.1],'Visible','off');

    H.alpha_txt = uicontrol(H.figure,...
          'Units','normalized','HorizontalAlignment','center',...
          'Style','text','BackgroundColor',[0.8 0.8 0.8],...
          'Position',[H.pos.aslider(1) H.pos.aslider(2)-0.005 0.2 0.02],...
          'String','Overlay Opacity of Differences to Sample Mean',...
          'FontSize',H.FS-2,'Visible','off');

    H.mm = uicontrol(H.figure,...
          'Units','normalized','position',H.pos.sslider,...
          'Min',(1 - Orig(3))*vx(3),'Max',(H.V(1).dim(3) - Orig(3))*vx(3),...
          'Style','slider','HorizontalAlignment','center',...
          'callback',@update_slices_array,...
          'ToolTipString','Select slice for display',...
          'SliderStep',[0.005 0.05],'Visible','off');

    H.mm_txt = uicontrol(H.figure,...
          'Units','normalized','HorizontalAlignment','center',...
          'Style','text','BackgroundColor',[0.8 0.8 0.8],...
          'Position',[H.pos.sslider(1) H.pos.sslider(2)-0.005 0.2 0.02],...
          'String','Slice [mm]','Visible','off','FontSize',H.FS-2);

    update_slices_array;
  end

  show_boxplot(H.mean_cov,'Mean correlation  ',1);

  if isfield(job,'save') && job.save
    %% filenames
    if ~isempty(job.fname)
      dpi = cat_get_defaults('print.dpi'); 
      if isempty(dpi), dpi = 150; end

      fignames   = {'matrix','boxplot'};
      figuresids = {figure(2),spm_figure('FindWin','Graphics')};
      if isempty(job.outdir{1}), job.outdir{1} = pwd; end

      % save
      warning('OFF','MATLAB:print:UIControlsScaled');
      for i=1:2
        fname = fullfile(job.outdir{1},[job.fname fignames{i} '.png']);
        print(figuresids{i}, '-dpng', '-opengl', sprintf('-r%d',dpi), fname);
      end
      warning('ON','MATLAB:print:UIControlsScaled');
    end


    %% close
    if job.save>1
      spm_figure('Clear','Graphics');
      for i=2:26, try, close(i); end; end
    end
  end
end

%-End
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj)
%-----------------------------------------------------------------------
% Old check worst function. The spm_input could be replaced by an popup 
% window. A specification of the data range would rather than the x worst 
% images would be useful. 
%-----------------------------------------------------------------------
global H

if isempty(spm_figure('FindWin','Interactive')), spm('createintwin'); end

n = length(H.V);
number = min([n 24]);
number = spm_input('How many files ?',1,'e',number);
number = min([number 24]);
number = min([number length(H.V)]);
  
ind_sorted_decreased = H.ind_sorted_display(n:-1:1);
list = char(H.V(ind_sorted_decreased).fname);
sample = H.sample(ind_sorted_decreased);
list2 = list(1:number,:);

if H.mesh_detected
  % display single meshes and correct colorscale of colorbar
  for i=1:number
    h = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.data(:,ind_sorted_decreased(i))));
    
    % shift each figure slightly
    if i==1
        pos = get(h.figure,'Position');
    else
        pos = pos - [20 20 0 0];
    end
    
    % remove menubar and toolbar, use filename as title
    set(h.figure,'MenuBar','none','Toolbar','none','Name',sprintf('Sample %d: %s',sample(i),list2(i,:)),...
         'NumberTitle','off','Position',pos);
    cat_surf_render2('ColourMap',h,jet);
    cat_surf_render2('ColourBar',h,'on');
    cat_surf_render2('CLim',h,H.range_data98);
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
global H

  H.show_name = get(H.fnambox,'Value');
  show_boxplot;
  
return
        
%-----------------------------------------------------------------------
function checkbox_plot(obj, event_obj)
%-----------------------------------------------------------------------
  global H
  
  H.show_violin = get(H.plotbox,'Value');
  show_boxplot;
return

%-----------------------------------------------------------------------
function show_IQRratio(X)
%-----------------------------------------------------------------------
global H

% close second render window which is not needed here
if isfield(H,'hy') && isgraphics(H.hy.figure)
  close(H.hy.figure)
end

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',H.pos.scat,'Parent',H.figure);

cmap = [jet(64); gray(64)];

% estimate ratio between weighted overall quality (IQR) and mean corr. 
H.IQRratio = (X(:,2)/std(X(:,2)))./(X(:,1)/std(X(:,1)));

% because we use a splitted colormap we have to set the color
% values explicitely
IQRratio_scaled = 63*(H.IQRratio-min(H.IQRratio))/(max(H.IQRratio)-min(H.IQRratio));
C = zeros(length(H.IQRratio),3);
for i=1:length(H.IQRratio)
  C(i,:) = cmap(round(IQRratio_scaled(i))+1,:);
end
scatter(X(:,1),X(:,2),30,C,'o','Linewidth',2);

xlabel('<----- Worst ---      Mean correlation      --- Best ------>  ','FontSize',H.FS-1,'FontWeight','Bold');
ylabel('<----- Best ---      Weighted overall image quality (IQR)      --- Worst ------>  ','FontSize',H.FS-1,'FontWeight','Bold');
title('<--- Smallest -- Norm. Ratio IQR/Mean Correlation -- Largest ---->  ','FontSize',H.FS+1,'FontWeight','Bold');

% add colorbar
H.cbar = axes('Position',H.pos.cbar+[0 0.9 0 0],'Parent',H.figure);
image((1:64));

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
  round(100*linspace(min(H.IQRratio),max(H.IQRratio),5))/100,'TickLength',[0 0]);

% update index of worst files
[tmp, H.ind_sorted_display] = sort(H.IQRratio,'ascend');

colormap(cmap)

H.isscatter = 1;

return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
global H

% get sorting order
H.sorted = order;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',H.pos.corr,'Parent',H.figure);
names = spm_file(char(H.filename.m(:)),'short40');
group = H.sample;
if H.sorted
  names = names(H.ind_sorted_display,:);
  group = group(H.ind_sorted_display);
end
cat_plot_cov(data,struct('ax',H.ax,'pos_cbar',H.pos.cbar,'name',names,'fontsize',H.FS,'group',group))
t = title('Sample Correlation');
set(t,'Fontsize',H.FS+2)
cmap = [jet(64); gray(64)];
colormap(cmap)

H.isscatter = 0;

% remove sliders and text
try
  if ~H.mesh_detected
    set(H.mm,'Visible','off'); 
    set(H.mm_txt,'Visible','off');
    set(H.alpha,'Visible','off');
    set(H.alpha_txt,'Visible','off');
  end
end

return

%-----------------------------------------------------------------------
function show_boxplot(data_boxp, name_boxp, quality_order)
%-----------------------------------------------------------------------
global H

if nargin == 0
  data_boxp = H.bp.data;
  name_boxp = H.bp.name;
  quality_order = H.bp.order;
end

H.Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',H.Fgraph);
set(H.Fgraph,'Renderer','OpenGL');

n_samples = max(H.sample);

xpos = cell(1,n_samples);
data = cell(1,n_samples);

hold on
allow_violin = 1;
for i=1:n_samples
  ind = find(H.sample == i);
  if length(ind)<10
    allow_violin = 0; 
    H.show_violin = 0;
  end
  data{i} = data_boxp(ind);
  
  if n_samples == 1
    xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
  else
    xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
  end

  for j=1:length(ind)
    if H.show_name
      text(xpos{i}(j),data{i}(j),H.filename.m{ind(j)},'FontSize',H.FS-2,'HorizontalAlignment','center')
    else
      plot(xpos{i}(j),data{i}(j),'k.');
    end
  end
end

H.fnambox = uicontrol(H.figure,...
    'string','Show filenames','Units','normalized',...
    'position',H.pos.fnambox,'callback',@checkbox_names,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Show filenames in boxplot','value',H.show_name,...
    'BackgroundColor',[0.8 0.8 0.8],...
    'Interruptible','on','Visible','on','FontSize',H.FS-2);

% allow violin plot onl if samples are all large enough
if allow_violin
  H.plotbox = uicontrol(H.figure,...
    'string','Violinplot','Units','normalized',...
    'position',H.pos.plotbox,'callback',@checkbox_plot,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Switch to Violinplot','value',H.show_violin,...
    'BackgroundColor',[0.8 0.8 0.8],...
    'Interruptible','on','Visible','on','FontSize',H.FS-2);
end

% colormap for samples
if exist('lines')
  cm = lines(n_samples);
else
  cm = jet(n_samples);
end
opt = struct('groupnum',0,'ygrid',0,'violin',2*H.show_violin,'median',2,'groupcolor',cm);
ylim_add = 0.075;

cat_plot_boxplot(data,opt);

set(gca,'XTick',[],'XLim',[-.25 n_samples+1.25]);
if max(data_boxp) > min(data_boxp)
  yamp = max(data_boxp) - min(data_boxp);
  ylim_min = min(data_boxp) - ylim_add*yamp;
  ylim_max = max(data_boxp) + ylim_add*yamp;
  set(gca,'YLim',[ylim_min ylim_max]);
end

% add colored labels and title
if n_samples > 1
  [tmp,  tmp2] = spm_str_manip(char(H.filename.s),'C');
  title_str = sprintf('Boxplot: %s  \n%s ',name_boxp, strrep(tmp,tmp2.s,''));
  fprintf('\nCommon filename: %s\n',tmp);
else
  title_str = sprintf('Boxplot: %s  \nCommon filename: %s*',name_boxp,spm_file(char(H.filename.s),'short25'));
end
title(title_str,'FontSize',H.FS-1,'FontWeight','Bold');
xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',H.FS+1,...
    'FontWeight','Bold');

xpos = -0.40 - n_samples*0.1;

if (length(data_boxp) > 2)
  if quality_order > 0
    text(xpos, ylim_min,'<----- Low rating (poor quality)  ','Color','red','Rotation',...
        90,'HorizontalAlignment','left','FontSize',H.FS,'FontWeight','Bold')
    text(xpos, ylim_max,'High rating (good quality) ------>  ','Color','green','Rotation',...
        90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
  else
      text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
          90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
      text(xpos, ylim_min,'<----- High rating (good quality) ','Color','green','Rotation',...
          90,'HorizontalAlignment','left','FontSize',H.FS,'FontWeight','Bold')
  end
  text(xpos, (ylim_max+ylim_min)/2,sprintf('%s',name_boxp),'Color','black','Rotation',...
        90,'HorizontalAlignment','center','FontSize',H.FS,'FontWeight','Bold')
end

hold off

% estimate sorted index new for displaying worst files
if quality_order > 0
  [tmp, H.ind_sorted_display] = sort(data_boxp,'descend');
else
  [tmp, H.ind_sorted_display] = sort(data_boxp,'ascend');
  
end

H.bp = struct('data',data_boxp,'name',name_boxp,'order',quality_order);

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H,'alpha')
  H.alphaval = get(H.alpha,'Value');
else
  H.alphaval = 0.5;
end

% display image with 2nd colorbar (gray)
image(65 + H.img);
if ~H.mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

% prepare alpha overlays for red and green colors
if H.alphaval > 0
  % get 2%/98% ranges of difference image
  range = cat_vol_iscaling(H.img_alpha(:),[0.02 0.98]);

  hold on
  alpha_g = cat(3, zeros(size(H.img_alpha)), H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)));
  alpha_r = cat(3, H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
  hg = image(alpha_g); set(hg, 'AlphaData', H.img_alpha.*(H.img_alpha>range(2)),'AlphaDataMapping','scaled')
  if ~H.mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-H.img_alpha.*(H.img_alpha<range(1)),'AlphaDataMapping','scaled')
  if ~H.mesh_detected, axis image; end
  hold off
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H,'mm')
  slice_mm = get(H.mm,'Value');
else
  slice_mm = 0;
end

if H.names_changed
  P = H.Vchanged;
else
  P = H.V;
end

vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
Orig = P(1).mat\[0 0 0 1]';
sl   = round(slice_mm/vx(3)+Orig(3));

% if slice is outside of image use middle slice
if (sl>P(1).dim(3)) || (sl<1)
  sl = round(P(1).dim(3)/2);
end

M  = spm_matrix([0 0 sl]);
H.data_diff = H.data;

for i = 1:length(H.V)
  H.img = single(spm_slice_vol(P(i),M,P(1).dim(1:2),[1 0]));
  H.img(isnan(H.img)) = 0;
  
  % rescue unscaled data
  H.data_diff(:,:,i) = H.img;

  % scale image according to mean
  H.data(:,:,i) = H.img/mean(H.img(H.img ~= 0));
end

% calculate individual difference to mean image
for i=1:size(H.data_diff,3)
  H.data_diff(:,:,i) = H.data_diff(:,:,i) - mean(H.data_diff,3);
end

% enhance contrast and scale image to 0..64
mn = min(H.data(:));
mx = max(H.data(:));
H.data = 64*((H.data - mn)/(mx-mn));

if H.sorted
  if isfield(H.pos,'x')
    x = H.ind_sorted(H.pos.x);
    if ~H.isscatter
      y = H.ind_sorted(H.pos.y);
    end
  end
else
  if isfield(H.pos,'x')
    x = H.pos.x;
    if ~H.isscatter
      y = H.pos.y;
    end
  end
end

% check whether mouse position is defined
if isfield(H.pos,'x')
  if H.isscatter
    H.img       = H.data(:,:,x)';
    H.img_alpha = H.data_diff(:,:,x)';
  else
    H.img       = [H.data(:,:,y) H.data(:,:,x)]';
    H.img_alpha = [H.data_diff(:,:,y) H.data_diff(:,:,x)]';
  end
  
  % correct orientation
  H.img = rot90(H.img,2);
  H.img_alpha = rot90(H.img_alpha,2);
  
  % use gray scale colormap for values > 64
  axes('Position',H.pos.slice);
  image(65 + H.img);
  axis image
  set(gca,'XTickLabel','','YTickLabel','');
  
  % prepare alpha overlays for red and green colors
  if H.alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(H.img_alpha(:),[0.02 0.98]);

    hold on
    alpha_g = cat(3, zeros(size(H.img_alpha)), H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)));
    alpha_r = cat(3, H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
    hg = image(alpha_g); set(hg, 'AlphaData', H.img_alpha.*(H.img_alpha>range(2)),'AlphaDataMapping','scaled')
    axis image
    hr = image(alpha_r); set(hr, 'AlphaData',-H.img_alpha.*(H.img_alpha<range(1)),'AlphaDataMapping','scaled')
    axis image
    hold off
  end
  
  if H.isscatter
    txt = {sprintf('%s',spm_file(H.filename.m{x},'short25')),[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  else
    txt = {sprintf('Correlation: %3.3f',H.YpY(x,y)),[],['Top: ',spm_file(H.filename.m{x},'short25')],...
      ['Bottom: ',spm_file(H.filename.m{y},'short25')],[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  end
  set(H.text,'String',txt,'FontSize',H.FS-2);
  set(H.mm_txt,'String',[num2str(round(get(H.mm,'Value'))),' mm'],...
      'FontSize',H.FS-2);
end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global H

if ~H.mesh_detected
  set(H.alpha,'Visible','on');
  set(H.alpha_txt,'Visible','on');
end

pos_mouse = get(event_obj, 'Position');

if H.isscatter
  H.pos.x = find(H.X(:,1) == pos_mouse(1));
  if isempty(H.pos.x)
    H.pos.x = find(H.X(:,2) == pos_mouse(2));
  end

  % text info for data cursor window
  txt = {sprintf('%s',H.filename.m{H.pos.x})};

  % text info for textbox
  txt2 = {[],sprintf('%s',spm_file(H.filename.m{H.pos.x},'short25')),[],'Difference to Sample Mean','(red: - green: +)'};

  set(H.text,'String',txt2,'FontSize',H.FS-2);

  if ~H.mesh_detected
    axes('Position',H.pos.slice);
  end
  x = H.pos.x;
else
  % check for valid mouse position
  if pos_mouse(1) > pos_mouse(2) || pos_mouse(1)>length(H.sample) || pos_mouse(2)>length(H.sample)
    txt = {''};
    return
  end

  % save position of mouse
  H.pos.x = pos_mouse(1);
  H.pos.y = pos_mouse(2);

  if H.sorted
    if isfield(H.pos,'x')
      x = H.ind_sorted(H.pos.x);
      y = H.ind_sorted(H.pos.y);
    end
  else
    if isfield(H.pos,'x')
      x = H.pos.x;
      y = H.pos.y;
    end
  end

  % text info for data cursor window
  if H.mesh_detected
    txt = {sprintf('Correlation: %3.3f',H.YpY(x,y)),['Left: ',H.filename.m{x}],...
      ['Right: ',H.filename.m{y}]};
  else
    txt = {sprintf('Correlation: %3.3f',H.YpY(x,y)),['Top: ',H.filename.m{x}],...
      ['Bottom: ',H.filename.m{y}]};
  end

  % text info for textbox
  if H.mesh_detected
    txt2 = {[],sprintf('Correlation: %3.3f',H.YpY(x,y)),[],'right (1st row) and left (2nd row) hemisphere',['Left: ',...
      spm_file(H.filename.m{x},'short25') '     Right: ',spm_file(H.filename.m{y},'short25')],...
      [],'Difference to Sample Mean', '(red: - green: +)'};
  else
    txt2 = {[],sprintf('Correlation: %3.3f',H.YpY(x,y)),['Top: ',...
      spm_file(H.filename.m{x},'short25')],['Bottom: ',spm_file(H.filename.m{y},'short25')],...
      [],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm'],...
      'Difference to Sample Mean', '(red: - green: +)'};
  end      

  set(H.text,'String',txt2,'FontSize',H.FS-2);
  
  if ~H.mesh_detected
    axes('Position',H.pos.slice);
  end
end

% show two render views for meshes
if H.mesh_detected 
  if isfield(H,'hx') && isgraphics(H.hx.figure)
    H.hx = cat_surf_render2('Overlay',H.hx,H.data(:,x));
  else
    H.hx = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.data(:,x)));
    H.hx = cat_surf_render2('Colourbar',H.hx);
  end
  H.hx = cat_surf_render2('clim',H.hx,H.range_data98);
  pos_Hy = get(H.hx.figure,'Position');
  pos_Hy = pos_Hy + [pos_Hy(3)+30 0 0 0];
  set(H.hx.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Name',sprintf('Sample %d: %s',H.sample(x),H.filename.m{x}))
  figure(H.hx.figure)
  
  if ~H.isscatter % don't show that for scatter plot
    if isfield(H,'hy') && isgraphics(H.hy.figure)
      H.hy = cat_surf_render2('Overlay',H.hy,H.data(:,y));
    else
      H.hy = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.data(:,y)));
      H.hy = cat_surf_render2('Colourbar',H.hy);
      set(H.hy.figure,'Position',pos_Hy);  
    end
    H.hy = cat_surf_render2('clim',H.hy,H.range_data98);
    set(H.hy.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Name',sprintf('Sample %d: %s',H.sample(y),H.filename.m{y}));  
    figure(H.hy.figure)
  end
else
  % add slider for colume data
  set(H.mm,'Visible','on');
  set(H.mm_txt,'Visible','on');
  if H.isscatter
    H.img = H.data(:,:,x)';
    % alpha overlay
    H.img_alpha = H.data_diff(:,:,x)';
  else
    H.img = [H.data(:,:,y) H.data(:,:,x)]';
    % alpha overlay
    H.img_alpha = [H.data_diff(:,:,y) H.data_diff(:,:,x)]';
  end
  
  % correct orientation
  H.img = rot90(H.img,2);
  H.img_alpha = rot90(H.img_alpha,2);

  % display image with 2nd colorbar (gray)
  image(65 + H.img);
  if ~H.mesh_detected, axis image; end
  set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

  % prepare alpha overlays for red and green colors
  if H.alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(H.img_alpha(:),[0.02 0.98]);

    hold on
    alpha_g = cat(3, zeros(size(H.img_alpha)), H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)));
    alpha_r   = cat(3, H.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
    hg = image(alpha_g); set(hg, 'AlphaData', H.img_alpha.*(H.img_alpha>range(2)),'AlphaDataMapping','scaled')
    if ~H.mesh_detected, axis image; end
    hr = image(alpha_r); set(hr, 'AlphaData',-H.img_alpha.*(H.img_alpha<range(1)),'AlphaDataMapping','scaled')
    if ~H.mesh_detected, axis image; end
    hold off
  end
end

return
