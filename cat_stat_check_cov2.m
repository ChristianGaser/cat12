function varargout = cat_stat_check_cov2(job)
%cat_stat_check_cov. To check Z-score across sample.
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%
% varargout = cat_stat_check_cov2(job)
%  
% job                .. SPM job structure
%  .data             .. volume and surface files
%  .globals          .. optionally correct TIV using global scaling (for VBM only)
%  .gSF              .. global scaling values
%  .c                .. confounds
%  .data_xml         .. optional xml QC data
%  .verb             .. print figures
%  .new_fig          .. use new window instead of SPM Fgraph
%  .SPM              .. SPM design matrix structure
%
% varargout          .. output structure 
%  .zscore           .. absolute mean Z-score
%  .table            .. Z-score table
%  .sorttable        .. sorted Z-score table
%  .threshold_zsc    .. lower threshold for Z-score (mean - 4*std)
%
% Example: 
%   cat_stat_check_cov2(struct('data',{{ files }} ,'gap',3,'c',[],'data_xml',{{}}));
%
% See also cat_stat_check_cov
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

clearvars -GLOBAL H;       % clear old
global H

% show data by fileorder
H.ui.sorted = false;

H.job = job;

if nargin == 0
  error('No argument given.');
end

H.sample   = [];
H.mouse.x  = 1;
G          = [];
n_subjects = 0;
H.ui.alphaval = 0.5;
H.names_changed = false;

if ~isfield(job,'verb')
  job.verb = true;
end

if isfield(job,'SPM')
  H.SPM = job.SPM;
end

% use this window for Fgraph
if isfield(job,'new_fig') && job.new_fig
  H.Fgraph = spm_figure('Create','Check','Check Z-score');
end

if isfield(job,'show_violin')
  H.ui.show_violin = job.show_violin;
else
  H.ui.show_violin = false;
end

if isfield(job,'show_name')
  H.ui.show_name = job.show_name;
else
  H.ui.show_name = false;
end

% consider new options, but keep compatibility
if isfield(job,'sel_xml')
  if isfield(job.sel_xml,'data_xml')
    job.data_xml = job.sel_xml.data_xml;
  else
    job.data_xml = '';
  end
end

% this function is also used for the longitudinal, thus we also have to
% test the thickness field
if ~spm_mesh_detect(char(job.data{1}(1,:))) && isempty(strfind(char(job.data{1}(1,:)),'thickness'))
  H.mesh_detected = false;
else
  H.mesh_detected = true;
end

% read filenames for each sample and indicate sample parameter
if H.mesh_detected
  n_samples = numel(job.data);
  sinfo = cat_surf_info(char(job.data{1}(1,:)));
  H.Pmesh = gifti(sinfo.Pmesh);
  for i=1:n_samples
    [pp,ff,ee] = spm_fileparts(char(job.data{i}(1,:))); 
    if any( ~isempty( strfind({'lh.thickness' },[ff ee]) ) ) && ~strcmp(ee,'.gii')
      %% native longitudinal surface
      sdata = gifti(fullfile(pp,[strrep(ff,'lh.thickness','lh.central') ee '.gii'])); 
      cdata = single(cat_io_FreeSurfer('read_surf_data',job.data{i})); 
      gdata = gifti(struct('vertices',sdata.vertices,'faces',sdata.faces,'cdata',cdata)); 
      V0 = struct('fname',job.data{i},'dim',size(cdata),'dt',[16 0], ...
             'pinfo',[1 0 0],'mat',eye(4),'n',[1 1],'descript','GMT'); 
      V0.private = gdata; 
    else
      V0 = spm_data_hdr_read(char(job.data{i}));
    end
    n_subjects = n_subjects + length(V0);
      
    if i==1, H.files.V = V0;
    else,    H.files.V = [H.files.V; V0]; end
    H.sample = [H.sample, i*ones(1,size(job.data{i},1))];
  end
  H.files.fname = cellstr({H.files.V.fname}'); 
  H.info = cat_surf_info(H.files.fname{1});
else
  n_samples = numel(job.data);
  for i=1:n_samples
    
    if size(job.data{i},1) == 1 % 4D data
      [pth,nam,ext] = spm_fileparts(char(job.data{i}));
      % remove ",1" at the end
      job.data{i} = fullfile(pth,[nam ext]);
    end
    
    V0 = nifti(char(job.data{i}));
    n_subjects = n_subjects + length(V0);
      
    if i==1, H.files.V = V0;
    else,    H.files.V = [H.files.V V0]; end

    H.sample = [H.sample, i*ones(1,length(V0))];
  end
  
  % we need that field to be comparable to V of mesh-structure
  H.files.fname = cellstr({H.files.V.dat.fname}');   
end

H.ind = true(1,n_subjects);

% give error because this is not yet tested
if isfield(job,'gSF') && numel(job.gSF) == n_subjects
  fprintf('Use global scaling from design matrix.\n');
  gSF = job.gSF;
end

% check for global scaling with TIV
if isfield(job,'globals') && job.globals
  if H.mesh_detected
    is_gSF = false;
    fprintf('Disabled global scaling with TIV, because this is not meaningful for surface data.\n');
  else
    is_gSF = true;
    gSF = ones(n_subjects,1);
  end
else
  is_gSF = false;
end

if isfield(job,'c') && ~isempty(job.c)
  for i=1:numel(job.c)
    G = [G job.c{i}];
  end
end

if isempty(char(job.data_xml))
  H.isxml = false;
  xml_defined = false;
  H.xml.QM_names = '';
  xml_files = [];
else
  xml_files = char(job.data_xml);
  if size(xml_files,1) ~= n_subjects
    fprintf('Only %d of %d report files were defined. Try to find xml-files for quality measures.\n',size(xml_files,1),n_subjects);
    H.isxml = false;
    xml_defined = false;
  else
    H.isxml = true;
    xml_defined = true;
  end
end

if H.mesh_detected
  H.xml.QM = ones(n_subjects,5);
  H.xml.QM_names = char('Noise','Bias','Weighted overall image quality (IQR)','Euler number','Size of topology defects');
  H.xml.QM_names_short = char('Noise & Mean absolute Z-score','Bias & Mean absolute Z-score','Weighted IQR & Mean absolute Z-score','Euler number & Mean absolute Z-score','Size of topology defects & Mean absolute Z-score');
else
  H.xml.QM = ones(n_subjects,3);
  H.xml.QM_names = char('Noise','Bias','Weighted overall image quality (IQR)');
  H.xml.QM_names_short = char('Noise & Mean absolute Z-score','Bias & Mean absolute Z-score','Weighted IQR & Mean absolute Z-score');
end

pth = spm_fileparts(H.files.fname{1});
if isfield(job,'sel_xml') && isfield(job.sel_xml,'select_dir')
  report_folder = char(job.sel_xml.select_dir);
else
  report_folder = fullfile(spm_fileparts(pth),'report');
end

% check whether report subfolder exists
if ~exist(report_folder,'dir')
  report_folder = pth;
end

% search xml report files if not defined
prep_str = '';
if ~xml_defined
  fprintf('Search xml-files ');
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
        Vfname = H.files.fname{i};
        ind = strfind(Vfname,fname);
        if ~isempty(ind)
          [pth, prep_str] = spm_fileparts(Vfname(1:ind-1));
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
  fprintf('\n');
end

n_xml_files = 0;
spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')

fprintf('Load xml-files ');
for i=1:n_subjects
  % get basename for data files
  [pth, data_name ee] = fileparts(H.files.fname{i});
  if ~strcmp(ee,'.nii') && ~strcmp(ee,'.gii'), data_name = [data_name ee]; end
  
  % remove ending for rigid or affine transformed files
  data_name = strrep(data_name,'_affine','');
  data_name = strrep(data_name,'_rigid','');
  
  % use xml-file if found by name
  if ~xml_defined

    % remove prep_str from name and use report folder and xml extension
    subjname = strrep(data_name,prep_str,'');
    % for meshes we also have to remove the additional "." from name
    if H.mesh_detected
      subjname = subjname(2:end);
    end
    xml_file = fullfile(report_folder,['cat_' subjname '.xml']);
  else % use defined xml-files

    [pth, subjname] = fileparts(deblank(xml_files(i,:)));
    % remove leading 'cat_'
    subjname = subjname(5:end);

    % check for filenames
    if isempty(strfind(data_name,subjname))
      fprintf('\nSkip use of xml-files for quality measures because of deviating subject names:\n%s vs. %s\n',H.files.fname{i},xml_files(i,:));
      H.isxml = false;
      break
    end
    xml_file = deblank(xml_files(i,:));
  end
  
  % get mri folder
  if strcmp(ee,'.gii'), mri_folder = fullfile(fileparts(pth),'mri');
  else mri_folder = pth; end
  
  % find raw/p0 files
  H.files.raw{i} = fullfile(fileparts(pth),[subjname ee]);
  H.files.p0{i}  = fullfile(mri_folder,['p0' subjname '.nii']);
  if ~exist(H.files.raw{i},'file'), H.files.raw{i} = ''; end
  if ~exist(H.files.p0{i}, 'file'), H.files.p0{i}  = ''; end

  if exist(xml_file,'file')
    H.job.data_xml{i} = xml_file;
    xml = cat_io_xml(xml_file);
    n_xml_files = n_xml_files + 1;
    H.isxml = true;
    
    % find jpg/pdf/log files
    H.files.jpg{i} = fullfile(report_folder,['catreportj_' subjname '.jpg']);
    H.files.pdf{i} = fullfile(report_folder,['catreport_' subjname '.pdf']);
    H.files.log{i} = fullfile(report_folder,['catlog_' subjname '.txt']);
    if ~exist(H.files.jpg{i},'file'), H.files.jpg{i} = ''; end
    if ~exist(H.files.pdf{i},'file'), H.files.pdf{i} = ''; end
    if ~exist(H.files.log{i},'file'), H.files.log{i} = ''; end
     
    % get TIV
    if is_gSF && isfield(xml,'subjectmeasures') && isfield(xml.subjectmeasures,'vol_TIV')
      gSF(i) = xml.subjectmeasures.vol_TIV;
    else
      is_gSF = false;
    end
    
  else
    if is_gSF
      fprintf('\nFile %s not found. Skip use of xml-files for quality measures and TIV.\n',xml_file);
    else
      fprintf('\nFile %s not found. Skip use of xml-files for quality measures.\n',xml_file);
    end
    fprintf('Please select xml-files manually.\n\n');
    H.isxml = false;
    is_gSF  = false;
    break
  end

  if ~isfield(xml,'qualityratings') && ~isfield(xml,'QAM')
    fprintf('\nQuality rating is not saved for %s. Report file %s is incomplete.\nPlease repeat preprocessing amd check for potential errors in the ''err'' folder.\n',H.files.fname{i},xml_files(i,:));  
    H.isxml = false;
    break
  end

  if H.mesh_detected
    if isfield(xml.qualityratings,'NCR')
    % check for newer available surface measures
      if isfield(xml.subjectmeasures,'EC_abs')
        H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR xml.subjectmeasures.EC_abs xml.subjectmeasures.defect_size];
      else
        H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR NaN NaN];
      end
    else % also try to use old version
      H.xml.QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
    end
  else
    if isfield(xml.qualityratings,'NCR')
      H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR];
    else % also try to use old version
      H.xml.QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
    end
  end
  spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');
fprintf('\n');

if H.isxml
  if n_xml_files ~= n_subjects
    fprintf('Only %d of %d report files found. Skip use of xml-files for quality measures.\n',n_xml_files,n_subjects);
    H.isxml = false;
  else
    fprintf('%d report files with quality measures were found.\n',n_xml_files);
  end
end

% remove last two columns if EC_abs and defect_size are not defined
if H.mesh_detected && all(isnan(H.xml.QM(:,4))) && all(isnan(H.xml.QM(:,5)))
  H.xml.QM = H.xml.QM(:,1:3);
  H.xml.QM_names = H.xml.QM_names(1:3,:);
end

H.data.Ymean = 0.0;
Yss   = 0.0; % sum of squares

% preload surface data for later render view
if H.mesh_detected
  % load surface texture data
  H.texture = single(spm_data_read(H.files.V));
end

spm_progress_bar('Init',n_subjects,'Load data','subjects completed')

fprintf('Load data ');
for i = 1:n_subjects
  fprintf('.');
  if H.mesh_detected
    tmp = spm_data_read(H.files.V(i));
  else
    tmp(:,:,:) = H.files.V(i).dat(:,:,:);
  end
  tmp(isnan(tmp)) = 0;
  
  % either global scaling was externally defined using job or values were
  % used from xml-file
  if is_gSF || isfield(job,'gSF')
    tmp = tmp*gSF(i)/mean(gSF);
  end
  
  if i>1 && numel(H.data.Ymean) ~= numel(tmp)
    fprintf('\n\nERROR: File %s has different data size: %d vs. %d\n\n',job.data{i},numel(H.data.Ymean),numel(tmp));
    return
  end
  
  H.data.Ymean = H.data.Ymean + tmp(:);
  Yss   = Yss + tmp(:).^2;
  spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');

% get mean and SD
H.data.Ymean = H.data.Ymean/n_subjects;

% get range 10..98%
H.data.range98 = cat_vol_iscaling(H.data.Ymean(H.data.Ymean~=0),[0.10 0.98]);
H.data.global = mean(H.data.Ymean(H.data.Ymean~=0))/8; H.data.global = mean(H.data.Ymean(H.data.Ymean>H.data.global));

H.data.global = 0.25*H.data.global;

% we have sometimes issues with number precision
Yvar   = 1.0/(n_subjects-1)*(Yss - n_subjects*H.data.Ymean.*H.data.Ymean);
Yvar(Yvar<0) = 0;
H.data.Ystd   = sqrt(Yvar);

% only consider non-zero areas for Ystd and threshold Ymean
ind = H.data.Ystd ~= 0 & H.data.Ymean > H.data.global;

spm_progress_bar('Init',n_subjects,'Calculate Z-score','subjects completed')

H.data.avg_abs_zscore = zeros(n_subjects,1);
for i = 1:n_subjects
  fprintf('.');
  
  if H.mesh_detected
    tmp = spm_data_read(H.files.V(i));
  else
    tmp(:,:,:) = H.files.V(i).dat(:,:,:);
  end
  tmp(isnan(tmp)) = 0;
  
  if is_gSF
    tmp = tmp*gSF(i)/mean(gSF);
  end
  
  % calculate Z-score  
  zscore = (tmp(ind) - H.data.Ymean(ind))./H.data.Ystd(ind);
  
  % and use mean of Z-score as overall measure
  H.data.avg_abs_zscore(i) = mean(abs(zscore));
  spm_progress_bar('Set',i);  
end
fprintf('\n');
spm_progress_bar('Clear');

% voxelsize and origin of volume data
if ~H.mesh_detected
  H.data.vx =  sqrt(sum(H.files.V(1).mat(1:3,1:3).^2));
  H.data.Orig = H.files.V(1).mat\[0 0 0 1]';
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

popb = [0.038 0.035];  % size of the small buttons
popm = 0.780;          % x-position of the control elements

H.pos = struct(...
    'fig',    [10 10 1.3*ws(3) 1.1*ws(3)],...% figure
    'cbar',   [0.045 0.050 0.700 0.020],...  % colorbar for figure
    'plot',   [0.050 0.050 0.700 0.825],...  % scatter plot
    ...
    'close',  [0.775 0.935 0.100 0.040],...  % close button
    'show',   [0.875 0.935 0.100 0.040],...  % button to show worst cases
    'scat',   [0.772 0.880 0.110 0.050],...  % button to enable ordered matrix
    'boxp',   [0.872 0.880 0.110 0.050],...  % button to display boxplot
    ...
    ... == navigation unit ==
    'scSelect',    [popm+popb(1)*0 0.835 popb],... % select (default) 
    'scZoomReset', [popm+popb(1)*1 0.835 popb],... % standard zoom
    'scZoomIn',    [popm+popb(1)*2 0.835 popb],... % zoom in 
    'scZoomOut',   [popm+popb(1)*3 0.835 popb],... % zoom out
    'scPan',       [popm+popb(1)*4 0.835 popb],... % pan (moving hand)
    ...
    ... == remove unit ==
    'rmDel',       [popm+popb(1)*0 0.775 popb],... % delete 
    'rmUndo',      [popm+popb(1)*1 0.775 popb],... % undo deletion
    'rmNew',       [popm+popb(1)*2 0.775 popb],... % calculate new
    'rmListNew',   [popm+popb(1)*3 0.775 popb],... % list remaining data
    'rmListDel',   [popm+popb(1)*4 0.775 popb],... % list deleted data
    ...
    ... == display unit ==
    'dpReport',    [popm+popb(1)*0 0.715 popb],... % report 
    'dpRaw',       [popm+popb(1)*1 0.715 popb],... % raw data
    'dpRawP0',     [popm+popb(1)*2 0.715 popb],... % raw data + p0
    'dpLog',       [popm+popb(1)*3 0.715 popb],... % log
    ...
    ... == check boxes ==
    'fnambox',[0.775 0.600 0.200 0.050],... % show filenames?
    'plotbox',[0.875 0.600 0.200 0.050],... % switch between boxplot and violin plot 
    ...
    ... == slice display ==
    'text',   [0.775 0.450 0.200 0.150],... % textbox with info
    'aslider',[0.775 0.405 0.200 0.040],... % slider for alpha overlay
    'slice',  [0.775 0.030 0.200 0.400],... % slice images according to position of mouse pointer
    'sslider',[0.775 0.020 0.200 0.040]);   % slider for z-slice   

if ~H.mesh_detected
  % correct filenames for 4D data
  if strcmp(H.files.fname{1}, H.files.fname{2})
    H.names_changed = true;
    H.files.Vchanged = H.files.V;
    for i=1:n_subjects
      [pth,nam,ext] = spm_fileparts(H.files.fname{i});
      H.files.fname{i} = fullfile(pth, [nam sprintf('%04d',i) ext]);
    end
  end  
end

if job.verb, fprintf('\n'); end
fname_m = [];
fname_tmp = cell(n_samples,1);
fname_s   = cell(n_samples,1);
fname_e   = cell(n_samples,1);
for i=1:n_samples
  [tmp, fname_tmp{i}] = spm_str_manip(char(H.files.fname{H.sample == i}),'C');
  if ~isempty(fname_tmp{i})
    fname_m    = [fname_m; fname_tmp{i}.m]; 
    fname_s{i} = fname_tmp{i}.s;
    fname_e{i} = fname_tmp{i}.e;
  else
    fname_s{i} = '';
    fname_e{i} = '';
  end
  if job.verb
    fprintf('Compressed filenames sample %d: %s  \n',i,tmp);
  end
end

H.filename = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});

% sort data
[H.data.avg_abs_zscore_sorted, H.ind_sorted] = sort(H.data.avg_abs_zscore,'ascend');

threshold_zsc = mean(H.data.avg_abs_zscore) + 2*std(H.data.avg_abs_zscore);
n_thresholded = min(find(H.data.avg_abs_zscore_sorted > threshold_zsc));

if ~isempty(n_thresholded) && job.verb
  fprintf('\nThese data have a mean absolute Z-score above 2 standard deviations.\n');
  fprintf('This does not necessarily mean that you have to exclude these data. However, these data have to be carefully checked:\n');
  for i=n_thresholded:n_subjects
    fprintf('%s: %3.3f\n',H.files.fname{H.ind_sorted(i)},H.data.avg_abs_zscore_sorted(i));
  end
  fprintf('\n');
end

if nargout>0
  varargout{1} = struct('table',{[H.files.fname,num2cell(H.data.avg_abs_zscore)]},...
                        'sorttable',{[H.files.fname(H.ind_sorted),num2cell(H.data.avg_abs_zscore_sorted)]},...
                        'zscore',H.data.avg_abs_zscore,...
                        'threshold_zsc',threshold_zsc);
end


if job.verb
  
  show_menu;
  show_boxplot(H.data.avg_abs_zscore,'Mean absolute Z-score  ',-1);

  if isfield(job,'save') && job.save
    %% filenames
    if ~isempty(job.fname)
      dpi = cat_get_defaults('print.dpi'); 
      if isempty(dpi), dpi = 150; end

      fignames   = {'matrix','boxplot'};
      figuresids = {figure(2),H.Fgraph};
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
function show_menu
%-----------------------------------------------------------------------
global H

% create figure
H.mainfig = figure;
clf(H.mainfig);

set(H.mainfig,...
   'MenuBar','none',...
   'Position',H.pos.fig,...
...   'DefaultTextFontSize',H.FS,...
...   'DefaultUicontrolFontSize',H.FS,...
   'NumberTitle','off');

if H.mesh_detected
  set(H.mainfig,'Name','Click in image to display surfaces');
else
  set(H.mainfig,'Name','Click in image to display slices');
end

cm = datacursormode(H.mainfig);
set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
try
  set(cm,'Interpreter','none','NewDataCursorOnClick',false);
end

% add colorbar
H.ui.cbar = axes('Position',H.pos.cbar,'Parent',H.mainfig);
try, image(H.ui.cbar);end
set(get(H.ui.cbar,'children'),'HitTest','off','Interruptible','off');
set(H.ui.cbar,'Ytick',[],'YTickLabel',''); 

% create two colormaps
cmap = [jet(64); gray(64)];
colormap(cmap)

% display YTick with 5 values (limit accuracy for floating numbers)
%  set(H.ui.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
%    round(100*linspace(min(H.data.avg_abs_zscore),max(H.data.avg_abs_zscore),5))/100,'TickLength',[0 0]);

% add button for closing all windows
H.ui.close = uicontrol(H.mainfig,...
        'String','Close','Units','normalized',...
        'Position',H.pos.close,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'Callback','for i=2:26, try close(i); end; end;',...
        'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on');

% check button
H.ui.show = uicontrol(H.mainfig,...
        'String','Check worst','Units','normalized',...
        'Position',H.pos.show,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'Callback',@check_worst_data,...
        'ToolTipString','Display most deviating files',...
        'Interruptible','on','Enable','on');

% create popoup menu 
if H.isxml

  % estimate product between weighted overall quality (IQR) and mean absolute Z-score 
  H.X = [H.data.avg_abs_zscore H.xml.QM];
  H.xml.QMzscore = H.X(:,1).*H.X(:,2);

  str  = { 'Boxplot','Mean absolute Z-score',H.xml.QM_names,'Weighted IQR x Mean absolute Z-score'};

  if size(H.xml.QM,2) == 5
    tmp  = { {@show_boxplot, H.data.avg_abs_zscore, 'Mean absolute Z-score  ', -1},...
             {@show_boxplot, H.xml.QM(:,1), H.xml.QM_names(1,:), -1},...
             {@show_boxplot, H.xml.QM(:,2), H.xml.QM_names(2,:), -1},...
             {@show_boxplot, H.xml.QM(:,3), H.xml.QM_names(3,:), -1},...
             {@show_boxplot, H.xml.QM(:,4), H.xml.QM_names(4,:), -1},...
             {@show_boxplot, H.xml.QM(:,5), H.xml.QM_names(5,:), -1},...
             {@show_boxplot, H.xml.QMzscore, 'Weighted IQR x Mean absolute Z-score  ', -1} };
  else
    tmp  = { {@show_boxplot, H.data.avg_abs_zscore, 'Mean absolute Z-score  ', -1},...
             {@show_boxplot, H.xml.QM(:,1), H.xml.QM_names(1,:), -1},...
             {@show_boxplot, H.xml.QM(:,2), H.xml.QM_names(2,:), -1},...
             {@show_boxplot, H.xml.QM(:,3), H.xml.QM_names(3,:), -1},...
             {@show_boxplot, H.xml.QMzscore, 'Weighted IQR x Mean absolute Z-score  ', -1} };
  end

  % show IQR x mean score as default
  show_QMzscore(H.X,4);

else
  str  = { 'Boxplot','Mean absolute Z-score'};
  tmp  = { {@show_boxplot, H.data.avg_abs_zscore, 'Mean absolute Z-score  ', -1} };

  H.X = [H.data.avg_abs_zscore (1:numel(H.data.avg_abs_zscore))'];
  show_QMzscore(H.X,0);
end

H.ui.boxp = uicontrol(H.mainfig,...
        'String',str,'Units','normalized',...
        'Position',H.pos.boxp,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'Callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Display boxplot',...
        'Interruptible','on','Visible','on');

% if QM values are available allow IQR and surface parameters, but skip
% noise and bias as first 2 entries
if H.isxml
  str  = { 'Scatterplot',H.xml.QM_names_short(3:end,:)};
  tmp  = {{@show_QMzscore, H.X,4}}; % IQR
  if size(H.xml.QM,2) == 5
  tmp  = {{@show_QMzscore, H.X,4}, ... % IQR
          {@show_QMzscore, H.X,5}, ... % Euler number
          {@show_QMzscore, H.X,6}};    % size of topology defect
  end    
else
  str  = { 'Scatterplot','Mean absolute Z-score'};
  tmp  = {{@show_QMzscore, H.X,0}}; % just mean absolute Z-score with file order
end

H.ui.sort = uicontrol(H.mainfig,...
        'String',str,'Units','normalized',...
        'Position',H.pos.scat,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'Callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Sort matrix',...
        'Interruptible','on','Visible','on');

H.text = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.text,...
        'String',{'','Click in image to display slices'},...
        'Style','text','HorizontalAlignment','center',...
        'ToolTipString','Select slice for display',...
        'FontSize',H.FS-2);

%% == zoom unit ==
H.naviui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.scSelect(1) H.pos.scSelect(2)+0.035 0.2 0.02],...
  'String','Zoom options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

H.naviui.select = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scSelect,'callback',...
  ['datacursormode(''on''); global H; H.isdel = false; set(H.delui.remove,''BackGroundColor'',[0.94 0.94 0.94]);' ...
  'set(H.naviui.select,''BackGroundColor'',[0.95 0.95 0.95]);'],...
  'Style','Pushbutton','enable','on','ToolTipString','Data selection','CData',load_icon('tool_pointer.png'),'BackGroundColor',[0.95 0.95 0.95]);

H.naviui.zoomReset = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomReset,'callback',...
  ['zoom out; zoom reset; datacursormode(''on''); global H; H.isdel = false;' ...
  'set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'],...
  'Style','Pushbutton','enable','on','ToolTipString','Reset view','CData',load_icon('tool_fit.png')); 

H.naviui.zoomIn = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomIn,'callback',...
  ['global H; ' ...
   'hz = zoom(H.ax); ' ...
   'set(hz,''enable'',''on'',''direction'',''in''); set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'], ... 
  'Style','Pushbutton','enable','on','ToolTipString','Zoom in','CData',load_icon('tool_zoom_in.png'));

H.naviui.zoomOut = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomOut,'callback',...
  ['global H; ' ...
   'hz = zoom(H.ax); ' ...
   'set(hz,''enable'',''on'',''direction'',''out''); set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'], ...
  'Style','Pushbutton','enable','on','ToolTipString','Zoom out','CData',load_icon('tool_zoom_out.png'));

H.naviui.pan = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scPan,'Enable','off','callback',...
  'pan on; set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);',...
  'Style','Pushbutton','enable','on','ToolTipString','Hand','CData',load_icon('tool_hand.png'));

%% == remove unit ==
H.delui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.rmDel(1) H.pos.rmDel(2)+0.035 0.2 0.02],...
  'String','Data remove options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

H.delui.remove = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmDel,'callback',...
  ['datacursormode(''on''); global H; H.isdel = true; set(H.delui.remove,''BackGroundColor'',[0.95 0.95 0.95]);' ...
  'set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'],...
  'Style','Pushbutton','enable','on','ToolTipString','Remove data','CData',load_icon('file_delete.png'));

H.delui.undo = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmUndo,'callback',{@do_rerun,1},...
  'Style','Pushbutton','enable','off','ToolTipString','Undo deletion','CData',load_icon('file_delete_restore.png')); 

H.delui.new = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmNew,'callback',{@get_new_list,1},...
  'Style','Pushbutton','enable','off','ToolTipString','Re-calculate without removed data','CData',load_icon('greenarrowicon.png')); 

H.delui.list_del = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmListDel,'callback',{@get_new_list,-1},...
  'Style','Pushbutton','enable','off','ToolTipString','List deleted file','CData',load_icon('list_del.png')); 

H.delui.list_new = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmListNew,'callback',{@get_new_list,0},...
  'Style','Pushbutton','enable','off','ToolTipString','List remaining files','CData',load_icon('list_new.png')); 

%% == display unit ==
H.dpui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.dpReport(1) H.pos.dpReport(2)+0.035 0.2 0.02],...
  'String','Data display options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

% enable some buttons only if respective files are available
if H.isxml, status_xml = 'on';
else status_xml = 'off'; end
if ~isempty(H.files.raw{numel(H.sample)}), status_raw = 'on';
else status_raw = 'off'; end
if ~isempty(H.files.p0{numel(H.sample)}) && ~isempty(H.files.raw{numel(H.sample)}), status_rawp0 = 'on';
else status_rawp0 = 'off'; end
if ~isempty(H.files.log{numel(H.sample)}), status_log = 'on';
else status_log = 'off'; end
if ~isempty(H.files.jpg{numel(H.sample)}), status_report = 'on';
else status_report = 'off'; end

H.dpui.report = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpReport,'callback',{@show_report},...
  'Style','Pushbutton','enable',status_report,'ToolTipString','Show report file','CData',load_icon('file_cat_report.png'));

H.dpui.raw = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpRaw,'callback',{@show_raw,false},...
  'Style','Pushbutton','enable',status_raw,'ToolTipString','Show raw data','CData',load_icon('file_spm_view.png')); 

H.dpui.rawp0 = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpRawP0,'callback',{@show_raw,true},...
  'Style','Pushbutton','enable',status_rawp0,'ToolTipString','Show overlayed label','CData',load_icon('file_spm_view_p0.png')); 

H.dpui.log = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpLog,'callback',{@show_log},...
  'Style','Pushbutton','enable',status_log,'ToolTipString','Show log file','CData',load_icon('file_cat_log.png')); 

% add slider and opacity control only for volume data
if ~H.mesh_detected
  H.ui.alpha = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.aslider,...
        'Min',0,'Max',1,...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_alpha,'Value',0.5,...
        'ToolTipString','Change Opacity of pos. (green colors) and neg. (red colors) Z-scores',...
        'SliderStep',[0.01 0.1],'Visible','off');

  H.ui.alpha_txt = uicontrol(H.mainfig,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[H.pos.aslider(1) H.pos.aslider(2)-0.005 0.2 0.02],...
        'String','Overlay Opacity of Z-score',...
        'FontSize',H.FS-2,'Visible','off');

  H.ui.mm = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.sslider,...
        'Min',(1 - H.data.Orig(3))*H.data.vx(3),'Max',(H.files.V(1).dat.dim(3) - H.data.Orig(3))*H.data.vx(3),...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_slices_array,...
        'ToolTipString','Select slice for display',...
        'SliderStep',[0.005 0.05],'Visible','off');

  H.ui.mm_txt = uicontrol(H.mainfig,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[H.pos.sslider(1) H.pos.sslider(2)-0.005 0.2 0.02],...
        'String','Slice [mm]','Visible','off','FontSize',H.FS-2);

  update_slices_array;
end
  
return

%-----------------------------------------------------------------------
function icon = load_icon(name)
%-----------------------------------------------------------------------

icon = imread(fullfile(spm('dir'),'toolbox','cat12','html','icons',name)); 
icon = double(icon)./double(max(icon(:))); icon(icon==0) = 0.94; 

return

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj)
%-----------------------------------------------------------------------
% Old check worst function. The spm_input could be replaced by an popup 
% window. A specification of the data range would rather than the x worst 
% images would be useful. 
%-----------------------------------------------------------------------
global H

if isempty(spm_figure('FindWin','Interactive')), spm('createintwin'); end

n = length(H.files.V);
number = min([n 24]);
number = spm_input('How many files ?',1,'e',number);
number = min([number 24]);
number = min([number length(H.files.V)]);
  
ind_sorted_decreased = flipud(H.ind_sorted_display);

list = char(H.files.fname{ind_sorted_decreased});
sample = H.sample(ind_sorted_decreased);
list2 = list(1:number,:);

if H.mesh_detected
  % display single meshes and correct colorscale of colorbar
  for i=1:number
    h = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.texture(:,ind_sorted_decreased(i))));
    
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
    cat_surf_render2('CLim',h,H.data.range98);
  end
else
  spm_check_registration(list2);
  spm_orthviews('Resolution',0.2);
  set(H.ui.boxp,'Visible','on');
end
return

%-----------------------------------------------------------------------
function checkbox_names(obj, event_obj)
%-----------------------------------------------------------------------
global H

  H.ui.show_name = get(H.ui.fnambox,'Value');
  show_boxplot;
  
return
        
%-----------------------------------------------------------------------
function checkbox_plot(obj, event_obj)
%-----------------------------------------------------------------------
  global H
  
  H.ui.show_violin = get(H.ui.plotbox,'Value');
  show_boxplot;
return

%-----------------------------------------------------------------------
function show_QMzscore(X,sel)
%-----------------------------------------------------------------------
global H

H.sel = sel;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.mainfig);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',H.pos.plot,'Parent',H.mainfig);
axes(H.ax);

cmap = [jet(64); gray(64)];

% estimate product between QM-value and mean absolute Z-score
if sel
  H.xml.QMzscore = X(:,1).*X(:,sel);
else
  H.xml.QMzscore = X(:,1);
end

% get min/max in 0.25 steps
min_QMzscore = min(H.xml.QMzscore(H.ind)); min_QMzscore = floor(4*min_QMzscore)/4;
max_QMzscore = max(H.xml.QMzscore(H.ind)); max_QMzscore = ceil(4*max_QMzscore)/4;
if min_QMzscore == max_QMzscore
  max_QMzscore = max_QMzscore + 0.25;
end

% because we use a splitted colormap we have to set the color
% values explicitely
if sel == 4 % IQR
  QMzscore_scaled = 63*(H.xml.QMzscore-1)/2; % scale 1..3
else
  QMzscore_scaled = 63*(H.xml.QMzscore-min_QMzscore)/(max_QMzscore-min_QMzscore); % scale min..max
end

H.C = zeros(length(H.xml.QMzscore),3);
for i=1:length(H.xml.QMzscore)
  H.C(i,:) = cmap(round(QMzscore_scaled(i))+1,:);
end

% create marker for different samples
marker = char('o','s','d','^','v','<','>','+','*','-','|');
marker = char('o','+','*','.','s','d','^','v','<','>','-','|');
while max(H.sample) > marker, marker = [marker; marker]; end

if sel
  for i=1:max(H.sample)
    ind  = H.sample.*H.ind == i;
    H.ui.scatter = scatter(X(ind,sel),X(ind,1),30,H.C(ind,:),marker(i),'Linewidth',2);
    hold on
  end
  hold off
  xstr = sprintf('<----- Best ---      %s      --- Worst ------>  ',deblank(H.xml.QM_names(sel-1,:)));
else
  % we have to consider removed points
  x = 1:numel(X(:,1));
  if isfield(H,'del'), x(H.del) = []; end
  
  for i=1:max(H.sample)
    ind  = H.sample.*H.ind == i;
    H.ui.scatter = scatter(x(ind),X(ind,1),30,H.C(ind,:),marker(i),'Linewidth',2);
    hold on
  end
  hold off
  xstr = sprintf('<----- First ---      File      --- Last ------>  ');
  xlim([0 numel(H.sample)+1]);
end
xlabel(xstr,'FontSize',H.FS-1,'FontWeight','Bold');
ylabel('<----- Best ---      Mean absolute Z-score      --- Worst ------>  ','FontSize',H.FS-1,'FontWeight','Bold');

% add colorbar
H.ui.cbar = axes('Position',H.pos.cbar+[0 0.9 0 0],'Parent',H.mainfig);
image((1:64));

if sel
  xstr = sprintf('<----- Best ---      %s x mean absolute Z-score     --- Worst ------>  ',deblank(H.xml.QM_names(sel-1,:)));
else
  xstr = sprintf('<----- Best ---      mean absolute Z-score     --- Worst ------>  ');
end
title(xstr,'FontSize',H.FS+1,'FontWeight','Bold');

% display YTick with 5 values (limit accuracy for floating numbers)
if sel == 4 % IQR
  set(H.ui.cbar,'YTickLabel','','YTick','', 'XTick',linspace(1,64,5),'XTickLabel',...
    round(100*linspace(1,3,5))/100,'TickLength',[0 0]);
else
  set(H.ui.cbar,'YTickLabel','','YTick','', 'XTick',linspace(1,64,5),'XTickLabel',...
    round(100*linspace(min_QMzscore,max_QMzscore,5))/100,'TickLength',[0 0]);
end

% update index of worst files
[tmp, H.ind_sorted_display] = sort(H.xml.QMzscore(H.ind),'ascend');

colormap(cmap)

return

%-----------------------------------------------------------------------
function show_boxplot(data_boxp, name_boxp, quality_order)
%-----------------------------------------------------------------------
global H

if nargin == 0
  data_boxp = H.ui.bp.data;
  name_boxp = H.ui.bp.name;
  quality_order = H.ui.bp.order;
end

H.show_sel = 1;
set(H.dpui.report,'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.log,   'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.raw,   'BackGroundColor',[0.94 0.94 0.94]);

% only use SPM window if not defined
if ~isfield(H,'Fgraph')
  H.Fgraph = spm_figure('GetWin','Graphics');
end
spm_figure('Clear',H.Fgraph);
set(H.Fgraph,'Renderer','OpenGL');
figure(H.Fgraph);
spm_figure('Select',H.Fgraph);

n_samples = max(H.sample);

xpos = cell(1,n_samples);
data = cell(1,n_samples);

hold on
allow_violin = true;
for i=1:n_samples
  ind = find(H.sample(H.ind) == i);
  if length(ind)<10
    allow_violin = false; 
    H.ui.show_violin = false;
  end
  data{i} = data_boxp(ind);
  
  if n_samples == 1
    xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
  else
    xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
  end
end

H.ui.fnambox = uicontrol(H.mainfig,...
    'String','Show filenames','Units','normalized',...
    'Position',H.pos.fnambox,'callback',@checkbox_names,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Show filenames in boxplot','value',H.ui.show_name,...
    'BackgroundColor',[0.8 0.8 0.8],...
    'Interruptible','on','Visible','on','FontSize',H.FS-2);

% allow violin plot onl if samples are all large enough
if allow_violin
  H.ui.plotbox = uicontrol(H.mainfig,...
    'String','Violinplot','Units','normalized',...
    'Position',H.pos.plotbox,'callback',@checkbox_plot,...
    'Style','CheckBox','HorizontalAlignment','center',...
    'ToolTipString','Switch to Violinplot','value',H.ui.show_violin,...
    'BackgroundColor',[0.8 0.8 0.8],...
    'Interruptible','on','Visible','on','FontSize',H.FS-2);
end

% colormap for samples
if exist('lines')
  cm = lines(n_samples);
else
  cm = jet(n_samples);
end
opt = struct('groupnum',0,'ygrid',0,'violin',2*H.ui.show_violin,'median',2,'groupcolor',cm);
ylim_add = 0.075;

cat_plot_boxplot(data,opt);

hold on
for i=1:n_samples
  ind = find(H.sample(H.ind) == i);

  for j=1:length(ind)
    if H.ui.show_name
      text(xpos{i}(j),data{i}(j),H.filename.m{ind(j)},'FontSize',H.FS-2,'HorizontalAlignment','center')
    else
      plot(xpos{i}(j),data{i}(j),'k.');
    end
  end
end

set(gca,'XTick',[],'XLim',[-.25 n_samples+1.25]);

yamp = max(data_boxp) - min(data_boxp) + 0.001;
ylim_min = min(data_boxp) - ylim_add*yamp;
ylim_max = max(data_boxp) + ylim_add*yamp;
set(gca,'YLim',[ylim_min ylim_max]);

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
    text(xpos, ylim_max,'High rating (good quality) ------>  ','Color','blue','Rotation',...
        90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
  else
      text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
          90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
      text(xpos, ylim_min,'<----- High rating (good quality) ','Color','blue','Rotation',...
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

H.ui.bp = struct('data',data_boxp,'name',name_boxp,'order',quality_order);

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H,'alpha')
  H.ui.alphaval = get(H.ui.alpha,'Value');
else
  H.ui.alphaval = 0.5;
end

% display image with 2nd colorbar (gray)
axes('Position',H.pos.slice);
image(65 + H.img);
if ~H.mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

% prepare alpha overlays for red and green colors
if H.ui.alphaval > 0

  hold on
  alpha_b = cat(3, zeros(size(H.img_alpha)), zeros(size(H.img_alpha)), H.ui.alphaval*ones(size(H.img_alpha)));
  alpha_r = cat(3, H.ui.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
  hg = image(alpha_b); set(hg, 'AlphaData', 0.25*H.img_alpha.*(H.img_alpha>=0),'AlphaDataMapping','none')
  if ~H.mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-0.25*H.img_alpha.*(H.img_alpha<0),'AlphaDataMapping','none')
  if ~H.mesh_detected, axis image; end
  hold off
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H,'mm')
  slice_mm = get(H.ui.mm,'Value');
else
  slice_mm = 0;
end

if H.names_changed
  P = H.files.Vchanged;
else
  P = H.files.V;
end

H.data.vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
H.data.Orig = P(1).mat\[0 0 0 1]';
sl   = round(slice_mm/H.data.vx(3)+H.data.Orig(3));

% if slice is outside of image use middle slice
if (sl>P(1).dat.dim(3)) || (sl<1)
  sl = round(P(1).dat.dim(3)/2);
end

M  = spm_matrix([0 0 sl]);
H.data.zscore = zeros([P(1).dat.dim(1:2) length(H.files.V)]);
Ymean = reshape(H.data.Ymean,P(1).dat.dim(1:3));
Ystd  = reshape(H.data.Ystd,P(1).dat.dim(1:3));
Ymean = Ymean(:,:,sl);
Ystd  = Ystd(:,:,sl);

for i = 1:length(H.files.V)
  img(:,:) = single(P(i).dat(:,:,sl));
  img(isnan(img)) = 0;
  
  % rescue unscaled data
  H.data.zscore(:,:,i) = img;

  % scale image according to mean
  H.data.vol(:,:,i) = img/mean(img(img ~= 0));
end

% calculate individual Z-score map
for i=1:size(H.data.zscore,3)
  img = H.data.zscore(:,:,i);
  ind = Ystd > 0 & (Ymean > H.data.global | img > H.data.global);
  img(ind) = (img(ind) - Ymean(ind))./Ystd(ind);
  img(~ind) = 0;
  H.data.zscore(:,:,i) = img;
end

% enhance contrast and scale image to 0..64
mn = min(H.data.vol(:));
mx = max(H.data.vol(:));
H.data.vol = 64*((H.data.vol - mn)/(mx-mn));

if isfield(H,'mouse') && isfield(H.mouse,'x')
  if H.ui.sorted
    x = H.ind_sorted(H.mouse.x);
  else
    x = H.mouse.x;
  end
  
  % check whether mouse position is defined
  H.img       = H.data.vol(:,:,x)';
  H.img_alpha = H.data.zscore(:,:,x)';
  
  % correct orientation
  H.img       = rot90(H.img,2);
  H.img_alpha = rot90(H.img_alpha,2);
  
  % use gray scale colormap for values > 64
  axes('Position',H.pos.slice);
  image(65 + H.img);
  axis image
  set(gca,'XTickLabel','','YTickLabel','');
  title('Z-score')
  
  % prepare alpha overlays for red and green colors
  if H.ui.alphaval > 0

    hold on
    alpha_b = cat(3, zeros(size(H.img_alpha)), zeros(size(H.img_alpha)), H.ui.alphaval*ones(size(H.img_alpha)));
    alpha_r = cat(3, H.ui.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
    hg = image(alpha_b); set(hg, 'AlphaData', 0.25*H.img_alpha.*(H.img_alpha>=0),'AlphaDataMapping','none')
    axis image
    hr = image(alpha_r); set(hr, 'AlphaData',-0.25*H.img_alpha.*(H.img_alpha<0),'AlphaDataMapping','none')
    axis image
    hold off
  end
  
  txt = {sprintf('%s',spm_file(H.filename.m{x},'short25')),[],['Displayed slice: ',num2str(round(get(H.ui.mm,'Value'))),' mm']};

  set(H.text,'String',txt,'FontSize',H.FS-2);
  set(H.ui.mm_txt,'String',[num2str(round(get(H.ui.mm,'Value'))),' mm'],...
      'FontSize',H.FS-2);
end

return

%-----------------------------------------------------------------------
function get_new_list(obj,event_obj, option)
%-----------------------------------------------------------------------
global H

if ~nargin
  option = false;
end

if isfield(H,'del')
  if isempty(H.del)
    fprintf('No data removed.\n');
    return
  end
else
  fprintf('No data removed.\n');
  return
end

Hdel = H.del;
job = H.job;

% create new list of xmlf-files if necessary
if ~isempty(job.data_xml{1})
  n = 0;
  for i=1:numel(job.data_xml)
    if any(Hdel == i), continue; end
    n = n + 1;
    data_xml{n,1} = job.data_xml{i};
  end
  job.data_xml = data_xml;  
end

% create new list without deleted data in each sample
data = job.data;
for i=1:numel(job.data)
  data_sel = data{i};
  if ~iscell(data_sel), data_sel = cellstr(data_sel); end
  ind = find(H.sample==i);
  n_subjects = numel(ind);
  del_list = ismember(ind,Hdel);
  Hdel(Hdel < n_subjects) = [];
  
  if exist('data_del','var')
    data_del = [data_del;data_sel( del_list)];
    data_new = [data_new;data_sel(~del_list)];
  else
    data_del = data_sel( del_list);
    data_new = data_sel(~del_list);
  end
  data{i}(del_list) = [];
  job.data{i} = data{i};
end

if isfield(H,'SPM')
  H.SPM.nscan = sum(H.ind);
  H.SPM.xY.P = H.SPM.xY.P(H.ind);
  H.SPM.xY.VY = H.SPM.xY.VY(H.ind);
  
end

switch option
  case -1,
    fprintf('Data that are removed from list:\n');
    for i=1:numel(data_del)
      fprintf('%s\n',data_del{i});
    end
  case 0,
    fprintf('Data that remain in list:\n');
    for i=1:numel(data_new)
      fprintf('%s\n',data_new{i});
    end
  case 1,
    do_rerun(obj,event_obj,false);
    set(H.delui.remove, 'BackGroundColor',[0.94 0.94 0.94]);
    set(H.naviui.select,'BackGroundColor',[0.95 0.95 0.95]);
    datacursormode('on');
    H.isdel = false;
end

return

%-----------------------------------------------------------------------
function do_rerun(obj, event_obj, undo)
%-----------------------------------------------------------------------
global H

set(H.dpui.rawp0,   'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.raw,     'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.log,     'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.report,  'BackGroundColor',[0.94 0.94 0.94]);
set(H.delui.remove, 'BackGroundColor',[0.94 0.94 0.94]);
set(H.naviui.select,'BackGroundColor',[0.95 0.95 0.95]);

  if undo
  H.ind = true(size(H.sample));
  set(H.delui.undo,    'enable','off');
  set(H.delui.new,     'enable','off');
  set(H.delui.list_del,'enable','off');
  set(H.delui.list_new,'enable','off');
  H.del = [];
  H.isdel = false;
  datacursormode('on');
else
  H.ind = ~ismember((1:numel(H.sample)),H.del);
end

show_boxplot(H.data.avg_abs_zscore(H.ind),'Mean absolute Z-score  ',-1);  
if H.isxml
  show_QMzscore(H.X,4);
else
  show_QMzscore(H.X,0);
end

return

%-----------------------------------------------------------------------
function show_glassbrain
%-----------------------------------------------------------------------
global H

vol = H.files.V(H.mouse.x).dat(:,:,:);

% get Z-score
Ymean  = reshape(H.data.Ymean,H.files.V(1).dat.dim);
Ystd   = reshape(H.data.Ystd,H.files.V(1).dat.dim);
zscore = (vol - Ymean)./Ystd;
ind = Ystd > 0 & (Ymean > H.data.global | vol > H.data.global);
zscore(~ind) = 0;

% glassbrain
d1 = squeeze(sum(zscore,1));
d2 = squeeze(sum(zscore,2));
d3 = squeeze(sum(zscore,3));

if ~isfield(H,'mipfig')
  H.mipfig = figure;
end

figure(H.mipfig);

cm = hot(64);
pos  = get(H.mipfig,'Position');
set(H.mipfig,'Menubar','none','NumberTitle','off','Name',sprintf('Sample %d: Z-score %s',...
    H.sample(H.mouse.x),H.filename.m{H.mouse.x}),'Position',[10 H.pos.fig(4)+pos(4) pos(3:4)]);
colormap([1-(cm); cm])

mx2 = 2*max(H.files.V(H.mouse.x).dat.dim);

subplot(2,2,1)
imagesc(rot90(-d1),[-mx2 mx2])
axis off image

subplot(2,2,2)
imagesc(rot90(-d2),[-mx2 mx2])
axis off image

subplot(2,2,3)
imagesc(-d3,[-mx2 mx2])
axis off image

subplot(2,2,4)
colorbar
set(gca,'CLim',[-3 3]);
axis off image
  
return

%-----------------------------------------------------------------------
function show_render_views
%-----------------------------------------------------------------------
global H

if isfield(H,'hx') && isgraphics(H.hx.figure)
  H.hx = cat_surf_render2('Overlay',H.hx,H.texture(:,H.mouse.x));
else
  H.hx = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.texture(:,H.mouse.x)));
  H.hx = cat_surf_render2('Colourbar',H.hx);
end
H.hx = cat_surf_render2('clim',H.hx,H.data.range98);
set(H.hx.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Name',...
  sprintf('Sample %d: %s %s',H.sample(H.mouse.x),H.info.texture,H.filename.m{H.mouse.x}))
figure(H.hx.figure)

% get Z-score
zscore = (H.texture(:,H.mouse.x) - H.data.Ymean)./H.data.Ystd;
zscore(H.data.Ystd == 0) = 0;

if isfield(H,'hy') && isgraphics(H.hy.figure)
  H.hy = cat_surf_render2('Overlay',H.hy,zscore);
else
  pos_Hy = get(H.hx.figure,'Position');
  pos_Hy = pos_Hy + [pos_Hy(3)+30 0 0 0];
  H.hy = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',zscore));
  H.hy = cat_surf_render2('Colourbar',H.hy);
  H.hy = cat_surf_render2('ColourMap',H.hy,cat_io_colormaps('BWR',64));
  set(H.hy.figure,'Position',pos_Hy);  
end
H.hy = cat_surf_render2('clim',H.hy,[-3 3]);
set(H.hy.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Name',sprintf('Sample %d: Z-score %s',H.sample(H.mouse.x),H.filename.m{H.mouse.x}));  
figure(H.hy.figure)
  
return

%-----------------------------------------------------------------------
function show_image_slice
%-----------------------------------------------------------------------
global H

% add sliders for volume data
set(H.ui.mm,'Visible','on');
set(H.ui.mm_txt,'Visible','on');
H.img = H.data.vol(:,:,H.mouse.x)';

% alpha overlay
H.img_alpha = H.data.zscore(:,:,H.mouse.x)';

% correct orientation
H.img = rot90(H.img,2);
H.img_alpha = rot90(H.img_alpha,2);

% display image with 2nd colorbar (gray)
image(65 + H.img);
if ~H.mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);
title('Z-score')

% prepare alpha overlays for red and green colors
if H.ui.alphaval > 0

  hold on
  alpha_b = cat(3, zeros(size(H.img_alpha)), zeros(size(H.img_alpha)), H.ui.alphaval*ones(size(H.img_alpha)));
  alpha_r = cat(3, H.ui.alphaval*ones(size(H.img_alpha)), zeros(size(H.img_alpha)), zeros(size(H.img_alpha)));
  hg = image(alpha_b); set(hg, 'AlphaData', 0.25*H.img_alpha.*(H.img_alpha>=0),'AlphaDataMapping','none')
  if ~H.mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-0.25*H.img_alpha.*(H.img_alpha<0),'AlphaDataMapping','none')
  if ~H.mesh_detected, axis image; end
  hold off
end
  
return

%-----------------------------------------------------------------------
function show_report(obj, event_obj)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  H.show_sel = 2;
  set(H.ui.plotbox, 'Visible', 'off');
  set(H.ui.fnambox, 'Visible', 'off');
  set(H.dpui.report, 'BackGroundColor',[0.95 0.95 0.95]);
  set(H.dpui.log,    'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.raw,    'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.rawp0,  'BackGroundColor',[0.94 0.94 0.94]);
end

jpg_file = H.files.jpg{H.mouse.x};
if ~isempty(jpg_file)

  spm_figure('Clear',H.Fgraph);
  figure(H.Fgraph);

  ppos = [0 0 1 1];
  jpg  = imread(jpg_file); 
  set(gca,'position',ppos(1,:));
  gpos = H.Fgraph.Position; 
  [Xq,Yq] = meshgrid(1:size(jpg,2)/gpos(3)/2:size(jpg,2),...
                     1:size(jpg,1)/gpos(4)/2:size(jpg,1));
  jpgi = zeros([size(Xq,1) size(Xq,2) 3],'uint8');
  for i=1:3, jpgi(:,:,i) = uint8(interp2(single(jpg(:,:,i)),Xq,Yq,'linear')); end
  image(H.Fgraph.CurrentAxes,jpgi);
  set(gca,'Visible','off'); 

end

return

%-----------------------------------------------------------------------
function show_raw(obj, event_obj, overlay)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  set(H.ui.plotbox, 'Visible', 'off');
  set(H.ui.fnambox, 'Visible', 'off');
  set(H.dpui.report,  'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.log,     'BackGroundColor',[0.94 0.94 0.94]);
  if overlay
    set(H.dpui.rawp0, 'BackGroundColor',[0.95 0.95 0.95]);
    set(H.dpui.raw,   'BackGroundColor',[0.94 0.94 0.94]);
    H.show_sel = 4;
  else
    set(H.dpui.rawp0, 'BackGroundColor',[0.94 0.94 0.94]);
    set(H.dpui.raw,   'BackGroundColor',[0.95 0.95 0.95]);
    H.show_sel = 3;
  end
end

raw_file = H.files.raw{H.mouse.x};
p0_file  = H.files.p0{H.mouse.x};
if ~isempty(raw_file)

%  spm_orthviews('Reset')

  job.colormapc = flipud(cat_io_colormaps('BCGWHcheckcov'));
  job.prop  = 0.2; 
 
  spm_check_registration(char(raw_file));

  if overlay && exist(p0_file,'file')
    spm_orthviews('addtruecolourimage',1,p0_file,...
      job.colormapc,job.prop,0,5);
  end

  spm_orthviews('Reposition',[0 0 0]); 
  spm_orthviews('redraw');  
  
end

return

%-----------------------------------------------------------------------
function show_log(obj, event_obj)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  H.show_sel = 5;
  set(H.ui.plotbox, 'Visible', 'off');
  set(H.ui.fnambox, 'Visible', 'off');
  set(H.dpui.report,  'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.log,     'BackGroundColor',[0.95 0.95 0.95]);
  set(H.dpui.raw,     'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.rawp0,   'BackGroundColor',[0.94 0.94 0.94]);
end

log_file = H.files.log{H.mouse.x};
if ~isempty(log_file)

  spm_figure('Clear',H.Fgraph);
  figure(H.Fgraph);
  axis off
  
  textbox = [0 0 1 1];  

  fid = fopen(log_file);
  ph  = uipanel(H.Fgraph,'Units','normalized','position',textbox, ...
    'BorderWidth',0,'title',spm_str_manip(log_file,'k100'),'ForegroundColor',[0 0 0.8]);
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
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global H

if H.sel
  sel = H.sel;
else
  sel = 2; % Z-score by default
end

pos_mouse = get(event_obj, 'Position');

x = find(H.X(:,1) == pos_mouse(2));
if isempty(x) || numel(x) > 1
  x = find(H.X(:,sel) == pos_mouse(1));
end

% text info for data cursor window
txt = {sprintf('%s',H.filename.m{x})};

% prevent that that function is called again if position has not changed
if isfield(H,'mouse') && isfield(H.mouse,'x') && x == H.mouse.x
  return
else
  H.mouse.x = x;
end

% build list of selected data points to remove and return
if isfield(H,'isdel') && H.isdel
  
  axes(H.ax)
  hold on
  % reconsider this data point if already in the list
  if isfield(H,'del') && any(H.del == x)
    plot(H.X(x,sel),H.X(x,1),'wx','MarkerSize',10,'Linewidth',2);
    plot(H.X(x,sel),H.X(x,1),'wo','MarkerSize',10,'Linewidth',2,'MarkerFaceColor','w');
    plot(H.X(x,sel),H.X(x,1),'ko','MarkerSize',5,'Linewidth',2);
    txt = {'Reconsider this data point'};
    H.del(H.del==x) = [];
    return
  else
    txt = {'Remove this data point from list'};
    plot(H.X(x,sel),H.X(x,1),'rx','MarkerSize',10,'Linewidth',2);
    
  end
  hold off
  
  % add point to the list
  if isfield(H,'del')
    H.del = unique([H.del x]);
  else
    H.del = x;
  end
    
  % also update index of considered data and enable icons
  H.ind = ~ismember((1:numel(H.sample)),H.del);
  set(H.delui.undo,'enable','on');
  set(H.delui.new, 'enable','on');
  set(H.delui.list_del,'enable','on');
  set(H.delui.list_new,'enable','on');

  return
end

if ~H.mesh_detected
  set(H.ui.alpha,'Visible','on');
  set(H.ui.alpha_txt,'Visible','on');
end

% text info for textbox
if H.mesh_detected
  txt2 = {[],sprintf('%s',spm_file(H.filename.m{H.mouse.x},'short25'))};
else
  txt2 = {[],sprintf('%s',spm_file(H.filename.m{H.mouse.x},'short25')),[],'Individual Z-score','(red: - blue: +)'};
end

set(H.text,'String',txt2,'FontSize',H.FS-2);

if ~H.mesh_detected
  axes('Position',H.pos.slice);
end

if H.mesh_detected 
  % show two render views for meshes: texture and Z-score
  show_render_views;
else
  % show image slice and glassbrain
  show_image_slice;
  show_glassbrain;
end

switch(H.show_sel)
  case 2, show_report;
  case 3, show_raw(obj, event_obj, false);
  case 4, show_raw(obj, event_obj, true);
  case 5, show_log;
end

return
