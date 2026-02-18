function varargout = cat_stat_homogeneity(job)
%cat_stat_homogeneity. To check Z-score across sample.
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%
% varargout = cat_stat_homogeneity(job)
%  
% job                .. SPM job structure
%  .data             .. volume and surface files
%  .globals          .. optionally correct TIV using global scaling (for VBM only)
%  .gSF              .. global scaling values
%  .c                .. confounds
%  .data_xml         .. optional xml QC data
%  .verb             .. print figures
%  .new_fig          .. use new window instead of SPM Fgraph
%  .xM               .. optional mask information from SPM.xM
%
% varargout          .. output structure 
%  .zscore           .. quartic mean Z-score
%  .table            .. Z-score table
%  .sorttable        .. sorted Z-score table
%  .threshold_zsc    .. lower threshold for Z-score (mean - 4*std)
%
% Example: 
%   cat_stat_homogeneity(struct('data',{{ files }} ,'c',[],'data_xml',{{}}));
%
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

if nargin == 0
  error('No argument given.');
end

H.sample        = [];
H.mouse.x       = 1;
H.mouse.xold    = 1;
H.del           = [];
H.ui.alphaval   = 0.5;
H.names_changed = false;
H.cmap          = [jet(64); gray(64)]; % create two colormaps
G               = [];
n_subjects      = 0;

% use new figure by default
if ~isfield(job,'new_fig')
  job.new_fig = true;
end

if ~isfield(job,'verb')
  job.verb = true;
end

H.job = job;

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

% check for repeated anova design with long. data
if isfield(job,'factorial_design') && isfield(job.factorial_design,'des') && isfield(job.factorial_design.des,'fblock')
  H.repeated_anova = true;
elseif isfield(job,'factorial_design') && isfield(job.factorial_design,'spmmat')
  load(job.factorial_design.spmmat{1})
  H.repeated_anova = ~isempty(SPM.xX.iB);
else
  H.repeated_anova = false;
end

% read filenames for each sample and indicate sample parameter
if H.mesh_detected
  n_samples = numel(job.data);
  sinfo = cat_surf_info(char(job.data{1}(1,:)));
  H.Pmesh = gifti(sinfo.Pmesh);
  for i=1:n_samples
    [pp,ff,ee] = spm_fileparts(char(deblank(job.data{i}(1,:)))); 
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

% use global scaling from design matrix
if isfield(job,'gSF') && numel(job.gSF) == n_subjects
  fprintf('Use global scaling from design matrix (i.e. with TIV).\n');
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

% prepare design matrix for adjusting nuisance parameters
if isfield(job,'c') && ~isempty(job.c) 
  for i=1:numel(job.c)
    G = [G job.c{i}];
  end
  if size(G,1) ~= n_subjects
    G = G';
  end
  % mean correction
  G = G - mean(G);
  iG = pinv(G);
end

if isempty(char(job.data_xml))
  H.isxml        = false;
  xml_defined    = false;
  H.xml.QM_names = '';
  xml_files      = [];
else
  xml_files = char(job.data_xml);
  if size(xml_files,1) ~= n_subjects
    fprintf('Only %d of %d report files were defined. Try to find xml-files for quality measures.\n',size(xml_files,1),n_subjects);
    H.isxml     = false;
    xml_defined = false;
  else
    H.isxml     = true;
    xml_defined = true;
  end
end

% select school marks (range best-low = 0.5 - 10.5) or percentage system (range low-best: 0-100%) 
if isfield(job,'userps'), H.userps = job.userps; else, H.userps = -1; end

H.xml.QM = ones(n_subjects,3);
H.xml.QM_names = char('Noise','Bias','Structural image quality rating (SIQR)');
H.xml.QM_names_multi = char('Noise & Quartic Mean Z-score','Bias & Quartic Mean Z-score','SIQR & Quartic Mean Z-score');
H.xml.QM_order = H.userps .* ones(1,3);

% add some more entries for surfaces
if H.mesh_detected
  H.xml.QM = ones(n_subjects,5);
  H.xml.QM_names = char(H.xml.QM_names,'Euler number','Size of topology defects');
  H.xml.QM_names_multi = char(H.xml.QM_names_multi,'Euler number & Quartic Mean Z-score','Size of topology defects & Quartic Mean Z-score');
  H.xml.QM_order = -ones(1,5);
end


% To find other related files of each subject, we first have to figure out
% what the general prefix of the data itself is (asuming that it is only
% useful to anlyse on data class (defined by the prefix) at once. 
% To find the prefix, we will try to find the XML file of the first subject
% that thould be in the report directory (if exist) or otherwise in the
% same directory in case of BIDS.
% The test for BIDS has to be done later as far BIDS and non BIDS can be
% mixed.

% try to detect report folder and check if this fits for all
pth = spm_fileparts( H.files.fname{1} );
if isfield(job,'sel_xml') && isfield(job.sel_xml,'select_dir')
  report_folder = char(job.sel_xml.select_dir);
else
  % if there is a report directory than use it
  ppth = spm_fileparts(pth); 
  if exist( fullfile( ppth, 'report' ), 'dir' )
    report_folder = fullfile( ppth, 'report' );
  else
    report_folder = pth; 
  end
end

% search xml report files if not defined
prep_str = '';
if ~xml_defined
  fprintf('Search xml-files ');

  % we now try to find all XML files in the report folder 
  xml_files = spm_select('List',report_folder,'^cat_.*\.xml$');
  if ~isempty(xml_files) || size(xml_files,1) == 1

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
if job.verb, cat_progress_bar('Init',n_subjects,'Load xml-files'); end
res_RMS = nan(n_subjects,1); 
for i=1:n_subjects
  % get basename for data files
  [pth, data_name, ee] = fileparts(H.files.fname{i});
  if ~strcmp(ee,'.nii') && ~strcmp(ee,'.gii'), data_name = [data_name ee]; end
  
  % remove ending for rigid or affine transformed files
  data_name = strrep(data_name,'_affine','');
  data_name = strrep(data_name,'_rigid','');
  
  % detect if BIDS file structure is used
  % (hope this is unique enough - would be possible to add CAT12 as subdirectory) 
  BIDSdir   = [filesep 'derivatives' filesep]; 
  isBIDS    = cat_io_contains( pth , BIDSdir );

  if isBIDS 
    % in case of BIDS CAT wrote all files into this directory 
    report_folder = pth; 
    [ppm,pps] = spm_fileparts(report_folder);
    % just in case that there are subdirs
    if cat_io_contains(pps,{'mri','surf','label'})
      pps = 'report';
      report_folder = fullfile(ppm,pps);
    end
  elseif ( isfield(job,'sel_xml') && isfield(job.sel_xml,'select_dir') ) 
    % in case of BIDS CAT wrote all files into this directory 
    report_folder = job.sel_xml.select_dir{1}; 
  else % use relative folder for autom. search
    report_folder = fullfile(pth,'..','report');
    if ~exist(report_folder,'dir'), report_folder = pth; end
  end
  
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
    if i > size(xml_files,1) 
      cat_io_cprintf('warn','\nSkip use of xml-files for quality measures because of not enough XML files were found.\n');
      H.isxml = false;
      break
    elseif isempty(strfind(data_name,subjname))
      cat_io_cprintf('warn','\nSkip use of xml-files for quality measures because of deviating subject names:\n%s vs. %s\n',H.files.fname{i},xml_files(i,:));
      H.isxml = false;
      break
    end
    xml_file = deblank(xml_files(i,:));
  end
  
  %% get mri folder
  [pth, data_name, ee] = fileparts(H.files.fname{i});
  if strcmp(ee,'.gii') && ~isBIDS, mri_folder = fullfile(fileparts(pth),'mri');
  else, mri_folder = pth; end
  
  % find raw/p0 files
  H.files.raw{i} = fullfile(fileparts(pth),[subjname ee]);
  H.files.p0{i}  = fullfile(mri_folder,['p0' subjname '.nii']);
  
  if isBIDS
  % get BIDS raw diretory of this subject by looking for the derivatives
  % directory. The parent path give us the BIDS main directory where the 
  % original files should be located in sub path given behind the derivative 
  % directory 
    BIDSfst         = strfind( pth , BIDSdir ) - 1;
    BIDSlst         = strfind( pth , BIDSdir ) + numel(BIDSdir); 
    BIDSlst         = BIDSlst + find(pth(BIDSlst:end)==filesep,1,'first');
    BIDSrawdir      = pth(1:BIDSfst);
    BIDSsubdirs     = pth(BIDSlst:end); 

    H.files.raw{i}  = fullfile(BIDSrawdir,BIDSsubdirs,[subjname ee]);
  end
  
  H.files.rawgz{i} = [H.files.raw{i} '.gz']; 
  if ~exist(H.files.raw{i},'file'), H.files.raw{i} = ''; end
  if isempty(H.files.raw{i}) && exist(H.files.rawgz{i},'file') && ~exist('foundrawgz','var')
    cat_io_cprintf('warn','Gzipped (original) files are not supported in SPM display yet.\n');
    foundrawgz = true; %#ok<NASGU> 
  end
  if ~exist(H.files.p0{i}, 'file'), H.files.p0{i}  = ''; end

  if exist(xml_file,'file')
    H.job.data_xml{i} = xml_file;
    xml = cat_io_xml(xml_file);
    n_xml_files = n_xml_files + 1;
    H.isxml = true;
    
    % find jpg/pdf/log files
    H.files.jpg{i} = fullfile(report_folder,['catreportj_' subjname '.jpg']);
    if H.repeated_anova
      H.files.jpg_long{i} = fullfile(report_folder,['catlongreportj_' data_name '.jpg']);
      if ~exist(H.files.jpg_long{i},'file'), H.files.jpg_long{i} = ''; end
    end
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
      cat_io_cprintf('warn','\nFile "%s" not found. \nSkip use of xml-files for quality measures and TIV. ',xml_file);
    else
      cat_io_cprintf('warn','\nFile "%s" not found. \nSkip use of xml-files for quality measures. ',xml_file);
    end
    cat_io_cprintf('warn','Please check if only one data type (e.g., mwp1) is used or select xml-files manually.\n\n');
    H.isxml = false;
    is_gSF  = false;
    break
  end

  if i > n_xml_files 
    cat_io_cprintf('warn','\nSkip use of xml-files for quality measures because of not enough XML files were found.\n');
    H.isxml = false;
    break
  elseif ~isfield(xml,'qualityratings') && ~isfield(xml,'QAM')
    cat_io_cprintf('warn',['\nQuality rating is not saved for %s. Report file %s is incomplete. ' ...
      '\nPlease repeat preprocessing and check for potential errors in the "err" folder.\n'],H.files.fname{i},xml_files(i,:));  
    H.isxml = false;
    break
  end

  if H.mesh_detected
    if isfield(xml.qualityratings,'NCR')
    % check for newer available surface measures
      if isfield(xml.subjectmeasures,'EC_abs') && isfinite(xml.subjectmeasures.EC_abs) && isfinite(xml.subjectmeasures.defect_size)
        H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.SIQR xml.subjectmeasures.EC_abs xml.subjectmeasures.defect_size];
      else
        H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.SIQR NaN NaN];
      end
    else % also try to use old version
      H.xml.QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
    end
  else
    if isfield(xml.qualityratings,'NCR')
      H.xml.QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.SIQR];
    else % also try to use old version
      H.xml.QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
    end
  end

  res_RMS(i,:) = xml.qualityratings.res_RMS; 
  if job.verb, cat_progress_bar('Set',i); end
end
if H.userps > 0
  mark2rps = @(mark) min(100,max(0,105 - mark*10));
  H.xml.QM = mark2rps( H.xml.QM ); 
end
if job.verb, cat_progress_bar('Clear'); end

if H.isxml
  if n_xml_files ~= n_subjects
    cat_io_cprintf('warn','Only %d of %d report files found. Skip use of xml-files for quality measures.\n',n_xml_files,n_subjects);
    H.isxml = false;
  else
    fprintf('%d report files with quality measures were found.\n',n_xml_files);
  end
end

% delete QM entries
if ~H.isxml
  H.xml.QM = [];
  H.xml.QM_order = [];
  H.xml.QM_names = '';
  H.xml.QM_names_multi = '';
end

% remove last two columns if EC_abs and defect_size are not defined
if H.isxml && H.mesh_detected && all(~isfinite(H.xml.QM(:,4))) && all(~isfinite(H.xml.QM(:,5)))
  H.xml.QM = H.xml.QM(:,1:3);
  H.xml.QM_order = H.xml.QM_order(1:3);
  H.xml.QM_names = H.xml.QM_names(1:3,:);
  H.xml.QM_names_multi = H.xml.QM_names_multi(1:3,:);
end

% add covariates to list if it's not from a statistical analysis where the
% covariates include all columns of the design matrix
if isfield(job,'c') && ~isempty(job.c) && ~isfield(job,'factorial_design')
  for i=1:numel(job.c)
    if ~isempty(H.xml.QM)
      H.xml.QM = [H.xml.QM job.c{i}];
      H.xml.QM_order = [H.xml.QM_order 0];
      H.xml.QM_names = char(H.xml.QM_names, sprintf('Covariate %d',i));
      H.xml.QM_names_multi = char(H.xml.QM_names_multi, sprintf('Covariate %d & Quartic Mean Z-score',i));
    else
      H.xml.QM = job.c{i};
      H.xml.QM_order = 0;
      H.xml.QM_names = sprintf('Covariate %d',i);
      H.xml.QM_names_multi = sprintf('Covariate %d & Quartic Mean Z-score',i);
    end
  end
end

H.data.Ymean = 0.0;
Yss = 0.0; % sum of squares

% preload surface data for later render view
if H.mesh_detected
  % load surface texture data
  H.texture = single(spm_data_read(H.files.V));
end

if job.verb, cat_progress_bar('Init',n_subjects,'Load data'); end

% prepare Beta
if ~isempty(G) 
  if H.mesh_detected
    dim = (H.files.V(1).dim);
  else
    dim = H.files.V(1).dat.dim;
  end
  Beta = zeros(prod(dim),size(G,2),'single');
end

% make initial mask
if H.mesh_detected
  dim = (H.files.V(1).dim);
else
  dim = H.files.V(1).dat.dim;
end
mask = true(dim);
M = H.files.V(1).mat;

%-Get explicit mask(s)
%==========================================================================
if isfield(job,'xM')
  for i = 1:numel(job.xM.VM)
    if ~H.mesh_detected
      C = spm_bsplinc(job.xM.VM(i), [0 0 0 0 0 0]');
      v = true(dim);
      [x1,x2] = ndgrid(1:dim(1),1:dim(2));
      for x3 = 1:dim(3)
        M2  = inv(M\job.xM.VM(i).mat);
        y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
        y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
        y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
        v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
      end
      mask = mask & v;
      clear C v x1 x2 x3 M2 y1 y2 y3
    else
      v = full(job.xM.VM(i).private.cdata) > 0;
      mask = mask & v(:);
      clear v
    end
  end
end

for i = 1:n_subjects
  if H.mesh_detected
    Ytmp = spm_data_read(H.files.V(i));
  else
    Ytmp(:,:,:) = H.files.V(i).dat(:,:,:);
  end
  
  % get mask
  if isfield(job,'xM')
    mask(mask) = Ytmp(mask) > job.xM.TH(i);      %-Threshold (& NaN) mask
  end
  
  Ytmp(~isfinite(Ytmp)) = 0;
  
  % either global scaling was externally defined using job or values were
  % used from xml-file
  if is_gSF || isfield(job,'gSF')
    Ytmp = Ytmp*gSF(i)/mean(gSF);
  end
  
  if i>1 && numel(H.data.Ymean) ~= numel(Ytmp)
    cat_io_cprintf('err','\n\nERROR: File %s has different data size: %d vs. %d\n\n',job.data{i},numel(H.data.Ymean),numel(Ytmp));
    return
  end

  % estimate Beta
  if ~isempty(G) 
    for j = 1:size(Beta,2)
      Beta(:,j) = Beta(:,j) + single(iG(j,i)*Ytmp(:));
    end
  end
  
  H.data.Ymean = H.data.Ymean + Ytmp(:);
  Yss   = Yss + Ytmp(:).^2;
  if job.verb, cat_progress_bar('Set',i); end
end
if job.verb, cat_progress_bar('Clear'); end

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
ind = H.data.Ystd ~= 0 & H.data.Ymean > H.data.global & mask(:);

if job.verb, cat_progress_bar('Init',n_subjects,'Calculate Z-score'); end

H.data.avg_abs_zscore = zeros(n_subjects,1);
for i = 1:n_subjects
  
  if H.mesh_detected
    Ytmp = spm_data_read(H.files.V(i));
  else
    Ytmp(:,:,:) = H.files.V(i).dat(:,:,:);
  end
  Ytmp(~isfinite(Ytmp)) = 0;
    
  if is_gSF
    Ytmp = Ytmp*gSF(i)/mean(gSF);
  end
  
  % correct for nuisance
  if ~isempty(G) 
    for j = 1:size(Beta,2)
      Ytmp(:) = Ytmp(:) - double(G(i,j)*Beta(:,j));
    end
  end
  
  % calculate Z-score
  zscore = (Ytmp(ind) - H.data.Ymean(ind))./H.data.Ystd(ind);
  
  % use mean of Z-score as overall measure, but emphasize outliers by
  % using power operation
  power_scale = 4;
  H.data.avg_abs_zscore(i) = mean((abs(zscore).^power_scale))^(1/power_scale);
  if job.verb, cat_progress_bar('Set',i); end
end

if job.verb, cat_progress_bar('Clear'); end

if isfield(job,'sites') && ~isempty(job.sites)
  if numel(job.sites{1}) == numel(H.files.fname)
    sites = job.sites{1}';
  else
    cat_io_cprintf('err', ...
      sprintf('  Warning number of site enties does not fit to the number of scans (%d/=%d).\n', ...
        numel(job.sites{1}), numel(H.files.fname) ) ); 
    sites = [];
  end
else
  sites = []; 
end

if H.isxml
% RD20250913: added normalized SIQR measure
% -------------------------------------------------------------------------
  if isfield(xml.qualityratings,'NCR'), SIQR = H.xml.QM(:,3); else, SIQR = H.xml.QM(:,1); end  % use SIQR or NCR
  if isempty(sites)
    % use site or resolution measures as approximation
    [nsites,~,sites] = unique(round(res_RMS,3)); 
    % in case of too many sites focus on one
    if numel(sites) / numel(nsites) < 5, sites = ones(size(res_RMS)); end
  end

  % estimate rating (use the "better" quantil for linear normalization for each site)
  NSIQR = cat_tst_qa_normer( SIQR , struct( 'sites' , sites , 'figure' , 0, 'model' , 0, 'cmodel' , 1));   

  % extend 
  H.xml.QM(:,end+1)     = NSIQR; 
  H.xml.QM_names        = char( [ cellstr(H.xml.QM_names); {'Normalized SIQR (nSIQR)'} ] );
  H.xml.QM_names_multi  = char( [ cellstr(H.xml.QM_names_multi); {'nSIQR & Quartic Mean Z-score'} ]);
  H.xml.QM_order(end+1) = H.userps; 
else
  if isempty(sites)
    sites = ones(size(H.files.fname));
  end
end
H.sites = sites; 


% save the rating
Pcsv = fullfile(pwd, sprintf('CheckSampleHomogeneity_%0.0fsubjects_%0.0fsites_%s.csv', ...
  numel(H.files.fname), numel(unique(sites)), char(datetime('now','format','yyyyMMdd-HHmm')) )); 
if H.isxml
  tab = [ 
    {'fname', 'site', 'res_RMS', 'SIQR','nSIQR','zscore'};
    H.files.fname, num2cell(sites), num2cell(res_RMS), num2cell(SIQR), num2cell(NSIQR), num2cell(H.data.avg_abs_zscore);
    ];
else
  tab = [ 
    {'fname', 'site', 'zscore'};
    H.files.fname, num2cell(sites), num2cell(H.data.avg_abs_zscore);
  ];
end
cat_io_csv(Pcsv,tab);
fprintf('Write csv-table:\n'); 
cat_io_cprintf('blue',sprintf('  %s\n',Pcsv)); 




% get data for each subject in longitudinal designs
if H.repeated_anova
  if isfield(job.factorial_design,'des')
    fsubject = job.factorial_design.des.fblock.fsuball.fsubject;
    n_subjects_long = numel(fsubject);
    H.ind_subjects_long = cell(numel(fsubject),1);
    n = 0;
    for i = 1:n_subjects_long
      n_scans = numel(fsubject(i).scans);
      % find time points in all subjects
      H.ind_subjects_long{i} = ismember(1:n_subjects,n + (1:n_scans));
      n = n + n_scans;
    end
  else
    [rw,cl] = find(SPM.xX.I == length(SPM.xX.iB)); % find column which codes subject factor (length(SPM.xX.iB) -> n_subj)
    subj_col = cl(1);
    n_subjects_long = max(SPM.xX.I(:,subj_col));

    H.ind_subjects_long = cell(n_subjects_long,1);
    n = 0;
    for i = 1:n_subjects_long
      ind_subj = find(SPM.xX.I(:,subj_col)==i);
      n_scans = numel(ind_subj);
      % find time points in all subjects
      H.ind_subjects_long{i} = ismember(1:n_subjects,n + (1:n_scans));
      n = n + n_scans;
    end
    
  end
end

% voxelsize and origin of volume data
if ~H.mesh_detected
  H.data.vx =  sqrt(sum(H.files.V(1).mat(1:3,1:3).^2));
  H.data.Orig = H.files.V(1).mat\[0 0 0 1]';
end

% positions & font size
ws = spm('Winsize','Graphics');
H.FS = cat_get_defaults('extopts.fontsize');

popb = [0.038 0.035];  % size of the small buttons
popm = 0.780;          % x-position of the control elements

H.pos = struct(...
    'fig',    [10 10 1.3*ws(3) 1.1*ws(3)],...% figure
    'cbar',   [0.045 0.035 0.700 0.020],...  % colorbar for figure
    'plot',   [0.050 0.050 0.700 0.825],...  % scatter plot
    ...
    'close',  [0.775 0.935 0.100 0.040],...  % close button
    'show',   [0.875 0.935 0.100 0.040],...  % button to show worst cases
    'scat',   [0.775 0.880 0.100 0.050],...  % button to enable ordered matrix
    'boxp',   [0.875 0.880 0.100 0.050],...  % button to display boxplot
    ...
    ... == navigation unit ==
    'scSelect',    [popm+popb(1)*0 0.820 popb],... % select (default) 
    'scZoomReset', [popm+popb(1)*1 0.820 popb],... % standard zoom
    'scZoomIn',    [popm+popb(1)*2 0.820 popb],... % zoom in 
    'scZoomOut',   [popm+popb(1)*3 0.820 popb],... % zoom out
    'scPan',       [popm+popb(1)*4 0.820 popb],... % pan (moving hand)
    ...
    ... == remove unit ==
    'rmDel',       [popm+popb(1)*0 0.750 popb],... % delete 
    'rmUndo',      [popm+popb(1)*1 0.750 popb],... % undo deletion
    'rmNew',       [popm+popb(1)*2 0.750 popb],... % calculate new
    'rmListNew',   [popm+popb(1)*3 0.750 popb],... % list remaining data
    'rmListDel',   [popm+popb(1)*4 0.750 popb],... % list removed data
    ...
    ... == display unit ==
    'dpReport',    [popm+popb(1)*0 0.680 popb],... % report 
    'dpReportLong',[popm+popb(1)*1 0.680 popb],... % report 
    'dpRaw',       [popm+popb(1)*2 0.680 popb],... % raw data
    'dpRawP0',     [popm+popb(1)*3 0.680 popb],... % raw data + p0
    'dpLog',       [popm+popb(1)*4 0.680 popb],... % log
    ...
    ... == check boxes ==
    'fnambox',[0.775 0.600 0.200 0.050],... % show filenames?
    'plotbox',[0.875 0.600 0.200 0.050],... % switch between boxplot and violin plot 
    ...
    ... == slice display ==
    'text',   [0.775 0.450 0.200 0.150],... % textbox with info
    'aslider',[0.775 0.405 0.200 0.040],... % slider for alpha overlay
    'slice',  [0.775 0.030 0.200 0.400],... % slice images according to position of mouse pointer
    'zslider',[0.775 0.020 0.200 0.040]);   % slider for z-slice   

% use this window for Fgraph
if isfield(job,'new_fig') && job.new_fig
  H.Fgraph = spm_figure('GetWin','Homogeneity');
else
  H.Fgraph = spm_figure('GetWin','Graphics');
end

% correct position so to prevent overlapping windows
pos =   get(H.Fgraph,'Position');
set(H.Fgraph,'Position',[H.pos.fig(1)+H.pos.fig(3)+10 pos(2:4)]);

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
  
  % get common filename (for repeated Anova use all data otherwise data for
  % that sample)
  if H.repeated_anova
    [tmp, fname_tmp{i}] = spm_str_manip(char(H.files.fname),'C');
  else
    [tmp, fname_tmp{i}] = spm_str_manip(char(H.files.fname{H.sample == i}),'C');
  end
  if ~isempty(fname_tmp{i})
    fname_m    = [fname_m; fname_tmp{i}.m]; 
    fname_s{i} = fname_tmp{i}.s;
    fname_e{i} = fname_tmp{i}.e;
  else
    fname_s{i} = '';
    fname_e{i} = '';
  end
  if job.verb
    try
      %% try some colorful output to make it easier to read
      % suppress too long outputs
      breaks(1) = find(tmp=='{',1,'first');
      breaks(2) = find(tmp=='}',1,'last');

      fprintf('Compressed filenames sample %d: ',i); 

      cat_io_cprintf([0.0 0.2 .8], '%s', tmp(1:breaks(1)-1));
      % to long cmd line output can cause java errors 
      if numel(tmp(breaks(1):breaks(2))) < 1000 
        cmdlinelim = 120; % 1.5 times as usual ?
        tmptmp = [' ...\n  ' tmp(breaks(1):breaks(2))];
        for tmpi = flip(cmdlinelim+8:cmdlinelim:numel(tmptmp))
          closekomma = find(tmptmp(1:tmpi)==',',1,'last');
          tmptmp = [tmptmp(1:closekomma) ' ...\n  ' tmptmp(closekomma+1:end)]; 
        end
        cat_io_cprintf([0.5 0.0 .5], sprintf('%s', tmptmp));
      else
        cat_io_cprintf([0.5 0.0 0], '%s', '{...TOO_LONG_SUPPRESSED...}');
      end
      cat_io_cprintf([0.0 0.2 .8], '%s', tmp(breaks(2)+1:end));
      fprintf('\n');
    catch
      fprintf('Compressed filenames sample %d: %s  \n',i,tmp);
    end
  end
end

H.filename = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});

% sort data
[H.data.avg_abs_zscore_sorted, H.ind_sorted] = sort(H.data.avg_abs_zscore,'ascend');

threshold_zsc = mean(H.data.avg_abs_zscore) + 2*std(H.data.avg_abs_zscore);
n_thresholded = find(H.data.avg_abs_zscore_sorted > threshold_zsc, 1 );

if ~isempty(n_thresholded) && job.verb
  fprintf('\nThese data have a quartic mean Z-score above 2 standard deviations.\n');
  fprintf('This does not necessarily mean that you have to exclude these data. However, these data have to be carefully checked:\n');
  
  fprintf('SIQR / nSIQR / Quartic mean Z-score / filename\n');
  for i = n_thresholded:n_subjects % just switch this improve readability in case of different fname length
    if isfield(H.xml,'QM') && ~isempty(H.xml.QM)
      cat_io_cprintf([0.5 0 0.5],'  %3.3f:', H.xml.QM(i,end-1:end) );
    end
    cat_io_cprintf([0.5 0 0.5],'  %3.3f: %s \n', H.data.avg_abs_zscore_sorted(i) ,H.files.fname{H.ind_sorted(i)});
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
  
  create_menu;
  show_boxplot(H.data.avg_abs_zscore,'Quartic Mean Z-score  ',-1);

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
  
  % we have to update slice array first if not defined
  if ~isfield(H.data,'vol') && ~H.mesh_detected
    preload_slice_data;
  end
end

%-----------------------------------------------------------------------
function create_menu
%-----------------------------------------------------------------------
global H

% create figure
H.mainfig = figure(22);
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

%colormap(H.cmap)
colormap(jet)

% add button for closing all windows
H.ui.close = uicontrol(H.mainfig,...
        'String','Close','Units','normalized',...
        'Position',H.pos.close,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'Callback','for i=2:26, try close(i); end; end;',...
        'ToolTipString','Close windows',...d 
        'Interruptible','on','Enable','on');

% check button
H.ui.show = uicontrol(H.mainfig,...
        'String','Check worst','Units','normalized',...
        'Position',H.pos.show,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'Callback',@check_worst_data,...
        'ToolTipString','Display most deviating files',...
        'Interruptible','on','Enable','on');

%% create popoup menu for boxplot

% check whether we have to add entries from quality measures or covariates
if isempty(H.xml.QM)
  str  = { 'Boxplot','Quartic Mean Z-score'};
  % average quartic Z-score vs. file order
  H.X = [H.data.avg_abs_zscore (1:numel(H.data.avg_abs_zscore))'];
  show_QMzscore(H.X,0, H.userps); % show file order on x-axis
else
  % average quartic Z-score vs. QM measures
  H.X = [H.data.avg_abs_zscore H.xml.QM];
  
  str  = { 'Boxplot','Quartic Mean Z-score'};
  for i = 1:size(H.xml.QM,2)
    str{i+2} = deblank(H.xml.QM_names(i,:));
  end
    
  if H.isxml
    % estimate product between structural quality rating (SIQR) and quartic mean Z-score 
    H.xml.QMzscore = H.X(:,1).*H.X(:,2);
    str{i+3} = 'SIQR (grad) x Quartic Mean Z-score';
    show_QMzscore(H.X,4, H.userps); % show SIQR on x-axis
  else
    show_QMzscore(H.X,0, H.userps); % show file order on x-axis
  end
end

tmp  = { {@show_boxplot, H.data.avg_abs_zscore, 'Quartic Mean Z-score', -1}};
for i = 1:size(H.xml.QM,2)
  tmp{i+1} = {@show_boxplot, H.xml.QM(:,i), deblank(H.xml.QM_names(i,:)), H.xml.QM_order(i)};
end

if H.isxml
  tmp{i+2} = {@show_boxplot, H.xml.QMzscore, 'SIQR x Quartic Mean Z-score  ', -1};
end

H.ui.boxp = uicontrol(H.mainfig,...
        'String',str,'Units','normalized',...
        'Position',H.pos.boxp,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'Callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Display boxplot',...
        'Interruptible','on','Visible','on');

%% create popoup menu for scatterplot

% if QM values are available allow SIQR and surface parameters, but skip
% noise and bias as first 2 entries
str  = { 'Scatterplot','Quartic Mean Z-score'};
if isempty(H.xml.QM)
  tmp  = {{@show_QMzscore, H.X, 0, H.userps}}; % just quartic mean Z-score with file order
else
  tmp  = {{@show_QMzscore, H.X, 0, H.userps}}; % file order
  if H.isxml
    for i = 1:size(H.xml.QM,2)-2
      str{i+2} = deblank(H.xml.QM_names_multi(i+2,:));
      tmp{i+1} = {@show_QMzscore, H.X, i+3, H.xml.QM_order(i+2)};
    end  
  else
    for i = 1:size(H.xml.QM,2)
      str{i+2} = deblank(H.xml.QM_names_multi(i,:));
      tmp{i+1} = {@show_QMzscore, H.X, i+1, H.xml.QM_order(i)};
    end    
  end
end

H.ui.scat = uicontrol(H.mainfig,...
        'String',str,'Units','normalized',...
        'Position',H.pos.scat,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'Callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Sort matrix',...
        'Interruptible','on','Visible','on');

H.ui.text = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.text,...
        'String',{'','Click in image to display slices'},...
        'Style','text','HorizontalAlignment','center',...
        'ToolTipString','Select slice for display',...
        'FontSize',H.FS-2,'Visible','off','BackgroundColor',[0.8 0.8 0.8]);

%% == zoom unit ==
H.naviui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.scSelect(1) H.pos.scSelect(2)+0.042 0.2 0.02],...
  'String','Zoom options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

H.naviui.select = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scSelect,'callback',...
  ['datacursormode(''on''); global H;' ...
  'set(H.naviui.select,''BackGroundColor'',[0.95 0.95 0.95]);'],...
  'Style','Pushbutton','enable','on','ToolTipString','Data selection','CData',load_icon('tool_pointer.png'),'BackGroundColor',[0.95 0.95 0.95]);

H.naviui.zoomReset = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomReset,'callback',...
  ['zoom out; zoom reset; datacursormode(''on''); global H;' ...
  'set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'],...
  'Style','Pushbutton','enable','on','ToolTipString','Reset view','CData',load_icon('tool_fit.png')); 

H.naviui.zoomIn = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomIn,'callback',...
  ['global H; ' ...
   'hz = zoom(H.ax); set(hz,''enable'',''on'',''direction'',''in'');' ...
   'set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'], ... 
  'Style','Pushbutton','enable','on','ToolTipString','Zoom in','CData',load_icon('tool_zoom_in.png'));

H.naviui.zoomOut = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scZoomOut,'callback',...
  ['global H; ' ...
   'hz = zoom(H.ax); set(hz,''enable'',''on'',''direction'',''out''); ' ...
   'set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);'], ...
  'Style','Pushbutton','enable','on','ToolTipString','Zoom out','CData',load_icon('tool_zoom_out.png'));

H.naviui.pan = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.scPan,'Enable','off','callback',...
  'pan on; set(H.naviui.select,''BackGroundColor'',[0.94 0.94 0.94]);',...
  'Style','Pushbutton','enable','on','ToolTipString','Hand','CData',load_icon('tool_hand.png'));

%% == remove unit ==
H.delui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.rmDel(1) H.pos.rmDel(2)+0.042 0.2 0.02],...
  'String','Data remove options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

H.delui.remove = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmDel,'callback',{@remove_point},...
  'Style','Pushbutton','enable','off','ToolTipString','Remove this data point','CData',load_icon('file_delete.png'));

H.delui.undo = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmUndo,'callback',{@do_rerun,1},...
  'Style','Pushbutton','enable','off','ToolTipString','Undo all deletions','CData',load_icon('file_delete_restore.png')); 

H.delui.new = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmNew,'callback',{@get_new_list,1},...
  'Style','Pushbutton','enable','off','ToolTipString','Refresh without removed data','CData',load_icon('refresh.png')); 

H.delui.list_del = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmListDel,'callback',{@get_new_list,-1},...
  'Style','Pushbutton','enable','off','ToolTipString','List removed data','CData',load_icon('list_del.png')); 

if isfield(H.job,'factorial_design')
  icon = load_icon('greenarrowicon.png');
  str  = 'Create new analysis without removed data';
else
  icon = load_icon('list_new.png');
  str  = ' List remaining data'; 
end

H.delui.analysis_new = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.rmListNew,'callback',{@get_new_list,0},...
  'Style','Pushbutton','enable','off','ToolTipString',str,'CData',icon); 

%% == display unit ==
H.dpui.text = uicontrol(H.mainfig,...
  'Units','normalized','Style','text',...
  'Position',[H.pos.dpReport(1) H.pos.dpReport(2)+0.042 0.2 0.02],...
  'String','Data display options','FontSize',H.FS,'BackgroundColor',[0.8 0.8 0.8]);

% enable some buttons only if respective files are available
if H.isxml, H.status.xml = true;
else H.status.xml = false; end

if isfield(H.files,'raw') && ~isempty(H.files.raw{1}), H.status.raw = true;
else H.status.raw = false; end

if isfield(H.files,'p0') && ~isempty(H.files.p0{1}), H.status.p0 = true;
else H.status.p0 = false; end %  && ~isempty(H.files.raw{numel(H.sample)})

H.status.rawp0 = H.status.raw && H.status.p0;

if isfield(H.files,'log') && ~isempty(H.files.log{1}), H.status.log = true;
else H.status.log = false; end

if H.repeated_anova && isfield(H.files,'jpg_long') && ~isempty(H.files.jpg_long{1}), H.status.reportlong = true;
else H.status.reportlong = false; end

if isfield(H.files,'jpg') && ~isempty(H.files.jpg{1}), H.status.report = true;
else H.status.report = false; end

H.dpui.report = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpReport,'callback',{@show_report,false},...
  'Style','Pushbutton','enable','off','ToolTipString','Show report file','CData',load_icon('file_cat_report.png'));

H.dpui.reportlong = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpReportLong,'callback',{@show_report,true},...
  'Style','Pushbutton','enable','off','ToolTipString','Show longitudinal report file','CData',load_icon('file_cat_reportlong.png'));

H.dpui.raw = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpRaw,'callback',{@show_raw,false},...
  'Style','Pushbutton','enable','off','ToolTipString','Show raw data','CData',load_icon('file_spm_view.png')); 

H.dpui.rawp0 = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpRawP0,'callback',{@show_raw,true},...
  'Style','Pushbutton','enable','off','ToolTipString','Show overlayed label','CData',load_icon('file_spm_view_p0.png')); 

H.dpui.log = uicontrol(H.mainfig,...
  'Units','normalized','position',H.pos.dpLog,'callback',{@show_log},...
  'Style','Pushbutton','enable','off','ToolTipString','Show log file','CData',load_icon('file_cat_log.png')); 

% add slider and opacity control only for volume data
if ~H.mesh_detected
  H.ui.alpha = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.aslider,...
        'Min',0,'Max',1,...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_alpha,'Value',0.5,...
        'ToolTipString','Change Opacity of pos. (blue colors) and neg. (red colors) Z-scores',...
        'SliderStep',[0.01 0.1],'Visible','off');

  H.ui.alpha_txt = uicontrol(H.mainfig,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[H.pos.aslider(1) H.pos.aslider(2)-0.005 0.2 0.02],...
        'String','Overlay Opacity of Z-score',...
        'FontSize',H.FS-2,'Visible','off');

  H.ui.mm = uicontrol(H.mainfig,...
        'Units','normalized','position',H.pos.zslider,...
        'Min',(1 - H.data.Orig(3))*H.data.vx(3),'Max',(H.files.V(1).dat.dim(3) - H.data.Orig(3))*H.data.vx(3),...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_slices_array,...
        'ToolTipString','Select slice for display',...
        'SliderStep',[0.005 0.05],'Visible','off');

  H.ui.mm_txt = uicontrol(H.mainfig,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[H.pos.zslider(1) H.pos.zslider(2)-0.005 0.2 0.02],...
        'String','Slice [mm]','Visible','off','FontSize',H.FS-2);

end
  
return

%-----------------------------------------------------------------------
function icon = load_icon(name)
%-----------------------------------------------------------------------

icon = imread(fullfile(fileparts(mfilename('fullpath')),'doc','icons',name)); 
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

if isfield(H,'delui') && isfield(H.delui,'remove') set(H.delui.remove,'Enable','off'); end
if isfield(H,'status')
  if H.status.report,     set(H.dpui.report,    'enable','off','BackgroundColor',[0.94 0.94 0.94]); end
  if H.status.raw,        set(H.dpui.raw,       'enable','off','BackgroundColor',[0.94 0.94 0.94]); end
  if H.status.rawp0,      set(H.dpui.rawp0,     'enable','off','BackgroundColor',[0.94 0.94 0.94]); end
  if H.status.log,        set(H.dpui.log,       'enable','off','BackgroundColor',[0.94 0.94 0.94]); end
  if H.status.reportlong, set(H.dpui.reportlong,'enable','off','BackgroundColor',[0.94 0.94 0.94]); end
end

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
  
  % add short name to caption
  for i=1:number
    txt = {{spm_str_manip(list2(i,:),'k30'),sprintf('Mean abs Z-score: %g',H.data.avg_abs_zscore(ind_sorted_decreased(i)))}};
    if H.isxml
      txt{1}{3} = sprintf('SIQR: %g',H.X(ind_sorted_decreased(i),4));
    end
    spm_orthviews('Caption',i,txt);
  end
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
function show_QMzscore(X, sel, quality_order)
%-----------------------------------------------------------------------
global H

if nargin < 3
  quality_order = -1; 
end

% delete old data tip
delete(findall(H.mainfig,'Type','hggroup'))

if ~H.mesh_detected && isfield(H.ui,'alpha')
  set(H.ui.alpha,    'Visible','off');
  set(H.ui.alpha_txt,'Visible','off');
  set(H.ui.mm,       'Visible','off');
  set(H.ui.mm_txt,   'Visible','off');
end

if isfield(H.ui,'alpha'), set(H.ui.text,'Visible','off'); end
if isfield(H,'delui') && isfield(H.delui,'remove') 
  set(H.delui.remove,'Enable','off'); 
end

if isfield(H,'status')
  if H.status.report,     set(H.dpui.report,    'enable','off'); end
  if H.status.raw,        set(H.dpui.raw,       'enable','off'); end
  if H.status.rawp0,      set(H.dpui.rawp0,     'enable','off'); end
  if H.status.log,        set(H.dpui.log,       'enable','off'); end
  if H.status.reportlong, set(H.dpui.reportlong,'enable','off'); end
end

H.sel = sel;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.mainfig);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',H.pos.plot,'Parent',H.mainfig,'Color',[.6 .6 .6]);
axes(H.ax);
grid on

% estimate product between QM-value and quartic mean Z-score
if sel
  if H.userps > 0
    % in case of the percentage scoring, we first have to go back to marks
    % in case of the normalized percentage scoring, no tranthe order changes the default is zero
    if sel == 5
      rps2mark = @(rps)  quality_order * (0 - rps/10); 
    else
      rps2mark = @(rps)  quality_order * (10.5 - rps/10); 
    end
  else
    rps2mark = @(rps)  -quality_order * rps; 
  end
  H.xml.QMzscore = X(:,1) .* max(0,rps2mark( X(:,sel) ));
else
  H.xml.QMzscore = X(:,1);
end

% get min/max in 0.25 steps
min_QMzscore = min(H.xml.QMzscore(H.ind)); min_QMzscore = floor(4*min_QMzscore)/4; 
max_QMzscore = max(H.xml.QMzscore(H.ind)); max_QMzscore = ceil(4*max_QMzscore)/4;  
if min_QMzscore == max_QMzscore
  max_QMzscore = max_QMzscore + 0.5;
end
if sel==5
  % RD202509: in principle the QMs are scaled with .5 for light and 1.0 for clear motion artifacts
  max_QMzscore = max_QMzscore + max(0,min_QMzscore + 1 - max_QMzscore); % grad system
end

% because we use a splitted colormap we have to set the color
% values explicitely
QMzscore_scaled = 63*(H.xml.QMzscore-min_QMzscore) / (max_QMzscore - min_QMzscore); % scale min..max

H.C = zeros(length(H.xml.QMzscore),3);
for i=1:length(H.xml.QMzscore)
  indc = min(128,round(QMzscore_scaled(i))+1);
  H.C(i,:) = H.cmap(indc,:);
end

% create marker for different samples
marker = char('o','s','d','^','v','<','>','.','+','*');
while max(H.sample) > numel(marker), marker = [marker; marker]; end

if sel % show QM measure on x-axis
  xx = X(:,sel);
  yy = X(:,1);
  if H.userps > 0, set(gca, 'XDir','reverse'); end
  if quality_order > 0 % reverse xdir !
    xstr = sprintf('<----- Best ---      %s      --- Worst ------>  ',deblank(H.xml.QM_names(sel-1,:)));
  elseif quality_order < 0 % reverse xdir !
    xstr = sprintf('<----- Worst ---      %s      --- Best ------>  ',deblank(H.xml.QM_names(sel-1,:)));
  else
    xstr = sprintf('%s',deblank(H.xml.QM_names(sel-1,:)));
  end
  
else % show file order on x-axis
  xx = 1:numel(X(:,1));
  yy = X(:,1);
  xstr = sprintf('<----- First ---      File      --- Last ------>  ');
end

% scatterplot
hold on
for i = 1:max(H.sample)
  ind = H.sample.*H.ind == i;
  H.ui.scatter = scatter(H.ax,xx(ind),yy(ind),30,H.C(ind,:),marker(i),'Linewidth',1);
  MarkerEdgeColor = get(H.ui.scatter,'MarkerEdgeColor');
  set(H.ui.scatter,'MarkerFaceColor',MarkerEdgeColor,'MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',0.5);
end

% connect points of each subject for long. designs
if H.repeated_anova
  
  % use SIQR*zscore if available
  if sel == 4
    measure = H.X(:,sel).*H.data.avg_abs_zscore;
  else
    measure = H.data.avg_abs_zscore;
  end
  
  cm = jet(64);
  ind = cell(numel(H.ind_subjects_long),1);
  diff_measure = zeros(numel(H.ind_subjects_long),1);
  
  % get difference between extreme values for each subject
  for i = 1:numel(H.ind_subjects_long)
    ind{i} = H.ind.*H.ind_subjects_long{i} > 0;
    if any(ind{i})
      diff_measure(i) = max(measure(ind{i}))-min(measure(ind{i}));
    end
  end
  
  % scale difference measure to a range 1..64 and use jet-colors for lines
  diff_measure = 1 + round(63*(diff_measure - min(diff_measure))/(max(diff_measure)-min(diff_measure)));
  for i = 1:numel(H.ind_subjects_long)
    hp = plot(xx(ind{i}),yy(ind{i}));
    set(hp, 'Color',cm(diff_measure(i),:),'LineWidth',diff_measure(i)/20);
  end
  title('The lines show the time points of each subject with color and thickness in relation to the magnitude of the differences');
end
hold off

if ~sel
  xlim([0 numel(H.sample)+1]);
end

xlabel(xstr,'FontSize',H.FS-1,'FontWeight','Bold');
ylabel('<----- Best ---      Quartic Mean Z-score      --- Worst ------>  ','FontSize',H.FS-1,'FontWeight','Bold');

% add colorbar
H.ui.cbar = axes('Position',H.pos.cbar+[0 0.9 0 0],'Parent',H.mainfig);
image((1:64));

if sel
  if ~quality_order
    xstr = sprintf('%s x quartic mean Z-score',deblank(H.xml.QM_names(sel-1,:)));
  else
    xstr = sprintf('<----- Best ---      %s (grad) x quartic mean Z-score     --- Worst ------>  ',deblank(H.xml.QM_names(sel-1,:)));
  end
else
  xstr = sprintf('<----- Best ---      quartic mean Z-score     --- Worst ------>  ');
end
title(xstr,'FontSize',H.FS+1,'FontWeight','Bold');

colormap(H.cmap)

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.ui.cbar,'YTickLabel','','YTick','', 'XTick',linspace(1,64,5),'XTickLabel',...
    round(100*linspace(min_QMzscore,max_QMzscore,5))/100,'TickLength',[0 0]);

% update index of worst files
[tmp, H.ind_sorted_display] = sort(H.xml.QMzscore(H.ind),'ascend');

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
set(H.dpui.report,    'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.reportlong,'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.log,       'BackGroundColor',[0.94 0.94 0.94]);
set(H.dpui.raw,       'BackGroundColor',[0.94 0.94 0.94]);

% only use SPM window if not defined
if ~isfield(H,'Fgraph')
  H.Fgraph = spm_figure('GetWin','Graphics');
end

set(H.Fgraph,'Renderer','OpenGL');
figure(H.Fgraph);
spm_figure('Select',H.Fgraph);
clf

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
set(gcf,'Color',[0.94 0.94 0.94])

hold on
for i=1:n_samples
  ind = find(H.sample(H.ind) == i);

  for j=1:length(ind)
    if H.ui.show_name
      text(xpos{i}(j),data{i}(j),H.filename.m{ind(j)},'FontSize',H.FS-2,'HorizontalAlignment','center')
    else
      p = plot(xpos{i}(j),data{i}(j),'.');
      set(p,'Color',0.8*cm(i,:));
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
  while ~isempty(strfind(tmp,',,')), tmp = strrep(tmp,',,',','); end
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
    text(xpos, ylim_min,'<----- Low rating  ','Color','red','Rotation',...
        90,'HorizontalAlignment','left','FontSize',H.FS,'FontWeight','Bold')
    text(xpos, ylim_max,'High rating ------>  ','Color','blue','Rotation',...
        90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
  elseif quality_order < 0
      text(xpos, ylim_max,'Low rating ------>  ','Color','red','Rotation',...
          90,'HorizontalAlignment','right','FontSize',H.FS,'FontWeight','Bold')
      text(xpos, ylim_min,'<----- High rating ','Color','blue','Rotation',...
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
set(H.delui.remove,'Enable','off');

figure(H.mainfig)

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

sz = spm('WinSize','0',1) - H.pos.fig; sz = sz*0.75; sz(3) = sz(4)*1.4;

if ~isfield(H,'mipfig')
  H.mipfig = figure(23);
end

figure(H.mipfig);

cm = hot(64);
set(H.mipfig,'Menubar','none','NumberTitle','off','Name',sprintf('Sample %d: Z-score %s',...
    H.sample(H.mouse.x),H.filename.m{H.mouse.x}),'Position',[10 H.pos.fig(4)+sz(4) sz(3:4)]);
colormap([1-(cm); cm]);

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
function show_mesh
%-----------------------------------------------------------------------
global H

if isfield(H,'hx') && isgraphics(H.hx.figure)
  H.hx = cat_surf_render2('Overlay',H.hx,H.texture(:,H.mouse.x(1)));
else
  H.hx = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',H.texture(:,H.mouse.x(1))));
  H.hx = cat_surf_render2('Colourbar',H.hx);
end
H.hx = cat_surf_render2('clim',H.hx,H.data.range98);

sz = spm('WinSize','0',1) - H.pos.fig; sz = sz*0.75; sz(3) = sz(4)*1.4;
pos_hx = [10 H.pos.fig(4)+sz(4) sz(3:4)];

set(H.hx.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Position',pos_hx,...
  'Name',sprintf('Sample %d: %s %s',H.sample(H.mouse.x),H.info.texture,H.filename.m{H.mouse.x(1)}))

figure(H.hx.figure)

% get Z-score
zscore = (H.texture(:,H.mouse.x(1)) - H.data.Ymean)./H.data.Ystd;
zscore(H.data.Ystd == 0) = 0;

if isfield(H,'hy') && isgraphics(H.hy.figure)
  H.hy = cat_surf_render2('Overlay',H.hy,zscore);
else
  pos_hy = pos_hx;
  pos_hy = pos_hy + [pos_hx(3)+5 0 0 0];
  H.hy = cat_surf_render2(struct('vertices',H.Pmesh.vertices,'faces',H.Pmesh.faces,'cdata',zscore));
  H.hy = cat_surf_render2('Colourbar',H.hy);
  H.hy = cat_surf_render2('ColourMap',H.hy,cat_io_colormaps('BWR',64));
  set(H.hy.figure,'Position',pos_hy);  
end

H.hy = cat_surf_render2('clim',H.hy,[-3 3]);
set(H.hy.figure,'Menubar','none','Toolbar','none','NumberTitle','off','Name',sprintf('Sample %d: Z-score %s',H.sample(H.mouse.x(1)),H.filename.m{H.mouse.x(1)}));  

figure(H.hy.figure)
  
return

%-----------------------------------------------------------------------
function show_image_slice
%-----------------------------------------------------------------------
global H

% add sliders for volume data
set(H.ui.mm,'Visible','on');
set(H.ui.mm_txt,'Visible','on');

% we have to update slice array first if not defined
if ~isfield(H.data,'vol')
  preload_slice_data;
end

H.img = H.data.vol(:,:,H.mouse.x(1))';

% alpha overlay
H.img_alpha = H.data.zscore(:,:,H.mouse.x(1))';

% correct orientation
H.img = rot90(H.img,2);
H.img_alpha = rot90(H.img_alpha,2);

if ~isfield(H,'ax_slice')
  H.ax_slice = axes('Position',H.pos.slice);
else
  axes(H.ax_slice);
end

% display image with 2nd colorbar (gray)
image(65 + H.img);
if ~H.mesh_detected, axis image; end
set(H.ax_slice,'XTickLabel','','YTickLabel','','TickLength',[0 0]);
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

figure(H.mainfig);
colormap(H.cmap)

return

%-----------------------------------------------------------------------
function show_report(obj, event_obj, long_report)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  if isfield(H.ui,'plotbox')
    set(H.ui.plotbox, 'Visible', 'off');
  end
  set(H.ui.fnambox, 'Visible', 'off');
  if long_report
    set(H.dpui.reportlong,'BackGroundColor',[0.95 0.95 0.95]);
    set(H.dpui.report,    'BackGroundColor',[0.94 0.94 0.94]);
    H.show_sel = 3;
  else
    set(H.dpui.report,    'BackGroundColor',[0.95 0.95 0.95]);
    set(H.dpui.reportlong,'BackGroundColor',[0.94 0.94 0.94]);
    H.show_sel = 2;
  end
  set(H.dpui.log,  'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.raw,  'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.rawp0,'BackGroundColor',[0.94 0.94 0.94]);
end

% select first time point for longitudinal report and selected data for
% normal report
if long_report
  jpg_file = H.files.jpg_long{min(H.mouse.x)};
else
  jpg_file = H.files.jpg{H.mouse.x(1)};
end

if ~isempty(jpg_file)
  
  figure(H.Fgraph);
  clf

  ppos = [0 0 1 1];
  jpg  = imread(jpg_file); 
  set(gca,'Position',ppos(1,:));
  gpos = H.Fgraph.Position; 
  [Xq,Yq] = meshgrid(1:size(jpg,2)/gpos(3)/2:size(jpg,2),...
                     1:size(jpg,1)/gpos(4)/2:size(jpg,1));
  jpgi = zeros([size(Xq,1) size(Xq,2) 3],'uint8');
  for i=1:3, jpgi(:,:,i) = uint8(interp2(single(jpg(:,:,i)),Xq,Yq,'linear')); end
  image(H.Fgraph.CurrentAxes,jpgi);
  set(gca,'Visible','off'); 

end

figure(H.mainfig)

return

%-----------------------------------------------------------------------
function show_raw(obj, event_obj, overlay)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  if isfield(H.ui,'plotbox')
    set(H.ui.plotbox, 'Visible', 'off');
  end
  set(H.ui.fnambox, 'Visible', 'off');
  set(H.dpui.report,    'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.reportlong,'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.log,      'BackGroundColor',[0.94 0.94 0.94]);
  if overlay
    set(H.dpui.rawp0,  'BackGroundColor',[0.95 0.95 0.95]);
    set(H.dpui.raw,    'BackGroundColor',[0.94 0.94 0.94]);
    H.show_sel = 5;
  else
    set(H.dpui.rawp0,  'BackGroundColor',[0.94 0.94 0.94]);
    set(H.dpui.raw,    'BackGroundColor',[0.95 0.95 0.95]);
    H.show_sel = 4;
  end
end

x_sort = sort(H.mouse.x);
raw_file = H.files.raw(x_sort);
p0_file  = H.files.p0(sort(H.mouse.x));
if (~isempty(raw_file) && ~isempty(raw_file{1})) ||  ...
   (~isempty(p0_file)  && ~isempty(p0_file{1}))

  job.colormapc = flipud(cat_io_colormaps('BCGWHcheckcov'));
  job.prop  = 0.2; 
 
  if isempty(char(raw_file)) && ~isempty(char(p0_file))
    spm_check_registration(char(p0_file));
  else
    spm_check_registration(char(raw_file));
  end

  % overlay p0image if available
  if ~isempty(char(raw_file)) && overlay
    for i = 1:numel(p0_file)
      if exist(p0_file{i},'file')
        spm_orthviews('addtruecolourimage',i,p0_file{i},...
          job.colormapc,job.prop,0,5);
      end
    end
  end

  % add short name to caption
  for i=1:numel(raw_file)
    txt = {{spm_str_manip(raw_file{i},'k30'),sprintf('Mean abs Z-score: %g',H.data.avg_abs_zscore(x_sort(i)))}};
    if H.isxml
      txt{1}{3} = sprintf('SIQR: %g',H.X(x_sort(i),4));
    end
    spm_orthviews('Caption',i,txt);
  end
  
  spm_orthviews('Reposition',[0 0 0]);
  if isempty(char(raw_file)) && ~isempty(char(p0_file))
    spm_orthviews('Zoom',0);
  end
  spm_orthviews('redraw');  
  
  % make annoying colorbars smaller
  if ~isempty(char(raw_file)) && overlay
    % find last axes that are the colorbars
    ax = findall(gcf,'type','Axes');
    for i = 1:numel(p0_file)
      pos = get(ax(i),'Position');
      set(ax(i),'Position',[pos(1:2) 0.5*pos(3:4)])
    end
  end
  
end

figure(H.mainfig)

return

%-----------------------------------------------------------------------
function show_log(obj, event_obj)
%-----------------------------------------------------------------------
global H

% change button status and checkboxes if button was pressed
if nargin
  H.show_sel = 6;
  if isfield(H.ui,'plotbox')
    set(H.ui.plotbox, 'Visible', 'off');
  end
  set(H.ui.fnambox, 'Visible', 'off');
  set(H.dpui.report,    'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.reportlong,'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.log,       'BackGroundColor',[0.95 0.95 0.95]);
  set(H.dpui.raw,       'BackGroundColor',[0.94 0.94 0.94]);
  set(H.dpui.rawp0,     'BackGroundColor',[0.94 0.94 0.94]);
end

log_file = H.files.log{H.mouse.x};
if ~isempty(log_file)

  figure(H.Fgraph);
  clf
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

figure(H.mainfig)

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H.ui,'alpha')
  H.ui.alphaval = get(H.ui.alpha,'Value');
else
  H.ui.alphaval = 0.5;
end

if ~isfield(H,'ax_slice') H.ax_slice = axes('Position',H.pos.slice); end
axes(H.ax_slice);

% display image with 2nd colorbar (gray)
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
function preload_slice_data
%-----------------------------------------------------------------------
global H

if isfield(H.ui,'mm')
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

M = spm_matrix([0 0 sl]);
H.data.zscore = zeros([P(1).dat.dim(1:2) length(H.files.V)]);
Ymean = reshape(H.data.Ymean,P(1).dat.dim(1:3));
Ystd  = reshape(H.data.Ystd,P(1).dat.dim(1:3));
Ymean = Ymean(:,:,sl);
Ystd  = Ystd(:,:,sl);

scl = max(Ymean(:));

% load slices and show progress bar only for large samples 
if length(H.files.V) > 500, cat_progress_bar('Init',length(H.files.V),'Load slices'); end
for i = 1:length(H.files.V)
  img(:,:) = single(P(i).dat(:,:,sl));
  img(~isfinite(img)) = 0;
  
  H.data.vol(:,:,i) = img;
  if length(H.files.V) > 500, cat_progress_bar('Set',i); end
end
if length(H.files.V) > 500, cat_progress_bar('Clear'); end

% calculate individual Z-score map
for i=1:size(H.data.zscore,3)
  img = H.data.vol(:,:,i);
  ind = Ystd > 0 & (Ymean > H.data.global | img > H.data.global);
  img(ind) = (img(ind) - Ymean(ind))./Ystd(ind);
  img(~ind) = 0;
  H.data.zscore(:,:,i) = img;
end

% enhance contrast and scale image to 0..64
H.data.vol = 64*((H.data.vol - H.data.range98(1))/(H.data.range98(2)-H.data.range98(1)));
H.data.vol(H.data.vol > 64) = 64;
H.data.vol(H.data.vol < 0)  = 0;

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global H

if isfield(H.ui,'mm')
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

M = spm_matrix([0 0 sl]);
H.data.zscore = zeros([P(1).dat.dim(1:2) length(H.files.V)]);
Ymean = reshape(H.data.Ymean,P(1).dat.dim(1:3));
Ystd  = reshape(H.data.Ystd,P(1).dat.dim(1:3));
Ymean = Ymean(:,:,sl);
Ystd  = Ystd(:,:,sl);

scl = max(Ymean(:));

% load slices and show progress bar only for large samples 
if length(H.files.V) > 500, cat_progress_bar('Init',length(H.files.V),'Load slices'); end
for i = 1:length(H.files.V)
  img(:,:) = single(P(i).dat(:,:,sl));
  img(~isfinite(img)) = 0;
  
  H.data.vol(:,:,i) = img;
  if length(H.files.V) > 500, cat_progress_bar('Set',i); end
end
if length(H.files.V) > 500, cat_progress_bar('Clear'); end

% calculate individual Z-score map
for i=1:size(H.data.zscore,3)
  img = H.data.vol(:,:,i);
  ind = Ystd > 0 & (Ymean > H.data.global | img > H.data.global);
  img(ind) = (img(ind) - Ymean(ind))./Ystd(ind);
  img(~ind) = 0;
  H.data.zscore(:,:,i) = img;
end

% enhance contrast and scale image to 0..64
H.data.vol = 64*((H.data.vol - H.data.range98(1))/(H.data.range98(2)-H.data.range98(1)));
H.data.vol(H.data.vol > 64) = 64;
H.data.vol(H.data.vol < 0)  = 0;

if isfield(H,'mouse') && isfield(H.mouse,'x')
  x = H.mouse.x(1);
  
  % check whether mouse position is defined
  H.img       = H.data.vol(:,:,x)';
  H.img_alpha = H.data.zscore(:,:,x)';
  
  % correct orientation
  H.img       = rot90(H.img,2);
  H.img_alpha = rot90(H.img_alpha,2);
  
  if ~isfield(H,'ax_slice')
    H.ax_slice = axes('Position',H.pos.slice);
  else
    axes(H.ax_slice);
  end
  
  % use gray scale colormap for values > 64
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

  set(H.ui.text,'String',txt,'FontSize',H.FS);
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

% create new list of xml-files if necessary
if H.isxml && ((iscell(job.data_xml) && ~isempty(job.data_xml{1})) || (~iscell(job.data_xml) && ~isempty(job.data_xml)))
  n = 0;
  for i=1:numel(job.data_xml)
    if any(Hdel == i), continue; end
    n = n + 1;
    data_xml{n,1} = job.data_xml{i};
  end
  job.data_xml = data_xml;  
end

% create new list without removed data in each sample
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
  data_sel(del_list) = [];
  job.data{i} = data_sel;
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
    
    if isfield(H.job,'factorial_design')
      modify_factorial_design(job.data);
    end

  case 1,
    do_rerun(obj,event_obj,false);
    set(H.naviui.select,'BackGroundColor',[0.95 0.95 0.95]);
    set(H.delui.new,'enable','off');
    datacursormode('on');
end

return

%-----------------------------------------------------------------------
function modify_factorial_design(data)
%-----------------------------------------------------------------------
global H

job = H.job.factorial_design;

% modify dir
[pth,name,ext] = fileparts(char(job.dir));
job.dir{1} = fullfile(pth,['wo_removed_data_' name ext]);
fprintf('\n------------------------------------------------------------------------------------------\n');
fprintf('Create new analysis without removed data in %s\n',job.dir{1});
fprintf('------------------------------------------------------------------------------------------\n');

% modify globals
if isfield(job,'globals') && isfield(job.globals,'g_user')
  job.globals.g_user.global_uval = job.globals.g_user.global_uval(H.ind);
end
if isfield(job,'globalc') && isfield(job.globalc,'g_user')
  job.globalc.g_user.global_uval = job.globalc.g_user.global_uval(H.ind);
end
if isfield(job,'globals') && isfield(job.globals,'g_ancova')
  job.globals.g_ancova.global_uval = job.globals.g_ancova.global_uval(H.ind);
end
if isfield(job,'globalc') && isfield(job.globalc,'g_ancova')
  job.globalc.g_ancova.global_uval = job.globalc.g_ancova.global_uval(H.ind);
end

% modify covariates
if isfield(job,'cov') 
  for i=1:numel(job.cov)
    job.cov(i).c = job.cov(i).c(H.ind);
  end
end

% modify files and factors for different designs
if isfield(job.des,'t2') % two-sample t-test
    job.des.t2.scans1 = data{1};
    job.des.t2.scans2 = data{2};
elseif isfield(job.des,'mreg') % multiple regression
    job.des.mreg.scans = data{1};
elseif isfield(job.des,'fd') % full factorial
  for i=1:numel(job.des.fd.icell)
    job.des.fd.icell(i).scans = data{i};
  end
elseif isfield(job.des,'fblock') % flexible factorial
  fsubject = job.des.fblock.fsuball.fsubject;
  
  % index of whole subjects is removed from list
  ind_remove_subject = [];
  
  % go through all subjects
  for i = 1:numel(H.ind_subjects_long)

    % index where time points are defined for this subject
    ind_subject = find(H.ind_subjects_long{i});
    
    % array where data are kept for this subject
    ind = H.ind.*H.ind_subjects_long{i} > 0;
    
    % we can keep that subject if we have at least 2 time points
    if sum(ind) > 1 % remove single time points and update scans and conditions
      job.des.fblock.fsuball.fsubject(i).scans = fsubject(i).scans(ind(ind_subject));
      if size(fsubject(i).conds,1) > 1
        job.des.fblock.fsuball.fsubject(i).conds = fsubject(i).conds(ind(ind_subject),:);
      else
        job.des.fblock.fsuball.fsubject(i).conds = fsubject(i).conds(ind(ind_subject));
      end
    elseif sum(ind) == 1 % indicate to remove whole subject because only one time point remains
      ind_remove_subject = [ind_remove_subject i];
      fprintf('Remove all time points of subject %d because only one time point remains.\n',i);
    elseif sum(ind) == 0 % indicate to remove whole subject
      ind_remove_subject = [ind_remove_subject i];
      fprintf('Remove all time points of subject %d\n',i);
    end
  end  
  fprintf('\n');
  
  % remove whole subject from list
  if ~isempty(ind_remove_subject)
    job.des.fblock.fsuball.fsubject(ind_remove_subject) = [];
  end
else
  fprintf('Other designs are not yet prepared.\n');
  return
end

% Sometimes window for DesRep is not accessable
try
  out = spm_run_factorial_design(job);
catch
  out = spm_run_factorial_design(job);
end
spm('alert!', sprintf('Create new analysis without removed data in %s\n',job.dir{1}), 0);

return

%-----------------------------------------------------------------------
function do_rerun(obj, event_obj, undo)
%-----------------------------------------------------------------------
global H

if H.status.report,     set(H.dpui.report,    'enable','on','BackGroundColor',[0.94 0.94 0.94]); end
if H.status.raw,        set(H.dpui.raw,       'enable','on','BackGroundColor',[0.94 0.94 0.94]); end
if H.status.rawp0,      set(H.dpui.rawp0,     'enable','on','BackGroundColor',[0.94 0.94 0.94]); end
if H.status.log,        set(H.dpui.log,       'enable','on','BackGroundColor',[0.94 0.94 0.94]); end
if H.status.reportlong, set(H.dpui.reportlong,'enable','on','BackGroundColor',[0.94 0.94 0.94]); end

set(H.naviui.select,  'BackGroundColor',[0.95 0.95 0.95]);

datacursormode('on');

if undo
  H.ind = true(size(H.sample));
  set(H.delui.undo,    'enable','off');
  set(H.delui.new,     'enable','off');
  set(H.delui.list_del,'enable','off');
  set(H.delui.analysis_new,'enable','off');
  H.del = [];
else
  
  % remove subjects where only one time point remains
  if H.repeated_anova
    for i = 1:numel(H.ind_subjects_long)

      % array where data are kept for this subject
      ind = H.ind.*H.ind_subjects_long{i} > 0;

      if sum(ind) == 1 % indicate to remove whole subject because only one time point remains
        H.del = unique([H.del find(H.ind_subjects_long{i})]);
        fprintf('Remove all time points of subject %d because only one time point remains.\n',i);
      end
    end
  end
  
  H.ind = ~ismember((1:numel(H.sample)),H.del);
end

show_boxplot(H.data.avg_abs_zscore(H.ind),'Quartic Mean Z-score  ',-1);  
if H.isxml
  show_QMzscore(H.X,4);
else
  show_QMzscore(H.X,0);
end

% delete old data tip
delete(findall(H.mainfig,'Type','hggroup'))
set(H.delui.remove,'enable','off');

return

%-----------------------------------------------------------------------
function remove_point(obj, event_obj)
%-----------------------------------------------------------------------
global H

if H.sel
  sel = H.sel;
else
  sel = 0; % file order by default
end

x = H.mouse.x(1);

% check whether we also have to exlude other time points for long. data
% if there are not enough timepoints anymore available
if H.repeated_anova
  for i = 1:numel(H.ind_subjects_long)
    ind = find(H.ind_subjects_long{i});
    if any(ismember(ind, x))
      if numel(ind) - 1 == sum(ismember(ind,unique([H.del x])))
        x = H.mouse.x;
      end
    end
  end
end

if sel % QM measure on x-axis
  xx = H.X(x,sel);
  yy = H.X(x,1);
else % file order on x-axis
  xx = x;
  yy = H.X(x,1);
end

axes(H.ax)
hold on

% reconsider this data point (and only this one) if already in the list
if ~isempty(H.del) && any(ismember(H.del,x(1)))
  plot(xx(1),yy(1),'wx','MarkerSize',15,'Linewidth',2);
  plot(xx(1),yy(1),'wo','MarkerSize',15,'Linewidth',2,'MarkerFaceColor','w');
  plot(xx(1),yy(1),'ko','MarkerSize',5,'Linewidth',2);
  H.del(ismember(H.del,x(1))) = [];
  return
else
  plot(xx,yy,'rx','MarkerSize',15,'Linewidth',2);
end
hold off

% add point to the list
H.del = unique([H.del x]);

% also update index of considered data and enable icons
H.ind = ~ismember((1:numel(H.sample)),H.del);
set(H.delui.undo,        'enable','on');
set(H.delui.new,         'enable','on');
set(H.delui.list_del,    'enable','on');
set(H.delui.analysis_new,'enable','on');

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global H

if H.sel
  sel = H.sel;
else
  sel = 0; % file order by default
end

% if it's the same point use old data tip and return
pos_mouse = get(event_obj, 'Position');
if H.mouse.xold == pos_mouse
  txt = get(findall(H.mainfig,'Type','hggroup'),'String');
  if isempty(txt), return; end
end

H.mouse.xold = pos_mouse;

x = find(H.X(:,1) == pos_mouse(2));
if isempty(x) || numel(x) > 1
  x = find(H.X(:,sel) == pos_mouse(1));
end

% empty data tip and return if no point was found
if isempty(x)
  txt = '';
  return
end

% text info for data cursor window
txt = {sprintf('%s',H.filename.m{x})};

% prevent that that function is called again if position has not changed or
% subject for long. data has not changed for showing raw data
if H.repeated_anova && (isfield(H,'show_sel') && (H.show_sel == 4 || H.show_sel == 5))
  if any(x == H.mouse.x)
    return
  else
    H.mouse.x = x;
  end
elseif x == H.mouse.x(1) % && (H.repeated_anova && (H.show_sel == 4 || H.show_sel == 5))
  return
else 
  H.mouse.x = x;
end

if H.mesh_detected 
  % show two render views for meshes: texture and Z-score
  show_mesh;
else
  % show image slice
  show_image_slice;
end

if ~H.mesh_detected
  set(H.ui.alpha,'Visible','on');
  set(H.ui.alpha_txt,'Visible','on');
end

% text info for textbox
txt2 = {[],sprintf('%s',spm_file(H.filename.m{H.mouse.x(1)},'short25'))};
if max(H.sample) > 1, txt2{end+1} = sprintf('Sample %d',H.sample(H.mouse.x(1))); end
if max(H.sites)>1, txt2{end} = sprintf('%s / Site %d', txt2{end}, H.sites(H.mouse.x(1))); end % just add the site variable if useful
txt2{end+1} = sprintf('Mean abs Z-score: %g',H.data.avg_abs_zscore(H.mouse.x(1)));

if H.isxml
  txt2{end+1} = sprintf('SIQR / nSIQR: %g / %g', H.X(H.mouse.x(1),4), H.X(H.mouse.x(1),5) );
end

if ~H.mesh_detected
  str = {[],'Individual Z-score','red: value < mean','blue: value > mean'};
  for i=1:numel(str)
    txt2{end+1} = str{i};
  end
end

set(H.ui.text,'String',txt2,'FontSize',H.FS);

% get list of time points for long. data
if H.repeated_anova
  for i = 1:numel(H.ind_subjects_long)
    if H.ind_subjects_long{i}(x)
      H.mouse.x = find(H.ind_subjects_long{i});
      
      % set actual position to 1st entry
      H.mouse.x(H.mouse.x == x) = [];
      H.mouse.x = [x H.mouse.x];
      break
    end
  end
end

% only enable select
% does not work and I am not sure whether this is necessary at all
%datacursormode('on');

% enable buttons
if strcmp(get(H.delui.remove,'enable'),'off')
  set(H.delui.remove,'enable','on');

  if H.status.report,     set(H.dpui.report,    'enable','on'); end
  if H.status.raw,        set(H.dpui.raw,       'enable','on'); end
  if H.status.rawp0,      set(H.dpui.rawp0,     'enable','on'); end
  if H.status.log,        set(H.dpui.log,       'enable','on'); end
  if H.status.reportlong, set(H.dpui.reportlong,'enable','on'); end
  
  if ~H.mesh_detected
    set(H.ui.alpha,    'Visible','on');
    set(H.ui.alpha_txt,'Visible','on');
    set(H.ui.mm,       'Visible','on');
    set(H.ui.mm_txt,   'Visible','on');
  end
  set(H.ui.text,       'Visible','on');

end


if isfield(H,'show_sel')
  switch(H.show_sel)
    case 2, show_report(obj, event_obj, false); % report
    case 3, show_report(obj, event_obj, true);  % report long
    case 4, show_raw(obj, event_obj, false);    % raw
    case 5, show_raw(obj, event_obj, true);     % raw + p0
    case 6, show_log;                           % log file
  end
end

return
