function cat_stat_analyze_ROIs
% Statistical analysis of ROI data using an existing SPM design saved in a SPM.mat file
%
%__________________________________________________________________________
% Christian Gaser
% $Id$

spmmat = spm_select(1,'SPM.mat','Select SPM.mat file to get design');
load(spmmat);

cwd = fileparts(spmmat);

%-Check that model has been estimated
if ~isfield(SPM, 'xVol')
    str = { 'This model has not been estimated.';...
              'Would you like to estimate it now?'};
    if spm_input(str,1,'bd','yes|no',[1,0],1)
        cd(cwd)
        SPM = spm_spm(SPM);
    else
        return
    end
end

% select contrast
[Ic,xCon] = spm_conman(SPM,'T&F',1,'Select contrast',' ',1);

% threshold for p-values
alpha = spm_input('p value',1,'r',0.05,1,[0,1]);

% mesh detected?
if isfield(SPM.xVol,'G')
  mesh_detected = 1;
  pattern = 'catROIs_';
  subfolder = 'surf';
else
  mesh_detected = 0;
  pattern = 'catROI_';
  subfolder = 'mri';
end

% filenames of existing design
P = SPM.xY.P;

% get names of ROI xml files using saved filename in SPM.mat
for i=1:numel(P)
  [pth,nam,ext] = fileparts(P{i});
  
  % check whether a label subfolder exists and replace the name with "label"
  if strcmp(pth(end-length(subfolder)+1:end),subfolder)
    pth(end-length(subfolder)+1:end) = [];
    pth_label = [pth 'label'];
  else
    pth_label = pth;
  end
  
  % check that name of ROI fits to SPM data for 1st file
  if i==1
    % check for catROI*-files
    files = cat_vol_findfiles(pth_label,[pattern '*']);
    if numel(files) == 0
      error(sprintf('No label files found in folder %s. Please check that you have not moved your data\n',pth_label));
    end
    [tmp, tmp_name, ext] = fileparts(files{1});
    
    % remove catROI[s]_ from first xml file
    tmp_name(1:length(pattern)) = [];

    % check wheter first filename in SPM.mat and xml-file are from the same subject
    ind = strfind(nam,tmp_name);
    if isempty(ind)
      error(sprintf('Label file %s does nit fit to analyzed file %s. Please check that you have not moved your data\n'),tmp_name,nam);
    end
    
    % get prepending pattern 
    if ~mesh_detected
      prepend = nam(1:ind-1);
      % check whether prepending pattern is not based on GM/WM
      if isempty(strfind(prepend,'mwp')) & isempty(strfind(prepend,'m0wp'))
        fprintf('\nWARNING: ROI analysis is only supported for VBM of GM/WM/CSF. No ROI values for DBM will be estimated.\n',prepend);
      end
    end

  end
  
  % get ROI name for all files and check whether the files are found
  roi_names{i} = fullfile(pth_label,[pattern nam(ind:end) '.xml']);
  if ~exist(roi_names{i},'file')
      error(sprintf('Label file %s not found. Please check that you have not moved your data\n',roi_names{i}));
  end
end

% use 1st xml file to get the available atlases
xml = convert(xmltree(roi_names{1}));
if ~isfield(xml,'ROI')
  error('XML file %s contains no ROI information.',roi_names{1});
end

atlases = fieldnames(xml.ROI);
n_atlases = numel(atlases);

% select one atlas
sel_atlas = spm_input('Select atlas','+1','m',atlases);
atlas = atlases{sel_atlas};

% measure names to search for
measure_names = char('Vgm','Vwm','Vcsf','mean_thickness','mean_fractaldimension','mean_amc','mean_gyrification','mean_sqrtsulc');
n_measure_names = size(measure_names,1);

% get header of selected atlas
measures = fieldnames(xml.ROI.(atlas));
if ~isfield(xml.ROI.(atlas),'tr')
  n_measures = numel(measures);
  if ~isfield(xml.ROI.(atlas).(measures{1}),'tr')
    error('Missing mandatory tr-field in XML file.');
  end
else n_measures = 1; end

% filed name for measure can be "tr" for surfaces
if strcmp(measures,'tr')
  hdr = xml.ROI.(atlas).tr{1}.td;

  found_measure_names = [];
  ind_measure_names = [];
  for k=1:numel(hdr)
    for l=1:n_measure_names
      if strcmp(hdr{k},deblank(measure_names(l,:)))
        found_measure_names = char(found_measure_names,hdr{k});
        ind_measure_names = [ind_measure_names k];
      end
    end
  end
  % remove 1st empty entry
  found_measure_names = found_measure_names(2:end,:);
else
  found_measure_names = measures;
end

% select a measure
if size(found_measure_names,1) > 1
  sel_measure = spm_input('Select measure','+1','m',found_measure_names);
else
  sel_measure = 1;
end
measure = found_measure_names(sel_measure,:);

% for surfaces add prepending "mean_" to field name
if mesh_detected, measure = ['mean_' measure{1}]; end

% remove spaces
measure = deblank(measure);

% get names IDs and values of selected atlas and measure inside ROIs
[ROInames ROIids ROIvalues] = get_ROI_measure(roi_names, sel_atlas, measure);

% get name of contrast
str_con = deblank(xCon(Ic).name);

% replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
str_con(strfind(str_con,' ')) = '_';
strpos = strfind(str_con,' > ');
if ~isempty(strpos), str_con = [str_con(1:strpos-1) '_gt_' str_con(strpos+1:end)]; end
strpos = strfind(str_con,' < ');
if ~isempty(strpos), str_con = [str_con(1:strpos-1) '_lt_' str_con(strpos+1:end)]; end
strpos = strfind(str_con,'>');
if ~isempty(strpos), str_con = [str_con(1:strpos-1) 'gt' str_con(strpos+1:end)]; end
strpos = strfind(str_con,'<');
if ~isempty(strpos), str_con = [str_con(1:strpos-1) 'lt' str_con(strpos+1:end)]; end
str_con = spm_str_manip(str_con,'v');

% build X and Y for GLM
Y = ROIvalues;
X = SPM.xX.X;

% get number of structures
n_structures = size(Y,2);

% estimate GLM and get p-value
p = estimate_GLM(Y,X,SPM,Ic);

% divide p-values and names into left and right hemisphere data
P{1}  = p(1:2:n_structures);
P{2}  = p(2:2:n_structures);
N{1}  = ROInames(1:2:n_structures);
N{2}  = ROInames(2:2:n_structures);
ID{1} = ROIids(1:2:n_structures);
ID{2} = ROIids(2:2:n_structures);

% select order according to label names to have left hemipshere always first
if strcmp(N{1}{1}(:,1),'l') % left hemisphere?
  order = [1 2];
else
  order = [2 1];
end

% prepare corrections for multiple comparisons
% corr = {'uncorrected','FDR corrected','Holm-Bonferroni corrected'};
% corr_short = {'','FDR','Holm'};
corr = {'uncorrected','FDR corrected'};
corr_short = {'','FDR'};
n_corr = numel(corr);
data = cell(size(corr));
Pcorr = cell(size(corr));

% go through left and right hemisphere with left first
for i = order

  % left hemisphere?
  if strcmp(N{i}{1}(:,1),'l')
    hemi = 'lh';
  else
    hemi = 'rh';
  end
  
  % select volume atlas only for the 1st hemisphere
  if ~mesh_detected & i==1
    V = spm_vol(fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm',[atlas '.nii']));
    data0 = round(spm_data_read(V));

    % create empty output data
    for c=1:n_corr
      data{c} = zeros(size(data0));
    end
  end

  % select surface atlas for each hemisphere
  if mesh_detected
    atlas_name = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces',...
        [hemi '.' atlas '.freesurfer.annot']);
    [vertices, rdata0, colortable, rcsv0] = cat_io_FreeSurfer('read_annotation',atlas_name);
    data0 = round(rdata0);

    % create empty output data
    for c=1:n_corr
      data{c} = zeros(size(data0));
    end
  end
  
  % sort p-values for FDR and sorted output
  [Psort, indP] = sort(P{i});

  % uncorrected p-values
  Pcorr{1} = P{i};
  
  % apply FDR correction
  Pfdr_sorted = spm_P_FDR(Psort);
  Pcorr{2}(indP) = Pfdr_sorted;

  % apply Holm-Bonferroni correction: correct lowest P by n, second lowest by n-1...
  if n_corr > 2
      n = length(Psort);
    Pcorr{3} = ones(size(P{i}));
    for k=1:n
      Pval = P{i}(indP(k))*(n+1-k);
      if Pval<alpha
        Pcorr{3}(indP(k)) = Pval;
      else
        % stop here if corrected p-value exceeds alpha
        break
      end
    end
  end
  
  output_name = [num2str(100*alpha) '_' str_con '_' atlas '_' measure];
  
  % display and save thresholded sorted p-values for each correction
  for c = 1:n_corr
    ind = find(Pcorr{c}(indP)<alpha);
    ind_corr{c} = ind;
    if ~isempty(ind)
      fprintf('\n%s (P<%g, %s):\n',hemi,alpha,corr{c});
      fprintf('P-value\t\t%s\n',atlas);
      for j=1:length(ind)
        data{c}(data0 == ID{i}(indP(ind(j)))) = -log10(Pcorr{c}(indP(ind(j))));
        fprintf('%g\t%s\n',Pcorr{c}(indP(ind(j))),N{i}{indP(ind(j))});
      end
    end
    
    % write label surface with thresholded p-values
    if mesh_detected
      % save P-alues as float32
      fname = [hemi '.logP' corr_short{c} output_name '.gii'];
      save(gifti(struct('cdata',data{c})),fname);
      fprintf('\nLabel file with thresholded logP values (%s) saved as %s.',corr{c},fname);
    end

  end
  fprintf('\n');
  
end

% prepare display ROI results according to found results
corr{n_corr+1} = 'Do not display';
ind_show = [];
for c=1:n_corr
  if ~isempty(ind_corr{c})
    ind_show = [ind_show c];
  end
end

% display ROI results
if ~isempty(ind_show)
  show_results = spm_input('Display ROI results?','+1','m',corr([ind_show n_corr+1]),[ind_show 0]);
else
  show_results = 0;
end

% write label volume with thresholded p-values
if ~mesh_detected

  if isempty(ind_show)
    fprintf('No results found.\n');
  end

  % go through all corrections and save label image if sign. results were found
  for c=1:n_corr 
    if ~isempty(ind_corr{c})
      V.fname = ['logP' corr_short{c} output_name '.nii'];
      V.dt(1) = 16;
      spm_write_vol(V,data{c});
      fprintf('\nLabel file with thresholded logP values (%s) was saved as %s.',corr{c},V.fname);
    end
  end
  fprintf('\n');
  
  % display ROI results for label image
  if show_results
    % display image as overlay
    OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152.nii');
    OV.reference_range = [0.2 1.0];                         % intensity range for reference image
    OV.opacity = Inf;                                      % transparence value for overlay (<1)
    OV.cmap    = jet;                                      % colormap for overlay
    OV.range   = [-log10(alpha) round(max(data{show_results}(:)))];
    OV.name = ['logP' corr_short{show_results} output_name '.nii'];
    OV.slices_str = char('-30:4:60');
    OV.transform = char('axial');
    cat_vol_slice_overlay(OV);
  end
  
end

% surface results display
if mesh_detected
  if ~isempty(ind_show)
    fprintf('\nLabels for both hemispheres are saved, thus it is not necessary to also estimate labels from the SPM.mat file of the opposite hemisphere.\n');
  else
    fprintf('No results found.\n');
  end
  
  % remove both hemisphere results if no results were significant
  for c=1:n_corr 
    if isempty(ind_corr{c})
      spm_unlink(['lh.logP' corr_short{c} output_name '.gii']);
      spm_unlink(['rh.logP' corr_short{c} output_name '.gii']);
    end
  end
    
  % display ROI surface results
  if show_results
    lh = ['lh.logP' corr_short{show_results} output_name '.gii'];
    rh = ['rh.logP' corr_short{show_results} output_name '.gii'];
    cat_surf_results('Disp',lh,rh,0);
  end
end

%_______________________________________________________________________
function p = estimate_GLM(Y,X,SPM,Ic);
% estimate GLM and return p-value for F- or T-test
%
% FORMAT p = estimate_GLM(Y,X,SPM,Ic);
% Y   - data matrix
% X   - design matrix
% SPM - SPM structure
% Ic  - selected contrast
% p   - returned p-value

df   = [SPM.xCon(Ic).eidf SPM.xX.erdf];
c = SPM.xCon(Ic).c;
n = size(Y,1);
n_structures = size(Y,2);

% estimate statistics for F- or T-test
if strcmp(SPM.xCon(Ic).STAT,'F')
  c0 = eye(size(X,2)) - c*pinv(c);
  Xc = X*c;
  X0 = X*c0;

  R  = eye(n) - X*pinv(X);
  R0 = eye(n) - X0*pinv(X0);
  M = R0 - R;

  pKX = pinv(X);
  trRV = n - rank(X);
  p = rank(X);
  p1 = rank(Xc);

  Beta = pKX * Y;

  yhat = X*Beta;
  F = zeros(n,1);
  for i=1:n_structures
    F(i) = (yhat(:,i)'*M*yhat(:,i))*(n-p)/((Y(:,i)'*R*Y(:,i))*p1);
  end
  F(find(isnan(F))) = 0;

  p = 1-spm_Fcdf(F,df);
else

  pKX = pinv(X);
  trRV = n - rank(X);

  Beta = pKX * Y;
  ResSS = sum((X*Beta - Y).^2);

  ResMS = ResSS/trRV;

  con = (c'*Beta);
  Bcov = pinv(X'*X);

  ResSD = sqrt(ResMS.*(c'*Bcov*c));
  t = con./(eps+ResSD);
  t(find(isnan(t))) = 0;

  p = 1-spm_Tcdf(t,df(2));
end

%_______________________________________________________________________
function [ROInames ROIids ROIvalues] = get_ROI_measure(roi_names, sel_atlas, measure_name)
% get names, IDs and values inside ROI for a selected atlas
%
% FORMAT [ROInames ROIids ROIvalues] = get_ROI_measure(roi_names, sel_atlas, measure_name);
% roi_names    - cell of ROI xml files
% sel_atlas    - index of selected atlas
% measure_name - name of selected measure
%
% ROInames     - cell rx1 of ROI names (r - # of ROIs)
% ROIids       - cell rx2 of ROI IDs for left and right hemipshere
% ROIvalues    - cell nxr of values inside ROI (n - # of data)

n_data = length(roi_names);

spm_progress_bar('Init',n_data,'Load xml-files','subjects completed')
for i=1:n_data        

  xml = convert(xmltree(deblank(roi_names{i})));

  if ~isfield(xml,'ROI')
    error('XML file contains no ROI information.');
  end

  % remove leading catROI*_ part from name
  [path2, ID] = fileparts(roi_names{i});
  ind = strfind(ID,'_');
  ID = ID(ind(1)+1:end);

  atlases = fieldnames(xml.ROI);
  
  for j=sel_atlas
    measures = fieldnames(xml.ROI.(atlases{j}));
    if ~isfield(xml.ROI.(atlases{j}),'tr')
      n_measures = numel(measures);
      if ~isfield(xml.ROI.(atlases{j}).(measures{1}),'tr')
        error('Missing mandatory tr-field in XML file.');
      end
    else n_measures = 1; end
    
    for m=1:n_measures
      
      try
        tr = xml.ROI.(atlases{j}).(measures{m}).tr;
      catch
        tr = xml.ROI.(atlases{j}).tr;
      end
      
      n_ROIs = numel(tr) - 1; % ignore header
      hdr = tr{1}.td;
    
      for k=1:numel(hdr)

        % check for field with ROI names
        if strcmp(hdr{k},'ROIappr') || strcmp(hdr{k},'ROIabbr') || strcmp(hdr{k},'lROIname') || strcmp(hdr{k},'rROIname')
          name_index = k;
        end

        if strcmp(hdr{k},'ROIid')
          ROIid_tmp = k;
        end

        % look for pre-defined ROI measures
        if strcmp(hdr{k},measure_name)
        
          for r=1:n_ROIs
            ROInames{r}    = char(tr{r+1}.td(name_index));
            ROIids(r,1)    = str2num(char(tr{r+1}.td(ROIid_tmp)));
            ROIvalues(i,r) = str2num(char(tr{r+1}.td(k)));
          end
          
        end
                
      end
    end
  end
  spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');
