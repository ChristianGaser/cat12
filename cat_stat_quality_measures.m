function cat_stat_quality_measures(job)
%To check Z-score across sample and save quality 
% measures in csv file.
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%
% varargout = cat_stat_quality_measures(job)
%  
% job                .. SPM job structure
%  .data             .. volume files
%  .globals          .. global scaling
%  .csv_name         .. csv output name
%
% Example: 
%   cat_stat_quality_measures(struct('data',{{ files }},'globals',1,'csv_name','test.csv'));
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

n_subjects = 0;
sample     = [];

% read header
n_subjects = numel(job.data);
mesh_detected = spm_mesh_detect(char(job.data{1}));

% use faster nifti function for reading data
if mesh_detected
  V = spm_data_hdr_read(char(job.data));
else
  V = nifti(char(job.data));
end

pth = spm_fileparts(job.data{1});
report_folder = fullfile(spm_fileparts(pth),'report');
subfolder = 1;
% check whether report subfolder exists
if ~exist(report_folder,'dir')
  report_folder = pth;
  subfolder = 0;
end

isxml = 0;
% search xml report files
xml_files = spm_select('List',report_folder,'^cat_.*\.xml$');
if ~isempty(xml_files)
  fprintf('Search xml-files\n');
  % find part of xml-filename in data files to get the prepending string
  % (e.g. mwp1)
  i = 1; j = 1;
  while i <= n_subjects
    while j <= size(xml_files,1)
      % remove "cat_" and ".xml" from name
      fname = deblank(xml_files(j,:));
      fname = fname(5:end-4);
      
      % and find that string in data filename
      ind = strfind(job.data{i},fname);
      if ~isempty(ind)
        [pth, prep_str] = spm_fileparts(job.data{1}(1:ind-1));
        isxml = 1;
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

% check for global scaling with TIV
if job.globals
  if mesh_detected
    is_gSF = false;
    fprintf('Disabled global scaling with TIV, because this is not meaningful for surface data.\n');
  else
    if isxml
      is_gSF = true;
      gSF = ones(n_subjects,1);
    else
      is_gSF = false;
      fprintf('No xml-files found. Disable global scaling with TIV.\n');
    end
  end
else
  is_gSF = false;
end

if isxml
  if mesh_detected
    QM = ones(n_subjects,5);
    QM_names = char('Noise','Bias','Weighted overall image quality (IQR)','Euler number','Size of topology defects');
  else
    QM = ones(n_subjects,3);
    QM_names = char('Noise','Bias','Weighted overall image quality (IQR)');
  end
  
  cat_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
    
    % get basename for data files
    [pth, data_name] = fileparts(job.data{i});
    
    % remove ending for rigid or affine transformed files
    data_name = strrep(data_name,'_affine','');
    data_name = strrep(data_name,'_rigid','');

    % get report folder
    if subfolder
      report_folder = fullfile(spm_fileparts(pth),'report');
    else
      report_folder = pth;
    end
    
    % remove prep_str from name and use report folder and xml extension
    if mesh_detected
      % for meshes we also have to remove the additional "." from name
      tmp_str = strrep(data_name,prep_str,'');
      xml_file = fullfile(report_folder,['cat_' tmp_str(2:end) '.xml']);
    else
      xml_file = fullfile(report_folder,['cat_' strrep(data_name,prep_str,'') '.xml']);
    end
    
    if ~exist(xml_file,'file')
      isxml = 0;
      fprintf('Cannot use quality ratings, because xml-file %s was not found\n',xml_file);
      break
    end
    
    if exist(xml_file,'file')
      xml = cat_io_xml(xml_file);
    else
      fprintf('File %s not found. Skip use of xml-files for quality measures.\n',xml_file);
      isxml = 0;
      break
    end
    
    % get TIV
    if is_gSF && isfield(xml,'subjectmeasures') && isfield(xml.subjectmeasures,'vol_TIV')
      gSF(i) = xml.subjectmeasures.vol_TIV;
    else
      is_gSF = false;
    end
    
    if ~isfield(xml,'qualityratings') && ~isfield(xml,'QAM')
      fprintf('Quality rating is not saved for %s. Report file %s is incomplete.\nPlease repeat preprocessing amd check for potential errors in the ''err'' folder.\n',job.data{i},xml_files(i,:));    
      return
    end
    if mesh_detected
      if isfield(xml.qualityratings,'NCR')
      % check for newer available surface measures
        if isfield(xml.subjectmeasures,'EC_abs') && isfinite(xml.subjectmeasures.EC_abs) && isfinite(xml.subjectmeasures.defect_size)
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
    cat_progress_bar('Set',i);  
  end
  cat_progress_bar('Clear');
  
  % remove last two columns if EC_abs and defect_size are not defined
  if mesh_detected && all(isnan(QM(:,4))) && all(isnan(QM(:,5)))
    QM = QM(:,1:3);
  end
  
end

[pth,nam] = spm_fileparts(job.data{1});

if ~mesh_detected  
  % voxelsize and origin
  vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
  Orig = V(1).mat\[0 0 0 1]';

  if length(V)>1 && any(any(diff(cat(1,V.dat.dim),1,1),1))
    error('images don''t all have same dimensions')
  end
  if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
    error('images don''t all have same orientation & voxel size')
  end
end

if is_gSF
  fprintf('Use global scaling with TIV\n');
end

Ymean = 0.0;
Yss   = 0.0; % sum of squares

fprintf('Load data ');
for i = 1:n_subjects
  fprintf('.');
  if mesh_detected
    tmp = spm_data_read(V(i));
  else
    tmp(:,:,:) = V(i).dat(:,:,:);
  end
  tmp(isnan(tmp)) = 0;
  if is_gSF
    tmp = tmp*gSF(i)/mean(gSF);
  end
  if i>1 && numel(Ymean) ~= numel(tmp)
    fprintf('\n\nERROR: File %s has different data size: %d vs. %d\n\n',job.data{i},numel(Ymean),numel(tmp));
    return
  end
  Ymean = Ymean + tmp(:);
  Yss   = Yss + tmp(:).^2;
end

% get mean and SD
Ymean = Ymean/n_subjects;

% we have sometimes issues with number precision
Yvar   = 1.0/(n_subjects-1)*(Yss - n_subjects*Ymean.*Ymean);
Yvar(Yvar<0) = 0;
Ystd   = sqrt(Yvar);

% only consider non-zero areas
ind = Ystd ~= 0;

% prepare glassbrain
Ytmp = zeros(size(tmp));
d1 = squeeze(sum(Ytmp,1));
d2 = squeeze(sum(Ytmp,2));
d3 = squeeze(sum(Ytmp,3));

mean_zscore = zeros(n_subjects,1);
for i = 1:n_subjects
  fprintf('.');
  if mesh_detected
    tmp = spm_data_read(V(i));
  else
    tmp(:,:,:) = V(i).dat(:,:,:);
  end
  tmp(isnan(tmp)) = 0;
  if is_gSF
    tmp = tmp*gSF(i)/mean(gSF);
  end
  % calculate Z-score  
  zscore = ((tmp(ind) - Ymean(ind)).^2)./Ystd(ind);

  % calculate glassbrain with emphasized Z-score
  Ytmp(ind) = zscore.^5;
  Ytmp = reshape(Ytmp,size(tmp));
  d1 = d1 + squeeze(sum(Ytmp,1));
  d2 = d2 + squeeze(sum(Ytmp,2));
  d3 = d3 + squeeze(sum(Ytmp,3));
  
  % and use mean of Z-score as overall measure
  mean_zscore(i) = mean(zscore);
end
fprintf('\n');

% not yet finished
if 0
  mx = max([d1(:); d2(:); d3(:)]);
  
  figure(11)
  colormap(hot(64))
  subplot(2,2,1)
  imagesc(rot90(d1),[0 mx])
  axis off image
  
  subplot(2,2,2)
  imagesc(rot90(d2),[0 mx])
  axis off image
  
  subplot(2,2,3)
  imagesc(d3,[0 mx])
  axis off image

end

if isxml
  % estimate product between weighted overall quality (IQR) and mean Z-score
  IQR = QM(:,3);
  IQRratio = (mean_zscore/std(mean_zscore)).*(IQR/std(IQR));
  if mesh_detected
    Euler_number = QM(:,4);
    Topo_defects = QM(:,5);
  end
end

figure
cat_plot_boxplot(mean_zscore,struct('style',2));

fid = fopen(job.csv_name,'w');

if fid < 0
  error('No write access for %s: check file permissions or disk space.',job.csv_name);
end

fprintf(fid,'Path;Name;Mean Z-score');
if isxml
  fprintf(fid,';Weighted overall image quality (IQR);Normalized product of IQR and Mean Z-score');
  if mesh_detected
    fprintf(fid,';Euler Number;Size of topology defects\n');
  else
    fprintf(fid,'\n');
  end
else
  fprintf(fid,'\n');
end
for i = 1:n_subjects
  [pth, data_name] = fileparts(job.data{i});
  fprintf(fid,'%s;%s;%g',pth,data_name,mean_zscore(i));
  if isxml
    fprintf(fid,';%g;%g',IQR(i),IQRratio(i));
    if mesh_detected
      fprintf(fid,';%d;%g\n',Euler_number(i),Topo_defects(i));
    else
      fprintf(fid,'\n');
    end
  else
    fprintf(fid,'\n');
  end
end

if fclose(fid)==0
  fprintf('\nValues saved in %s.\n',job.csv_name);
end
