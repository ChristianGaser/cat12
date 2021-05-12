function ROI_summarize_all(sel)
% tool to select the best threshold for the MPM atlases by comparing volumes of gray matter inside the labels

if nargin == 0
  sel = spm_input('Which Atlas?', 1, 'm',{'LPBA40','Cobra','IBSR2','Neuromorphometrics','Hammers'});
end

switch sel

case 1
  atlas_file = cellstr(spm_select('FPList','/Volumes/UltraMax/LPBA40/delineation_space','^MPM_th.*S.delineation'));
  def_file   = cellstr(spm_select('FPList','/Volumes/UltraMax/LPBA40/delineation_space/mri','^y_'));
  value_file = cellstr(spm_select('FPList','/Volumes/UltraMax/LPBA40/delineation_space/mri','^p1'));
  label_file = cellstr(spm_select('FPList','/Volumes/UltraMax/LPBA40/delineation_space','^S.*delineation.structure.label'));
case 2
  atlas_file = cellstr(spm_select('FPList','/Volumes/UltraMax/WinterburnHippocampalAtlas_NIfTI/labels','^MPM_th*'));
  def_file   = cellstr(spm_select('FPList','/Volumes/UltraMax/WinterburnHippocampalAtlas_NIfTI/resampled_registered/mri','^y_'));
  value_file = cellstr(spm_select('FPList','/Volumes/UltraMax/WinterburnHippocampalAtlas_NIfTI/resampled_registered/mri','^p1lowresR2x2x2'));
  label_file = cellstr(spm_select('FPList','/Volumes/UltraMax/WinterburnHippocampalAtlas_NIfTI/resampled_registered','^lowresR2x2x2_subject.*labels.nii'));
case 3
  atlas_file = cellstr(spm_select('FPList','/Volumes/UltraMax/IBSR2','^MPM_th.*IBSR'));
  def_file   = cellstr(spm_select('FPList','/Volumes/UltraMax/IBSR2/mri','^y_'));
  value_file = cellstr(spm_select('FPList','/Volumes/UltraMax/IBSR2/mri','^p1'));
  label_file = cellstr(spm_select('FPList','/Volumes/UltraMax/IBSR2','^IBSR.*seg_ana'));
case 4
  atlas_file = cellstr(spm_select('FPList','/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/','^MPM_th.*glm'));
  def_file   = cellstr(spm_select('FPList',{'/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/testing-images/mri','/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/training-images/mri'},'^y_'));
  value_file = cellstr(spm_select('FPList',{'/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/testing-images/mri','/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/training-images/mri'},'^p1'));
  label_file = cellstr(spm_select('FPList',{'/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/testing-labels','/Volumes/UltraMax/Neuromorphometrics_MICCAI2012/training-labels'},'^1.*glm'));
case 5
  atlas_file = cellstr(spm_select('FPList','/Volumes/UltraMax/HammersAtlas/Hammers_mith-n30r95','^MPM_th.*'));
  def_file   = cellstr(spm_select('FPList','/Volumes/UltraMax/HammersAtlas/Hammers_mith-n30r95/mri','^y_'));
  value_file = cellstr(spm_select('FPList','/Volumes/UltraMax/HammersAtlas/Hammers_mith-n30r95/mri','^p1'));
  label_file = cellstr(spm_select('FPList','/Volumes/UltraMax/HammersAtlas/Hammers_mith-n30r95','^a.*-seg'));
end

if numel(def_file) ~= numel(value_file) | numel(def_file) ~= numel(label_file)
  error('Number of files for deformations (n=%d), values (n=%d) or label (n=%d) differ\n',numel(def_file),numel(value_file),numel(label_file));
end

n_subjects = numel(def_file);
n_atlas    = numel(atlas_file);

fprintf('Number of atlases:\t%d\nNumber of subjects:\t%d\n',n_atlas,n_subjects);

% use relative error or absolute error
calc_relative_error = 1;

err_arr = NaN(n_subjects,n_atlas);
mean_arr = NaN(n_subjects,n_atlas);

V = spm_vol(atlas_file{1});
atlas = round(spm_read_vols(V));
structures = sort(unique(atlas(atlas > 0)));
n_structures = numel(structures);

for si=1:n_subjects

  fprintf('%s:\n',value_file{si});
  V = spm_vol(value_file{si});
  value = spm_read_vols(V);
  label = round(spm_read_vols(spm_vol(label_file{si})));
  
  % IBSR2 atlas should be permuted because of wrong orientation
  if sel == 3
    label = permute(label,[1 3 2]);
  end
  vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  
  ref_val = zeros(n_structures,1);
  for ri = 1:n_structures
    ref_val(ri) = prod(vx_vol)*sum(value(label==structures(ri)))/1000;
  end
  
  Y  = cat_vol_ROI_summarize(struct('atlases',{atlas_file},'field1',{{def_file{si}}},'images',{{value_file{si}}},'fhandle','volume'));
  
  ind = (ref_val ~= 0);

  for ai=1:n_atlas
    [pp,atlas_name] = spm_fileparts(atlas_file{ai});
    if calc_relative_error
        perc_err = 100*mean((Y{ai,1}(ind)-ref_val(ind))./ref_val(ind));
    else
        perc_err = mean((Y{ai,1}(ind)-ref_val(ind)));
    end
    mean_Y = mean(Y{ai,1}(ind)./ref_val(ind));
    fprintf('%3.6f\t%3.6f\t%s\n',perc_err,mean_Y,atlas_name);
    err_arr(si,ai) = perc_err;
    mean_arr(si,ai) = mean_Y;
  end
end

[tmp, name] = spm_str_manip(atlas_file,'C');
fprintf('Atlas names:%s\n',tmp)
fprintf('%20s\t%30s\t%30s\n','Atlas','median error (should be low)','median ratio (should be close to one)');
med_err = median(err_arr);
med_mean = median(mean_arr);
for ai=1:n_atlas
  fprintf('%20s\t%30.5f\t%30.5f\n',name.m{ai},med_err(ai), med_mean(ai));
end

figure(11)
cat_plot_boxplot(err_arr(:,:,1),struct('violin',0,'showdata',1,'names',{name.m}));
title('median error (should be low)')

figure(12)
cat_plot_boxplot(mean_arr(:,:,1),struct('violin',0,'showdata',1,'names',{name.m}));
title('median ratio (should be close to one)')

