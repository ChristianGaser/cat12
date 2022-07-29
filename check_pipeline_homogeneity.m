function check_pipeline_homogeneity
%check_pipeline_homogeneity to sample homogeneity for different
% cat12 releases using check_pipeline.sh
%
% _________________________________________________________________________
% $Id$

min_release = 1800;
max_release = 3000;

% list nifti files and exclude longitudinal ADNI-data
data = spm_select('List','.','.nii');
ind = [];
for i=1:size(data,1)
  if ~isempty(strfind(data(i,:),'ADNI'))
    ind = [ind i];
  end
end
data(ind,:) = [];
data = char('ADNI-longitudinal',data);
sel = spm_input('Data',1,'m',data);

job = struct('c',[],'data_xml',[],'gap',3, 'verb',true,'show_name',1,'show_violin',0);

% longitudinal data
if sel == 1
  dirs = spm_select('List','.','dir','^check_r');
  
  % exclude older releases
  ind = [];
  for i = 1:size(dirs,1)
    ind_r = strfind(dirs(i,:),'check_r');
    release = str2num(dirs(i,ind_r+7:ind_r+10));
    if release < min_release || release > max_release
      ind = [ind i];
    end
  end
  dirs(ind,:) = [];
  
  j = 0;
  for i = 1:size(dirs,1)
    folder = deblank(dirs(i,:));
    if ~isempty(strfind(folder,'check_r'))
      files = spm_select('FPListRec',[folder '/long'],'^mwmwp1.*\.nii$');
      if size(files,1) >= 2
        j = j + 1;        
        job.data_vol{j} = files;
      end
    end
  end
else % selected cross-sectional data
  files = spm_select('FPListRec','.',['^mwp1' deblank(data(sel,:))]);
  ind = [];
  for i=1:size(files,1)
    ind_r = strfind(files(i,:),'check_r');
    release = str2num(files(i,ind_r+7:ind_r+10));
    if (release < min_release) || release > max_release || ~isempty(strfind(files(i,:),'not_used'))
      ind = [ind i];
    end
  end
  files(ind,:) = [];
  job.data_vol{1} = files;
end

cat_stat_homogeneity(job);
[pth filename] = spm_fileparts(deblank(data(sel,:)));
 
name = ['check_cov' filename '.png'];

saveas(1, name);
