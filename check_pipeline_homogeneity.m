function check_pipeline_homogeneity
%check_pipeline_homogeneity to sample homogeneity for different
% cat12 releases using check_pipeline.sh
%
% _________________________________________________________________________
% $Id$

long = spm_input('Longitudinal data?','+1','y/n',[1,0],2);

[files,dirs] = spm_select('FPList','.','^check_r');

job = struct('c',[],'data_xml',[],'gap',3, 'verb',true);

j = 0;
for i = 1:size(dirs,1)
  folder = deblank(dirs(i,:));
  if ~isempty(strfind(folder,'check_r'))
    if long
      files = spm_select('FPListRec',[folder '/long'],'^mw.*p1.*\.nii$');
      mn_files = 2;
    else
      files = spm_select('FPListRec',folder,'^mwp1.*\.nii$');
      mn_files = 5;
    end
    if size(files,1) >= mn_files
      j = j + 1;
      
      % remove phantom data because there are only sparse entries for that data
      for k = 1:size(files)
        if ~isempty(strfind(files(k,:),'dilate'))
          files(k,:) = [];
          break
        end
      end
      
      job.data_vol{j} = files;
    end
  end
end

cat_stat_check_cov(job)
