function cat_surf_map_atlas(central_files, atlas_files)
% Map Freesurfer atlas annotation labels to individual surface data for left
% and right hemisphere
%
% central - array of file names for left central surface
% atlas   - array of file names for Freesurfer annot atlas files
%_______________________________________________________________________
% $Id$

fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
atlasDir  = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces');

if nargin < 1
  central_files = spm_select([1 Inf],'^lh.central.(?!nofix).*','Select left central surfaces');
end
nc = size(central_files,1);

if nargin < 2
  atlas_files = spm_select([1 Inf],'^lh.*\.annot$','Select left Freesurfer atlas labels',{},atlasDir);
end
na = size(atlas_files,1);

hemi = char('lh','rh');
for k = 1:size(hemi,1)
  for i = 1:nc
    central    = cat_surf_rename(deblank(central_files(i,:)),'side',hemi(k,:));
    spherereg  = cat_surf_rename(central,'dataname','sphere.reg');
    fs_central = fullfile(fsavgDir,[hemi(k,:) '.central.freesurfer.gii']);
    fs_sphere  = fullfile(fsavgDir,[hemi(k,:) '.sphere.freesurfer.gii']);
    
    for j = 1:na
      fs_annot  = cat_surf_rename(deblank(atlas_files(j,:)),'side',hemi(k,:));
      tmp       = cat_surf_info(fs_annot);
      
      % get new annot name and save that file in label instead of surf folder
      annot_tmp = cat_surf_rename(central,'dataname',tmp.dataname,'ee','.annot');
      [pth0,nam0,ext] = spm_fileparts(annot_tmp{1});
      [pth1,nam1]     = spm_fileparts(pth0);
      if strcmp(nam1,'surf')
        annot = fullfile(pth1,'label',[nam0 ext]);
      else
        annot = annot_tmp{1};
      end

      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "NULL" "%s" "%s"',...
        fs_central,fs_sphere,spherereg{1},fs_annot{1},annot);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS);
      if ~ST
        fprintf('Save %s\n',annot);
      end
    end
  
  end
end
