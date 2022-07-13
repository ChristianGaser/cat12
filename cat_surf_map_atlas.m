function cat_surf_map_atlas(central_files, atlas_files)
% Map Freesurfer atlas annotation labels to individual surface data for left
% and right hemisphere
%
% central - array of file names for left central surface
% atlas   - array of file names for Freesurfer annot atlas files
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
atlasDir  = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces');

if nargin < 1
  central_files = spm_select([1 Inf],'^lh.central.(?!nofix).*','Select left central surfaces. Right side will be automatically processed.');
end
nc = size(central_files,1);

if nargin < 2
  atlas_files = spm_select([1 Inf],'^lh.*\.annot$','Select left Freesurfer atlas labels. Right side will be automatically processed.',{},atlasDir);
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
            
      % get new annot name
      annot_tmp = cat_surf_rename(central,'dataname',tmp.dataname,'ee','.annot');

      pth0 = spm_fileparts(annot_tmp{1});

      % temporarily save annot file as txt file 
      annot_txt = cat_surf_rename(annot_tmp{1},'ee','.txt');
      [vertices, label, colortable] = cat_io_FreeSurfer('read_annotation',fs_annot{1});
      fs_annot_txt = fullfile(pth0,[hemi(k,:) '.' tmp.dataname '.txt']);
      fp = fopen(fs_annot_txt,'w');
      if fp
        fprintf(fp,'%d\n',label);
        fclose(fp);
      else
        error('Cannot write %s. Please check file permissions.\n',fs_annot_txt);
      end
        
      cmd = sprintf('CAT_ResampleSurf -label "%s" "%s" "%s" "NULL" "%s" "%s"',...
        fs_central,fs_sphere,spherereg{1},fs_annot_txt,annot_txt{1});
      ST = cat_system(cmd);
      delete(fs_annot_txt)
      if ~ST
        fprintf('Save %s\n',annot_tmp{1});
        label = load(annot_txt{1});
        vertices = ((1:numel(label))-1)';
        cat_io_FreeSurfer('write_annotation',annot_tmp{1}, vertices, label, colortable);
        delete(annot_txt{1})
      end
    end
  
  end
end
