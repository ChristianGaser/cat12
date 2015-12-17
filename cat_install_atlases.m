function cat_install_atlases
% Add Dartel atlas labels to spm atlas folder
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

atlas = cat_get_defaults('extopts.atlas');

spm_dir = spm('dir');
[ST, RS] = mkdir(spm_dir,'atlas');
if ST
  atlas_dir = fullfile(spm_dir,'atlas');
  for i = 1:size(atlas,1)
    atlas_file = atlas{i,1};
    [pth,nam,ext] = spm_fileparts(atlas_file);
    xml_file = fullfile(pth,['label_dartel_' nam '.xml']);
    try
      copyfile(atlas_file,atlas_dir);
      copyfile(xml_file,atlas_dir);
      fprintf('Copy %s\n',atlas_file);
    catch
      disp('Writing error: Please check file permissions.');
    end
  end
else
  error(RS);
end

fprintf('Use atlas function in SPM Results or context menu in orthogonal view (via right mouse button): Display|Labels\n');