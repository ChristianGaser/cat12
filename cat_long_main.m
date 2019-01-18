%-----------------------------------------------------------------------
% Job for longitudinal batch
% Christian Gaser
% $Id$
%-----------------------------------------------------------------------

global opts extopts output modulate dartel warps delete_temp ROImenu surfaces

warning('off','MATLAB:DELETE:FileNotFound');
matlabbatch{1}.spm.tools.cat.tools.series.data = '<UNDEFINED>';

% use some options from gui or default file
for j=2:3
  if exist('opts','var')
    matlabbatch{j}.spm.tools.cat.estwrite.opts = opts;
  end
  if exist('extopts','var')
    matlabbatch{j}.spm.tools.cat.estwrite.extopts = extopts;
  end
  if exist('output','var')
    matlabbatch{j}.spm.tools.cat.estwrite.output = output;
  end
  matlabbatch{j}.spm.tools.cat.estwrite.nproc = 0;
end

% ROI options
if exist('ROImenu','var')
  matlabbatch{3}.spm.tools.cat.estwrite.output.ROImenu = ROImenu;
end

% modulation option for applying deformations
if modulate
  matlabbatch{4}.spm.tools.cat.tools.defs.modulate = modulate;
end

% surface estimation
if surfaces
  matlabbatch{3}.spm.tools.cat.estwrite.output.surface = 1;
else
  matlabbatch{3}.spm.tools.cat.estwrite.output.surface = 0;
end

% dartel export option
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.dartel = dartel;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.dartel = dartel;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.dartel = dartel;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.dartel = dartel;

% longitudinal rigid registration with final masking
matlabbatch{1}.spm.tools.cat.tools.series.bparam = 1000000;
matlabbatch{1}.spm.tools.cat.tools.series.use_brainmask = 1;

% cat12 segmentation of midpoint average and saving deformation field
matlabbatch{2}.spm.tools.cat.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Midpoint Average', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{2}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.ROImenu.noROI = struct([]);
matlabbatch{2}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.biaswarped = 1;
matlabbatch{2}.spm.tools.cat.estwrite.output.warps = [1 0];

% cat12 segmentation of realigned images
matlabbatch{3}.spm.tools.cat.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{3}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.biaswarped = 0;

% applying deformations to native segmentations
matlabbatch{4}.spm.tools.cat.tools.defs.field1(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.images(1) = cfg_dep('CAT12: Segmentation: p1 Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.images(2) = cfg_dep('CAT12: Segmentation: p2 Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.interp = 1;

if delete_temp
  matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('CAT12: Segmentation: p1 Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
  matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('CAT12: Segmentation: p2 Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
end

% save deformations
if ~warps & delete_temp
  matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
end

if delete_temp
  matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
end


