%-----------------------------------------------------------------------
% Job for longitudinal batch
% Christian Gaser
% $Id$
%-----------------------------------------------------------------------

global opts extopts output modulate

warning('off','MATLAB:DELETE:FileNotFound');
matlabbatch{1}.spm.tools.vbm.tools.series.data = '<UNDEFINED>';

% use some options from gui or default file
for j=2:3
  if exist('opts','var')
    matlabbatch{j}.spm.tools.vbm.estwrite.opts = opts;
  end
  if exist('extopts','var')
    matlabbatch{j}.spm.tools.vbm.estwrite.extopts = extopts;
  end
  if exist('output','var')
    matlabbatch{j}.spm.tools.vbm.estwrite.output = output;
  end
end

% modulation option for applying deformations
if exist('modulate','var')
  matlabbatch{4}.spm.tools.vbm.tools.defs.modulate = modulate;
end

matlabbatch{1}.spm.tools.vbm.tools.series.bparam = 1000000;
matlabbatch{2}.spm.tools.vbm.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Midpoint Average', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{2}.spm.tools.vbm.estwrite.output.surface = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.ROI = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.GM.native = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.GM.warped = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.GM.modulated = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.GM.dartel = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.WM.native = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.WM.warped = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.WM.modulated = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.WM.dartel = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.CSF.native = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.CSF.warped = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.CSF.modulated = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.CSF.dartel = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.label.native = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.label.warped = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.label.dartel = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.bias.native = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.bias.warped = 1;
matlabbatch{2}.spm.tools.vbm.estwrite.output.bias.dartel = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.jacobian.warped = 0;
matlabbatch{2}.spm.tools.vbm.estwrite.output.warps = [1 0];
matlabbatch{3}.spm.tools.vbm.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{3}.spm.tools.vbm.estwrite.output.GM.native = 1;
matlabbatch{3}.spm.tools.vbm.estwrite.output.GM.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.GM.modulated = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.GM.dartel = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.WM.native = 1;
matlabbatch{3}.spm.tools.vbm.estwrite.output.WM.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.WM.modulated = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.WM.dartel = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.CSF.native = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.CSF.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.CSF.modulated = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.CSF.dartel = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.label.native = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.label.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.label.dartel = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.bias.native = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.bias.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.bias.dartel = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.jacobian.warped = 0;
matlabbatch{3}.spm.tools.vbm.estwrite.output.warps = [0 0];
matlabbatch{4}.spm.tools.vbm.tools.defs.field1(1) = cfg_dep('VBM12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{4}.spm.tools.vbm.tools.defs.images(1) = cfg_dep('VBM12: Segmentation: p1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{4}.spm.tools.vbm.tools.defs.images(2) = cfg_dep('VBM12: Segmentation: p2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{4}.spm.tools.vbm.tools.defs.interp = 4;
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('VBM12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('VBM12: Segmentation: p1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('VBM12: Segmentation: p2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;