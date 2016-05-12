%-----------------------------------------------------------------------
% Job for fmri batch
% Christian Gaser
% $Id: cat_fmri_main.m 932 2016-04-22 11:11:51Z gaser $
%-----------------------------------------------------------------------
global opts extopts output

warning('off','MATLAB:DELETE:FileNotFound');

% use some options from gui or default file
if exist('opts','var')
  matlabbatch{2}.spm.tools.cat.estwrite.opts = opts;
end
if exist('extopts','var')
  matlabbatch{2}.spm.tools.cat.estwrite.extopts = extopts;
end
if exist('output','var')
  matlabbatch{2}.spm.tools.cat.estwrite.output = output;
end
matlabbatch{2}.spm.tools.cat.estwrite.nproc = 0;

matlabbatch{1}.spm.spatial.coreg.estimate.ref = '<UNDEFINED>';
matlabbatch{1}.spm.spatial.coreg.estimate.source = '<UNDEFINED>';
matlabbatch{1}.spm.spatial.coreg.estimate.other = '<UNDEFINED>';
matlabbatch{2}.spm.tools.cat.estwrite.data = '<UNDEFINED>';
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.native = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{2}.spm.tools.cat.estwrite.output.jacobian.warped = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.warps = [1 0];
matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';
