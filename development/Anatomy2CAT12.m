function Anatomy2CAT12

load JuBrain_Data_public_v30.mat

V = JuBrain.Vo;
V.fname ='Anatomy3_refined.nii';
V.pinfo(1) = 1;
atlas = zeros(V.dim);
n_structures = numel(JuBrain.Namen);

% flip hemisphere coding
JuBrain.lr = 2 - JuBrain.lr;

fid = fopen('Anatomy3.csv','w');
fprintf(fid,'ROIid;ROIabbr;ROIname;ROIoname\n');

% code left/right hemispheres
tmp = 2*JuBrain.mpm - JuBrain.lr;

% correct atlas regions that were zero before
tmp(JuBrain.mpm == 0) = 0;
atlas(JuBrain.idx) = tmp;

% replace remaining holes with median value
holes = atlas > 0;
holes = (cat_vol_morph(holes,'c') - holes) > 0;
tmp = cat_vol_median3c(single(atlas));
atlas(holes) = double(tmp(holes));

spm_write_vol(V,atlas);

% write csv file for CAT12
for i=1:n_structures  
  fprintf(fid,'%d;l%s;Left %s;%s\n',2*i-1,deblank(JuBrain.Files{i}),deblank(JuBrain.Namen{i}),deblank(JuBrain.Namen{i}));
  fprintf(fid,'%d;r%s;Right %s;%s\n',2*i-1,deblank(JuBrain.Files{i}),deblank(JuBrain.Namen{i}),deblank(JuBrain.Namen{i}));
end
fclose(fid);

return

% not needed with new template anymore...
matlabbatch{1}.spm.tools.cat.tools.defs.field1 = {'/Users/gaser/matlab/cat12/development/y_mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
matlabbatch{1}.spm.tools.cat.tools.defs.images = {'Anatomy3.nii,1'};
matlabbatch{1}.spm.tools.cat.tools.defs.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{1}.spm.tools.cat.tools.defs.vox = [1 1 1];
matlabbatch{1}.spm.tools.cat.tools.defs.interp = -1;
matlabbatch{1}.spm.tools.cat.tools.defs.modulate = 0;

spm_jobman('run',matlabbatch)