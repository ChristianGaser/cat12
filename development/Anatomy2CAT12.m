function convert2CAT12

load JuBrain_Data_public_v30.mat

V = JuBrain.Vo;
V.fname ='Anatomy3.nii';
V.pinfo(1) = 1;
atlas = zeros(V.dim);

JuBrain.lr = 2 - JuBrain.lr;

fid = fopen('Anatomy3.csv','w');
fprintf(fid,'ROIid;ROIabbr;ROIname;ROIoname\n');
for i=1:93
  % code left/right hemispheres
  tmp = 2*JuBrain.mpm - JuBrain.lr;
  
  % correct atlas regions zhat were zero before
  tmp(JuBrain.mpm == 0) = 0;
  atlas(JuBrain.idx) = tmp;
  
  % write csv file for CAT12
  fprintf(fid,'%d;l%s;Left %s;%s\n',2*i-1,deblank(JuBrain.Files{i}),deblank(JuBrain.Namen{i}),deblank(JuBrain.Namen{i}));
  fprintf(fid,'%d;r%s;Right %s;%s\n',2*i-1,deblank(JuBrain.Files{i}),deblank(JuBrain.Namen{i}),deblank(JuBrain.Namen{i}));
end
fclose(fid);

spm_write_vol(V,atlas);

matlabbatch{1}.spm.tools.cat.tools.defs.field1 = {'/Users/gaser/matlab/cat12/development/y_mni_icbm152_t1_tal_nlin_asym_09c.nii,1'};
matlabbatch{1}.spm.tools.cat.tools.defs.images = {'Anatomy3.nii,1'};
matlabbatch{1}.spm.tools.cat.tools.defs.bb = [NaN NaN NaN
                                              NaN NaN NaN];
matlabbatch{1}.spm.tools.cat.tools.defs.vox = [1 1 1];
matlabbatch{1}.spm.tools.cat.tools.defs.interp = -1;
matlabbatch{1}.spm.tools.cat.tools.defs.modulate = 0;

spm_jobman('run',matlabbatch)