function Atlas_Anatomy2CAT12
% Convert Anatomy3 atlas to CAT12 atlas
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

load JuBrain_Data_public_v30.mat

V = JuBrain.Vo;
V.fname ='Anatomy3_refined.nii';
V.pinfo(1) = 1;
atlas = zeros(V.dim);
n_structures = numel(JuBrain.Namen);

% flip hemisphere coding
JuBrain.lr = 2 - JuBrain.lr;

fid = fopen('anatomy3.csv','w');
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

matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(cat_get_defaults('extopts.pth_templates'),'aal3.nii')
                                        'Anatomy3_refined.nii'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'anatomy3.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;

spm_jobman('run',matlabbatch)

V = spm_vol('anatomy3.nii');
vol = round(spm_read_vols(V));
V.pinfo(1) = 1;
spm_write_vol(V,vol);