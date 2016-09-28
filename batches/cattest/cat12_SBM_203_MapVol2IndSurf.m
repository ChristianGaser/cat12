% ---------------------------------------------------------------------
% Test batch for volume to surface mapping of cat_tst_cattest.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_io_volctype.m 945 2016-06-06 06:54:56Z gaser $

% prepare filename
% ----------------------------------------------------------------------
vols = files; surfs = files; 
if exist('files','var') 
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    vols{fi,1}  = fullfile( pp , mridir  , ['m' ff '.nii'] );
    surfs{fi,1} = fullfile( pp , surfdir , ['lh.central.' ff '.gii'] ); 
  end
end  

% batch
% ----------------------------------------------------------------------

% cat default value - 3 average points in the GM (CSF/GM - GM - GM/WM)
matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_vol                           = vols; 
matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
matlabbatch{1}.spm.tools.cat.stools.vol2surf.interp                             = {'linear'};
matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'T1GMavg3';
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 0;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0.5;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 1;

% L1 intensity of equal distance layer
matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{2}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
matlabbatch{2}.spm.tools.cat.stools.vol2surf.interp                             = {'linear'};
matlabbatch{2}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'T1EL0';
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 1/13;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
% L4 intensity of equal distance layer
matlabbatch{3}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{3}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
matlabbatch{3}.spm.tools.cat.stools.vol2surf.interp                             = {'linear'};
matlabbatch{3}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'T1EL4';
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 7/13;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
% L6 intensity of equal distance layer
matlabbatch{4}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{4}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
matlabbatch{4}.spm.tools.cat.stools.vol2surf.interp                             = {'linear'};
matlabbatch{4}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'T1L6';
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 12/13;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;

% WM value at 1.5 times of the local thickness (half thickness into the WM)
matlabbatch{5}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{5}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
matlabbatch{5}.spm.tools.cat.stools.vol2surf.interp                             = {'linear'};
matlabbatch{5}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'T1WM';
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 1.5;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
