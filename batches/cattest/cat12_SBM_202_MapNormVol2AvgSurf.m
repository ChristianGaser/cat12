% ---------------------------------------------------------------------
% Test batch to test the Mapping of Normalized Volumes to the average 
% surface. Used in cat_tst_cattest.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_run_job.m 1013 2016-09-22 11:49:13Z dahnke $

%#ok<*SAGROW>
 
% prepare filename
%-----------------------------------------------------------------------
if exist('files','var') 
  for fi = 1:numel(files)
    [pp,ff]     = spm_fileparts(files{fi}); 
    vols{fi,1}  = fullfile( pp , mridir  , ['wm' ff '.nii'] );
  end
else
  vols = {''};
  exp = cat_get_defaults('extopts.expertgui'); 
end  
surfs = { fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152.gii') };

% batch
% ----------------------------------------------------------------------

% cat default value - 3 average points in the GM (CSF/GM - GM - GM/WM)
matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_vol                           = vols; 
matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
if exp, matlabbatch{1}.spm.tools.cat.stools.vol2surf.interp                     = {'linear'}; end
matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'GM-T1-intensity';
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 0;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0.5;
matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 1;

% L1 intensity of equal distance layer
matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{2}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
if exp, matlabbatch{2}.spm.tools.cat.stools.vol2surf.interp                     = {'linear'}; end
matlabbatch{2}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'EquiDist-Layer0-Intensity';
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 1/13;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{2}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
% L4 intensity of equal distance layer
matlabbatch{3}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{3}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
if exp, matlabbatch{3}.spm.tools.cat.stools.vol2surf.interp                     = {'linear'}; end
matlabbatch{3}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'EquiDist-Layer4-Intensity';
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 7/13;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{3}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
% L6 intensity of equal distance layer
matlabbatch{4}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{4}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
if exp, matlabbatch{4}.spm.tools.cat.stools.vol2surf.interp                     = {'linear'}; end
matlabbatch{4}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'EquiDist-Layer6-Intensity';
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 12/13;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{4}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;

% WM value at 1.5 times of the local thickness (half thickness into the WM)
matlabbatch{5}.spm.tools.cat.stools.vol2surf.data_vol                           = vols;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.data_mesh_lh                       = surfs; 
matlabbatch{5}.spm.tools.cat.stools.vol2surf.sample                             = {'avg'};
if exp, matlabbatch{5}.spm.tools.cat.stools.vol2surf.interp                     = {'linear'}; end
matlabbatch{5}.spm.tools.cat.stools.vol2surf.datafieldname                      = 'WM-Intensity-at150percentCT';
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class          = 'GM';
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint     = 1.5;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.stepsize       = 0;
matlabbatch{5}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint       = 0;
