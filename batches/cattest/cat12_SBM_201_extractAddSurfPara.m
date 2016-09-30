% ---------------------------------------------------------------------
% Test batch for surface measures of cat_tst_cattest.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_run_job.m 1013 2016-09-22 11:49:13Z dahnke $

%#ok<*SAGROW>

% prepare filename
% ----------------------------------------------------------------------
if exist('files','var') 
  for fi = 1:numel(files)
    [pp,ff]  = spm_fileparts(files{fi}); 
    surfs{fi*2-1,1} = fullfile( pp , surfdir , ['lh.central.' ff '.gii'] );  
    surfs{fi*2  ,1} = fullfile( pp , surfdir , ['rh.central.' ff '.gii'] );  
  end
else
  surfs = {''}; 
  exp = cat_get_defaults('extopts.expertgui'); 
end  

% batch
% ----------------------------------------------------------------------
matlabbatch{1}.spm.tools.cat.stools.surfextract.data_surf = surfs; % required the lh and rh surface!
matlabbatch{1}.spm.tools.cat.stools.surfextract.GI        = 1; % abs mean curvature
matlabbatch{1}.spm.tools.cat.stools.surfextract.FD        = 1; % fractal dimention
matlabbatch{1}.spm.tools.cat.stools.surfextract.SD        = 1; % sulcal depth
matlabbatch{1}.spm.tools.cat.stools.surfextract.nproc     = 0; % multi processes (not for this script)
% expert options that are stil in development
if exp
  matlabbatch{1}.spm.tools.cat.stools.surfextract.area      = 1; % EXPERT
  matlabbatch{1}.spm.tools.cat.stools.surfextract.GIA       = 0; % EXPERT 
  matlabbatch{1}.spm.tools.cat.stools.surfextract.GII       = 0; % EXPERT 
  matlabbatch{1}.spm.tools.cat.stools.surfextract.GIL       = 1; % EXPERT laplacican GI
  matlabbatch{1}.spm.tools.cat.stools.surfextract.GIS       = 0; % EXPERT sperical GI 
  matlabbatch{1}.spm.tools.cat.stools.surfextract.lazy      = 0; % EXPERT
end