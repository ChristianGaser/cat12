% ---------------------------------------------------------------------
% Test batch for surface data resampling and smoothing of cat_tst_cattest.
% This batch exports all ROI measures contained in the XML files.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat_run_job.m 1013 2016-09-22 11:49:13Z dahnke $

%#ok<*SAGROW>

% prepare filenames
% ---------------------------------------------------------------------
if exist('files','var') 
  xmls  = {}; sxmls = {};
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    xmls{end+1,1}  = fullfile( pp , roidir , ['catROI_'  ff '.xml'] ); 
    sxmls{end+1,1} = fullfile( pp , roidir , ['catROIs_' ff '.xml'] ); 
  end
end  

outdir = {fullfile( job.resdir , 'RBMexport' )};

% batch
% ---------------------------------------------------------------------

% VBM atlases
matlabbatch{1}.spm.tools.cat.tools.calcroi.roi_xml      = xmls;
matlabbatch{1}.spm.tools.cat.tools.calcroi.folder       = 0;               
matlabbatch{1}.spm.tools.cat.tools.calcroi.point        = '.';             
matlabbatch{1}.spm.tools.cat.tools.calcroi.outdir       = outdir;      
matlabbatch{1}.spm.tools.cat.tools.calcroi.calcroi_name = 'VBMatlases';

% SBM atlases with full path in subject name
matlabbatch{2}.spm.tools.cat.tools.calcroi.roi_xml      = sxmls;
matlabbatch{2}.spm.tools.cat.tools.calcroi.folder       = 1;               
matlabbatch{2}.spm.tools.cat.tools.calcroi.point        = ',';             
matlabbatch{2}.spm.tools.cat.tools.calcroi.outdir       = outdir;      
matlabbatch{2}.spm.tools.cat.tools.calcroi.calcroi_name = 'SBMatlases';

