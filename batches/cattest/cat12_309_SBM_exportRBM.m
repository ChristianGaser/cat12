% ---------------------------------------------------------------------
% Test batch for surface data resampling and smoothing of cat_tst_cattest.
% This batch exports all ROI measures contained in the XML files.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat12_SBM_302_exportRBM.m 1072 2016-11-09 15:01:16Z dahnke $

%#ok<*SAGROW>

% prepare filenames
% ---------------------------------------------------------------------
if exist('files','var') 
  sxmls = {};
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    sxmls{end+1,1} = fullfile( pp , roidir , ['catROIs_' ff '.xml'] ); 
  end
  outdir = {fullfile( pp , 'RBMexport' )};
else
  sxmls  = {''};   
  outdir = {'<UNDEFINED>'};  
end  




% batch
% ---------------------------------------------------------------------
 
% SBM atlases with full path in subject name
matlabbatch{1}.spm.tools.cat.tools.calcroi.roi_xml      = sxmls;
matlabbatch{1}.spm.tools.cat.tools.calcroi.folder       = 1;               
matlabbatch{1}.spm.tools.cat.tools.calcroi.point        = ',';             
matlabbatch{1}.spm.tools.cat.tools.calcroi.outdir       = outdir;      
matlabbatch{1}.spm.tools.cat.tools.calcroi.calcroi_name = 'SBMatlases';

