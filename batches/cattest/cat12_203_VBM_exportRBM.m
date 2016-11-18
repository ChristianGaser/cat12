% ---------------------------------------------------------------------
% Test batch for surface data resampling and smoothing of cat_tst_cattest.
% This batch exports all ROI measures contained in the XML files.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id$

%#ok<*SAGROW>

% prepare filenames
% ---------------------------------------------------------------------
if exist('files','var') 
  xmls  = {}; 
  for fi = 1:numel(files)
    [pp,ff] = spm_fileparts(files{fi}); 
    xmls{end+1,1}  = fullfile( pp , roidir , ['catROI_'  ff '.xml'] ); 
  end
  outdir = {fullfile( pp , 'RBMexport' )};
else
  xmls   = {''}; 
  outdir = {'<UNDEFINED>'};  
end  




% batch
% ---------------------------------------------------------------------

% VBM atlases
matlabbatch{1}.spm.tools.cat.tools.calcroi.roi_xml      = xmls;
matlabbatch{1}.spm.tools.cat.tools.calcroi.folder       = 0;               
matlabbatch{1}.spm.tools.cat.tools.calcroi.point        = '.';             
matlabbatch{1}.spm.tools.cat.tools.calcroi.outdir       = outdir;      
matlabbatch{1}.spm.tools.cat.tools.calcroi.calcroi_name = 'VBMatlases';
