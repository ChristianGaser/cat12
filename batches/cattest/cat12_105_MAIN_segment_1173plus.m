% ---------------------------------------------------------------------
% Test batch "Segment Data" for VBM, SBM, and RBM preprocessing of 
% cat_tst_cat.est.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id: cat12_101_MAIN_segment.m 1083 2016-11-21 14:11:02Z dahnke $

if ~exist('files_human1173plus','var')
  files_human1173plus = {'<UNDEFINED>'};
  exp = cat_get_defaults1173plus('extopts1173plus.expertgui'); 
elseif isempty(files_human1173plus)
  return;
end


% batch
% -- opts --------------------------------------------------------------
matlabbatch{1}.spm.tools.cat.estwrite1173plus.data                       = files_human1173plus;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.nproc                      = 0;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.tpm                   = { fullfile( spm('dir') , 'tpm' , 'TPM.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.affreg                = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.biasstr               = 0.5;
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.ngaus               = [1 1 2 3 4 2];
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.biasreg             = 0.001;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.biasfwhm            = 60;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.warpreg             = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.opts.samp                = 3;
end
% -- extopts -----------------------------------------------------------
% registration
matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.darteltpm          = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_1_IXI555_MNI152.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.shootingtpm        = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_0_IXI555_MNI152_GS.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.regstr             = 0;   matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.vox                             = 1.5;
if exp==0
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.APP                           = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.LASstr                        = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.gcutstr                       = 2;          
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.restypes.fixed                = [1 0.1]; % use 2 mm for faster default pp  
elseif exp>0 % EXPERT
  % segmentation
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.APP              = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.LASstr           = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.gcutstr          = 2;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.regstr           = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.cleanupstr       = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.NCstr            = -Inf;     % noise filter strength 
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.WMHCstr          = 0.5;     
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.WMHC             = 1;  
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.restypes.fixed   = [1 0.1]; % use 2 mm for faster default pp  
  % surfaces     
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.surface.pbtres                = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.surface.scale_cortex          = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.surface.add_parahipp          = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.surface.close_parahipp        = 0;
  % admin      
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.admin.experimental            = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.admin.ignoreErrors            = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.admin.verb                    = 2;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.admin.print                   = 2;
end
if exp>1 % DEVELOPER
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.segmentation.BVCstr              = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.cat12atlas          = { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'cat.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.brainmask           = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' , 'brainmask.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.registration.T1                  = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' ,  'T1.nii' ) };
% matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.mrf              = 1;  
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.extopts.admin.lazy          = 0;
end
% -- output ------------------------------------------------------------
% surfaces
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.surface                              = 1;  % SBM preprocessing
% ROIs
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.neuromorphometrics   = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.lpba40               = 0;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.cobra                = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.hammers              = 0;
if exp
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.ibsr               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.aal                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.mori               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.ROImenu.atlases.anatomy            = 0;
end
% GM
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.GM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.GM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.GM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.GM.dartel           = 1;
% WM
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WM.dartel           = 1;
% CSF
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.CSF.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.CSF.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.CSF.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.CSF.dartel        = 0;
  % WMH (for WMHC>0)
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WMH.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WMH.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WMH.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.WMH.dartel        = 0;
  % label map
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.label.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.label.warped      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.label.dartel      = 1;
  % global intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.bias.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.bias.warped       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.bias.dartel       = 1;
  % local intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.las.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.las.warped        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.las.dartel        = 1;
end
if exp>1 % DEVELOPER
  % TPMC (for TPM creation, e.g., to create an animal template)
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.TPMC.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.TPMC.warped       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.TPMC.mod          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.TPMC.dartel       = 0;
  % atlas maps (for atlas creation, e.g., in an animal tempate)
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.atlas.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.atlas.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.atlas.dartel      = 0;
end
% jacobian
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.jacobianwarped      = 1;
% deformation 
matlabbatch{1}.spm.tools.cat.estwrite1173plus.output.warps               = [1 1];
