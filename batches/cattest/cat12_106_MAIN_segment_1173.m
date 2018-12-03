% ---------------------------------------------------------------------
% Test batch "Segment Data" for VBM, SBM, and RBM preprocessing of 
% cat_tst_cat.est.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id$

if ~exist('files_human','var')
  files_human1173 = {'<UNDEFINED>'};
  exp = cat_get_defaults1173('extopts1173.expertgui'); 
elseif isempty(files_human1173)
  return;
end

% batch
% -- opts --------------------------------------------------------------
matlabbatch{1}.spm.tools.cat.estwrite1173.data                       = files_human1173;
matlabbatch{1}.spm.tools.cat.estwrite1173.nproc                      = 0;
matlabbatch{1}.spm.tools.cat.estwrite1173.opts.tpm                   = { fullfile( spm('dir') , 'tpm' , 'TPM.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173.opts.affreg                = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite1173.opts.biasstr               = 0.5;
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite1173.opts.ngaus               = [1 1 2 3 4 2];
  matlabbatch{1}.spm.tools.cat.estwrite1173.opts.biasreg             = 0.001;
  matlabbatch{1}.spm.tools.cat.estwrite1173.opts.biasfwhm            = 60;
  matlabbatch{1}.spm.tools.cat.estwrite1173.opts.warpreg             = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.cat.estwrite1173.opts.samp                = 3;
end
% -- extopts -----------------------------------------------------------
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.APP                             = 1070;
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.LASstr                          = 0.5;       
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.gcutstr                         = 0.5;       
% registration
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.darteltpm          = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_1_IXI555_MNI152.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.shootingtpm        = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_0_IXI555_MNI152_GS.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.regstr             = 0;   matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.vox                             = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.restypes.fixed                  = [1 0.1]; % use 2 mm for faster 
if exp>0 % EXPERT
  % segmentation
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.APP              = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.LASstr           = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.gcutstr          = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.regstr           = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.cleanupstr       = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.NCstr            = -Inf;     % noise filter strength 
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.WMHCstr          = 0.5;     
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.WMHC             = 1;  
  % surfaces     
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.surface.pbtres                = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.surface.scale_cortex          = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.surface.add_parahipp          = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.surface.close_parahipp        = 0;
  % admin      
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.admin.experimental  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.admin.ignoreErrors  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.admin.verb          = 2;
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.admin.print         = 2;
end
if exp>1 % DEVELOPER
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.segmentation.BVCstr              = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.cat12atlas          = { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'cat.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.brainmask           = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' , 'brainmask.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.registration.T1                  = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' ,  'T1.nii' ) };
% matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.mrf              = 1;  
  matlabbatch{1}.spm.tools.cat.estwrite1173.extopts.admin.lazy          = 0;
end
% -- output ------------------------------------------------------------
% surfaces
matlabbatch{1}.spm.tools.cat.estwrite1173.output.surface                              = 1;  % SBM preprocessing
% ROIs
matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.neuromorphometrics   = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.lpba40               = 0;
matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.cobra                = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.hammers              = 0;
if exp
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.ibsr               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.aal                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.mori               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.ROImenu.atlases.anatomy            = 0;
end
% GM
matlabbatch{1}.spm.tools.cat.estwrite1173.output.GM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite1173.output.GM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite1173.output.GM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173.output.GM.dartel           = 1;
% WM
matlabbatch{1}.spm.tools.cat.estwrite1173.output.WM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite1173.output.WM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite1173.output.WM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite1173.output.WM.dartel           = 1;
% CSF
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.CSF.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.CSF.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.CSF.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.CSF.dartel        = 0;
  % WMH (for WMHC>0)
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.WMH.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.WMH.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.WMH.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.WMH.dartel        = 0;
  % label map
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.label.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.label.warped      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.label.dartel      = 1;
  % global intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.bias.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.bias.warped       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.bias.dartel       = 1;
  % local intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.las.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.las.warped        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.las.dartel        = 1;
end
if exp>1 % DEVELOPER
  % TPMC (for TPM creation, e.g., to create an animal template)
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.TPMC.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.TPMC.warped       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.TPMC.mod          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.TPMC.dartel       = 0;
  % atlas maps (for atlas creation, e.g., in an animal tempate)
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.atlas.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.atlas.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite1173.output.atlas.dartel      = 0;
end
% jacobian
matlabbatch{1}.spm.tools.cat.estwrite1173.output.jacobianwarped      = 1;
% deformation 
matlabbatch{1}.spm.tools.cat.estwrite1173.output.warps               = [1 1];
