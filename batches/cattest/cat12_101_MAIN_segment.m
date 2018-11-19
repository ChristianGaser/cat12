% ---------------------------------------------------------------------
% Test batch "Segment Data" for VBM, SBM, and RBM preprocessing of 
% cat_tst_cat.est.
% ---------------------------------------------------------------------
% Robert Dahnke
% $Id$

if ~exist('files_human','var')
  files_human = {'<UNDEFINED>'};
  exp = cat_get_defaults('extopts.expertgui'); 
elseif isempty(files_human)
  return;
end

% batch
% -- opts --------------------------------------------------------------
matlabbatch{1}.spm.tools.cat.estwrite.data                       = files_human;
matlabbatch{1}.spm.tools.cat.estwrite.nproc                      = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                   = { fullfile( spm('dir') , 'tpm' , 'TPM.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg                = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr               = 0.5;
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus               = [1 1 2 3 4 2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.biasreg             = 0.001;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.biasfwhm            = 60;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg             = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.samp                = 3;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.tol                 = 1e-4;
end
% -- extopts -----------------------------------------------------------
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP                             = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr                          = 0.5;       
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr                         = 0.5;       
% registration
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm          = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_1_IXI555_MNI152.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm        = ...
  { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'Template_0_IXI555_MNI152_GS.nii' ) };
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr             = 0;   matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox                             = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed                  = [1 0.1]; % use 2 mm for faster   
if exp>0 % EXPERT
  % segmentation
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP              = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr           = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr          = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.regstr           = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr       = 0.5;       
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr            = -Inf;     % noise filter strength 
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHCstr          = 0.5;     
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC             = 1;  
  % surfaces     
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres                = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex          = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp          = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp        = 0;
  % admin      
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors  = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb          = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print         = 2;
end
if exp>1 % DEVELOPER
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr              = 0;       
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas          = { fullfile( spm('dir') , 'toolbox' , 'cat12' , 'templates_1.50mm' , 'cat.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask           = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' , 'brainmask.nii' ) };
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1                  = { fullfile( spm('dir') , 'toolbox' , 'FieldMap' ,  'T1.nii' ) };
% matlabbatch{1}.spm.tools.cat.estwrite.extopts.mrf              = 1;  
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy          = 0;
end
% -- output ------------------------------------------------------------
% surfaces
matlabbatch{1}.spm.tools.cat.estwrite.output.surface                              = 1;  % SBM preprocessing
% ROIs
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics   = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40               = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra                = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers              = 0;
if exp
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal                = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori               = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy            = 0;
end
% GM
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel           = 1;
% WM
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native           = 1;
if exp>0, matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1; end
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod              = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel           = 1;
% CSF
if exp>0 % EXPERT
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel        = 0;
  % WMH (for WMHC>0)
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped        = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod           = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel        = 0;
  % label map
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel      = 1;
  % global intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel       = 1;
  % local intensity normalized
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.native        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped        = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel        = 1;
end
if exp>1 % DEVELOPER
  % thickness
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel         = 0;
  % Stroke 
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native         = 0; 
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped         = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod            = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel         = 0;
  % TPMC (for TPM creation, e.g., to create an animal template)
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native       = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped       = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod          = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel       = 0;
  % atlas maps (for atlas creation, e.g., in an animal tempate)
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native      = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped      = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel      = 0;
end
% jacobian
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped      = 1;
% deformation 
matlabbatch{1}.spm.tools.cat.estwrite.output.warps               = [1 1];
