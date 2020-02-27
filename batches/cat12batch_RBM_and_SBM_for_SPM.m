%-----------------------------------------------------------------------
% Job saved on 26-Oct-2016 14:08:52 by cfg_util (rev $Rev$)
% spm SPM - SPM12 (6685)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
% sanlm denoising
matlabbatch{1}.spm.tools.cat.tools.sanlm.data                     = '<UNDEFINED>';
matlabbatch{1}.spm.tools.cat.tools.sanlm.prefix                   = 'sanlm_';
matlabbatch{1}.spm.tools.cat.tools.sanlm.NCstr                    = Inf;
matlabbatch{1}.spm.tools.cat.tools.sanlm.rician                   = 0;
% SPM segment
matlabbatch{2}.spm.spatial.preproc.channel.vols(1)                = cfg_dep('Spatially adaptive non-local means denoising filter: All Output Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','vfiles'));
matlabbatch{2}.spm.spatial.preproc.channel.biasreg                = 0.001;
matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm               = 60;
matlabbatch{2}.spm.spatial.preproc.channel.write                  = [0 1];
matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,1')};
matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus                = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(1).native               = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,2')};
matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus                = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(2).native               = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,3')};
matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus                = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(3).native               = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,4')};
matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus                = 3;
matlabbatch{2}.spm.spatial.preproc.tissue(4).native               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,5')};
matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus                = 4;
matlabbatch{2}.spm.spatial.preproc.tissue(5).native               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm                  = {fullfile(spm('dir'),'tpm','TPM.nii,6')};
matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus                = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(6).native               = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).warped               = [0 0];
matlabbatch{2}.spm.spatial.preproc.warp.mrf                       = 1;
matlabbatch{2}.spm.spatial.preproc.warp.cleanup                   = 1;
matlabbatch{2}.spm.spatial.preproc.warp.reg                       = [0 0.001 0.5 0.05 0.2];
matlabbatch{2}.spm.spatial.preproc.warp.affreg                    = 'mni';
matlabbatch{2}.spm.spatial.preproc.warp.fwhm                      = 0;
matlabbatch{2}.spm.spatial.preproc.warp.samp                      = 3;
matlabbatch{2}.spm.spatial.preproc.warp.write                     = [0 0];
% CAT SPM segment
matlabbatch{3}.spm.tools.cat.estwrite_spm.data(1)                 = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{3}.spm.tools.cat.estwrite_spm.nproc                   = max(0,round(feature('numcores') ./ (1+ispc)));
matlabbatch{3}.spm.tools.cat.estwrite_spm.extopts.darteltpm       = {fullfile(spm('dir'),'toolbox','cat12','templates_volumes','Template_1_IXI555_MNI152.nii')};
matlabbatch{3}.spm.tools.cat.estwrite_spm.extopts.vox             = 1.5;
matlabbatch{3}.spm.tools.cat.estwrite_spm.extopts.ignoreErrors    = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.ROI              = 1;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.surface          = 1;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.GM.warped        = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.GM.mod           = 1;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.GM.dartel        = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.WM.warped        = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.WM.mod           = 1;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.WM.dartel        = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.CSF.warped       = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.CSF.mod          = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.CSF.dartel       = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.jacobian.warped  = 0;
matlabbatch{3}.spm.tools.cat.estwrite_spm.output.warps            = [0 0];
