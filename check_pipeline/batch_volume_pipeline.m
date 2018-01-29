function batch_volume_pipeline

pth=pwd;

spm_jobman('initcfg');

matlabbatch{1}.spm.tools.cat.tools.calcvol.data_xml = {
  fullfile(pth,'report/cat_IXI002-Guys-0828-T1.xml')
  fullfile(pth,'report/cat_IXI016-Guys-0697-T1.xml')
  fullfile(pth,'report/cat_IXI017-Guys-0698-T1.xml')
  fullfile(pth,'report/cat_IXI019-Guys-0702-T1.xml')
                                                       };
matlabbatch{1}.spm.tools.cat.tools.calcvol.calcvol_TIV = 1;
matlabbatch{1}.spm.tools.cat.tools.calcvol.calcvol_name = 'TIV.txt';

%__________________________________________________________________________
matlabbatch{2}.spm.tools.cat.tools.check_cov.data_vol = {
                                                         {
  fullfile(pth,'mri/mwp1IXI002-Guys-0828-T1.nii,1')
  fullfile(pth,'mri/mwp1IXI016-Guys-0697-T1.nii,1')
  fullfile(pth,'mri/mwp1IXI017-Guys-0698-T1.nii,1')
  fullfile(pth,'mri/mwp1IXI019-Guys-0702-T1.nii,1')
                                                         }
                                                         }';
matlabbatch{2}.spm.tools.cat.tools.check_cov.data_xml = {
  fullfile(pth,'report/cat_IXI002-Guys-0828-T1.xml')
  fullfile(pth,'report/cat_IXI016-Guys-0697-T1.xml')
  fullfile(pth,'report/cat_IXI017-Guys-0698-T1.xml')
  fullfile(pth,'report/cat_IXI019-Guys-0702-T1.xml')
                                                         };
matlabbatch{2}.spm.tools.cat.tools.check_cov.gap = 3;
matlabbatch{2}.spm.tools.cat.tools.check_cov.c{1}(1) = cfg_dep('Estimate TIV and global tissue volumes: TIV', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','calcvol'));

%__________________________________________________________________________
matlabbatch{3}.spm.tools.cat.tools.T2x.data_T2x = {fullfile(pth,'analysis/volume/spmT_0002.nii,1')};
matlabbatch{3}.spm.tools.cat.tools.T2x.conversion.sel = 2;
matlabbatch{3}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{3}.spm.tools.cat.tools.T2x.conversion.inverse = 1;
matlabbatch{3}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{3}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{3}.spm.tools.cat.tools.T2x.atlas = 'Neuromorphometrics';

%__________________________________________________________________________
matlabbatch{4}.spm.tools.cat.tools.F2x.data_F2x = {fullfile(pth,'analysis/volume/spmF_0001.nii,1')};
matlabbatch{4}.spm.tools.cat.tools.F2x.conversion.sel = 2;
matlabbatch{4}.spm.tools.cat.tools.F2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{4}.spm.tools.cat.tools.F2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{4}.spm.tools.cat.tools.F2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{4}.spm.tools.cat.tools.F2x.atlas = 'Neuromorphometrics';

%__________________________________________________________________________
spm_jobman('run',matlabbatch);
%__________________________________________________________________________

OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152_GS.nii');
OV.reference_range = [0.2 1.0];                         % intensity range for reference image
OV.opacity = Inf;                                      % transparence value for overlay (<1)
OV.cmap    = jet;                                      % colormap for overlay
OV.name = char(fullfile(pth,'analysis/volume/logP_neg_p0.1_pkFWE5_k590_bi.nii'));
OV.range   =[[3 6]];
OV.slices_str = char('-20:5:45');
OV.transform = char('coronal');
OV.xy = [3 5];
OV.save = 'result.png';
OV.labels.format = '%3.1f';
cat_vol_slice_overlay(OV)

%__________________________________________________________________________
cat_stat_analyze_ROIs(fullfile(pth,'analysis/volume/SPM.mat'), 0.05, 1);

%__________________________________________________________________________
clear matlabbatch
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {pth};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = 'logP|png|gyrification';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
spm_jobman('run',matlabbatch);

