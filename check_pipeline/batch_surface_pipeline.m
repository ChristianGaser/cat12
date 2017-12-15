function batch_surface_pipeline

pth=pwd;

spm_jobman('initcfg');

matlabbatch{1}.spm.tools.cat.stools.surfextract.data_surf = {
  fullfile(pth,'surf/lh.central.IXI002-Guys-0828-T1.gii')
  fullfile(pth,'surf/lh.central.IXI016-Guys-0697-T1.gii')
                                                             };
matlabbatch{1}.spm.tools.cat.stools.surfextract.GI = 1;
matlabbatch{1}.spm.tools.cat.stools.surfextract.FD = 0;
matlabbatch{1}.spm.tools.cat.stools.surfextract.SD = 0;
matlabbatch{1}.spm.tools.cat.stools.surfextract.nproc = 0;

%__________________________________________________________________________
matlabbatch{2}.spm.tools.cat.stools.surfresamp.data_surf = {
  fullfile(pth,'surf/lh.thickness.IXI002-Guys-0828-T1')
  fullfile(pth,'surf/lh.thickness.IXI016-Guys-0697-T1')
  fullfile(pth,'surf/lh.thickness.IXI017-Guys-0698-T1')
  fullfile(pth,'surf/lh.thickness.IXI019-Guys-0702-T1')
                                                            };
matlabbatch{2}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.fwhm_surf = 15;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.nproc = 0;

%__________________________________________________________________________
matlabbatch{3}.spm.tools.cat.stools.check_mesh_cov.data_surf{1}(1) = cfg_dep('Resample and Smooth Surface Data: Merged Resample & Smooth', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','Psdata', '()',{':'}));
matlabbatch{3}.spm.tools.cat.stools.check_mesh_cov.data_xml = {
  fullfile(pth,'report/cat_IXI002-Guys-0828-T1.xml')
  fullfile(pth,'report/cat_IXI016-Guys-0697-T1.xml')
  fullfile(pth,'report/cat_IXI017-Guys-0698-T1.xml')
  fullfile(pth,'report/cat_IXI019-Guys-0702-T1.xml')
                                                               };
matlabbatch{3}.spm.tools.cat.stools.check_mesh_cov.c = cell(1, 0);

%__________________________________________________________________________
matlabbatch{4}.spm.tools.cat.stools.surf2roi.cdata{1}(1) = cfg_dep('Extract additional surface parameters: Left gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));

%__________________________________________________________________________
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.data_T2x = {fullfile(pth,'analysis/surface/spmT_0002.gii')};
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.sel = 2;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.inverse = 0;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.cluster.fwe2.noniso = 1;

%__________________________________________________________________________
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.data_F2x = {fullfile(pth,'analysis/surface/spmF_0001.gii')};
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.sel = 2;
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.cluster.En.noniso = 1;

%__________________________________________________________________________
matlabbatch{7}.spm.tools.cat.stools.surfresamp.data_surf(1) = cfg_dep('Extract additional surface parameters: Left gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));
matlabbatch{7}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{7}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{7}.spm.tools.cat.stools.surfresamp.fwhm_surf = 15;
matlabbatch{7}.spm.tools.cat.stools.surfresamp.nproc = 0;

%__________________________________________________________________________
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Extract additional surface parameters: Left gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Extract additional surface parameters: Right gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rPGI', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Resample and Smooth Surface Data: Merged Resample & Smooth', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','Psdata', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Threshold and transform spmT surfaces: Transform & Threshold spm surfaces', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','Pname', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Threshold and transform spmF surfaces: Transform & Threshold spm surfaces', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','Pname', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Resample and Smooth Surface Data: Merged Resample & Smooth', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','Psdata', '()',{':'}));
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

%__________________________________________________________________________
spm_jobman('run',matlabbatch);
%__________________________________________________________________________

%__________________________________________________________________________
cat_stat_analyze_ROIs(fullfile(pth,'analysis/surface/SPM.mat'), 0.05, 1);