function batch_surface_pipeline

pth = pwd;

spm_jobman('initcfg');

%__________________________________________________________________________
matlabbatch{1}.spm.tools.cat.stools.surfextract.data_surf = {
                                                             fullfile(pth,'surf/lh.central.IXI002-Guys-0828-T1.gii')
                                                             fullfile(pth,'surf/lh.central.IXI016-Guys-0697-T1.gii')
                                                             fullfile(pth,'surf/lh.central.IXI017-Guys-0698-T1.gii')
                                                             fullfile(pth,'surf/lh.central.IXI019-Guys-0702-T1.gii')
                                                             };
matlabbatch{1}.spm.tools.cat.stools.surfextract.GI = 1;
matlabbatch{1}.spm.tools.cat.stools.surfextract.FD = 0;
matlabbatch{1}.spm.tools.cat.stools.surfextract.SD = 0;
matlabbatch{1}.spm.tools.cat.stools.surfextract.nproc = 0;
%__________________________________________________________________________
matlabbatch{2}.spm.tools.cat.stools.surfresamp.data_surf(1) = cfg_dep('Extract additional surface parameters: Left MNI gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));
matlabbatch{2}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.fwhm_surf = 12;
matlabbatch{2}.spm.tools.cat.stools.surfresamp.nproc = 0;
%__________________________________________________________________________
matlabbatch{3}.spm.tools.cat.tools.check_homogeneity.data{1}(1) = cfg_dep('Resample and Smooth Surface Data: Merged resampled MNIgyrification', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sample', '()',{1}, '.','Psdata'));
matlabbatch{3}.spm.tools.cat.tools.check_homogeneity.sel_xml.data_xml = {
                                                               fullfile(pth,'report/cat_IXI002-Guys-0828-T1.xml')
                                                               fullfile(pth,'report/cat_IXI016-Guys-0697-T1.xml')
                                                               fullfile(pth,'report/cat_IXI017-Guys-0698-T1.xml')
                                                               fullfile(pth,'report/cat_IXI019-Guys-0702-T1.xml')
                                                               };
matlabbatch{3}.spm.tools.cat.tools.check_homogeneity.c = cell(1, 0);
%__________________________________________________________________________
matlabbatch{4}.spm.tools.cat.stools.surf2roi.cdata{1}(1) = cfg_dep('Extract additional surface parameters: Left MNI gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));
matlabbatch{4}.spm.tools.cat.stools.surf2roi.rdata = {fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces/lh.aparc_DK40.freesurfer.annot')};
%__________________________________________________________________________
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.data_T2x = {fullfile(pth,'analysis/surface/spmT_0002.gii')};
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.sel = 2;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{5}.spm.tools.cat.tools.T2x_surf.conversion.inverse = 0;
matlabbatch{1}.spm.tools.cat.tools.T2x_surf.conversion.cluster.none = 1;
%__________________________________________________________________________
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.data_F2x = {fullfile(pth,'analysis/surface/spmF_0001.gii')};
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.sel = 2;
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{6}.spm.tools.cat.tools.F2x_surf.conversion.cluster.none = 1;
%__________________________________________________________________________
matlabbatch{7}.spm.tools.cat.tools.calcvol.data_xml = {
                                                       fullfile(pth,'report/cat_IXI002-Guys-0828-T1.xml')
                                                       fullfile(pth,'report/cat_IXI016-Guys-0697-T1.xml')
                                                       fullfile(pth,'report/cat_IXI017-Guys-0698-T1.xml')
                                                       fullfile(pth,'report/cat_IXI019-Guys-0702-T1.xml')
                                                       };
matlabbatch{7}.spm.tools.cat.tools.calcvol.calcvol_TIV = 1;
matlabbatch{7}.spm.tools.cat.tools.calcvol.calcvol_name = 'TIV.txt';
%__________________________________________________________________________

matlabbatch{8}.spm.tools.cat.stools.vol2surf.data_vol = {fullfile(pth,'IXI002-Guys-0828-T1.nii,1')};
matlabbatch{8}.spm.tools.cat.stools.vol2surf.data_mesh_lh = {fullfile(pth,'surf/lh.central.IXI002-Guys-0828-T1.gii')};
matlabbatch{8}.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
matlabbatch{8}.spm.tools.cat.stools.vol2surf.interp = {'linear'};
matlabbatch{8}.spm.tools.cat.stools.vol2surf.datafieldname = 'intensity';
matlabbatch{8}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.class = 'GM';
matlabbatch{8}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.startpoint = -0.6;
matlabbatch{8}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.steps = 7;
matlabbatch{8}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.endpoint = 0.6;
%_____________________________________________________________
matlabbatch{9}.spm.tools.cat.stools.renderresults.cdata(1) = cfg_dep('Threshold and transform spmT surfaces: Transform & Threshold spm surfaces', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Pname'));
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.surface = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.view = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.texture = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.transparency = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.colormap = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.invcolormap = 0;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.background = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.render.showfilename = 1;
matlabbatch{9}.spm.tools.cat.stools.renderresults.stat.threshold = 0;
matlabbatch{9}.spm.tools.cat.stools.renderresults.stat.hide_neg = 0;
matlabbatch{9}.spm.tools.cat.stools.renderresults.fparts.outdir = {''};
matlabbatch{9}.spm.tools.cat.stools.renderresults.fparts.prefix = 'render_';
matlabbatch{9}.spm.tools.cat.stools.renderresults.fparts.suffix = '';
%__________________________________________________________________________
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Extract additional surface parameters: Left MNI gyrification', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lPGI', '()',{':'}));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Resample and Smooth Surface Data: Merged resampled MNIgyrification', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sample', '()',{1}, '.','Psdata'));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Extract ROI-based surface values: Extracted Surface ROIs', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','xmlname', '()',{':'}));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Threshold and transform spmT surfaces: Transform & Threshold spm surfaces', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Pname'));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Threshold and transform spmF surfaces: Transform & Threshold spm surfaces', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Pname'));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.files(7) = cfg_dep('Map Volume (Native Space) to Individual Surface: Right mapped values', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rh'));
matlabbatch{10}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
%__________________________________________________________________________
spm_jobman('run',matlabbatch);
%__________________________________________________________________________

%__________________________________________________________________________
cat_stat_analyze_ROIs(fullfile(pth,'analysis/surface/SPM.mat'), 0.05, 1);

%__________________________________________________________________________
clear matlabbatch
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {pth};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = 'logP|png|gyrification';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
spm_jobman('run',matlabbatch);

