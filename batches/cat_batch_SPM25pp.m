function matlabbatch = cat_batch_SPM25pp()

ngaus = [1 1 2 3 4 2]; 
Ptpm  = fullfile( spm('path'), 'tpm', 'TPM.nii'); 

%% copy files to BIDS dir


%% SPM preprocessing
mi = mi + 1; 
matlabbatch{mi}.spm.spatial.preproc.channel.vols = { files }; % '<UNDEFINED>';
matlabbatch{mi}.spm.spatial.preproc.channel.biasreg = 0.0001;
matlabbatch{mi}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{mi}.spm.spatial.preproc.channel.write = [0 0];
for ti = 1:6
  matlabbatch{mi}.spm.spatial.preproc.tissue(ti).tpm    = {fprintf('%s,%d',Ptpm,ti)};
  matlabbatch{mi}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
  matlabbatch{mi}.spm.spatial.preproc.tissue(ti).native = [ti<6 0];
  matlabbatch{mi}.spm.spatial.preproc.tissue(ti).warped = [0 0];
end
matlabbatch{mi}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{mi}.spm.spatial.preproc.warp.reg = [0 0 0.1 0.01 0.04];
matlabbatch{mi}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{mi}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{mi}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{mi}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{mi}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{mi}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];

%% CAT12 SPM surface processing
mi = mi + 1; 
matlabbatch{mi}.spm.tools.cat.estwrite_spm.data(1) = ...
  cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{mi}.spm.tools.cat.estwrite_spm.nproc = 6;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.segmentation.spmAMAP = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.segmentation.WMHtpm = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.segmentation.BVtpm = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.segmentation.SLtpm = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.T1 = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.brainmask = {''};
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.cat12atlas = cat_get_defaults('extopts.cat12atlas');
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.darteltpm = cat_get_defaults('extopts.darteltpm');
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.shootingtpm = cat_get_defaults('extopts.shootingtpm');
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.regstr = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.bb = 12;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.registration.vox = 1.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.pbtres = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.pbtmethod = 'pbtsimple';
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.SRP = 24;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.vdist = 2;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.scale_cortex = 0.5;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.add_parahipp = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.surface.close_parahipp = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.experimental = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.new_release = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.lazy = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.ignoreErrors = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.verb = 2;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.extopts.admin.print = 2;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.BIDS.BIDSyes2.BIDSfolder = '../../CAT12.9_2701_SPMamap_CS24';
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.surface = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.ROImenu.noROI = struct([]);
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.mod = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.GM.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.mod = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.WM.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.mod = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.CSF.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.TPMC.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.TPMC.mod = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.TPMC.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.native = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.label.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.bias.warped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.bias.dartel = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.labelnative = 1;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.jacobianwarped = 0;
matlabbatch{mi}.spm.tools.cat.estwrite_spm.output.warps = [0 0];

