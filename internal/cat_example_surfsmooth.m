%% function cat_example_volsurfsmooth

resdir = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/tmp'; 
%S = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_surfaces/lh.central.freesurfer.gii';
%S = '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cat_defaults_rd/BO/surf/lh.central.Collins.gii'; 
%S = '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cat_defaults_rd/BO/surf/s15mm.lh.thickness.resampled.Collins.gii'; 
S = '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cat_defaults_rd/BO/surf/s15mm.rh.thickness.resampled.Collins.gii'; 
P = '/Users/dahnke/Neuroimaging/spm12/toolbox/cat12/templates_1.50mm/Template_T1_IXI555_MNI152.nii';

if ~exist(resdir,'dir'), mkdir(resdir); end

SG = gifti(S);

s = 32; 
ROI.vertex = 59595; % motorcortex
ROI.vertex = 34111; % superior temporal lobe
ROI.vertex = [59595,34111]; % superior temporal lobe
ROI.vertex = [142556,138036,149965,65103]; % superior temporal lobe
ROI.vertex = [138036,96286,62104]; % superior temporal lobe
ROI.xyz    = SG.vertices(ROI.vertex,:);
sinfo = cat_surf_info(S);

% -- Volume smoothing ----------------------------------------------------
% smooth volume
V = spm_vol(P); 
Y = zeros(V.dim,'single');
mati = spm_imatrix(V.mat); % .*sign(mati(7:9))
for pi=1:numel(ROI.vertex)
  ROI.xyz2 = round((ROI.xyz(pi,:) - mati(1:3))./mati(7:9)); 
  Y(sub2ind(V.dim,ROI.xyz2(1),ROI.xyz2(2),ROI.xyz2(3))) = 1000; 
end
spm_smooth(Y,Y,repmat(s/1.5,1,3));
Y = Y + 1;
V.fname = fullfile(resdir,sprintf('cat_example_volsurfsmooth_v%0.0fs%0.0f_vol.%s.nii',ROI.vertex(1),s,sinfo.name));
V.dt(1) = 16;
spm_write_vol(V,Y);

% project volume to surface
job = struct(); 
job.data_mesh_lh  = S;
job.data_vol      = {V.fname}; 
job.mapping.abs_mapping.startpoint = -0.5; 
job.mapping.abs_mapping.endpoint   = +0.5; 
job.mapping.abs_mapping.stepsize   = 1; 
if sinfo.resampled, res = 'resampled.'; else res = ''; end
job.datafieldname = sprintf('cat_example_volsurfsmooth_v%0.0fs%0.0f_vol',ROI.vertex(1),s);
cat_surf_vol2surf(job);

% -- Surface smoothing ---------------------------------------------------
% create mapping
copyfile(S,resdir);
copyfile(sinfo.Psphere,resdir);
SG.cdata = ones(size(SG.vertices,1),1); 
for pi=1:numel(ROI.vertex)
  SG.cdata(ROI.vertex(pi)) = 1000; 
end
if sinfo.resampled, res = 'resampled.'; else res = ''; end
SS1 = fullfile(resdir,sprintf('lh.cat_example_volsurfsmooth_v%ds%0.0f_catblur.%s%s.gii',ROI.vertex(1),s,res,sinfo.name)); 
SS2 = fullfile(resdir,sprintf('lh.cat_example_volsurfsmooth_v%ds%0.0f_spmblur.%s%s.gii',ROI.vertex(1),s,res,sinfo.name)); 
save(gifti(SG),SS1);
save(gifti(SG),SS2);
% smoothing 
job = struct(); 
job.data = {SS1};
job.fwhm = s*1; 
job.datafieldname = sprintf('surfsmooth_%d',ROI.vertex(1));
job.catblur = 1; 
cat_surf_smooth(job);
% smoothing 2
job.data = {SS2};
job.catblur = 0; 
cat_surf_smooth(job);
