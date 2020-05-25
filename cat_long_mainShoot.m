%-----------------------------------------------------------------------
% Job for longitudinal batch
% Christian Gaser
% $Id$
%-----------------------------------------------------------------------

global opts extopts output modulate dartel delete_temp ROImenu sROImenu surfaces cat

if isempty(dartel),       dartel      = 0; end
if isempty(modulate),     modulate    = 1; end
if isempty(surfaces),     surfaces    = cat_get_defaults('output.surfaces'); end
if isempty(delete_temp),  delete_temp = 1; end

write_CSF = cat_get_defaults('output.CSF.mod') > 0;

warning('off','MATLAB:DELETE:FileNotFound');

% display start
if 0 %~isempty(extopts)  
  % The idea of simply repeat the input is not optimal.  
  % You have to use the DEP output otherwise it will result in more problems. 
  mbi = 1;
  matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.inputs{1}.images = '<UNDEFINED>';
  matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.outputs          = {}; % @(x) cat_io_depin2depout;
  matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.fun              = @(x)cat_io_cprintf('blue',sprintf([...
    '================================================================================================================================================\n' ...
    'Start CAT12 longitudinal processing of \n  %s\b\b\b\n' ...
    '================================================================================================================================================\n'],...
    sprintf('%s',char( cellfun(@(s) ([s(1:end-2) '\n  '])',x,'UniformOutput',0) )) ));
else
  mbi = 0;
end


% 1) longitudinal rigid registration with final masking (SAVG space)
% -----------------------------------------------------------------------
% Here we bring all time points to the same rigid orientation, reslice the
% images, correct roughly for inhomogeneities between time points and create 
% an average image.
% ########
% RD202005: In case of strong developmental differences due to head size
%           an affine registration or John's longitudinal average is maybe 
%           required.
% ########
mbi = mbi + 1; mb_rigid = mbi; 
matlabbatch{mbi}.spm.tools.cat.tools.series.bparam          = 1e6;
matlabbatch{mbi}.spm.tools.cat.tools.series.use_brainmask   = 1;
matlabbatch{mbi}.spm.tools.cat.tools.series.reduce          = 1;
matlabbatch{mbi}.spm.tools.cat.tools.series.data            = '<UNDEFINED>';


% 2) cat12 segmentation of average image 
% -----------------------------------------------------------------------
% The average image is used for a general segmentation and registration
% to the MNI template.  The rigid segmentation is used to create an 
% individual TPM in step 3.  
mbi = mbi + 1; mb_catavg = mbi;
matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)             = cfg_dep('Longitudinal Registration: Midpoint Average', substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.estwrite.nproc               = 0;
if exist('opts','var') && ~isempty(opts)
	matlabbatch{mbi}.spm.tools.cat.estwrite.opts              = opts;
end
if exist('extopts','var') && ~isempty(extopts)
	matlabbatch{mbi}.spm.tools.cat.estwrite.extopts           = extopts;
end
if exist('output','var') && ~isempty(output)
	matlabbatch{mbi}.spm.tools.cat.estwrite.output            = output;
end
% surface estimation
matlabbatch{mbi}.spm.tools.cat.estwrite.output.surface      = surfaces;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.ROImenu.noROI= struct([]);
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.native    = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.dartel    = 2; 
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.mod       = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.native    = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.dartel    = 2; 
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.mod       = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.dartel   = 2; 
matlabbatch{mbi}.spm.tools.cat.estwrite.output.TPMC.dartel  = 2;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.bias.warped  = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.warps        = [0 0];


% 3) creating longitudinal TPM
% -----------------------------------------------------------------------
% Using a subject-specific TPM allows to stabilize the preprocessing of the
% indivividual time points, mostly of the initial affine registration and 
% the Unified segmentation that also compensates for slight structural 
% changes between the time points.  However the effects on the final AMAP
% segmenation are relatively small. 
mbi = mbi + 1; mb_tpm = mbi;
matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.files(1) = cfg_dep('CAT12: Segmentation (current release): rp1 affine Image', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','rpa', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.fstrength = 2;
matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.writeBM = 0;
matlabbatch{mbi}.spm.tools.cat.tools.createTPMlong.verb = 1;


% 4) cat12 segmentation of realigned images with prior from step 2  
% -----------------------------------------------------------------------
% In this step each time point is estimated separatelly but uses the prior
% from the SAVG - the TPM from step 3 for segmentation (and the individual 
% surface from step 2)
mbi = mbi + 1; mb_cat = mbi;
% use average image as prior for affine transformation and surface extraction
matlabbatch{mbi}.spm.tools.cat.estwrite.data(1)             = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.estwrite.nproc               = 0;
if exist('opts','var') && ~isempty(opts)
	matlabbatch{mbi}.spm.tools.cat.estwrite.opts              = opts;
end
if exist('extopts','var') && ~isempty(extopts)
	matlabbatch{mbi}.spm.tools.cat.estwrite.extopts           = extopts;
end
if exist('output','var') && ~isempty(output)
	matlabbatch{mbi}.spm.tools.cat.estwrite.output            = output;
end
% surface estimation
matlabbatch{mbi}.spm.tools.cat.estwrite.output.surface      = surfaces;
if exist('ROImenu','var') && ~isempty(ROImenu)
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.ROImenu    = ROImenu;
end
if exist('sROImenu','var') && ~isempty(sROImenu)
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.sROImenu   = sROImenu;
end
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.native    = 1;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.dartel    = dartel;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.GM.mod       = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.native    = 1;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.dartel    = dartel;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.WM.mod       = 0;
if write_CSF
  matlabbatch{mbi}.spm.tools.cat.estwrite.output.CSF.native = 1; % also write CSF?
end
matlabbatch{mbi}.spm.tools.cat.estwrite.output.bias.warped  = 0;
matlabbatch{mbi}.spm.tools.cat.estwrite.output.warps        = [1 0];
matlabbatch{mbi}.spm.tools.cat.estwrite.useprior(1)         = cfg_dep('Longitudinal Registration: Midpoint Average', substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.estwrite.opts.tpm            = cfg_dep('Longitudinal TPM creation: Longitudinal TPMs', substruct('.','val', '{}',{mb_tpm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tpm', '()',{':'}));

% 5) averaging deformations
% -----------------------------------------------------------------------
% To map the data to the MNI space, the time point specific deformations
% were averaged. 
% #######
% RD202005: In case of developemental data, we may need to use the 
%           deformation from the SAVG to deal with larger affine changes
%           due to different head size.
% #######
mbi = mbi + 1; mb_avgdef = mbi;
matlabbatch{mbi}.spm.tools.cat.tools.avg_img.data(1)  = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.tools.avg_img.output   = '';
matlabbatch{mbi}.spm.tools.cat.tools.avg_img.outdir   = {''};


% 6) creating time point specific deformation 
% -----------------------------------------------------------------------
% To reduce longitudinal changes of moving structures between time points 
% a longitudinal Shooting template is estimated.
% #######
% RD202005: In case of developemental data, we may need to use different
%           Shooting parameters (e.g., more iterations, more low-freq.
%           changes to adapt for head size changes.
% #######

lowres = 2; % define resolution in mm
if lowres
  % reduce resolution 
  % It would be also possible to use the rigid output from the time points 
  % but those depend on user definition of extopts.vox and we are more
  % flexible and probably faster and more robust this way.
  mb_lr = zeros(1,2);
  for ci = 1:2 % only GM and WM are reqired for Shooting
    mbi = mbi + 1; mb_lr(ci) = mbi; % have to do this for all shooting tissues to get the dependencies
    matlabbatch{mbi}.spm.tools.cat.tools.resize.data(1)  = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci), substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{ci}, '.','p', '()',{':'}));
    matlabbatch{mbi}.spm.tools.cat.tools.resize.res      = lowres;
    matlabbatch{mbi}.spm.tools.cat.tools.resize.method   = 'cubic';
    matlabbatch{mbi}.spm.tools.cat.tools.resize.prefix   = 'l'; % need to be another file
  end
  % Shooting low res
  mbi = mbi + 1; mb_GS = mbi;
  matlabbatch{mbi}.spm.tools.cat.tools.warp.images{1}(1) = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{mb_lr(1)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.warp.images{2}(1) = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{mb_lr(2)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.warp.dfile        = {fullfile(spm('dir'),'toolbox','cat12','cat_long_shoot_defaults.m')};

  % interpolate to acceptable resolution (higher than vox) 
  mbi = mbi + 1; mb_GSI = mbi; % have to do this for all shooting tissues to get the dependencies
  matlabbatch{mbi}.spm.tools.cat.tools.resize.data(1)    = cfg_dep('Run Shooting (create Templates): Deformation Fields', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def', '()',{':'}));
  try
    matlabbatch{mbi}.spm.tools.cat.tools.resize.res      = matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.vox;
  catch
    try 
      matlabbatch{mbi}.spm.tools.cat.tools.resize.res    = min(matlabbatch{mbi}.spm.tools.cat.estwrite.extopts.registration.vox);
    catch
      matlabbatch{mbi}.spm.tools.cat.tools.resize.res    = 1; 
    end
  end
  matlabbatch{mbi}.spm.tools.cat.tools.resize.method     = 'cubic';
  matlabbatch{mbi}.spm.tools.cat.tools.resize.prefix     = ''; % has to be another name
else
  % Shooting full res
  mbi = mbi + 1; mb_GS = mbi;
  matlabbatch{mbi}.spm.tools.cat.tools.warp.images{1}(1) = cfg_dep('CAT12: Segmentation (current release): p1 Image', substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.warp.images{2}(1) = cfg_dep('CAT12: Segmentation (current release): p2 Image', substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.warp.dfile        = {fullfile(spm('dir'),'toolbox','cat12','cat_long_shoot_defaults.m')};
end




% 7) applying time point deformations to rigid native segmentations
% -----------------------------------------------------------------------
% this is the first simple approach with full resolution
mb_aGS = zeros(1,2 + single(write_CSF));
for ci = 1:2 + single(write_CSF)
  mbi = mbi + 1; mb_aGS(ci) = mbi;
  if lowres
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.field(1)   = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{mb_GSI}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'}));
  else 
    matlabbatch{mbi}.spm.tools.cat.tools.defs2.field(1)   = cfg_dep('Run Shooting (create Templates): Deformation Fields', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def', '()',{':'}));
  end
  matlabbatch{mbi}.spm.tools.cat.tools.defs2.images{1}(1) = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci), substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{ci}, '.','p', '()',{':'}));
  matlabbatch{mbi}.spm.tools.cat.tools.defs2.interp       = 1;
  if modulate, matlabbatch{mbi}.spm.tools.cat.tools.defs2.modulate  = modulate; end  % modulation option for applying deformations
end


% 8) applying deformations to time point optimized native segmentations
% -----------------------------------------------------------------------
% applying deformations to tissues
mbi = mbi + 1; 
matlabbatch{mbi}.spm.tools.cat.tools.defs.field1(1)       = cfg_dep('Image Average: Average Image: ', substruct('.','val', '{}',{mb_avgdef}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
for ci = 1:2 + single(write_CSF) % fill image sets
  matlabbatch{mbi}.spm.tools.cat.tools.defs.images(ci)    = cfg_dep('Apply deformations (many subjects): All Output Files', substruct('.','val', '{}',{mb_aGS(ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','vfiles'));
end
matlabbatch{mbi}.spm.tools.cat.tools.defs.interp          = 1;
if modulate, matlabbatch{mbi}.spm.tools.cat.tools.defs.modulate = modulate; end  % modulation option for applying deformations


% 9) applying deformations to average T1 image
% -----------------------------------------------------------------------
mbi = mbi + 1; 
matlabbatch{mbi}.spm.tools.cat.tools.defs.field1(1)       = cfg_dep('Image Average: Average Image: ', substruct('.','val', '{}',{mb_avgdef}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{mbi}.spm.tools.cat.tools.defs.images(1)       = cfg_dep('Longitudinal Registration: Midpoint Average', substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{mbi}.spm.tools.cat.tools.defs.interp          = 1;
matlabbatch{mbi}.spm.tools.cat.tools.defs.modulate        = 0;




% 10) delete temporary files
% -----------------------------------------------------------------------
if delete_temp
  mbi = mbi + 1; 
  c = 1;
  % time point specific preprocessing data
  for ci = 1:2 + single(write_CSF)
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep(sprintf('CAT12: Segmentation (current release): p%d Image',ci), substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{ci}, '.','p', '()',{':'})); c = c+1;
  end
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{mb_cat}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'})); c = c+1;
  % average preprocessing data
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT Report JGP', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','catreportjpg', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT Report', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','catxml', '()',{':'}));
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): CAT log-file', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','catlog', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): rp1 affine Image', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','rpa', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): rp2 affine Image', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','rpa', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): rp3 affine Image', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','rpa', '()',{':'})); c = c+1;  
  % surfaces
  if surfaces
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Central Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhcentral', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Sphere Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhsphere', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Spherereg Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhspherereg', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Thickness', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhthickness', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Left Pbt', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhpbt', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Central Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rhcentral', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Sphere Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rhsphere', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Spherereg Surface', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rhspherereg', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Thickness', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rhthickness', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('CAT12: Segmentation (current release): Right Pbt', substruct('.','val', '{}',{mb_catavg}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','rhpbt', '()',{':'})); c = c+1;
  end
  % timepoint deformations
  if lowres
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{mb_lr(1)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'})); c = c+1;
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Resize images: Resized', substruct('.','val', '{}',{mb_lr(2)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res', '()',{':'})); c = c+1;
  end
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Template (0)', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Velocity Fields', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','vel', '()',{':'})); c = c+1;
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Deformation Fields', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def', '()',{':'})); c = c+1; 
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Run Shooting (create Templates): Jacobian Fields', substruct('.','val', '{}',{mb_GS}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','jac', '()',{':'})); c = c+1;
  for ci = 1:2 + single(write_CSF)
    matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.files(c) = cfg_dep('Apply deformations (many subjects): All Output Files', substruct('.','val', '{}',{mb_aGS(ci)}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','vfiles')); c = c+1;
  end
  % final command of this batch 
  matlabbatch{mbi}.cfg_basicio.file_dir.file_ops.file_move.action.delete  = false;
end


% 9) display finishing
%-----------------------------------------------------------------------
mbi = mbi + 1; 
matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.inputs{1}.images(1)  = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{mb_rigid}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.outputs              = {};
matlabbatch{mbi}.cfg_basicio.run_ops.call_matlab.fun                  = @(x)cat_io_cprintf('blue',sprintf([...
  '================================================================================================================================================\n' ...
  'Finished CAT12 longitudinal processing of \n  %s\b\b\b\n' ...
  '================================================================================================================================================\n'],...
  sprintf('%s',char( cellfun(@(s) ([s(1:end-2) '\n  '])',x,'UniformOutput',0) )) ));
