function cg_vbm_defaults
% Sets the defaults for VBM
% FORMAT cg_vbm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

global vbm

% Estimation options
%=======================================================================
vbm.opts.tpm       = {fullfile(spm('dir'),'tpm','TPM.nii')}; % TPM.nii
vbm.opts.ngaus     = [1 1 2 3 4 2];  % Gaussians per class
vbm.opts.affreg    = 'mni';    % Affine regularisation
vbm.opts.warpreg   = [0 0.001 0.5 0.05 0.2];      % Warping regularisation
vbm.opts.biasreg   = 0.001; % Bias regularisation
vbm.opts.biasfwhm  = 60;   % Bias FWHM
vbm.opts.samp      = 3;      % Sampling distance


% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

% bias corrected
vbm.output.bias.native = 0;
vbm.output.bias.warped = 0;
vbm.output.bias.affine = 0;

% bias and noise corrected
vbm.output.mnT.native = 0;
vbm.output.mnT.warped = 1;
vbm.output.mnT.affine = 0;

% global intensity, bias and noise corrected
vbm.output.mgT.native = 0;
vbm.output.mgT.warped = 0;
vbm.output.mgT.mod    = 0;
vbm.output.mgT.dartel = 0;

% local intensity, bias and noise corrected
% WARNING: This is the resulting T1 image for the local adaptive segmentation
% pipeline (vbm.extopts.LAS == 1). Do not use this map for any statistical analysis!
vbm.output.mlT.native = 0;
vbm.output.mlT.warped = 0;
vbm.output.mlT.mod    = 0;
vbm.output.mlT.dartel = 0;

% preprocessing changes map
% WARNING: This map describes the changes of the preprocessing and is under development. 
vbm.output.pcT.native = 0;
vbm.output.pcT.warped = 0;
vbm.output.pcT.mod    = 0;
vbm.output.pcT.dartel = 0;

% GM tissue maps
vbm.output.GM.native = 0;
vbm.output.GM.warped = 0;
vbm.output.GM.mod    = 2;
vbm.output.GM.dartel = 0;

% WM tissue maps
vbm.output.WM.native = 0;
vbm.output.WM.warped = 0;
vbm.output.WM.mod    = 2;
vbm.output.WM.dartel = 0;

% CSF tissue maps
vbm.output.CSF.native = 0;
vbm.output.CSF.warped = 0;
vbm.output.CSF.mod    = 0;
vbm.output.CSF.dartel = 0;

% label
vbm.output.label.native = 0; 
vbm.output.label.warped = 0;
vbm.output.label.dartel = 0;

% GM thickness maps
vbm.output.th1T.native = 0;
vbm.output.th1T.warped = 0;
vbm.output.th1T.dartel = 0;

% WARNING: WM thickness is under development! 
vbm.output.th2T.native = 0;
vbm.output.th2T.warped = 0;
vbm.output.th2T.dartel = 0;

% partitioning atlas maps
vbm.output.l1T.native = 0; 
vbm.output.l1T.warped = 0; 
vbm.output.l1T.dartel = 0; 

vbm.output.l2T.native = 0;
vbm.output.l2T.warped = 0;
vbm.output.l2T.dartel = 0;

% jacobian determinant 0/1 (none/yes)
vbm.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
vbm.output.warps = [0 0];

% Extended writing options
%=======================================================================
vbm.extopts.dartelwarp  = 1; % dartel normalization: 0 - spm default; 1 - yes
vbm.extopts.darteltpm   = {fullfile(spm('dir'),'toolbox','vbm12','Template_1_IXI550_MNI152.nii')}; % Indicate first Dartel template
vbm.extopts.print       = 1; % Display and print results

% bias correction options
%=======================================================================
vbm.bias.nits_bias    = 8;
vbm.bias.biasfwhm     = 60;
vbm.bias.biasreg      = 1e-6;
vbm.bias.lmreg        = 1e-6;

% realign options
%=======================================================================
vbm.realign.halfway   = 1; % use halfway registration: 0 - no; 1 - yes
vbm.realign.weight    = 1; % weight registration with inverse std: 0 - no; 1 - yes
vbm.realign.ignore_mat= 0; % ignore exisiting positional information: 0 - no; 1 - yes

% apply deformations options
%=======================================================================
vbm.defs.interp    = 5;  % 5th degree B-spline

% expert options
%=======================================================================
vbm.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
vbm.extopts.finalmask   = 1;    % Final masking: 0 - no; 1 - yes
vbm.extopts.gcut        = 1;    % Skull-stripping with graph-cut: 0 - no; 1 - yes
vbm.extopts.kmeans      = 0;    % segmentation initialization: 0 - new segment; 1 - Kmeans
vbm.extopts.mrf         = 0.15; % MRF weighting
vbm.extopts.sanlm       = 2;    % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading
vbm.extopts.bias_fwhm   = 60;   % FWHM of Kmeans internal bias correction
vbm.extopts.vox         = 1.5;  % voxel size for normalized data
vbm.extopts.bb          = [[-90 -126 -72];[90 90 108]];   % bounding box for normalized data
vbm.extopts.LAS         = 1;    % Local Adaptive Segmentation (VMB12i)
vbm.extopts.pbt.interpV = 1;    % resolution for thickness estimation: 1 - default; 0.5 high res
vbm.extopts.pbt.atlas   = 0;    % use atlas map (thickness only for cortical regions:  0 - no;  1 - yes (default)

% experimental (not yet working)
%=======================================================================
vbm.extopts.mask        = {fullfile(spm('dir'),'toolbox','vbm','submask.nii')}; % mask for subcortical areas + ventricles
vbm.output.surf.dartel  = 0; % WM-surface 
