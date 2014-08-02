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
vbm.opts.tpm       = {fullfile(spm('dir'),'tpm','TPM.nii')};
vbm.opts.ngaus     = [1 1 2 3 4 2];           % Gaussians per class
vbm.opts.affreg    = 'mni';                   % Affine regularisation
vbm.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation
vbm.opts.biasreg   = 0.0001;                  % Bias regularisation - smaller for stronger bias fields
vbm.opts.biasfwhm  = 60;                      % Bias FWHM - lower for stronger bias fieds, but look for overfitting in subcortical GM
vbm.opts.samp      = 3;                       % Sampling distance - smaller 'better' and slower


% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

% bias and noise corrected, (localy) intensity normalized
vbm.output.bias.native = 0;
vbm.output.bias.warped = 1;
vbm.output.bias.affine = 0;

% GM tissue maps
vbm.output.GM.native  = 0;
vbm.output.GM.warped  = 0;
vbm.output.GM.mod     = 2;
vbm.output.GM.dartel  = 0;

% WM tissue maps
vbm.output.WM.native  = 0;
vbm.output.WM.warped  = 0;
vbm.output.WM.mod     = 2;
vbm.output.WM.dartel  = 0;
 
% CSF tissue maps
vbm.output.CSF.native = 0;
vbm.output.CSF.warped = 0;
vbm.output.CSF.mod    = 0;
vbm.output.CSF.dartel = 0;

% label
vbm.output.label.native = 0; 
vbm.output.label.warped = 0;
vbm.output.label.dartel = 0;

% jacobian determinant 0/1 (none/yes)
vbm.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
vbm.output.warps        = [0 0];

% experimental maps
%=======================================================================
% WARNING: This map describes the changes of the preprocessing and is under development. 

% WMH tissue maps
vbm.output.WMH.native  = 0;
vbm.output.WMH.warped  = 0;
vbm.output.WMH.mod     = 0;
vbm.output.WMH.dartel  = 0;

% partitioning atlas maps
vbm.output.atlas.native = 0; 
vbm.output.atlas.warped = 0; 
vbm.output.atlas.dartel = 0; 

% preprocessing changes map
vbm.output.pc.native = 0;
vbm.output.pc.warped = 0;
vbm.output.pc.mod    = 0;
vbm.output.pc.dartel = 0;

% tissue expectation map
vbm.output.te.native = 0;
vbm.output.te.warped = 0;
vbm.output.te.mod    = 0;
vbm.output.te.dartel = 0;

% Extended writing options
%=======================================================================
vbm.extopts.dartelwarp  = 1;  % dartel normalization: 0 - spm default; 1 - yes
vbm.extopts.darteltpm   = ... % Indicate first Dartel template
  {fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','Template_1_IXI550_MNI152.nii')}; 
vbm.extopts.print       = 1; % Display and print results

% bias correction options
%=======================================================================
vbm.bias.nits_bias      = 8;
vbm.bias.biasfwhm       = 60;
vbm.bias.biasreg        = 1e-6;
vbm.bias.lmreg          = 1e-6;

% realign options
%=======================================================================
vbm.realign.halfway     = 1; % use halfway registration: 0 - no; 1 - yes
vbm.realign.weight      = 1; % weight registration with inverse std: 0 - no; 1 - yes
vbm.realign.ignore_mat  = 0; % ignore exisiting positional information: 0 - no; 1 - yes

% apply deformations options
%=======================================================================
vbm.defs.interp         = 5;  % 5th degree B-spline

% expert options
%=======================================================================
% skull-stripping options
vbm.extopts.gcut         = 1;     % Skull-stripping with graph-cut:  0 - no; 1 - yes (default)
vbm.extopts.gcutstr      = 0.5;   % Strengh of skull-stripping with 0 for softer and wider, and 1 for harder and closer (default = 0.5)
vbm.extopts.cleanup      = 3;     % Cleanup of meninges: 0 - no; 1 - light; 2 - thorough; 3 - new improved cleanup (default)
% segmenation options
vbm.extopts.LAS          = 1;     % Local adaptive segmentation (VMB12i):  0 - no adaption; 1 - adaption (default); 2 - adaption & sharpening (this is just a test - do not use!)
vbm.extopts.LASstr       = 0.5;   % Strength of the local adaption:  0 - lower adaption; 1 - strong adaption (default = 0.5)
vbm.extopts.BVCstr       = 0.5;   % Strength of the Blood Vessel Correction (in development):  0 - no correction; 1 - strong correction
vbm.extopts.WMHC         = 1;     % correction for WM hyperintensities (in development):  0 - no; 1 - only for Dartel (default); 2 - also for segmentation 
vbm.extopts.WMHCstr      = 0.5;   % strength of WM hyperintensity correction (in development):  0 for lower, 1 for stronger corrections (default = 0.5)
vbm.extopts.INV          = 1;     % Invert PD/T2 images for standard preprocessing:  0 - no processing, 1 - try invertation (default), 2 - synthesize T1 image
vbm.extopts.mrf          = 1;     % MRF weighting:  0-1 - manuell setting; 1 - auto (default)
vbm.extopts.sanlm        = 3;     % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading (not stable!); 
                                  %                   3 - SANLM with single-threading + ORNLM filter; 2 - SANLM with multi-threading (not stable!) + ORNLM filter; 
vbm.extopts.bias_fwhm    = 60;    % REMOVE ME - because I had to turn it to 0 segmenation due to error % FWHM of Kmeans internal bias correction
vbm.extopts.kmeans       = 0;     % REMOVE ME - because only 0 works correct % segmentation initialization: 0 - new segment; 1 - Kmeans
% normalization options
vbm.extopts.vox          = 1.5;   % voxel size for normalized data
vbm.extopts.bb           = [[-90 -126 -72];[90 90 108]];   % bounding box for normalized data; 
% surface options
vbm.extopts.surface      = 1;     % surface and thickness creation
vbm.extopts.pbtres       = 0.5;   % resolution for thickness estimation in mm: 1 - normal res (default); 0.5 high res 
% visualisation, print and debugging options
vbm.extopts.colormap     = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
vbm.extopts.ROI          = 2;     % write csv-files with ROI data: 1 - subject space; 2 - normalized space; 3 - both (default 2)
vbm.extopts.debug        = 1;     % debuging option: 0 - default; 1 - write debuging files 
vbm.extopts.verb         = 2;     % Verbose: 1 - default; 2 - details
vbm.extopts.ignoreErrors = 1;     % catching preprocessing errors: 1 - catch errors (default); 0 - stop with error 
% QA options
vbm.extopts.QAcleanup    = 1;     % NOT IMPLEMENTED % move images with questionable or bad quality (see QAcleanupth) to subdirectories
vbm.extopts.QAcleanupth  = [3 5]; % NOT IMPLEMENTED % mark threshold for questionable and bad quality for QAcleanup


% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to VBM-space [IXI550]
%  { filename , refinement , tissue }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','gm','none']                                  - refinement of ROIs in subject space
%  tissue      = {['csf','gm','wm','brain','none','']}                  - tissue classes for volume estimation
vbm.extopts.atlas       = { ... 
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','vbm12.nii')      'none'  {''}            ; ... % VBM atlas with major regions for SBM & ROIs
 %fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','aal.nii')      'gm'    {'gm'}            ; ... 
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','anatomy.nii')  'none'  {'gm','wm'}       ; ...
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','hammers.nii')  'gm'    {'csf','gm','wm'} ; ...
 %fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','ibsr.nii')     'brain' {'gm'}            ; ...
 %fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','mori.nii')     'brain' {'gm'}            ; ...
  }; 

% IDs of the ROIs in the vbm12 atlas map (vbm12.nii). 
% Do not change this!
vbm.extopts.LAB.CT =  1; % cortex
vbm.extopts.LAB.MB = 13; % MidBrain
vbm.extopts.LAB.BS = 13; % BrainStem
vbm.extopts.LAB.CB =  3; % Cerebellum
vbm.extopts.LAB.ON = 11; % Optical Nerv
vbm.extopts.LAB.BG =  5; % BasalGanglia 
vbm.extopts.LAB.TH =  9; % Hypothalamus 
vbm.extopts.LAB.HC = 19; % Hippocampus 
vbm.extopts.LAB.VT = 15; % Ventricle
vbm.extopts.LAB.NV = 17; % no Ventricle
vbm.extopts.LAB.BV =  7; % Blood Vessels
vbm.extopts.LAB.NB =  0; % no brain 
vbm.extopts.LAB.HD = 21; % head
vbm.extopts.LAB.HI = 23; % WM hyperintensities
