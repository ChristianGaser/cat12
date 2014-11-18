function monkey_oldworld_cg_vbm_defaults
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

% important fields of the animal version
%=======================================================================
% - vbm.opts.tpm 
% - vbm.extopts.darteltpm
% - vbm.extopts.vbm12atlas
% - vbm.extopts.brainmask
% - vbm.extopts.bb         > [-inf -inf -inf; inf inf inf] 
% - vbm.extopts.vox        > inf
% - vbm.opts.affreg        > subj
% - vbm.opts.biasreg       > 0.00001
% - vbm.opts.biasfwhm      > 40
% - vbm.opts.samp          > 2 mm
%=======================================================================


% Estimation options
%=======================================================================
vbm.opts.tpm       = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_TPM.nii')};
vbm.opts.ngaus     = [3 3 2 3 4 2];           % Gaussians per class - 3 GM and 3 WM classes for robustness
vbm.opts.affreg    = 'subj';                  % Affine regularisation - '';'mni';'eastern';'subj';'none';'rigid';
vbm.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation - see Dartel instructions
vbm.opts.biasreg   = 0.00001;                 % Bias regularisation - smaller values for stronger bias fields
vbm.opts.biasfwhm  = 40;                      % Bias FWHM - lower values for stronger bias fields, but look for overfitting in subcortical GM (values <50 mm)
vbm.opts.samp      = 2;                       % Sampling distance - smaller 'better', but slower - maybe usefull for >6 Tesla 

                                              
                                              
                                              
                                              
% Writing options
%=======================================================================

% options:
%   native    0/1     (none/yes)
%   warped    0/1     (none/yes)
%   mod       0/1/2   (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2   (none/rigid/affine)
%   affine    0/1     (none/affine)

% bias and noise corrected, (localy - if LAS>0) intensity normalized
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

% WMH tissue maps (only for opt.extopts.WMHC==3) - in development
% no modulation available, due to the high spatial variation of WMHs
vbm.output.WMH.native  = 0;
vbm.output.WMH.warped  = 0;
vbm.output.WMH.dartel  = 0;

% label 
% background=0, CSF=1, GM=2, WM=3, WMH=4 (if opt.extropts.WMHC==3)
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

% partitioning atlas maps (vbm12 atlas)
vbm.output.atlas.native = 0; 
vbm.output.atlas.warped = 0; 
vbm.output.atlas.dartel = 0; 

% preprocessing changes map
% this is the map of the MPC QA measure   
vbm.output.pc.native = 0;
vbm.output.pc.warped = 0;
vbm.output.pc.dartel = 0;

% tissue expectation map
vbm.output.te.native = 0;
vbm.output.te.warped = 0;
vbm.output.te.dartel = 0;


% Longitudinal pipeline
%=======================================================================
% bias correction options
vbm.bias.nits_bias      = 8;
vbm.bias.biasfwhm       = 60;
vbm.bias.biasreg        = 1e-6;
vbm.bias.lmreg          = 1e-6;
% realign options 
vbm.realign.halfway     = 1;      % use halfway registration: 0 - no; 1 - yes
vbm.realign.weight      = 1;      % weight registration with inverse std: 0 - no; 1 - yes
vbm.realign.ignore_mat  = 0;      % ignore exisiting positional information: 0 - no; 1 - yes
% apply deformations options
vbm.defs.interp         = 5;      % 5th degree B-spline


% expert options
%=======================================================================

% Subject species: - 'human';'ape_greater';'ape_lesser';'monkey_oldworld';'monkey_newwold' (in development)
vbm.extopts.species      = 'monkey_oldworld';  

% skull-stripping options
vbm.extopts.gcutstr      = 0.5;   % Strengh of skull-stripping:               0 - no gcut; eps - softer and wider; 1 - harder and closer (default = 0.5)
vbm.extopts.cleanupstr   = 0.5;   % Strength of the cleanup process:          0 - no cleanup; eps - soft cleanup; 1 - strong cleanup (default = 0.5) 

% segmentation options
vbm.extopts.LASstr       = 0.5;   % Strength of the local adaption:           0 - no adaption; eps - lower adaption; 1 - strong adaption (default = 0.5)
vbm.extopts.BVCstr       = 0.5;   % Strength of the Blood Vessel Correction:  0 - no correction; eps - low correction; 1 - strong correction (default = 0.5)
vbm.extopts.WMHC         = 1;     % Correction of WM hyperintensities:        0 - no (VBM8); 1 - only for Dartel (default); 
                                  %                                           2 - also for segmentation (corred to WM like SPM); 3 - separate class
vbm.extopts.WMHCstr      = 0.5;   % Strength of WM hyperintensity correction: 0 - no correction; eps - for lower, 1 for stronger corrections (default = 0.5)
vbm.extopts.mrf          = 1;     % MRF weighting:                            0-1 - manuell setting; 1 - auto (default)
vbm.extopts.NCstr        = 0.5;   % Strength of the noise correction:         0 - no noise correction; eps - low correction; 1 - strong corrections (default = 0.5)
vbm.extopts.sanlm        = 3;     % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading (not stable!); 
                                  %                   3 - SANLM with single-threading + ORNLM filter; 4 - SANLM with multi-threading (not stable!) + ORNLM filter;
                                  %                   5 - only ORNLM filter for the final result
vbm.extopts.INV          = 1;     % Invert PD/T2 images for standard preprocessing:  0 - no processing, 1 - try invertation (default), 2 - synthesize T1 image

% resolution options:
vbm.extopts.restype      = 'best';        % resolution handling: 'native','fixed','best'
vbm.extopts.resval       = [1.00 0.10];   % resolution value and its variance for the 'fixed' and 'best' restype

% registration and normalization options 
vbm.extopts.vox          = inf;                                % voxel size for normalized data (not yet working):  inf - use Tempate values
vbm.extopts.bb           = [[-inf -inf -inf];[inf inf inf]];   % bounding box for normalized data (not yet working): inf - use Tempate values
vbm.extopts.dartelwarp   = 1;                                  % dartel normalization: 0 - spm default; 1 - yes
vbm.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_Template_1.nii')}; % Indicate first Dartel template
vbm.extopts.vbm12atlas   = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_vbm12.nii')};      % VBM atlas with major regions for VBM, SBM & ROIs
vbm.extopts.brainmask    = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_brainmask.nii')};  % brainmask for affine registration
vbm.extopts.T1           = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_T1.nii')};         % T1 for affine registration

% surface options
vbm.extopts.surface      = 0;     % surface and thickness creation
vbm.extopts.pbtres       = 0.5;   % resolution for thickness estimation in mm: 1 - normal res (default); 0.5 high res 

% visualisation, print and debugging options
vbm.extopts.colormap     = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
vbm.extopts.ROI          = 2;     % write csv-files with ROI data: 1 - subject space; 2 - normalized space; 3 - both (default 2)
vbm.extopts.print        = 1;     % Display and print results
vbm.extopts.verb         = 2;     % Verbose: 1 - default; 2 - details
vbm.extopts.debug        = 1;     % debuging option: 0 - default; 1 - write debuging files 
vbm.extopts.ignoreErrors = 1;     % catching preprocessing errors: 1 - catch errors (default); 0 - stop with error 

% QA options -  NOT IMPLEMENTED - just the idea
%vbm.extopts.QAcleanup    = 1;     % NOT IMPLEMENTED % move images with questionable or bad quality (see QAcleanupth) to subdirectories
%vbm.extopts.QAcleanupth  = [3 5]; % NOT IMPLEMENTED % mark threshold for questionable and bad quality for QAcleanup

vbm.extopts.gui           = 1;     % use GUI 

% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to VBM-space [IXI555]
%  { filename , refinement , tissue }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','gm','none']                                  - refinement of ROIs in subject space
%  tissue      = {['csf','gm','wm','brain','none','']}                  - tissue classes for volume estimation
vbm.extopts.atlas       = { ... 
  fullfile(spm('dir'),'toolbox','vbm12','templates_animals','monkey_oldworld_dartel_atlas_inia19NeuroMaps.nii') 'gm'    {'csf','gm','wm'} ; ... % good atlas based on 19 subjects
  }; 


% IDs of the ROIs in the vbm12 atlas map (vbm12.nii). Do not change this!
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
