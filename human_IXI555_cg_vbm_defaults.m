function human_IXI555_cat_defaults
% Sets the defaults for VBM
% FORMAT cat_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

global cat

% important fields of the animal version
%=======================================================================
% - cat12.opts.tpm 
% - cat12.extopts.darteltpm
% - cat12.extopts.cat12atlas
% - cat12.extopts.brainmask
% - cat12.extopts.bb         > [-inf -inf -inf; inf inf inf] 
% - cat12.extopts.vox        > inf
% - cat12.opts.affreg        > subj
% - cat12.opts.biasreg       > 0.00001
% - cat12.opts.biasfwhm      > 40
% - cat12.opts.samp          > 2 mm
%=======================================================================


% Estimation options
%=======================================================================
cat12.opts.tpm       = {fullfile(spm('dir'),'tpm','TPM.nii')};
cat12.opts.ngaus     = [3 3 2 3 4 2];           % Gaussians per class - 3 GM and 3 WM classes for robustness
cat12.opts.affreg    = 'mni';                   % Affine regularisation - '';'mni';'eastern';'subj';'none';'rigid';
cat12.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation - see Dartel instructions
cat12.opts.biasreg   = 0.0001;                  % Bias regularisation - smaller values for stronger bias fields
cat12.opts.biasfwhm  = 60;                      % Bias FWHM - lower values for stronger bias fields, but look for overfitting in subcortical GM (values <50 mm)
cat12.opts.samp      = 3;                       % Sampling distance - smaller 'better', but slower - maybe usefull for >6 Tesla 

                                              
                                              
                                              
                                              
% Writing options
%=======================================================================

% options:
%   native    0/1     (none/yes)
%   warped    0/1     (none/yes)
%   mod       0/1/2   (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2   (none/rigid/affine)
%   affine    0/1     (none/affine)

% bias and noise corrected, (localy - if LAS>0) intensity normalized
cat12.output.bias.native = 0;
cat12.output.bias.warped = 1;
cat12.output.bias.affine = 0;

% GM tissue maps
cat12.output.GM.native  = 0;
cat12.output.GM.warped  = 0;
cat12.output.GM.mod     = 2;
cat12.output.GM.dartel  = 0;

% WM tissue maps
cat12.output.WM.native  = 0;
cat12.output.WM.warped  = 0;
cat12.output.WM.mod     = 2;
cat12.output.WM.dartel  = 0;
 
% CSF tissue maps
cat12.output.CSF.native = 0;
cat12.output.CSF.warped = 0;
cat12.output.CSF.mod    = 0;
cat12.output.CSF.dartel = 0;

% WMH tissue maps (only for opt.extopts.WMHC==3) - in development
% no modulation available, due to the high spatial variation of WMHs
cat12.output.WMH.native  = 0;
cat12.output.WMH.warped  = 0;
cat12.output.WMH.dartel  = 0;

% label 
% background=0, CSF=1, GM=2, WM=3, WMH=4 (if opt.extropts.WMHC==3)
cat12.output.label.native = 0; 
cat12.output.label.warped = 0;
cat12.output.label.dartel = 0;

% jacobian determinant 0/1 (none/yes)
cat12.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
cat12.output.warps        = [0 0];


% experimental maps
%=======================================================================

% partitioning atlas maps (cat12 atlas)
cat12.output.atlas.native = 0; 
cat12.output.atlas.warped = 0; 
cat12.output.atlas.dartel = 0; 

% preprocessing changes map
% this is the map of the MPC QA measure   
cat12.output.pc.native = 0;
cat12.output.pc.warped = 0;
cat12.output.pc.dartel = 0;

% tissue expectation map
cat12.output.te.native = 0;
cat12.output.te.warped = 0;
cat12.output.te.dartel = 0;


% Longitudinal pipeline
%=======================================================================
% bias correction options
cat12.bias.nits_bias      = 8;
cat12.bias.biasfwhm       = 60;
cat12.bias.biasreg        = 1e-6;
cat12.bias.lmreg          = 1e-6;
% realign options 
cat12.realign.halfway     = 1;      % use halfway registration: 0 - no; 1 - yes
cat12.realign.weight      = 1;      % weight registration with inverse std: 0 - no; 1 - yes
cat12.realign.ignore_mat  = 0;      % ignore exisiting positional information: 0 - no; 1 - yes
% apply deformations options
cat12.defs.interp         = 5;      % 5th degree B-spline


% expert options
%=======================================================================

% Subject species: - 'human';'ape_greater';'ape_lesser';'monkey_oldworld';'monkey_newwold' (in development)
cat12.extopts.species      = 'human';  

% skull-stripping options
cat12.extopts.gcutstr      = 0.5;   % Strengh of skull-stripping:               0 - no gcut; eps - softer and wider; 1 - harder and closer (default = 0.5)
cat12.extopts.cleanupstr   = 0.5;   % Strength of the cleanup process:          0 - no cleanup; eps - soft cleanup; 1 - strong cleanup (default = 0.5) 

% segmentation options
cat12.extopts.LASstr       = 0.5;   % Strength of the local adaption:           0 - no adaption; eps - lower adaption; 1 - strong adaption (default = 0.5)
cat12.extopts.BVCstr       = 0.5;   % Strength of the Blood Vessel Correction:  0 - no correction; eps - low correction; 1 - strong correction (default = 0.5)
cat12.extopts.WMHC         = 1;     % Correction of WM hyperintensities:        0 - no (VBM8); 1 - only for Dartel (default); 
                                  %                                           2 - also for segmentation (corred to WM like SPM); 3 - separate class
cat12.extopts.WMHCstr      = 0.5;   % Strength of WM hyperintensity correction: 0 - no correction; eps - for lower, 1 for stronger corrections (default = 0.5)
cat12.extopts.mrf          = 1;     % MRF weighting:                            0-1 - manuell setting; 1 - auto (default)
cat12.extopts.NCstr        = 0.5;   % Strength of the noise correction:         0 - no noise correction; eps - low correction; 1 - strong corrections (default = 0.5)
cat12.extopts.sanlm        = 3;     % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading (not stable!); 
                                  %                   3 - SANLM with single-threading + ORNLM filter; 4 - SANLM with multi-threading (not stable!) + ORNLM filter;
                                  %                   5 - only ORNLM filter for the final result
cat12.extopts.INV          = 1;     % Invert PD/T2 images for standard preprocessing:  0 - no processing, 1 - try invertation (default), 2 - synthesize T1 image

% resolution options:
cat12.extopts.restype      = 'best';        % resolution handling: 'native','fixed','best'
cat12.extopts.resval       = [1.00 0.10];   % resolution value and its variance for the 'fixed' and 'best' restype

% registration and normalization options 
cat12.extopts.vox          = 1.5;                                % voxel size for normalized data (not yet working):  inf - use Tempate values
cat12.extopts.bb           = [[-90 -126 -72];[90 90 108]];       % bounding box for normalized data (not yet working): inf - use Tempate values
cat12.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_1_IXI555_MNI152.nii')};  % Indicate first Dartel template
cat12.extopts.cat12atlas   = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','cat12.nii')};                     % VBM atlas with major regions for VBM, SBM & ROIs
cat12.extopts.brainmask    = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','brainmask.nii')};                 % brainmask for affine registration
cat12.extopts.T1           = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152.nii')}; % T1 for affine registration

% surface options
cat12.extopts.surface      = 0;     % surface and thickness creation
cat12.extopts.pbtres       = 0.5;   % resolution for thickness estimation in mm: 1 - normal res (default); 0.5 high res 

% visualisation, print and debugging options
cat12.extopts.colormap     = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
cat12.extopts.ROI          = 2;     % write csv-files with ROI data: 1 - subject space; 2 - normalized space; 3 - both (default 2)
cat12.extopts.print        = 1;     % Display and print results
cat12.extopts.verb         = 2;     % Verbose: 1 - default; 2 - details
cat12.extopts.debug        = 1;     % debuging option: 0 - default; 1 - write debuging files 
cat12.extopts.ignoreErrors = 1;     % catching preprocessing errors: 1 - catch errors (default); 0 - stop with error 

% QA options -  NOT IMPLEMENTED - just the idea
%cat12.extopts.QAcleanup    = 1;     % NOT IMPLEMENTED % move images with questionable or bad quality (see QAcleanupth) to subdirectories
%cat12.extopts.QAcleanupth  = [3 5]; % NOT IMPLEMENTED % mark threshold for questionable and bad quality for QAcleanup

cat12.extopts.gui           = 1;     % use GUI 

% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to VBM-space [IXI555]
%  { filename , refinement , tissue }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','gm','none']                                  - refinement of ROIs in subject space
%  tissue      = {['csf','gm','wm','brain','none','']}                  - tissue classes for volume estimation
cat12.extopts.atlas       = { ... 
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','hammers.nii')             'gm'    {'csf','gm','wm'} ; ... % good atlas based on 20 subjects
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','neuromorphometrics.nii')  'gm'    {'csf','gm'};       ... % good atlas based on 35 subjects
 %fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','ibsr.nii')     'brain' {'gm'}            ; ... % less regions than hammer, 18 subjects, low T1 image quality
 %fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','anatomy.nii')  'none'  {'gm','wm'}       ; ... % ROIs requires further work >> use Anatomy toolbox
 %fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','aal.nii')      'gm'    {'gm'}            ; ... % only one subject 
 %fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','mori.nii')     'brain' {'gm'}            ; ... % only one subject, but with WM regions
  }; 


% IDs of the ROIs in the cat12 atlas map (cat12.nii). Do not change this!
cat12.extopts.LAB.CT =  1; % cortex
cat12.extopts.LAB.MB = 13; % MidBrain
cat12.extopts.LAB.BS = 13; % BrainStem
cat12.extopts.LAB.CB =  3; % Cerebellum
cat12.extopts.LAB.ON = 11; % Optical Nerv
cat12.extopts.LAB.BG =  5; % BasalGanglia 
cat12.extopts.LAB.TH =  9; % Hypothalamus 
cat12.extopts.LAB.HC = 19; % Hippocampus 
cat12.extopts.LAB.VT = 15; % Ventricle
cat12.extopts.LAB.NV = 17; % no Ventricle
cat12.extopts.LAB.BV =  7; % Blood Vessels
cat12.extopts.LAB.NB =  0; % no brain 
cat12.extopts.LAB.HD = 21; % head
cat12.extopts.LAB.HI = 23; % WM hyperintensities
