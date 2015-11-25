function cat_defaults
% Sets the defaults for CAT
% FORMAT cat_defaults_monkey_newworld
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

if exist('cat12','var'), clear cat12; end 
global cat12

% Important fields for the processing of animal data
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


% Options for inital SPM12 segmentation that is used as starting point for CAT12
%=======================================================================
cat12.opts.tpm       = {fullfile(spm('dir'),'toolbox','vbm12','templates_animals','ape_greater_TPM.nii')};
cat12.opts.ngaus     = [3 3 2 3 4 2];           % Gaussians per class    - 3 GM and 3 WM classes for robustness
cat12.opts.affreg    = 'none';                  % Affine regularisation  - '';'mni';'eastern';'subj';'none';'rigid';
cat12.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation - see Dartel instructions
cat12.opts.biasreg   = 0.001;                   % Bias regularisation    - smaller values for stronger bias fields
cat12.opts.biasfwhm  = 60;                      % Bias FWHM              - lower values for stronger bias fields, but check for overfitting in subcortical GM (values <50 mm)
cat12.opts.samp      = 3;                       % Sampling distance      - smaller 'better', but slower - maybe useful for >= 7 Tesla 

                                              
% Writing options
%=======================================================================

% options:
%   native    0/1     (none/yes)
%   warped    0/1     (none/yes)
%   mod       0/1/2/3 (none/affine+nonlinear/nonlinear only/both)
%   dartel    0/1/2/3 (none/rigid/affine/both)

% save surface and thickness
cat12.output.surface     = 0;     % surface and thickness creation

% save ROI values
cat12.output.ROI         = 0;     % write csv-files with ROI data: 1 - subject space; 2 - normalized space; 3 - both (default 2)

% bias and noise corrected, (locally - if LAS>0) intensity normalized
cat12.output.bias.native = 1;
cat12.output.bias.warped = 0;
cat12.output.bias.dartel = 0;

% GM tissue maps
cat12.output.GM.native  = 0;
cat12.output.GM.warped  = 0;
cat12.output.GM.mod     = 0;
cat12.output.GM.dartel  = 3;

% WM tissue maps
cat12.output.WM.native  = 0;
cat12.output.WM.warped  = 0;
cat12.output.WM.mod     = 0;
cat12.output.WM.dartel  = 3;
 
% CSF tissue maps
cat12.output.CSF.native = 0;
cat12.output.CSF.warped = 0;
cat12.output.CSF.mod    = 0;
cat12.output.CSF.dartel = 3;

% WMH tissue maps (only for opt.extopts.WMHC==3) - in development
cat12.output.WMH.native  = 0;
cat12.output.WMH.warped  = 0;
cat12.output.WMH.mod     = 0;
cat12.output.WMH.dartel  = 0;

% label 
% background=0, CSF=1, GM=2, WM=3, WMH=4 (if opt.extropts.WMHC==3)
cat12.output.label.native = 1; 
cat12.output.label.warped = 0;
cat12.output.label.dartel = 0;

% jacobian determinant 0/1 (none/yes)
cat12.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
cat12.output.warps        = [0 0];


% Expert options
%=======================================================================

% skull-stripping options
cat12.extopts.gcutstr      = 0.5;   % Strengh of skull-stripping:               0 - no gcut; eps - softer and wider; 1 - harder and closer (default = 0.5)
cat12.extopts.cleanupstr   = 0.5;   % Strength of the cleanup process:          0 - no cleanup; eps - soft cleanup; 1 - strong cleanup (default = 0.5) 

% segmentation options
cat12.extopts.sanlm        = 3;     % use SANLM filter: 0 - no SANLM; 1 - SANLM; 3 - SANLM + ORNLM filter; 5 - only ORNLM filter for the final result
cat12.extopts.NCstr        = 0.5;   % Strength of the noise correction:         0 - no noise correction; eps - low correction; 1 - strong corrections (default = 0.5)
cat12.extopts.LASstr       = 0.5;   % Strength of the local adaption:           0 - no adaption; eps - lower adaption; 1 - strong adaption (default = 0.5)
cat12.extopts.BVCstr       = 0.5;   % Strength of the Blood Vessel Correction:  0 - no correction; eps - low correction; 1 - strong correction (default = 0.5)
cat12.extopts.WMHC         = 1;     % Correction of WM hyperintensities:        0 - no (VBM8); 1 - only for Dartel (default); 
                                  %                                           2 - also for segmentation (corred to WM like SPM); 3 - separate class
cat12.extopts.WMHCstr      = 0.5;   % Strength of WM hyperintensity correction: 0 - no correction; eps - for lower, 1 for stronger corrections (default = 0.5)
cat12.extopts.mrf          = 1;     % MRF weighting:                            0 - no MRF; 0 > mrf < 1 - manual setting; 1 - auto (default)
cat12.extopts.INV          = 1;     % Invert PD/T2 images for standard preprocessing:  0 - no processing, 1 - try intensity inversion (default), 2 - synthesize T1 image

% resolution options:
cat12.extopts.restype      = 'best';        % resolution handling: 'native','fixed','best'
cat12.extopts.resval       = [0.70 0.10];   % resolution value and its variance for the 'fixed' and 'best' restype

%{
native:
    Preprocessing with native resolution.
    In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). 

    Examples:
      native resolution       internal resolution 
       0.95 0.95 1.05     >     0.95 0.95 1.05
       0.45 0.45 1.70     >     0.45 0.45 1.50 (if voxel size for normalized images is 1.50 mm)

best:
    Preprocessing with the best (minimal) voxel dimension of the native image.'
    The first parameters defines the lowest spatial resolution for every dimension, while the second is used to avoid tiny interpolations for almost correct resolutions.
    In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). 

    Examples:
      Parameters    native resolution       internal resolution
      [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.00 1.00
      [1.00 0.10]    0.45 0.45 1.50     >     0.45 0.45 1.00
      [0.75 0.10]    0.45 0.45 1.50     >     0.45 0.45 0.75  
      [0.75 0.10]    0.45 0.45 0.80     >     0.45 0.45 0.80  
      [0.00 0.10]    0.45 0.45 1.50     >     0.45 0.45 0.45  

fix:
    This options prefers an isotropic voxel size that is controled by the first parameters.  
    The second parameter is used to avoid tiny interpolations for almost correct resolutions. 
    In order to avoid interpolation artifacts in the Dartel output the lowest spatial resolution is always limited to the voxel size of the normalized images (default 1.5mm). 
    There is no upper limit, but we recommend to avoid unnecessary interpolation.

    Examples: 
      Parameters     native resolution       internal resolution
      [1.00 0.10]     0.45 0.45 1.70     >     1.00 1.00 1.00
      [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00
      [1.00 0.02]     0.95 1.05 1.25     >     1.00 1.00 1.00
      [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00
      [0.75 0.10]     0.75 0.95 1.25     >     0.75 0.75 0.75

%}


% registration and normalization options 
% Subject species: - 'human';'ape_greater';'ape_lesser';'monkey_oldworld';'monkey_newwold' (in development)
cat12.extopts.species      = 'monkey_newworld';  
% Affine PreProcessing (APP) with rough bias correction and brain extraction for special anatomies (nonhuman/neonates) - EXPERIMENTAL  
cat12.extopts.APP          = 3;   % 0 - none (default); 1 - APP with init. affreg; 2 - APP without init. affreg (standard in non human); 
cat12.extopts.vox          = 0.70; % voxel size for normalized data (EXPERIMENTAL:  inf - use Tempate values
cat12.extopts.bb           = [[-inf -inf -inf];[inf inf inf]];; % bounding box for normalized data (not yet working): inf - use Tempate values
cat12.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','cat12','templates_animals','monkey_newworld_Template_1.nii')}; % Indicate first Dartel template
cat12.extopts.cat12atlas   = {fullfile(spm('dir'),'toolbox','cat12','templates_animals','monkey_newworld_vbm12.nii')};      % VBM atlas with major regions for VBM, SBM & ROIs
cat12.extopts.brainmask    = {fullfile(spm('dir'),'toolbox','cat12','templates_animals','monkey_newworld_brainmask.nii')};  % brainmask for affine registration
cat12.extopts.T1           = {fullfile(spm('dir'),'toolbox','cat12','templates_animals','monkey_newworld_T1.nii')};         % T1 for affine registration

% surface options
cat12.extopts.pbtres       = 0.35;   % internal resolution for thickness estimation in mm: 
                                     % 1   - normal resolution
                                     % 0.5 - high res (default) 

% visualisation, print and debugging options
cat12.extopts.colormap     = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
cat12.extopts.print        = 1;     % Display and print results
cat12.extopts.verb         = 2;     % Verbose: 1 - default; 2 - details
cat12.extopts.debug        = 0;     % debuging option: 0 - default; 1 - write debugging files 
cat12.extopts.ignoreErrors = 1;     % catching preprocessing errors: 1 - catch errors (default); 0 - stop with error 
cat12.extopts.gui          = 1;     % use GUI 
cat12.extopts.expertgui    = 2;     % 0 - common user modus; 1 - expert modus with full GUI; 2 - experimental modus with experimental, unsafe functions!


% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to Dartel VBM-space of IXI-template
%  { filename , refinement , tissue }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','gm','none']                                  - refinement of ROIs in subject space
%  tissue      = {['csf','gm','wm','brain','none','']}                  - tissue classes for volume estimation
cat12.extopts.atlas       = { ... 
  }; 








%=======================================================================
% PRIVATE PARAMETER (NOT FOR GENERAL USE)
%=======================================================================


% further maps
%=======================================================================
% Tissue classes 4-6 to create own TPMs
cat12.output.TPMC.native = 0; 
cat12.output.TPMC.warped = 0;
cat12.output.TPMC.mod    = 0;
cat12.output.TPMC.dartel = 0;

% partitioning atlas maps (cat12 atlas)
cat12.output.atlas.native = 0; 
cat12.output.atlas.warped = 0; 
cat12.output.atlas.dartel = 0; 

% preprocessing changes map
% this is the map that include local changes by preprocessing   
cat12.output.pc.native = 0;
cat12.output.pc.warped = 0;
cat12.output.pc.mod    = 0;
cat12.output.pc.dartel = 0;

% tissue expectation map
% this is a map that describes that difference to the TPM
cat12.output.te.native = 0;
cat12.output.te.warped = 0;
cat12.output.te.mod    = 0; % meaningfull?
cat12.output.te.dartel = 0;

% IDs of the ROIs in the cat12 atlas map (cat12.nii). Do not change this!
cat12.extopts.LAB.NB =  0; % no brain 
cat12.extopts.LAB.CT =  1; % cortex
cat12.extopts.LAB.CB =  3; % Cerebellum
cat12.extopts.LAB.BG =  5; % BasalGanglia 
cat12.extopts.LAB.BV =  7; % Blood Vessels
cat12.extopts.LAB.TH =  9; % Hypothalamus 
cat12.extopts.LAB.ON = 11; % Optical Nerve
cat12.extopts.LAB.MB = 13; % MidBrain
cat12.extopts.LAB.BS = 13; % BrainStem
cat12.extopts.LAB.VT = 15; % Ventricle
cat12.extopts.LAB.NV = 17; % no Ventricle
cat12.extopts.LAB.HC = 19; % Hippocampus 
cat12.extopts.LAB.HD = 21; % Head
cat12.extopts.LAB.HI = 23; % WM hyperintensities
