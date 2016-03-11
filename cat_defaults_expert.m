function cat_defaults_expert
% Sets the defaults for CAT
% FORMAT cat_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

if exist('cat','var'), clear cat; end 
global cat

% Important fields for the processing of animal data
%=======================================================================
% - cat.opts.tpm 
% - cat.extopts.darteltpm
% - cat.extopts.cat12atlas
% - cat.extopts.brainmask
% - cat.extopts.bb         > [-inf -inf -inf; inf inf inf] 
% - cat.extopts.vox        > inf
% - cat.opts.affreg        > subj
% - cat.opts.biasreg       > 0.00001
% - cat.opts.biasfwhm      > 40
% - cat.opts.samp          > 2 mm
%=======================================================================


% Options for inital SPM12 segmentation that is used as starting point for CAT12
%=======================================================================
cat.opts.tpm       = {fullfile(spm('dir'),'tpm','TPM.nii')};
cat.opts.ngaus     = [3 3 2 3 4 2];           % Gaussians per class    - 3 GM and 3 WM classes for robustness
cat.opts.affreg    = 'mni';                   % Affine regularisation  - '';'mni';'eastern';'subj';'none';'rigid';
cat.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation - see Dartel instructions
cat.opts.biasreg   = 0.001;                   % Bias regularisation    - smaller values for stronger bias fields
cat.opts.biasfwhm  = 60;                      % Bias FWHM              - lower values for stronger bias fields, but check for overfitting in subcortical GM (values <50 mm)
cat.opts.samp      = 3;                       % Sampling distance      - smaller 'better', but slower - maybe useful for >= 7 Tesla 

                                              
% Writing options
%=======================================================================

% options:
%   native    0/1     (none/yes)
%   warped    0/1     (none/yes)
%   mod       0/1/2/3 (none/affine+nonlinear/nonlinear only/both)
%   dartel    0/1/2/3 (none/rigid/affine/both)

% save surface and thickness
cat.output.surface     = 1;     % surface and thickness creation

% save ROI values
cat.output.ROI         = 1;     % write xml-file and csv-files with ROI data (0 - no, 1 - yes (default))

% bias and noise corrected, (locally - if LAS>0) intensity normalized
cat.output.bias.native = 1;
cat.output.bias.warped = 0;
cat.output.bias.dartel = 0;

% GM tissue maps
cat.output.GM.native  = 0;
cat.output.GM.warped  = 0;
cat.output.GM.mod     = 1;
cat.output.GM.dartel  = 0;

% WM tissue maps
cat.output.WM.native  = 0;
cat.output.WM.warped  = 0;
cat.output.WM.mod     = 0;
cat.output.WM.dartel  = 0;
 
% CSF tissue maps
cat.output.CSF.native = 0;
cat.output.CSF.warped = 0;
cat.output.CSF.mod    = 0;
cat.output.CSF.dartel = 0;

% WMH tissue maps (only for opt.extopts.WMHC==3) - in development
cat.output.WMH.native  = 0;
cat.output.WMH.warped  = 0;
cat.output.WMH.mod     = 0;
cat.output.WMH.dartel  = 0;

% label 
% background=0, CSF=1, GM=2, WM=3, WMH=4 (if opt.extropts.WMHC==3)
cat.output.label.native = 1; 
cat.output.label.warped = 0;
cat.output.label.dartel = 0;

% jacobian determinant 0/1 (none/yes)
cat.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
cat.output.warps        = [0 0];


% Expert options
%=======================================================================

% skull-stripping options
cat.extopts.gcutstr      = 0.5;   % Strengh of skull-stripping:               0 - no gcut; eps - softer and wider; 1 - harder and closer (default = 0.5)
cat.extopts.cleanupstr   = 0.5;   % Strength of the cleanup process:          0 - no cleanup; eps - soft cleanup; 1 - strong cleanup (default = 0.5) 

% segmentation options
cat.extopts.sanlm        = 2;     % use SANLM filter: 0 - no SANLM; 1 - SANLM; 2 - ISARNLM
cat.extopts.NCstr        = 0.5;   % Strength of the noise correction:         0 - no noise correction; eps - low correction; 1 - strong corrections (default = 0.5)
cat.extopts.LASstr       = 0.5;   % Strength of the local adaption:           0 - no adaption; eps - lower adaption; 1 - strong adaption (default = 0.5)
cat.extopts.BVCstr       = 0.5;   % Strength of the Blood Vessel Correction:  0 - no correction; eps - low correction; 1 - strong correction (default = 0.5)
cat.extopts.WMHC         = 3;     % Correction of WM hyperintensities:        0 - no (VBM8); 1 - only for Dartel (default); 
                                    %                                           2 - also for segmentation (corred to WM like SPM); 3 - separate class
cat.extopts.WMHCstr      = 0.5;   % Strength of WM hyperintensity correction: 0 - no correction; eps - for lower, 1 for stronger corrections (default = 0.5)
cat.extopts.mrf          = 1;     % MRF weighting:                            0 - no MRF; 0 > mrf < 1 - manual setting; 1 - auto (default)
cat.extopts.INV          = 1;     % Invert PD/T2 images for standard preprocessing:  0 - no processing, 1 - try intensity inversion (default), 2 - synthesize T1 image

% resolution options:
cat.extopts.restype      = 'best';        % resolution handling: 'native','fixed','best'
cat.extopts.resval       = [1.00 0.10];   % resolution value and its variance for the 'fixed' and 'best' restype

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
cat.extopts.species      = 'human';  
% Affine PreProcessing (APP) with rough bias correction and brain extraction for special anatomies (nonhuman/neonates) - EXPERIMENTAL  
cat.extopts.APP          = 1;   % 0 - none; 1 - light; 2 - medium; 3 - strong; 4 - heavy
cat.extopts.vox          = 1.5; % voxel size for normalized data (EXPERIMENTAL:  inf - use Tempate values
cat.extopts.bb           = [[-90 -126 -72];[90 90 108]]; % bounding box for normalized data (not yet working): inf - use Tempate values
cat.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_1_IXI555_MNI152.nii')};     % Indicate first Dartel template (Tempalte_1)
%cat.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_0_NKI174_MNI152_GS.nii')};  % Indicate first Shooting template (Template 0)
cat.extopts.cat12atlas   = {fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','cat.nii')};                     % CAT atlas with major regions for VBM, SBM & ROIs
cat.extopts.brainmask    = {fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii')};                                 % Brainmask for affine registration
%cat.extopts.T1           = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};                                        % T1 for affine registration
cat.extopts.T1           = {fullfile(spm('Dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152.nii')};  % T1 for affine registration

% surface options
cat.extopts.pbtres       = 0.5;   % internal resolution for thickness estimation in mm: 
                                    % 1   - normal resolution
                                    % 0.5 - high res (default) 

% visualisation, print and debugging options
cat.extopts.colormap     = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
cat.extopts.verb         = 2;     % Verbose: 1 - default; 2 - details
cat.extopts.debug        = 0;     % debuging option: 0 - default; 1 - write debugging files 
cat.extopts.ignoreErrors = 1;     % catching preprocessing errors: 1 - catch errors (default); 0 - stop with error 
cat.extopts.gui          = 1;     % use GUI 
cat.extopts.expertgui    = 1;     % 0 - common user modus; 1 - expert modus with full GUI; 2 - experimental modus with experimental, unsafe functions!
cat.extopts.subfolders   = 1;     % use subfolders such as mri, surf, report and label to organize your data 


% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to Dartel CAT-space of IXI-template
%  { filename , refinement , tissue }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','tissue','gm','none']                         - not working                  
%  tissue      = {['csf','gm','wm','brain','none']}                     - tissue classes for volume estimation
cat.extopts.atlas       = { ... 
  ... filename                                                                        refinement  defined tissues
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','hammers.nii')             'gm'        {'csf','gm','wm'}; ... % atlas based on 20 subjects
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','neuromorphometrics.nii')  'gm'        {'csf','gm'};      ... % atlas based on 35 subjects
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','lpba40.nii')              'gm'        {'gm'};            ... % atlas based on 40 subjects
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','ibsr.nii')                'brain'     {'csf','gm'};      ... % less regions than hammers, 18 subjects, low T1 image quality
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','aal.nii')                 'gm'        {'gm'};            ... % only one subject 
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','cobra.nii')               'gm'        {'gm'};            ... % only one subject 
  fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','mori.nii')                'brain'     {'gm','wm'};       ... % only one subject, but with WM regions
  %fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','anatomy.nii')             'none'      {'gm','wm'};      ... % ROIs requires further work >> use Anatomy toolbox
  }; 







%=======================================================================
% PRIVATE PARAMETER (NOT FOR GENERAL USE)
%=======================================================================


% further maps
%=======================================================================
% Tissue classes 4-6 to create own TPMs
cat.output.TPMC.native = 0; 
cat.output.TPMC.warped = 0;
cat.output.TPMC.mod    = 0;
cat.output.TPMC.dartel = 0;

% partitioning atlas maps (for evaluation)
cat.output.atlas.native = 0; 
cat.output.atlas.warped = 0; 
cat.output.atlas.dartel = 0; 

% preprocessing changes map
% this is the map that include local changes by preprocessing   
cat.output.pc.native = 0;
cat.output.pc.warped = 0;
cat.output.pc.mod    = 0;
cat.output.pc.dartel = 0;

% tissue expectation map
% this is a map that describes that difference to the TPM
cat.output.te.native = 0;
cat.output.te.warped = 0;
cat.output.te.mod    = 0; % meaningfull?
cat.output.te.dartel = 0;

% IDs of the ROIs in the cat atlas map (cat.nii). Do not change this!
cat.extopts.LAB.NB =  0; % no brain 
cat.extopts.LAB.CT =  1; % cortex
cat.extopts.LAB.CB =  3; % Cerebellum
cat.extopts.LAB.BG =  5; % BasalGanglia 
cat.extopts.LAB.BV =  7; % Blood Vessels
cat.extopts.LAB.TH =  9; % Hypothalamus 
cat.extopts.LAB.ON = 11; % Optical Nerve
cat.extopts.LAB.MB = 13; % MidBrain
cat.extopts.LAB.BS = 13; % BrainStem
cat.extopts.LAB.VT = 15; % Ventricle
cat.extopts.LAB.NV = 17; % no Ventricle
cat.extopts.LAB.HC = 19; % Hippocampus 
cat.extopts.LAB.HD = 21; % Head
cat.extopts.LAB.HI = 23; % WM hyperintensities
