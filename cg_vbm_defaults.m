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
vbm.opts.ngaus     = [2 2 2 3 4 2];           % Gaussians per class
vbm.opts.affreg    = 'mni';                   % Affine regularisation
vbm.opts.warpreg   = [0 0.001 0.5 0.05 0.2];  % Warping regularisation
vbm.opts.biasreg   = 0.001;                   % Bias regularisation
vbm.opts.biasfwhm  = 60;                      % Bias FWHM
vbm.opts.samp      = 3;                       % Sampling distance


% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

% bias and noise corrected, (localy) intensity normalized
vbm.output.bias.native = 1;
vbm.output.bias.warped = 0;
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

% partitioning atlas maps
vbm.output.atlas.native = 0; 
vbm.output.atlas.warped = 0; 
vbm.output.atlas.dartel = 0; 

% jacobian determinant 0/1 (none/yes)
vbm.output.jacobian.warped = 0;

% deformations
% order is [forward inverse]
vbm.output.warps        = [0 0];

% experimental maps
%=======================================================================
% WARNING: This map describes the changes of the preprocessing and is under development. 

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
vbm.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
vbm.extopts.finalmask   = 1;    % Final masking: 0 - no; 1 - yes
vbm.extopts.gcut        = 2;    % Skull-stripping with graph-cut: 0 - no; 1 - old; 2 - new
vbm.extopts.gcutstr     = 0.5;  % gcut+: strengh of skull-stripping with 0 for softer and wider and 1 for harder and closer (default 0.5)
vbm.extopts.kmeans      = 0;    % segmentation initialization: 0 - new segment; 1 - Kmeans
vbm.extopts.mrf         = 1;    % MRF weighting: 0-1 - manuell setting; 1 - auto
vbm.extopts.sanlm       = 2;    % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading
vbm.extopts.bias_fwhm   = 60;   % FWHM of Kmeans internal bias correction
vbm.extopts.vox         = 1.5;  % voxel size for normalized data
vbm.extopts.bb          = [[-90 -126 -72];[90 90 108]];   % bounding box for normalized data; 
vbm.extopts.LAS         = 1;    % Local Adaptive Segmentation (VMB12i)
vbm.extopts.BVC         = 1;    % Blood Vessel Correction: 
vbm.extopts.INV         = 1;    % Invert PD/T2 images for standard preprocessing  
vbm.extopts.pbtres      = 1;    % resolution for thickness estimation in mm: 1 - normal res (default); 0.5 high res 
vbm.extopts.colormap    = 'BCGWHw'; % {'BCGWHw','BCGWHn'} and matlab colormaps {'jet','gray','bone',...};
vbm.extopts.ROI         = 3;    % write csv-files with ROI data: 1 - subject space; 2 - normalized space; 3 - both (default 1)

% expert options - ROIs
%=======================================================================
% ROI maps from different sources mapped to VBM-space [IXI550]
%  { filename , refinement , tissue , space }
%  filename    = ''                                                     - path to the ROI-file
%  refinement  = ['brain','gm','none']                                  - refinemnt of ROIs in subject space
%  tissue      = {['csf','gm','wm','brain','none','']}                  - tissue classes for volume estimation
%  space       = [0|1|2|3] = [no ROIs|subject space|group space|both];  - space for value estimation
vbm.extopts.atlas       = { ... 
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','l1a.nii')      'none'  {''}              0; ... % VBM atlas with major regions for SBM & ROIs
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','aal.nii')      'gm'    {'gm'}            3; ... 
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','anatomy.nii')  'none'  {'gm','wm'}       3; ...
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','hammers.nii')  'gm'    {'csf','gm','wm'} 3; ...
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','ibsr.nii')     'brain' {'gm'}            3; ...
  fullfile(spm('dir'),'toolbox','vbm12','templates_1.50mm','mori.nii')     'brain' {'gm'}            3; ...
  }; 

% experimental (not yet working)
%=======================================================================
vbm.extopts.mask        = {fullfile(spm('dir'),'toolbox','vbm','submask.nii')}; % mask for subcortical areas + ventricles
vbm.output.surf.dartel  = 0; % WM-surface 
