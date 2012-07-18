function cg_vbm12_defaults
% Sets the defaults for VBM
% FORMAT cg_vbm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

global vbm12

% Estimation options
%=======================================================================
vbm12.opts.tpm       = {fullfile(spm('dir'),'toolbox','Seg','TPM.nii')}; % TPM.nii
vbm12.opts.ngaus     = [2 2 2 3 4 2];  % Gaussians per class
vbm12.opts.affreg    = 'mni';    % Affine regularisation
vbm12.opts.warpreg   = 4;      % Warping regularisation
vbm12.opts.biasreg   = 0.0001; % Bias regularisation
vbm12.opts.biasfwhm  = 60;   % Bias FWHM
vbm12.opts.samp      = 3;      % Sampling distance

% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

vbm12.output.bias.native  = 1;
vbm12.output.bias.warped  = 0;
vbm12.output.bias.affine  = 0;

vbm12.output.label.native = 1;
vbm12.output.label.warped = 0;
vbm12.output.label.dartel = 0;

% order is [native normalised modulated dartel]
vbm12.output.GM.native = 1;  % GM
vbm12.output.GM.warped = 0;  % GM
vbm12.output.GM.mod    = 0;  % GM
vbm12.output.GM.dartel = 0;  % GM

vbm12.output.WM.native = 0;  % WM
vbm12.output.WM.warped = 0;  % WM
vbm12.output.WM.mod    = 0;  % WM
vbm12.output.WM.dartel = 0;  % WM

vbm12.output.CSF.native = 0; % CSF
vbm12.output.CSF.warped = 0; % CSF
vbm12.output.CSF.mod    = 0; % CSF
vbm12.output.CSF.dartel = 0; % CSF

% jacobian determinant 0/1 (none/yes)
vbm12.output.jacobian.warped = 0;

% order is [forward inverse]
vbm12.output.warps = [0 0];

% Extended writing options
%=======================================================================
vbm12.extopts.dartelwarp  = 0; % dartel normalization: 0 - spm default; 1 - yes
vbm12.extopts.darteltpm   = {fullfile(spm('dir'),'toolbox','vbm12','Template_1_IXI550_MNI152.nii')}; % Indicate first Dartel template
vbm12.extopts.print       = 1; % Display and print results

% bias correction options
%=======================================================================
vbm12.bias.nits_bias    = 8;
vbm12.bias.biasfwhm     = 60;
vbm12.bias.biasreg      = 1e-6;
vbm12.bias.lmreg        = 1e-6;

% realign options
%=======================================================================
vbm12.realign.halfway   = 1; % use halfway registration: 0 - no; 1 - yes
vbm12.realign.weight    = 1; % weight registration with inverse std: 0 - no; 1 - yes
vbm12.realign.ignore_mat= 0; % ignore exisiting positional information: 0 - no; 1 - yes

% apply deformations options
%=======================================================================
vbm12.defs.interp    = 5;  % 5th degree B-spline

% expert options
%=======================================================================
vbm12.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
vbm12.extopts.finalmask   = 1;    % Final masking: 0 - no; 1 - yes
vbm12.extopts.gcut        = 1;    % Skull-stripping with graph-cut: 0 - no; 1 - yes
vbm12.extopts.kmeans      = 1;    % segmentation initialization: 0 - new segment; 1 - Kmeans
vbm12.extopts.mrf         = 0.15; % MRF weighting
vbm12.extopts.sanlm       = 2;    % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading
vbm12.extopts.bias_fwhm   = 60;   % FWHM of Kmeans internal bias correction
vbm12.extopts.histeq_deep = 0;    % weighting of local histogram equalization: 0 - no; 1 - full weighting (not recommended)
vbm12.extopts.histeq_mask = {fullfile(spm('dir'),'toolbox','vbm12','histeq_mask.nii')};

% experimental (not yet working)
%=======================================================================
vbm12.extopts.mask      = {fullfile(spm('dir'),'toolbox','vbm12','submask.nii')}; % mask for subcortical areas + ventricles
vbm12.output.surf.dartel= 0; % WM-surface 
