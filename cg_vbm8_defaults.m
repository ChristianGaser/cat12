function cg_vbm8_defaults
% Sets the defaults for VBM
% FORMAT cg_vbm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id$

global defaults

% Estimation options
%=======================================================================
defaults.vbm8.opts.ngaus     = [2 2 2 3 4 2];	% Gaussians per class
defaults.vbm8.opts.affreg    = 'mni';		% Affine regularisation
defaults.vbm8.opts.warpreg   = 4;			% Warping regularisation
defaults.vbm8.opts.affmethod = 1;			% Affine registration method
defaults.vbm8.opts.biasreg   = 0.0001;	% Bias regularisation
defaults.vbm8.opts.biasfwhm  = 60;		% Bias FWHM
defaults.vbm8.opts.samp      = 3;			% Sampling distance

% Writing options
%=======================================================================

% segmentations:
%   native		0/1   (none/yes)
%   warped		0/1   (none/yes)
%   modulated	0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

defaults.vbm8.output.bias.native  = 0;
defaults.vbm8.output.bias.warped  = 1;
defaults.vbm8.output.bias.affine  = 0;

defaults.vbm8.output.label.native = 0;
defaults.vbm8.output.label.warped = 0;
defaults.vbm8.output.label.dartel = 0;

% order is [native normalised modulated dartel]
defaults.vbm8.output.GM.native = 0;	% GM
defaults.vbm8.output.GM.warped = 1;	% GM
defaults.vbm8.output.GM.mod    = 2;	% GM
defaults.vbm8.output.GM.dartel = 0;	% GM

defaults.vbm8.output.WM.native = 0;	% WM
defaults.vbm8.output.WM.warped = 1;	% WM
defaults.vbm8.output.WM.mod    = 2;	% WM
defaults.vbm8.output.WM.dartel = 0;	% WM

defaults.vbm8.output.CSF.native = 0;	% CSF
defaults.vbm8.output.CSF.warped = 0;	% CSF
defaults.vbm8.output.CSF.mod    = 0;	% CSF
defaults.vbm8.output.CSF.dartel = 0;	% CSF

% jacobian determinant 0/1 (none/yes)
defaults.vbm8.output.jacobian.warped = 0;

% order is [inverse forward]
defaults.vbm8.output.warps = [0 0];

% Extended writing options
%=======================================================================
defaults.vbm8.extopts.dartelwarp  = 1; % dartel normalization: 0 - spm default; 1 - yes
defaults.vbm8.extopts.vox          = NaN;	% Voxel size to write
defaults.vbm8.extopts.bb           = [[-78 78]' [-112 76]' [-70 85]'];	% bounding box
defaults.vbm8.extopts.print        = 1;	% Display and print results
