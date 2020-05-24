function d = cat_long_shoot_defaults
%cat_long_shoot_defaults: Defaults settings for intra subject deformations. 
% Shoooting settings for tiny deformations between time points in the CAT 
% longidtudinal preprocessing.  These settings should correct for tiny 
% disaplacements in aging.
% 
% Defaults file
% FORMAT d = spm_shoot_defaults
% This file contains settings that are intended to be customised
% according to taste.  Some of them will influence the speed/accuracy
% tradeoff, whereas others are various regularisation settings
% (registration and template blurring)...

%_______________________________________________________________________
% Copyright (C) Wellcome Trust Centre for Neuroimaging (2009)

% John Ashburner
% $Id: spm_shoot_defaults.m 6489 2015-06-26 11:46:29Z john $


%_______________________________________________________________________
% The following settings are intended to be customised according to
% taste.  Some of them will influence the speed/accuracy tradeoff,
% whereas others are various regularisation settings (registration
% and template blurring)...

%
% shooting defaults for cross-sectional processing (just for comparison)
if 0
  d.tname   = 'Template'; % Base file name for templates
  d.issym   = false;      % Use a symmetric template?

  d.cyc_its = [2 2];      % No. multigrid cycles and iterations

  nits      = 24;         % No. iterations of Gauss-Newton

  % Schedule for coarse to fine
  lam     = 0.5;          % Decay of coarse to fine schedule
  inter   = 32;           % Scaling of parameters at first iteration
  d.sched = (inter-1)*exp(-lam*((1:(nits+1))-1))+1;
  d.sched = d.sched/d.sched(end);

  maxoil    = 8;                          % Maximum number of time steps for integration
  d.eul_its = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps

  d.rparam  = [1e-4 0.001 0.2 0.05 0.2];  % Regularisation parameters for deformation
  d.sparam  = [0.0001 0.08 0.8];          % Regularisation parameters for blurring
  d.smits   = 16;                         % No. smoothing iterations
  %d.smits   = 0; d.sparam = [];
  d.scale   = 0.8;                        % Fraction of Gauss-Newton update step to use

  d.bs_args = [2 2 2  1 1 1];             % B-spline settings for interpolation
else
  %% RD202005: definition for individual aging deformation 
  
  d.tname     = 'Template_FNAME'; % Base file name for templates ... FNAME is replaced by filename
  d.issym     = false;            % Use a symmetric template?
  d.cyc_its   = [1 1];            % No. multigrid cycles and iterations
  nits        = 8;                % No. iterations of Gauss-Newton (4-8)
  % Schedule for coarse to fine
  lam         = 0.25;             % Decay of coarse to fine schedule
  inter       = 64;               % Scaling of parameters at first iteration (32)
  d.sched     = (inter-1)*exp(-lam*((1:(nits+1))-1))+1;
  d.sched     = d.sched/d.sched(end) * 32;  % 32 as gener offset to have low-frequency deformations
  maxoil      = 8;                          % Maximum number of time steps for integration
  d.eul_its   = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps
  d.eul_its   = flip(d.eul_its);
  d.rparam    = [1e-4 0.001 0.2 0.05 0.2];  % Regularisation parameters for deformation (just the defaults)
  d.sparam    = [0 0 0.0000001];  % Regularisation parameters for blurring - strongly reduced as far as it create incorrect results in the individual cases
  d.smits     = 1;                % Minimum smoothing iterations (turn it off will produce an error at the end of the Shooting iteration)
  d.scale     = 0.8;              % Fraction of Gauss-Newton update step to use
  d.bs_args   = [2 2 2  1 1 1];   % B-spline settings for interpolation
end
