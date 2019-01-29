function out = cat_vol_series_align(job)
% Longitudinal rigid registration of image series
% FORMAT out = cat_vol_series_align(job)
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
%
% modified version of
% John Ashburner
% spm_series_align.m 5044 2012-11-09 13:40:35Z john
%
% $Id cat_vol_series_align.m $

N = numel(job.data);

if isfield(job.reg,'nonlin')
	tim = job.reg.nonlin.times(:);
	if all(isfinite(tim))
		if numel(tim) ~= N,
				error('Incompatible numbers of times and scans.');
		end
		if any(abs(diff(tim)) > 50),
				error('Time differences should be in years.');
		end;
		wparam0   = job.reg.nonlin.wparam;
		
		midtim = median(tim);
		tim    = tim - midtim;
		w_settings = kron(wparam0,1./(abs(tim)+1/365));
		s_settings = round(3*abs(tim)+2);
	else % use default regularization if tim is set to NAN
		w_settings = job.reg.nonlin.wparam;
		s_settings = 6;
	end
else
  w_settings = [Inf Inf Inf Inf Inf];
  s_settings = Inf;
end

for i=1:N,
    % Make an estimate of the scanner noise
    noise(i,1) = spm_noise_estimate(job.data{i});
    fprintf('Estimated noise sd for "%s" = %g\n', job.data{i}, noise(i,1));
end
prec   = noise.^(-2);

b_settings = [0 0 job.bparam];
Nii = nifti(strvcat(job.data));
ord = [3 3 3 0 0 0];

output = {};
if job.write_avg,  output = [output, {'wavg'}]; end
if job.write_rimg, output = [output, {'wimg'}]; end

if isfield(job.reg,'nonlin')
  if job.reg.nonlin.write_jac, output = [output, {'wjac'} ]; end
  if job.reg.nonlin.write_def, output = [output, {'wdef'} ]; end
end

if ~isfield(job,'use_brainmask')
  use_brainmask = 1;
else
  use_brainmask = job.use_brainmask;
end

if use_brainmask
  fprintf('Use initial brainmask for final rigid registration\n');
end

out    = cat_vol_groupwise_ls(Nii, output, prec, w_settings, b_settings, s_settings, ord, use_brainmask);

return

