function out = cat_vol_tfce(in, dh, E, H, calc_neg, invert)
% Apply threshold-free cluster enhancement (TFCE) and scale the image to
% have the same intensity distribution as the input to strengthen large 
% connected regions for clean up or removing vessels in WM
%
% FORMAT out = cat_vol_tfce(in, dh, E, H, calc_neg, invert)
% in        - input image 
% dh        - step size (e.g. dh = max(abs(in))/100)
%             leave empty or set to Inf or NaN to automatically estimate dh
% E         - TFCE parameter for extent (default 2)
% H         - TFCE parameter for height  (default 1)
% calc_neg  - also calc neg. TFCE values (not by default)
% invert    - apply TFCE to inverted image to strengthen fine structures
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% call interactive mode if no arguments are given
if ~nargin
  P = spm_select(Inf,'image','Select images for TFCE filtering');
  EH = spm_input('Weighting for extent and height',1,'r',[2 1],[2,1]);
  invert = spm_input('Invert to strengthen fine structures?','+1','yes|no',[1 0],2);
  V = spm_vol(P);
  for i = 1:numel(V)
    [pth,nam,ext] = spm_fileparts(V(i).fname);
    in = spm_read_vols(V(i));
    out = cat_vol_tfce(in, Inf, EH(1), EH(2), 0, invert);
    if invert
      V(i).fname = fullfile(pth,sprintf('tfce_invE%gH%g_%s%s',EH(1),EH(2),nam,ext));
    else
      V(i).fname = fullfile(pth,sprintf('tfceE%gH%g_%s%s',EH(1),EH(2),nam,ext));
    end
    spm_write_vol(V(i),out);
  end
  clear out
  return
end

% check minimum arguments
if nargin < 4
  error('You have to at least define "in, dh, E, H" as input parameters.');
end

% use 100 steps from 0..max
if isempty(dh) || ~isfinite(dh)
  dh = max(abs(in(:)))/100;
end

% set defaults for E and H
if nargin < 3
  E = 2;
end
if nargin < 4
  H = 1;
end

% also calc neg. TFCE values
if nargin < 5
  calc_neg = 0;
end

% use multithreading
if nargin < 6
  invert = 0;
end

if invert
  mx = max(abs(in(:)));
  in = mx - in;
end

% use multithreading
single_threaded = 0;

out = tfceMex_pthread(in, dh, E, H, calc_neg, single_threaded);

% estimate the scaling factor from the threshold values and scale the output 
% to have the same intensity distribution as the input
ind = in > 0.01*max(abs(in(:)));
scl = median(in(ind)./out(ind));
out = scl*out;

if invert
  out = mx - out;
end
