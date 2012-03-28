function cg_ornlm(vargin)
% 
% Optimized Blockwise Non Local Means Denoising Filter
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

if nargin == 1
	P = char(vargin.data);
	ornlm_weight = vargin.weight;
else
  P = spm_select(Inf,'image','Select images to filter');
  % get ORNLM weight
  ornlm_weight = spm_input('ORNLM weighting (0.7 to segment)?',1,'e',0.7);
end

V = spm_vol(P);
n = size(P,1);

spm_progress_bar('Init',n,'Filtering','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = fileparts(deblank(V(i).fname));
	in = spm_read_vols(V(i));
	
  try
    h = cg_noise_estimation(in);
  catch
    % use additionally reduced size to prevent memory issues
    fprintf('\nUse reduced images size for noise estimation.');
    h = cg_noise_estimation(in(1:2:V(i).dim(1),1:2:V(i).dim(2)/2,1:2:V(i).dim(3)/2));
  end
  fprintf('Noise estimate for %s: %3.2f\n',nm,h);

  % ORNLM weighting
  h = ornlm_weight*h;
  
  out = ornlmMex(in,3,1,h);
  V(i).fname = fullfile(pth,['ornlm_' nm xt vr]);
  V(i).descrip = sprintf('ORNLM filtered h=%3.2f',h);
  spm_write_vol(V(i), out);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

return
