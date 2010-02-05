function cg_ornlm(vargin)
% 
% Optimized Blockwise Non Local Means Denoising Filter
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_ornlm.m 224 2009-12-02 23:39:15Z gaser $

if nargin == 1
	P = [];
	for i=1:numel(vargin.data)
		P = strvcat(P,deblank(vargin.data{i}));
	end
	ornlm_weight = vargin.weight;
else
  P = spm_select(Inf,'image','Select images to filter');
  % get ORNLM weight
  ornlm_weight = spm_input('ORNLM weighting (0.7 to segment)?',1,'e',cg_vbm8_get_defaults('extopts.ornlm'));
end

V = spm_vol(P);
n = size(P,1);

spm_progress_bar('Init',n,'Filtering','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = fileparts(deblank(V(i).fname));
	in = spm_read_vols(V(i));
	
	h = cg_noise_estimation(in);
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
