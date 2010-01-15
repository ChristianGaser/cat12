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
else
  P = spm_select(Inf,'image','Select images to filter');
end

V = spm_vol(P);
n = size(P,1);

spm_progress_bar('Init',n,'Filtering','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = fileparts(deblank(V(i).fname));
	in = spm_read_vols(V(i));
	h = rician_noise_estimation(in);

	if h>0
	  fprintf('Rician noise estimate for %s: %3.2f\n',nm,h);
  else
    h = gaussian_noise_estimation(in);
	  fprintf('Gaussian noise estimate for %s: %3.2f\n',nm,h);
  end

  h = 0.65*h;
  
  out = ornlmMex(in,3,1,h);
  V(i).fname = fullfile(pth,['ornlm_' nm xt vr]);
  V(i).descrip = sprintf('ORNLM filtered h=%3.2f',h);
  spm_write_vol(V(i), out);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

return
