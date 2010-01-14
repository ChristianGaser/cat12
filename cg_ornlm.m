function cg_ornlm
% 
% Optimized Blockwise Non Local Means Denoising Filter
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_ornlm.m 224 2009-12-02 23:39:15Z gaser $

P = spm_select(Inf,'image','Select images to filter');
V = spm_vol(P);
n = size(P,1);

spm_progress_bar('Init',n,'Smoothing','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = fileparts(deblank(V(i).fname));
	in = spm_read_vols(V(i));
	h = rician_local_estimate(in);
	if h==0
	  fprintf('Image %s has no background noise (probably skull stripped).\n',nm);
	else
	  fprintf('Rician noise estimate for %s: %3.2f\n',nm,h);
    out = ornlmMex(in,3,1,h);
    V(i).fname = fullfile(pth,['ornlm_' nm xt vr]);
    spm_write_vol(V(i), out);
  end
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

return
