function cg_sanlm(vargin)
% 
% Spatial Adaptive Non Local Means Denoising Filter
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_sanlm.m 224 2009-12-02 23:39:15Z gaser $

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
	
	in = single(spm_read_vols(V(i)));
  sanlmMex(in,3,1);
  
  V(i).fname = fullfile(pth,['sanlm_' nm xt vr]);
  V(i).descrip = sprintf('SANLM filtered');
  spm_write_vol(V(i), in);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

return
