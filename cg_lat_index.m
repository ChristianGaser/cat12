function cg_lat_index(P)
% calculation of lateralization index
% LI = 2*(R-L)/(R+L)
% new filenames are prefixed with 'li_'
% FORMAT cg_lat_index(P)
%__________________________________________________________________________
% @(#)cg_lat_index.m	cg 1.0 02/07/13

% get filenames
%----------------------------------------------------------------------------
if nargin == 0
	P = spm_select(Inf,'image','select scans');
end

V = spm_vol(P);
n = prod(size(V));

spm('Pointer','Watch');
spm_progress_bar('Init',n,'Calculating','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
	VO = V(i);
	VO.fname    = fullfile(pth,['li_' nm xt vr]);
	VO.descrip  = [VO.descrip ' LI=2*(R-L)/(R+L)'];
	li = '2*(i1-flipud(i1))./(i1+flipud(i1)+eps)';
	spm_imcalc_ui(V(i).fname,VO.fname,li,{0 0 4 1});
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear',i);
spm('Pointer');

return;
