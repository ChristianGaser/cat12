function cg_vbm8_brainmask_ui

% get filenames
%----------------------------------------------------------------------------
spm_figure('Clear','Interactive');
set(spm_figure('FindWin','Interactive'),'Name','Brainmask')

P     = spm_select(Inf,'image','select scans');
V = spm_vol(P);

n     = size(P,1);
% implement the computing
%---------------------------------------------------------------------------
set(spm_figure('FindWin','Interactive'),'Name','executing','Pointer','watch');
spm_progress_bar('Init',n,'Brainmask','Volumes Complete');
for i = 1:n
	mask = cg_vbm8_brainmask(V(i),1);
	[pth,nam,ext] = spm_fileparts(V(i).fname);
	V(i).fname = fullfile(pth,['brain_' nam ext]);
	V(i).dt(1) = 2;
	V(i).pinfo(1:2) = [1 0]';
	spm_write_vol(V(i),mask);
	spm_progress_bar('Set',i);
end
spm_figure('Clear','Interactive');
