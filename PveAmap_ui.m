function PveAmap_ui

!rm *.mexmaci

pve = 1;		% 0 - No PVE; 1 - marginalized; 2 - Kmeans
method = 3;	% 2- Kmeans; 3 - Bayes
warp = 1;

P = spm_select(Inf,'image','Select normalized images to segment');
n = size(P,1);
V = spm_vol(P);	

Vmask = spm_vol('brainmask_LPBA40.nii');

Vpriors = spm_vol(str2mat('/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,1',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,2',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,3',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,4',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,5',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,6'));

for i=1:n

    mask = zeros(V(i).dim(1:3),'uint8');
    priors = zeros([V(i).dim(1:3) 6],'uint8');
	
	vol = spm_read_vols(V(i));


	for j=1:V(i).dim(3),

		Mi  = spm_matrix([0 0 j]);

    % Load slice j from all images
    M1  = V(i).mat\Vpriors(1).mat\Mi;
    mask(:,:,j) = uint8(round(255*spm_slice_vol(Vmask,M1,V(i).dim(1:2),[1 0])));
    for k=1:6
        priors(:,:,j,k) = uint8(round(255*spm_slice_vol(Vpriors(k),M1,V(i).dim(1:2),[1 0])));
    end
	end

	% dilate the mask by convolving
	k = [1 1 1 1 1];
	spm_conv_vol(mask,mask,k,k,k,-[1 1 1]);

	vx = sqrt(sum(V(1).mat(1:3,1:3).^2));

	[pth,nm,xt,vr] = fileparts(deblank(V(i).fname));
for pve=1
  for method=3
	tic;prob = PveAmapMex(vol, priors, mask, vx, pve, method, warp);toc
	[maxi,label] = max(prob,[],4);
	label(find(sum(prob,4)==0)) = 0;

	V(i).fname    = fullfile(pth,['p0' nm '_' num2str(pve) '_' num2str(method) '.nii']);
	V(i).pinfo = [1 0 0]';
	spm_write_vol(V(i),label);

	V(i).fname    = fullfile(pth,['p1' nm '_' num2str(pve) '_' num2str(method) '.nii']);
	V(i).pinfo = [1/255 0 0]';
	spm_write_vol(V(i),prob(:,:,:,2));
end
end

end