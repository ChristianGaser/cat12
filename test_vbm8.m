function test_vbm8

if 1
V = spm_vol('/Users/gaser/Desktop/SVE.LPBA40.testdata/S01.native.mri.nii');
vol = spm_read_vols(V);

Vmask = spm_vol('brainmask_LPBA40.nii');
mask = zeros(V(1).dim(1:3),'uint8');

Vpriors = spm_vol(str2mat('/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,1',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,2',...
    '/Users/gaser/spm/spm8b/toolbox/Seg/TPM.nii,3'));
priors = zeros([V(1).dim(1:3) 3],'uint8');

for j=1:V(1).dim(3),

    Mi  = spm_matrix([0 0 j]);

    % Load slice j from all images
    M1  = V(1).mat\Vpriors(1).mat\Mi;
    mask(:,:,j) = uint8(round(255*spm_slice_vol(Vmask,M1,V(1).dim(1:2),[1 0])));
    for i=1:3
        priors(:,:,j,i) = uint8(round(255*spm_slice_vol(Vpriors(i),M1,V(1).dim(1:2),[1 0])));
    end
end

vx = sqrt(sum(V(1).mat(1:3,1:3).^2));

save all
else

!rm *.mexmaci
load all
tic;prob = PveAmapMex(vol, priors, mask, vx);toc

end