function mask = cg_vbm8_brainmask(VF)

open_th = 0.25;

if ~isstruct(VF)
  VF = spm_vol(VF);
end

% get registered brainmask
priors = get_registered_priors(VF);

% get src
mex -O BayesMex.c Bayes.c vollib.c WarpPriors.c optimizer3d.c diffeo3d.c
src = spm_read_vols(VF);
vx = sqrt(sum((VF.mat(1:3,1:3)).^2))';

% Bayesian segmentation with 6 classes
correct_nu = 1;
[label, cls] = BayesMex(src, priors, vx, correct_nu);

mask = single(cls(:,:,:,1));
mask = mask + single(cls(:,:,:,2));

% keep largest connected component after 2 its of opening
mask = cg_morph_vol(mask,'open',1,open_th);
mask = mask_largest_cluster(mask,0.5);

% dilate and close to fill ventricles
mask = cg_morph_vol(mask,'dilate',1,0.5);
mask = cg_morph_vol(mask,'close',10,0.5);

% remove sinus
mask = mask & single((cls(:,:,:,5)<cls(:,:,:,1)) | ...
                     (cls(:,:,:,5)<cls(:,:,:,2)) | ...
                     (cls(:,:,:,5)<cls(:,:,:,3)));                

% and fill holes that may remain
mask = cg_morph_vol(mask,'close',2,0.5);

return

%=======================================================================
function priors = get_registered_priors(VF);

addpath(fullfile(spm('dir'),'toolbox','Seg'));

affreg = 'mni';
fudge  = 5;
samp   = 4;
tpm    = {fullfile(spm('dir'),'toolbox','Seg','TPM.nii')};

Affine  = eye(4);
VG = spm_vol(fullfile(spm('Dir'),'templates','T1.nii'));
                    
% smooth source with 8mm
VF1 = spm_smoothto8bit(VF,8);

% Rescale images so that globals are better conditioned
VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
VG.pinfo(1:2,:)  = VG.pinfo(1:2,:)/spm_global(VG);

fprintf('Initial Coarse Affine Registration..\n');
aflags    = struct('sep',8, 'regtype',affreg,'WG',[],'WF',[],'globnorm',0);
aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

M = eye(4);
spm_chi2_plot('Init','Affine Registration','Mean squared difference','Iteration');
Affine  = spm_affreg(VG, VF1, aflags, M);

% load TPMs
[pth,nam,ext,num] = spm_fileparts(tpm{1});
for i=1:6
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
end
tpm = spm_load_priors8(strvcat(cat(1,tissue(:).tpm)));

fprintf('Fine Affine Registration..\n');
Affine  = spm_maff8(VF,samp,fudge,tpm,Affine,affreg);

% affine parameters
Ma      = tpm.M\Affine*VF.mat;

% load priors
n_priors = length(tpm.V);
priors = zeros([VF.dim(1:3) n_priors],'uint8');
for i=1:VF.dim(3),
    for j=1:n_priors
        img = spm_slice_vol(tpm.V(j),Ma*spm_matrix([0 0 i]),VF.dim(1:2),[1,NaN]);
        priors(:,:,i,j) = uint8(round(255*img));
    end
end

return
%=======================================================================

%=======================================================================
function y = mask_largest_cluster(y, th)

if nargin < 2
	th = 0.5;
end

sz = size(y);

th = th*max(single(y(:)));

mask = y > th;
Q = find(mask);

Qth = find(y <= th & y>0);
yth = y(Qth);

% save memory by using bounding box where y > th
[indx, indy, indz] = ind2sub(size(mask),Q);
indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

[A,num] = spm_bwlabel(double(mask(indx,indy,indz)),26);

clear mask

if isempty(A)
  error('No cluster found!');
  return
end

% interrupt if cluster was > 7.5% of whole image to save time
max_A = max(A(:));
sz_cluster = zeros(max_A,1);
for i=1:max_A
	QA = find(A == i);
	ind = i;
	if length(QA)/prod(size(A)) > 0.075
		break
	end
	sz_cluster(i) = length(QA);
end

if length(QA)/prod(size(A)) <= 0.075
	[mx, ind] = max(sz_cluster);
	QA = find(A == ind);
end

QA0 = find(A ~= ind);
A = y(indx,indy,indz);
A(QA0) = 0;

y(indx,indy,indz) = A;
y(Qth) = yth;

return;

