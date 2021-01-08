function correct_tpm2one
% correct last class of TPM or GS-Template to ensure sum=1

P = spm_select(Inf,'image','Select TPMs or GS-Templates');

for i=1:size(P,1)
  N = nifti(deblank(P(i,:)));
  dat = zeros(N.dat.dim);
  dat(:,:,:,:) = N.dat(:,:,:,:);
  nc = N.dat.dim(4);  
  sum_dat = sum(dat(:,:,:,1:(nc-1)),4);
  sum_dat(sum_dat>1) = 1;
  
  % correct last class and ensure sum=1
  dat(:,:,:,nc) = 1 - sum_dat;
  
  [pth,nam,ext] = spm_fileparts(N.dat.fname);

  % write corrected TPM
  No = nifti;
  No.dat  = file_array(fullfile(pth,['c' nam ext]),size(dat),N.dat.dtype,...
      N.dat.offset,N.dat.scl_slope,N.dat.scl_inter);
  No.mat  = N.mat;
  No.mat0 = N.mat;
  No.descrip = N.descrip;
  create(No);
  No.dat(:,:,:,:) = dat(:,:,:,:);
end