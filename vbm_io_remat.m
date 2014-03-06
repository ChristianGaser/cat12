function vbm_io_remat(P1,Pn)
% set orientation of Pn by orientation of P1

  if ~exist('P1','var') || isempty(P1)
    P1 = spm_select([1 1],'image','Select correct orientated image');
  else
    P1 = char(P1);
  end
  if ~exist('Pn','var') || isempty(Pn)
    Pn = spm_select(inf,'image','Select uncorrect orientated images');
  else
    Pn = char(Pn);
  end

  V1 = spm_vol(P1);
  Vn = spm_vol(Pn);

  for i=1:numel(Vn)
    Y = spm_read_vols(Vn(i));
    Vn(i).mat = V1.mat;
    spm_write_vol(Vn(i),Y);
  end
end