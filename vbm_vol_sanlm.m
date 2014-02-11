function vbm_vol_sanlm(varargin)
% 
% Spatial Adaptive Non Local Means Denoising Filter
%
%_______________________________________________________________________
% Christian Gaser/Users/gaser/matlab/vbm8/cg_sanlm.m
% $Id$

if nargin > 1 && ischar(varargin{2})
  prefix = varargin{2};
else
  prefix = 'sanlm_';
end
  
if nargin > 0 && ~isempty(varargin{1})
  if isstruct(varargin{1})
  	P = char(varargin{1}.data);
  else
    P = char(varargin{1});
  end
else
  P = spm_select(Inf,'image','Select images to filter');
end

V = spm_vol(P);
n = size(P,1);

spm_progress_bar('Init',n,'Filtering','Volumes Complete');
for i = 1:n
	[pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));
	
	src = single(spm_read_vols(V(i)));
	% prevent NaN
  src(isnan(src)) = 0;
  try
      sanlmMex(src,3,1);
  catch
      sanlmMex_noopenmp(src,3,1);
  end
  
  V(i).fname = fullfile(pth,[prefix nm '.nii' vr]);
  V(i).descrip = sprintf('%s SANLM filtered',V(i).descrip);
  % use at least float precision
  if  V(i).dt(1)<16 V(i).dt(1) = 16; end 
  spm_write_vol(V(i), src);
	spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

return
