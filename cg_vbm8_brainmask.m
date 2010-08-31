function mask = cg_vbm8_brainmask(in, th)

if nargin < 2
  th = 1;
end

if ~isstruct(in)
  in = spm_vol(in);
end

src = spm_read_vols(in);

h = cg_noise_estimation(src);

src = ornlmMex(src, 3, 1, 5*h);

% initial label
label = KmeansMex(src,3);

% refine label by removing background and increasing numbers of classes
label = KmeansMex(src.*(label>0),6);
%label = KmeansMex(src.*(label<5),3);
mask = (label > 1) & (label < 5);
D = bwdist(~mask);

% keep largest connected component after 2 its of opening
mask = cg_morph_vol(D>th,'open',3,0.5);
mask = mask_largest_cluster(mask,0.5);

% dilate and close to fill ventricles
mask = cg_morph_vol(mask,'dilate',3,0.5);
mask = cg_morph_vol(mask,'close',10,0.5);

return

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
