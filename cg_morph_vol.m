function vol = cg_morph_vol(in,action,n,th);
% morphological operations to 3D data
%__________________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin < 4, th = 0.5; end
if nargin < 3, n = 1; end
if nargin < 2, action = 'open'; end
if nargin < 1
	error('No arguments given.');
end

th = th*max(double(in(:)));

vol = uint8(in>th);

kx = [1 1 1];
ky = [1 1 1];
kz = [1 1 1];

order = sum(kx(:) ~= 0)*sum(ky(:) ~= 0);


switch lower(action)
	case 'dilate'
	%=======================================================================
	% enlarge image according to number of dilations
	sz = size(vol);
	vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
	vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2~=0);
	end
	vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
	clear vol2

	case 'erode'
	%=======================================================================
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end

	case 'close'
	%=======================================================================
	% enlarge image according to number of dilations
	sz = size(vol);
	vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
	vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2~=0);
	end
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2>=order);
	end
	vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
	clear vol2
	
	case 'open'
	%=======================================================================
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol~=0);
	end

	otherwise
		error('Unknown action');
end

if isa(in,'double')
	vol = double(vol);
end

if isa(in,'single')
	vol = single(vol);
end

return;
