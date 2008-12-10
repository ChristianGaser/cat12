function vol = cg_morph_vol(in,action,n,th);
% morphological operations to 3D data
%__________________________________________________________________________
% @(#)cg_morph_vol.m 1.01 2007/01/17 Christian Gaser

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
	case 'erode'
	%=======================================================================
	% dilate
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol~=0);
	end

	case 'dilate'
	%=======================================================================
	% erose
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end

	case 'close'
	%=======================================================================
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol~=0);
	end
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end
	
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
