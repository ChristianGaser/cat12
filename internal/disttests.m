% Image size
N = 128;
img = zeros(N, N);

% Spiral parameters
theta = linspace(0, 6*pi, 5000);
a = 0;
b = 3;

% Thickness parameter (radius in pixels)
thickness = 2;

% Convert spiral to Cartesian coordinates
r = a + b * theta;
x = r .* cos(theta);
y = r .* sin(theta);

% Normalize to image grid
x = round(x + N/2);
y = round(y + N/2);

% Keep valid indices
valid = x > 0 & x <= N & y > 0 & y <= N;
x = x(valid);
y = y(valid);

% Create a disk kernel for thickness
[dx, dy] = meshgrid(-thickness:thickness, -thickness:thickness);
mask = (dx.^2 + dy.^2) <= thickness^2;

% Draw thick spiral
for k = 1:length(x)
    xi = x(k);
    yi = y(k);
    
    % Neighborhood indices
    x_idx = xi + dx;
    y_idx = yi + dy;
    
    % Valid neighborhood mask
    valid_nb = x_idx > 0 & x_idx <= N & y_idx > 0 & y_idx <= N;
    
    % Apply mask
    idx = sub2ind(size(img), y_idx(valid_nb), x_idx(valid_nb));
    mask_valid = mask(valid_nb);
    
    img(idx(mask_valid)) = 1;
end

% Display
figure
imshow(img);
colormap(gray);
axis image;
title(['Spiral with Thickness = ', num2str(thickness)]);

%%
F = repmat(single(1-img),1,1,10); 
O = zeros(size(F),'single');  O(:,[1 size(O,2)],:) = 1;  O([1 size(O,2)],:,:) = 1;   
D1 = cat_vbdist3(O,F>0); 
D2 = cat_vbdist3(F<1,O<1); 
ds('d2sm','',1,D1/200,D2/200,5);
%%
mex cat_vol_pbtp.c
T1 = cat_vol_pbtp(1+F+O,D1,D2); mean(T1(T1(:)>0))
ds('d2sm','',1,D/100,T1/100,5)