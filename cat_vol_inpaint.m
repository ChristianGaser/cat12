function out = cat_vol_inpaint(vol,niter,smooth,reduce,init,verb)
% ----------------------------------------------------------------------
% This function uses the inpaintn function from Damien Garcia to replaces
% missing values (indicated by NaN or Inf) with interpolated/extrapolated
% values using discrete cosine transformation.
% It uses an iterative process baased on DCT and IDCT with 10 iterations 
% as default. Optionally the output can be smoothed with a Gaussian kernel.
% Missing areas are initialized using Laplace method.
%
%   out = cat_vol_inpaint(vol,niter,smooth,reduce,init,verb)
%
%   vol      .. input image
%   niter    .. number of iteratiosn for inpainting 
%   smooth   .. size for Gaussian smoothing
%   reduce   .. increase speed by using a reduced image for inpainting
%   init     .. use either euclidean distance (init=1) or Laplace method
%               (init = 2, default) for initialization of missing values
%   verb     .. show progress bar (default=0)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin < 2
  niter = 10;
end

if nargin < 3
  smooth = 2;
end

if nargin < 4
  reduce = 4;
end

% use Laplace method for initialization
if nargin < 5
  init = 2;
end

if nargin < 6
  verb = 0;
end


% check whether NaN or Inf exist
if 0 && sum(~isfinite(vol)) == 0
  error('Your image does not contain any areas for inpainting.')
end

if reduce
  [volr,resTr] = cat_vol_resize(vol,'reduceV',[1 1 1],reduce,32,'max');
  % set zero areas to NaN
  volr(volr == 0) = NaN;
  out = inpaintn(volr,niter,init,[],verb);
  out = cat_vol_resize(out,'dereduceV',resTr); 
else
  % set zero areas to NaN
  vol(vol == 0) = NaN;
  out = inpaintn(vol,niter,init,[],verb);
end

% optional smoothing
if smooth
  out = cat_vol_smooth3X(out,smooth);
end

end

function y = inpaintn(x,n,init,m,verb)

% INPAINTN Inpaint over missing data in N-D array
%   Y = INPAINTN(X) replaces the missing data in X by extra/interpolating
%   the non-missing elements. The non finite values (NaN or Inf) in X are
%   considered as missing data. X can be any N-D array.
%
%   INPAINTN (no input/output argument) runs the following 3-D example.
%
%   Important note:
%   --------------
%   INPAINTN uses an iterative process baased on DCT and IDCT.
%   Y = INPAINTN(X,N) uses N iterations. By default, N = 100. If you
%   estimate that INPAINTN did not totally converge, increase N:
%   Y = INPAINTN(X,1000)
%
%   Y = INPAINTN(X,N,INIT) allows to define the method to initialize
%   missing values:
%   1 - euclidean distance
%   2 - Laplace method
%
%   References (please refer to the two following references)
%   ---------- 
%   1) Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010;54:1167-1178. 
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/csda10.pdf')">download PDF</a>
%   2) Wang G, Garcia D et al. A three-dimensional gap filling method for
%   large geophysical datasets: Application to global satellite soil
%   moisture observations. Environ Modell Softw, 2012;30:139-142.
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/envirmodellsoftw12.pdf')">download PDF</a>
%
%   See also SMOOTHN, GRIDDATAN
%
%   -- Damien Garcia -- 2010/06, last update 2017/08
%   website: <a
%   href="matlab:web('http://www.biomecardio.com/en')">www.BiomeCardio.com</a>
%
%   Copyright (c) 2014, Damien Garcia
%   All rights reserved.


class0 = class(x);
x = double(x);

if nargin==1 || isempty(n), n = 100; end

sizx = size(x);
d = ndims(x);
Lambda = zeros(sizx);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizx(i);
    Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizx(i),siz0)-1)/sizx(i)));
end
Lambda = 2*(d-Lambda);

% Initial condition
if nargin<3 || isempty(init), init = 2; end
W = isfinite(x);
if any(~W(:))
    [y,s0] = InitialGuess(x,isfinite(x),init);
else
    y = x;
    return
end
x(~W) = 0;

% Smoothness parameters: from high to negligible values
s = logspace(s0,-6,n); 

RF = 2; % relaxation factor

if nargin<4 || isempty(m), m = 2; end
Lambda = Lambda.^m;

if verb, h = waitbar(0,'Inpainting...'); end
for i = 1:n
        Gamma = 1./(1+s(i)*Lambda);
        y = RF*idctn(Gamma.*dctn(W.*(x-y)+y)) + (1-RF)*y;
        if verb, waitbar(i/n,h); end
end
if verb, close(h); end

y(W) = x(W);
y = cast(y,class0);

end

%% Initial Guess
function [z,s0] = InitialGuess(y,I,init)

% Euclidean distance
if init == 1
  [D,L] = cat_vbdist(single(I));
  z = y;
  z(~I) = y(L(~I));
  s0 = 3; % note: s = 10^s0
else % Laplace method
  y(~I) = 0;
  z = cat_vol_laplace3R(single(y),I,0.00001);
  s0 = 6; % note: s = 10^s0

end

end

%% DCTN
function y = dctn(y)

%DCTN N-D discrete cosine transform.
%   Y = DCTN(X) returns the discrete cosine transform of X. The array Y is
%   the same size as X and contains the discrete cosine transform
%   coefficients. This transform can be inverted using IDCTN.
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   -- Damien Garcia -- 2008/06, revised 2011/11
%   -- www.BiomeCardio.com --

y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1, y = y.'; end
end

% Weighting vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end

% --- DCT algorithm ---
if ~isreal(y)
    y = complex(dctn(real(y)),dctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = y([1:2:n 2*floor(n/2):-2:2],:);
        y = reshape(y,n,[]);
        y = y*sqrt(2*n);
        y = ifft(y,[],1);
        y = bsxfun(@times,y,w{dim});
        y = real(y);
        y(1,:) = y(1,:)/sqrt(2);
        y = reshape(y,siz);
        y = shiftdim(y,1);
    end
end
        
y = reshape(y,sizy);

end

%% IDCTN
function y = idctn(y)

%IDCTN N-D inverse discrete cosine transform.
%   X = IDCTN(Y) inverts the N-D DCT transform, returning the original
%   array if Y was obtained using Y = DCTN(X).
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also DCTN, IDSTN, IDCT, IDCT2, IDCT3.
%
%   -- Damien Garcia -- 2009/04, revised 2011/11
%   -- www.BiomeCardio.com --

y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        y = y.';
    end
end

% Weighing vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end

% --- IDCT algorithm ---
if ~isreal(y)
    y = complex(idctn(real(y)),idctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = reshape(y,n,[]);
        y = bsxfun(@times,y,w{dim});
        y(1,:) = y(1,:)/sqrt(2);
        y = ifft(y,[],1);
        y = real(y*sqrt(2*n));
        I = (1:n)*0.5+0.5;
        I(2:2:end) = n-I(1:2:end-1)+1;
        y = y(I,:);
        y = reshape(y,siz);
        y = shiftdim(y,1);            
    end
end
        
y = reshape(y,sizy);

end
