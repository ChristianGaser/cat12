function y = cat_stat_polynomial(x,p)
% Polynomial expansion and orthogonalization of function x
% FORMAT y = cat_stat_polynomial(x,p)
% x   - data matrix
% p   - order of polynomial [default: 1]
% 
% y   - orthogonalized data matrix
%__________________________________________________________________________
%
% cat_stat_polynomial orthogonalizes a polynomial function of order p
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $ 

if nargin < 2, p = 1; end

if size(x,1) < size(x,2)
    x = x';
end

y = spm_detrend(x(:));
v = zeros(size(y,1),p + 1);

for j = 0:p
    v(:,(j + 1)) = (y.^j) - v*(pinv(v)*(y.^j));
end

for j = 2:p
    y = [y v(:,(j + 1))];
end
