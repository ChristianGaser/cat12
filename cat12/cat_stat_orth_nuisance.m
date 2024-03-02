function G_orth = cat_stat_orth_nuisance(C,G)
% Orthogonalization of nuisance parameters G w.r.t C
% FORMAT G_orth = cat_stat_orth_nuisance(C,G)
% C      - vector of covariate parameter
% G      - matrix of nuisance parameter
% 
% G_orth - orthogonalized nuisance parameter
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $ 

if nargin < 2
  fprintf('Syntax: G = cat_stat_orth_nuisance(C,G)\n');
  G = [];
  return
end

if size(C,2) > size(C,1)
  C = C';
  G = G';
  trans = 1;
else
  trans = 0;
end

if size(C,1) ~= size(G,1)
  fprintf('Parameters C and G have different length: %d vs %d\n',size(C,1),size(G,1));
  G = [];
  return
end

if size(C,2) ~= 1
  fprintf('Parameters C must have only one column: %d \n',size(C,2));
  G = [];
  return
end

C2 = [ones(size(C)) C];
G_orth = G - C2*(pinv(C2)*G);

if trans
  G_orth = G_orth';
end
