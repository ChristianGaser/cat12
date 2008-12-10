function prob = AmapMex(src, label, nc, BG, niters, nflips, sub, weight_MRF)
%
% Christian Gaser
% $Id$

disp('Compiling AmapMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O AmapMex.c Amap.c MrfPrior.c
cd(p_path);

prob = AmapMex(src, label, nc, BG, niters, nflips, sub, weight_MRF);

return
