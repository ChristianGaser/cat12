function label = BayesMex(src, priors, separations, iters_nu)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling BayesMex.c')
mex -O BayesMex.c Bayes.c vollib.c
label = BayesMex(src, priors, separations, iters_nu);

return
