function label = BayesMex(src, priors, separations, iters_nu)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling BayesMex.c')
mex -O BayesMex.c Bayes.c splineSmooth.cc -lebtks
% -L/usr/local/lib -I/usr/local/include
label = BayesMex(src, priors, separations, iters_nu);

return
