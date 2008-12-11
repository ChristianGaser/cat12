function label = KmeansMex(src, mask, separations, iters_nu)

disp('Compiling KmeansMex.c')
mex KmeansMex.c Kmeans.c splineSmooth.cc -lebtks
% -L/usr/local/lib -I/usr/local/include
label = KmeansMex(src, mask, separations, iters_nu);

return
