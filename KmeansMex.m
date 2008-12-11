function label = KmeansMex(src, mask, separations, iters_nu)

disp('Compiling KmeansMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex KmeansMex.c Kmeans.c splineSmooth.cc -lebtks
cd(p_path);

% -L/usr/local/lib -I/usr/local/include
label = KmeansMex(src, mask, separations, iters_nu);

return
