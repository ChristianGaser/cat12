function label = KmeansMex(src, n_classes)
%
% Christian Gaser
% $Id$

disp('Compiling KmeansMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O KmeansMex.c KmeansProper.c
cd(p_path);

label = KmeansMex(src, n_classes);

return
