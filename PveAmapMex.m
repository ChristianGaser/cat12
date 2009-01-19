function [prob, mean] = PveAmapMex(src, priors, mask, vx, pve, method)
%
% Christian Gaser
% $Id$

rev = '$Rev$';

disp('Compiling PveAmapMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
ext = mexext;
ext = ext(4:end);
matlabinclude=fullfile(matlabroot,'extern','include');
try
    eval(['mex -O PveAmapMex.c PveAmap.c Amap.c MrfPrior.c Pve5.c Kmeans.c WarpPriors.c Bayes.c optimizer3d.c diffeo3d.c splineSmooth.cc -lEBTKS -L' ext ' -I' matlabinclude ' -I.' fullefile(filesep,pwd)]);
catch
    mex -O PveAmapMex.c PveAmap.c Amap.c MrfPrior.c Pve5.c Kmeans.c WarpPriors.c Bayes.c optimizer3d.c diffeo3d.c splineSmooth.cc -lEBTKS

end

cd(p_path);

[prob, mean] = PveAmapMex(src, priors, mask, vx, pve, method);

return
