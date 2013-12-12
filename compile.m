function compile

  rand('state',0);
  d = single(rand(10,10,10));
  
  mex -O AmapMex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
  mex -O vbm_vol_median3.c
  mex -O vbm_vol_median3c.c
  mex -O vbm_vol_eikonal3.c
  mex -O vbm_vol_downcut.c
  mex -O vbm_vol_laplace3.c
  mex -O vbm_vol_laplace3R.c
  mex -O vbm_vol_gradient3.c
  mex -O vbm_vol_simgrow.c
  mex -O vbm_vol_localstat.c
  mex -O vbm_vol_pbtp.cpp
  mex -O vbm_vol_interp3f.cpp
  mex -O vbm_vol_eidist.c
  mex -O vbdist.c
  
  d2 = vbm_vol_median3(d);             disp('Compilation of vbm_vol_median3 successful')
  d2 = vbm_vol_median3c(d);            disp('Compilation of vbm_vol_median3c successful')
  d2 = vbm_vol_eikonal3(d);            disp('Compilation of vbm_vol_eikonal3 successful')
  d2 = vbm_vol_laplace3(d,0,0,0.001);  disp('Compilation of vbm_vol_laplace3 successful')
  d2 = vbm_vol_laplace3R(d,d>0.5,0.2); disp('Compilation of vbm_vol_laplace3R successful')
  [d2,d3,d4] = vbm_vol_gradient3(d);   disp('Compilation of vbm_vol_gradient3 successful')
  d2 = vbm_vol_downcut(d,d.^1.5,1);    disp('Compilation of vbm_vol_down_cut successful')
  d2 = vbdist(d);                      disp('Compilation of vbdist successful')
  d2 = vbm_vol_interp3f(d,d,d,d);      disp('Compilation of vbm_vol_interp3f successful')
  d2 = vbm_vol_localstat(d,d>0);       disp('Compilation of vbm_vol_localstat successful')
  d2 = vbm_vol_simgrow(d,d,0.01);      disp('Compilation of vbm_vol_simgrow successful')
  d2 = vbm_vol_eidist(d,d);            disp('Compilation of vbm_vol_eidist successful')
  d2 = vbm_vol_pbtp(3*d,d,d);          disp('Compilation of vbm_vol_pbtp successful')

  try % try OpenMP support
      if strcmp(mexext,'mexmaci64')
          mex CC='gcc-4.0' CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex CC='gcc-4.0' CFLAGS='-m64 -fPIC -O3' -O /usr/local/lib/libgomp.a sanlmMex.c sanlm_float.c
      elseif strcmp(mexext,'mexmaci')
          mex CC='gcc-4.0' CFLAGS='-U_OPENMP -m32 -fPIC -O3' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex CC='gcc-4.0' CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/local/lib/libgomp.a sanlmMex.c sanlm_float.c
      elseif strcmp(mexext,'mexa64')
          mex CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex CFLAGS='-fopenmp -m64 -fPIC -O3' -O -lgomp sanlmMex.c sanlm_float.c
      elseif strcmp(mexext,'mexglx')
          mex CFLAGS='-U_OPENMP -m32 -fPIC -O3' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/lib/i386-linux-gnu/gcc/i686-linux-gnu/4.4/libgomp.a sanlmMex.c sanlm_float.c
      elseif strcmp(mexext,'mexw64')
          mex CFLAGS='-U_OPENMP -m64' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex -O sanlmMex.c sanlm_float.c
      elseif strcmp(mexext,'mexw32')
          mex CFLAGS='-U_OPENMP -m32' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex -O sanlmMex.c sanlm_float.c
      end
      disp('Compiling sanlmMex with OpenMP')
  catch 
      disp('Compiling sanlmMex without OpenMP')
      mex CFLAGS='-fPIC -O3' -O sanlmMex.c sanlm_float.c 
  end
  
  sanlmMex(d,3,1);
  disp('Compilation of sanlmMex successful')
  
end
