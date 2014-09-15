function compile

  rand('state',0);
  d0  = single(rand(10,10,10));
  d0(5,5,5) = NaN;
    
  if strcmp(mexext,'mexmaci64')
    mexflag='-Dchar16_t=UINT16_T';
  else
    mexflag='';
  end

  eval(['mex ' mexflag ' -O AmapMex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c'])
  eval(['mex ' mexflag ' -O vbm_vol_median3.c'])
  eval(['mex ' mexflag ' -O vbm_vol_median3c.c'])
  eval(['mex ' mexflag ' -O vbm_vol_downcut.c'])
  eval(['mex ' mexflag ' -O vbm_vol_laplace3.c'])
  eval(['mex ' mexflag ' -O vbm_vol_laplace3R.c'])
  eval(['mex ' mexflag ' -O vbm_vol_gradient3.c'])
  eval(['mex ' mexflag ' -O vbm_vol_simgrow.c'])
  eval(['mex ' mexflag ' -O vbm_vol_localstat.c'])
  eval(['mex ' mexflag ' -O vbm_vol_pbtp.cpp'])
  eval(['mex ' mexflag ' -O vbm_vol_interp3f.cpp'])
  eval(['mex ' mexflag ' -O vbm_vol_eidist.c'])
  eval(['mex ' mexflag ' -O vbm_vol_genus0.c genus0.c'])
  eval(['mex ' mexflag ' -O vbdist.c'])
  eval(['mex ' mexflag ' -O ornlmMex.c ornlm_float.c'])
  
  %%
  
  d = cell(17, 1);
  d{1} = vbm_vol_pbtp(3*d0,d0,d0);         disp('Compilation of vbm_vol_pbtp successful')
  d{2} = vbm_vol_median3(d0);              disp('Compilation of vbm_vol_median3 successful')
  d{3} = vbm_vol_median3c(d0,d0==0);       disp('Compilation of vbm_vol_median3c successful')
  d{4} = ornlmMex(d0,3,1,0.1);             disp('Compilation of ornlmMex successful')
  d{5} = vbm_vol_laplace3(d0,0,0,0.001);   disp('Compilation of vbm_vol_laplace3 successful')
  d{6} = vbm_vol_laplace3R(d0,d0>0.5,0.2); disp('Compilation of vbm_vol_laplace3R successful')
  [d{7},d{8},d{9}] = vbm_vol_gradient3(d0);disp('Compilation of vbm_vol_gradient3 successful')
  d{10} = vbm_vol_downcut(d0,d0.^1.5,1);   disp('Compilation of vbm_vol_down_cut successful')
  d{11} = vbdist(d0);                      disp('Compilation of vbdist successful')
  d{12} = vbm_vol_interp3f(d0,d0,d0,d0);   disp('Compilation of vbm_vol_interp3f successful')
  d{13} = vbm_vol_localstat(d0,d0>0);      disp('Compilation of vbm_vol_localstat successful')
  d{14} = vbm_vol_simgrow(d0,d0,0.01);     disp('Compilation of vbm_vol_simgrow successful')
  d{15} = vbm_vol_eidist(d0,d0);           disp('Compilation of vbm_vol_eidist successful')
  [tmp,CS.faces,CS.vertices] = vbm_vol_genus0(d0,0.5); disp('Compilation of vbm_vol_genus0')
  
  %%
  try % try OpenMP support
      if strcmp(mexext,'mexmaci64')
          mex -Dchar16_t=UINT16_T CC='gcc' CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O sanlmMex.c sanlm_float.c
          movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
          mex -Dchar16_t=UINT16_T CC='gcc' CFLAGS='-m64 -fPIC -O3' -O /usr/local/lib/libgomp.a sanlmMex.c sanlm_float.c
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
  
  sanlmMex(d0,3,1);
  d{16} = d0;

  rand('state',0);
  d0  = single(rand(10,10,10));
  d0(5,5,5) = NaN;
  
  sanlmMex_noopenmp(d0,3,1);
  d{17} = d0;

  debugname = ['debug_' mexext '.mat'];
  disp(['save ' debugname]);
  save(debugname,'d','CS');
  
  disp('Compilation of sanlmMex successful')
  
end
