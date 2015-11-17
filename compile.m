function compile

  rand('state',0);
  d0  = single(rand(10,10,10));
  d0(5,5,5) = NaN;
    
  if strcmp(mexext,'mexmaci64')
    mexflag='-Dchar16_t=UINT16_T CFLAGS=''$CFLAGS -Wall -ansi -pedantic -Wextra'' CPPLAGS=''$CPPFLAGS -Wall -ansi -pedantic -Wextra''';
  else
    mexflag='';
  end

  eval(['mex ' mexflag ' -O AmapMex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c'])
  eval(['mex ' mexflag ' -O cat_vol_median3.c'])
  eval(['mex ' mexflag ' -O cat_vol_median3c.c'])
  eval(['mex ' mexflag ' -O cat_vol_downcut.c'])
  eval(['mex ' mexflag ' -O cat_vol_laplace3.c'])
  eval(['mex ' mexflag ' -O cat_vol_laplace3R.c'])
  eval(['mex ' mexflag ' -O cat_vol_gradient3.c'])
  eval(['mex ' mexflag ' -O cat_vol_simgrow.c'])
  eval(['mex ' mexflag ' -O cat_vol_localstat.c'])
  eval(['mex ' mexflag ' -O cat_vol_pbtp.cpp'])
  eval(['mex ' mexflag ' -O cat_vol_interp3f.cpp'])
  eval(['mex ' mexflag ' -O cat_vol_eidist.c'])
  eval(['mex ' mexflag ' -O cat_vol_eulernumber.c'])
  eval(['mex ' mexflag ' -O cat_vol_genus0.c genus0.c'])
  eval(['mex ' mexflag ' -O vbdist.c'])
  eval(['mex ' mexflag ' -O ornlmMex.c ornlm_float.c'])
  eval(['mex ' mexflag ' -O sanlmMex.c sanlm_float.c'])
  
  %%
  
  d = cell(17, 1);
  d{1} = cat_vol_pbtp(3*d0,d0,d0);          disp('Compilation of cat_vol_pbtp successful')
  d{2} = cat_vol_median3(d0);               disp('Compilation of cat_vol_median3 successful')
  d{3} = cat_vol_median3c(d0,d0==0);        disp('Compilation of cat_vol_median3c successful')
  d{4} = ornlmMex(d0,3,1,0.1);              disp('Compilation of ornlmMex successful')
  d{5} = cat_vol_laplace3(d0,0,0,0.001);    disp('Compilation of cat_vol_laplace3 successful')
  d{6} = cat_vol_laplace3R(d0,d0>0.5,0.2);  disp('Compilation of cat_vol_laplace3R successful')
  [d{7},d{8},d{9}] = cat_vol_gradient3(d0); disp('Compilation of cat_vol_gradient3 successful')
  d{10} = cat_vol_downcut(d0,d0.^1.5,1);    disp('Compilation of cat_vol_down_cut successful')
  d{11} = vbdist(d0);                       disp('Compilation of vbdist successful')
  d{12} = cat_vol_interp3f(d0,d0,d0,d0);    disp('Compilation of cat_vol_interp3f successful')
  d{13} = cat_vol_localstat(d0,d0>0);       disp('Compilation of cat_vol_localstat successful')
  d{14} = cat_vol_simgrow(d0,d0,0.01);      disp('Compilation of cat_vol_simgrow successful')
  d{15} = cat_vol_eidist(d0,d0);            disp('Compilation of cat_vol_eidist successful')
  d{16} = cat_vol_eulernumber(double(d0>0));disp('Compilation of cat_vol_eulernumber successful')
  [tmp,CS.faces,CS.vertices] = cat_vol_genus0(d0,0.5); disp('Compilation of cat_vol_genus0 successful')
    
  sanlmMex(d0,3,1);
  d{17} = d0;
  disp('Compilation of sanlmMex successful')

  debugname = ['debug_' mexext '.mat'];
  disp(['save ' debugname]);
  save(debugname,'d','CS');
  
  
end
