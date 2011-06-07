function compile

rand('state',0);
d = single(rand(50,50,50));

mex -O AmapMex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c
mex -O median3.c
d2 = median3(d);
disp('Compilation of median3 successful')
mex -O eikonal3.c
d2 = eikonal3(d);
disp('Compilation of eikonal3 successful')
mex -O down_cut.c
d2 = down_cut(d,d.^1.5,1);
disp('Compilation of down_cut successful')
mex -O vbdist.c
d2 = vbdist(d);
disp('Compilation of vbdist successful')

try % try OpenMP support
    if strcmp(mexext,'mexmaci64')
        mex CC='gcc-4.2' CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O sanlmMex.c sanlm_float.c
        movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.2' CFLAGS='-fopenmp -m64 -fPIC -O3' -O /usr/local/lib/x86_64/libgomp.a sanlmMex.c sanlm_float.c
    elseif strcmp(mexext,'mexmaci')
        mex CC='gcc-4.2' CFLAGS='-U_OPENMP -m32 -fPIC -O3' -O sanlmMex.c sanlm_float.c
        movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
        mex CC='gcc-4.2' CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/local/lib/x86/libgomp.a sanlmMex.c sanlm_float.c
    elseif strcmp(mexext,'mexa64')
        mex CFLAGS='-U_OPENMP -m64 -fPIC -O3' -O sanlmMex.c sanlm_float.c
        movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m64 -fPIC -O3' -O -lgomp sanlmMex.c sanlm_float.c
    elseif strcmp(mexext,'mexglx')
        mex CFLAGS='-U_OPENMP -m32 -fPIC -O3' -O sanlmMex.c sanlm_float.c
        movefile(['sanlmMex.' mexext], ['sanlmMex_noopenmp.' mexext],'f');
        mex CFLAGS='-fopenmp -m32 -fPIC -O3' -O /usr/lib/gcc/i686-linux-gnu/4.4.4/libgomp.a sanlmMex.c sanlm_float.c
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