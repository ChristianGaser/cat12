function compile

mex -O upfirdn2dMex.c 
mex -O ornlmMex.c ornlm.c 
mex -O KmeansMex.c Kmeans.c vollib.c
mex -O AmapMex.c Kmeans.c Amap.c MrfPrior.c Pve.c vollib.c