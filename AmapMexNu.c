/*
 * Christian Gaser
 * $Id: AmapMex.c 223 2009-12-01 16:03:56Z gaser $ 
 *
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"
#include "Amap.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label, *prob, *mask;
  double *src, *mean, *voxelsize;
  double max_vol, mrf;
  const int *dims;
  int dims2[4];
  int i, n_classes, pve, nvox;
  int niters, sub, init;
    
  if (nrhs!=9)
    mexErrMsgTxt("9 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  src       = (double*)mxGetPr(prhs[0]);
  label     = (unsigned char*)mxGetPr(prhs[1]);
  n_classes = (int)mxGetScalar(prhs[2]);
  niters    = (int)mxGetScalar(prhs[3]);
  sub       = (int)mxGetScalar(prhs[4]);
  pve       = (int)mxGetScalar(prhs[5]);
  init      = (int)mxGetScalar(prhs[6]);
  mrf       = (double)mxGetScalar(prhs[7]);
  voxelsize = (double*)mxGetPr(prhs[8]);

  if ( mxGetM(prhs[8])*mxGetN(prhs[8]) != 3) 
    mexErrMsgTxt("Voxelsize should have 3 values.");

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = dims[0]; dims2[1] = dims[1]; dims2[2] = dims[2];
  dims2[3] = n_classes;
  
  /* for PVE we need more classes */
  if(pve == 6) dims2[3] += 3;
  if(pve == 5) dims2[3] += 2;

  plhs[0] = mxCreateNumericArray(4, dims2, mxUINT8_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, n_classes+3, mxDOUBLE_CLASS, mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[0]);
  mean  = (double *)mxGetPr(plhs[1]);

  /* initial labeling using Kmeans */
  if (init) {
    nvox = dims[0]*dims[1]*dims[2];
    mask = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    for (i=0; i<nvox; i++)
      mask[i] = (src[i]>0) ? 255 : 0;
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, 0, 128, 50, KMEANS, 500.0);
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, 0, 128, 50, NOPVE,  500.0);
    free(mask);
  }
  
  Amap(src, label, prob, mean, n_classes, niters, sub, dims2, pve, mrf, voxelsize);
  if(pve==6) Pve6(src, prob, label, mean, dims2);
  if(pve==5) Pve5(src, prob, label, mean, dims2);

}

