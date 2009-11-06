/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label, *prob, *mask;
  double *src, *mean, separations[3];
  double weight_MRF, max_vol, mrf;
  const int *dims;
  int dims2[4];
  int i, nc, pve, update_label, nvox;
  int niters, sub, niters_nu, init;
    
  if (nrhs!=8)
    mexErrMsgTxt("8 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  src    = (double*)mxGetPr(prhs[0]);
  label  = (unsigned char*)mxGetPr(prhs[1]);
  nc     = (int)mxGetScalar(prhs[2]);
  niters = (int)mxGetScalar(prhs[3]);
  sub    = (int)mxGetScalar(prhs[4]);
  pve    = (double)mxGetScalar(prhs[5]);
  init   = (int)mxGetScalar(prhs[6]);
  mrf    = (double)mxGetScalar(prhs[7]);

  /* output label as PVE label */
  update_label = 2;

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = dims[0];
  dims2[1] = dims[1];
  dims2[2] = dims[2];
  dims2[3] = nc;
  
  /* for PVE we need tow more classes */
  if(pve) dims2[3] += 3;

  plhs[0] = mxCreateNumericArray(4,dims2,mxUINT8_CLASS,mxREAL);
  plhs[1] = mxCreateNumericMatrix(1,nc+3,mxDOUBLE_CLASS,mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[0]);
  mean  = (double *)mxGetPr(plhs[1]);

  /* initial labeling using Kmeans */
  if (init) {
    nvox = dims[0]*dims[1]*dims[2];
    mask = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    for (i=0; i<nvox; i++)
      mask[i] = (src[i]>0) ? 255 : 0;
    separations[0] = 1.0; separations[1] = 1.0; separations[2] = 1.0;
    max_vol = Kmeans( src, label, mask, 25, nc, separations, dims, 128, 128, 10, 2, 50.0);
    max_vol = Kmeans( src, label, mask, 25, nc, separations, dims, 0, 128, 10, 0, 50.0);
    free(mask);
  }
  
  Amap(src, label, prob, mean, nc, niters, sub, dims, pve, mrf);
  if(pve) Pve6(src, prob, label, mean, dims, update_label);

}

