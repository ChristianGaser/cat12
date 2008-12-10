#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label;
  unsigned char *prob;
  double *src, weight_MRF, mean[100];
  const int *dims;
  int dims2[4];
  int nc;
  int BG, niters, nflips, sub, niters_nu;
    
  if (nrhs!=8)
    mexErrMsgTxt("8 inputs required.");
  else if (nlhs>1)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  src    = (double*)mxGetPr(prhs[0]);
  label  = (unsigned char*)mxGetPr(prhs[1]);
  nc     = (int)mxGetScalar(prhs[2]);
  BG     = (int)mxGetScalar(prhs[3]);
  niters = (int)mxGetScalar(prhs[4]);
  nflips = (int)mxGetScalar(prhs[5]);
  sub    = (int)mxGetScalar(prhs[6]);
  weight_MRF = (double)mxGetScalar(prhs[7]);

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = dims[0];
  dims2[1] = dims[1];
  dims2[2] = dims[2];
  dims2[3] = 3;

  plhs[0] = mxCreateNumericArray(4,dims2,mxUINT8_CLASS,mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[0]);
  Amap(src, label, prob, mean, nc, BG, niters, nflips, sub, dims, weight_MRF);

}

