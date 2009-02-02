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
  unsigned char *label;
  unsigned char *prob;
  double *src, weight_MRF, *mean;
  const int *dims;
  int dims2[4];
  int nc, pve, update_label;
  int niters, sub, niters_nu;
    
  if (nrhs!=6)
    mexErrMsgTxt("6 inputs required.");
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
  Amap(src, label, prob, mean, nc, niters, sub, dims, pve);
  if(pve) Pve6(src, prob, label, mean, dims, update_label);

}

