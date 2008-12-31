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
  double *src, *mean, *vx;
  const int *dims;
  int dims2[4];
    
  if (nrhs!=4)
    mexErrMsgTxt("4 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  if (!mxIsUint8(prhs[2]))
	mexErrMsgTxt("Third argument must be uint8.");

  src   = (double*)mxGetPr(prhs[0]);
  label = (unsigned char*)mxGetPr(prhs[1]);
  mask  = (unsigned char*)mxGetPr(prhs[2]);
  vx    = (double*)mxGetPr(prhs[3]);

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = dims[0];
  dims2[1] = dims[1];
  dims2[2] = dims[2];
  dims2[3] = 5;
  
  plhs[0] = mxCreateNumericArray(4,dims2,mxUINT8_CLASS,mxREAL);
  plhs[1] = mxCreateNumericMatrix(1,5,mxDOUBLE_CLASS,mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[0]);
  mean  = (double *)mxGetPr(plhs[1]);
  
  PveAmap(src, label, mask, prob, mean, vx, dims);

}

