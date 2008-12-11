#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label;
  unsigned char *mask;
  double *src, *separations, mx;
  const int *dims;
  int niters, niters_nu;
    
  if (nrhs!=4)
    mexErrMsgTxt("4 inputs required.");
  else if (nlhs>1)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  src    = (double*)mxGetPr(prhs[0]);
  mask  = (unsigned char*)mxGetPr(prhs[1]);
  separations  = (double*)mxGetPr(prhs[2]);
  niters_nu = (int)mxGetScalar(prhs[3]);

  dims = mxGetDimensions(prhs[1]);

  plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
  label  = (unsigned char *)mxGetPr(plhs[0]);
  
  mx = Kmeans(src, label, mask, 25, 3, 1, separations, dims, 0, 50, 128, niters_nu);

}

