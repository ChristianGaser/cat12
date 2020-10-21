#include "mex.h"
#include "math.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label, *priors, *prob;
  double *src, *separations;
  const int *dims;
  int niters, niters_nu;
    
  if (nrhs!=4)
    mexErrMsgTxt("4 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  src    = (double*)mxGetPr(prhs[0]);
  priors  = (unsigned char*)mxGetPr(prhs[1]);
  separations  = (double*)mxGetPr(prhs[2]);
  niters_nu = (int)mxGetScalar(prhs[3]);

  dims = mxGetDimensions(prhs[1]);

  plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
  label  = (unsigned char *)mxGetPr(plhs[0]);

  plhs[1] = mxCreateNumericArray(4,dims,mxUINT8_CLASS,mxREAL);
  prob  = (unsigned char *)mxGetPr(plhs[1]);
  
  Bayes(src, label, priors, prob, separations, dims, niters_nu);

}

