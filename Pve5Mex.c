/*
 * Christian Gaser
 * $Id: PveMex.c 8 2008-12-10 20:21:27Z gaser $ 
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
  double *src, *mean;
  const int *dims;
  int nc;
  int BG;
    
  if (nrhs!=4)
    mexErrMsgTxt("4 inputs required.");
  else if (nlhs>0)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsUint8(prhs[1]))
	mexErrMsgTxt("Second argument must be uint8.");

  if (!mxIsUint8(prhs[2]))
	mexErrMsgTxt("Second argument must be uint8.");

  src   = (double*)mxGetPr(prhs[0]);
  prob  = (unsigned char*)mxGetPr(prhs[1]);
  label = (unsigned char*)mxGetPr(prhs[2]);
  mean  = (double*)mxGetPr(prhs[3]);

  dims = mxGetDimensions(prhs[0]);

  Pve5(src, prob, label, mean, dims);

}

