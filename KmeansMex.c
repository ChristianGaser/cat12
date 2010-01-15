#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label;
  unsigned char *mask;
  double *src, voxelsize[3], mx;
  int nvox, i, n_classes;
  const int *dims;
    
  if (nrhs!=2)
    mexErrMsgTxt("2 inputs required.");
  else if (nlhs>1)
    mexErrMsgTxt("Too many output arguments.");
  
  src = (double*)mxGetPr(prhs[0]);
  n_classes = (int)mxGetScalar(prhs[1]);
  
  for (i=0; i<3; i++)
    voxelsize[i] = 1.0;

  dims = mxGetDimensions(prhs[0]);

  plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
  label  = (unsigned char *)mxGetPr(plhs[0]);
  
  nvox = dims[0]*dims[1]*dims[2];
  mask = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
  for (i=0; i<nvox; i++)
    mask[i] = (src[i]>0) ? 255 : 0;
  mx = Kmeans(src, label, mask, 25, n_classes, voxelsize, dims, 0, 128, 0, 0,  500.0);
  free(mask);
  
}

