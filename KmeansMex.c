/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdio.h"

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *voxelsize, const int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label;
  unsigned char *mask;
  double *src, voxelsize[3], mx, thresh, thresh_kmeans_int, bias_fwhm;
  double *mean, n[100];
  int nvox, i, n_classes, iters_nu;
  const int *dims;
    
  if (nrhs!=2)
    mexErrMsgTxt("2 inputs required.");
  else if (nlhs>2)
    mexErrMsgTxt("Too many output arguments.");
  
  src = (double*)mxGetPr(prhs[0]);
  n_classes = (int)mxGetScalar(prhs[1]);
  
  for (i=0; i<3; i++)
    voxelsize[i] = 1.0;

  dims = mxGetDimensions(prhs[0]);

  plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, n_classes, mxDOUBLE_CLASS, mxREAL);
  label  = (unsigned char *)mxGetPr(plhs[0]);
  mean   = (double *)mxGetPr(plhs[1]);
  
  thresh = 0;
  thresh_kmeans_int = 128;
  iters_nu = 0;
  bias_fwhm = 60.0;

  nvox = dims[0]*dims[1]*dims[2];
  mask = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
  for (i=0; i<nvox; i++)
    mask[i] = (src[i]>0) ? 255 : 0;
  mx = Kmeans(src, label, mask, 25, n_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, 0,  bias_fwhm);

  /* calculate mean */
  for(i = 0; i < n_classes; i++) {
    n[i] = 0;
    mean[i] = 0.0;
  }
  for(i = 0; i < nvox; i++) {
    if(label[i] == 0) continue;
    n[label[i]-1]++;
    mean[label[i]-1] += src[i];
  }
  for(i = 0; i < n_classes; i++) mean[i] /= n[i];

  free(mask);
  
}

