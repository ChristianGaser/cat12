/* ______________________________________________________________________
 *
 * Christian Gaser, Robert Dahnke
 * Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
 * Departments of Neurology and Psychiatry
 * Jena University Hospital
 * ______________________________________________________________________
 * $Id$ 
 *
 */

#include "mex.h"
#include "math.h"
#include "stdio.h"

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *voxelsize, const int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label;
  unsigned char *mask;
  double *src, *bias, voxelsize[3], mx, thresh, thresh_kmeans_int, bias_fwhm;
  double *mean, n[100];
  unsigned long nvox, i;
  int n_classes, iters_nu;
  const unsigned long *dims0;
  int dims[3], ind;
  
  if (nrhs<2)
    mexErrMsgTxt("At least 2 inputs required.");
  else if (nlhs>3)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("Image must be double.");

  src = (double*)mxGetPr(prhs[0]);
  n_classes = (int)mxGetScalar(prhs[1]);
  if (nrhs>2)
    iters_nu = (int)mxGetScalar(prhs[2]);
  else iters_nu = 0;
        
  dims0 = mxGetDimensions(prhs[0]);

  for (i=0; i<3; i++) {
    voxelsize[i] = 1.0;
    dims[i] = dims0[i];
  }

  plhs[0] = mxCreateNumericArray(3,dims0,mxUINT8_CLASS,mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, n_classes, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericArray(3,dims0,mxDOUBLE_CLASS,mxREAL);
  label  = (unsigned char *)mxGetPr(plhs[0]);
  mean   = (double *)mxGetPr(plhs[1]);
  bias  =  (double *)mxGetPr(plhs[2]);
  
  thresh = 0;
  thresh_kmeans_int = 128;
  bias_fwhm = 50.0;

  nvox = dims[0]*dims[1]*dims[2];
  mask = (unsigned char *)mxMalloc(sizeof(unsigned char)*nvox);
  if(mask == NULL) {
    mexErrMsgTxt("Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  
  for (i=0; i<nvox; i++)
    bias[i] = src[i];

  for (i=0; i<nvox; i++)
    mask[i] = (src[i]!=0) ? 255 : 0;
  
  /* initial Kmeans estimation with pve classes */
/*  if (iters_nu > 0)
    mx = Kmeans(src, label, mask, 25, n_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, 0, bias_fwhm);
*/  
  mx = Kmeans(src, label, mask, 25, n_classes, voxelsize, dims, thresh, thresh_kmeans_int, iters_nu, 0, bias_fwhm);

  for (i=0; i<nvox; i++)
    bias[i] -= src[i];

  /* calculate mean */
  for(i = 0; i < n_classes; i++) {
    n[i] = 0;
    mean[i] = 0.0;
  }
  
  for(i = 0; i < nvox; i++) {
    if(label[i] == 0) continue;
    ind = label[i]-1;
    n[ind]++;
    mean[ind] += src[i];
  }
  for(i = 0; i < n_classes; i++) mean[i] /= n[i];

  mxFree(mask);
  
}

