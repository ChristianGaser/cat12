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

/*
 * TODO: 
 *  - use structure with defaults for input parameter
 *  - use long rather to indexing ultra-high-resolution data 
 */

#include "mex.h"
#include "math.h"
#include "stdio.h"
#include "Amap.h"
/* #include "matrix.h" */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  unsigned char *label0, *label, *prob, *mask;
  double *src0, *src, *srco, *mean, *fmeans, *fstds, *voxelsize;
  double max_vol = -1e15, weight_MRF, bias_fwhm, offset;
  const mwSize *dims;
  mwSize dims3[4];
  int dims2[4];
  int n_classes, pve, nvox, iters_icm, verb;
  int niters, iters_nu, sub, init, thresh, thresh_kmeans_int;
    
  if (nrhs<11 | nrhs > 12)
    mexErrMsgTxt("11 inputs required: \n  [prob, means, stds, srcb] = cat_amap(src, label, n_classes, n_iters, sub, pve, init, mrf_weight, voxelsize, iters_icm, bias_fwhm, verb)");
  else if (nlhs>4)
    mexErrMsgTxt("Too many output arguments.");
  
  if (!mxIsDouble(prhs[0]))  /* src */
    mexErrMsgTxt("First argument must be double.");
  if (!mxIsUint8(prhs[1]))   /* label */
    mexErrMsgTxt("Second argument must be uint8.");
  if (!mxIsDouble(prhs[2]))  /* n_classes */
    mexErrMsgTxt("Third argument must be double.");
  if (!mxIsDouble(prhs[3]))  /* n_iters */
    mexErrMsgTxt("4th argument must be double.");
  if (!mxIsDouble(prhs[4]))  /* sub */
    mexErrMsgTxt("5th argument must be double.");
  if (!mxIsDouble(prhs[5]))  /* pve */
    mexErrMsgTxt("6th argument must be double.");
  if (!mxIsDouble(prhs[6]))  /* init */
    mexErrMsgTxt("7th argument must be double.");
  if (!mxIsDouble(prhs[7]))  /* mrf_weight */
    mexErrMsgTxt("8th argument must be double.");
  if (!mxIsDouble(prhs[8]))  /* voxelsize */
    mexErrMsgTxt("9th argument must be double.");
  if (nrhs>9  && !mxIsDouble(prhs[9]))  /* iters_icm */
    mexErrMsgTxt("10th argument must be double.");
  if (nrhs>10 && !mxIsDouble(prhs[10])) /* bias_fwhm */
    mexErrMsgTxt("11th argument must be double.");
  if (nrhs>11 && !mxIsDouble(prhs[11])) /* verb */
    mexErrMsgTxt("12th argument must be double.");
  
  src0       = (double*)mxGetPr(prhs[0]);
  label0     = (unsigned char*)mxGetPr(prhs[1]);
  n_classes  = (int)mxGetScalar(prhs[2]);
  niters     = (int)mxGetScalar(prhs[3]);
  sub        = (int)mxGetScalar(prhs[4]);
  pve        = (int)mxGetScalar(prhs[5]);
  init       = (int)mxGetScalar(prhs[6]);
  weight_MRF = (double)mxGetScalar(prhs[7]);
  voxelsize  = (double*)mxGetPr(prhs[8]);
  iters_icm  = (int)mxGetScalar(prhs[9]);
  if (nrhs>10) bias_fwhm  = (double)mxGetScalar(prhs[10]); else bias_fwhm = 60.0;
  if (nrhs>11) verb       = (int)mxGetScalar(prhs[11]);    else verb = 0;
  
  if ( mxGetM(prhs[8])*mxGetN(prhs[8]) != 3) 
    mexErrMsgTxt("Voxelsize should have 3 values.");

  dims = mxGetDimensions(prhs[0]);
  dims2[0] = (int)dims[0]; dims2[1] = (int)dims[1]; dims2[2] = (int)dims[2]; dims2[3] = n_classes;
  
  /* for PVE we need more classes */
  if(pve == 6) dims2[3] += 3;
  if(pve == 5) dims2[3] += 2;

  /* mxCreateNumericArray expects mwSize data type */
  for(int i = 0; i < 4; i++) dims3[i] = (mwSize)dims2[i]; 

  /* final segmentation */
  plhs[0] = mxCreateNumericArray(4, dims3, mxUINT8_CLASS, mxREAL);
  prob    = (unsigned char *)mxGetPr(plhs[0]);

  /* internal mean and std values */
  mxArray *hlps[3]; 
  hlps[0] = mxCreateNumericMatrix(1, n_classes+3, mxDOUBLE_CLASS, mxREAL);  /* old segmentation mean values */
  hlps[1] = mxCreateNumericMatrix(1, n_classes+3, mxDOUBLE_CLASS, mxREAL);  /* new corrected mean values (may equal to old) */
  hlps[2] = mxCreateNumericMatrix(1, n_classes+3, mxDOUBLE_CLASS, mxREAL);  /* new std values */
  mean    = (double *)mxGetPr(hlps[0]); 
  fmeans  = (double *)mxGetPr(hlps[1]);
  fstds   = (double *)mxGetPr(hlps[2]);
  for (int i=0; i<n_classes+3; i++) mean[i]   = 0.0;
  for (int i=0; i<n_classes+3; i++) fmeans[i] = 0.0;
  for (int i=0; i<n_classes+3; i++) fstds[i]  = 0.0;
   /* new dynamic output for segmentation mean and std values */
  double *fmeanso; double *fstdso; 
  if ( nlhs>1 ) { /* means */
    plhs[1] = mxCreateNumericMatrix(1, n_classes, mxDOUBLE_CLASS, mxREAL);
    fmeanso = (double *)mxGetPr(plhs[1]);
    for (int i=0; i<n_classes; i++) fmeanso[i] = fmeans[i];
  }
  if ( nlhs>2 ) { /* stds */
    plhs[2] = mxCreateNumericMatrix(1, n_classes, mxDOUBLE_CLASS, mxREAL);
    fstdso  = (double *)mxGetPr(plhs[2]);
    for (int i=0; i<n_classes; i++) fstdso[i] = fstds[i];
  }
  if ( nlhs>3 ) { /* bias corrected */
    plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    srco  = (double *)mxGetPr(plhs[3]);
    for (int i=0; i<nvox; i++) srco[i] = src[i];
  }
  
  /* internal dublicat of the input - NOT YET WORKING - creating a bug somewhere else but should be correct in general */
  mxArray *prhsi[2]; 
  const mwSize *sL = mxGetDimensions(prhs[0]);
  int     nL = mxGetNumberOfElements(prhs[0]);
  int     dL = mxGetNumberOfDimensions(prhs[0]);
  prhsi[0]   = mxCreateNumericArray(dL, sL, mxDOUBLE_CLASS, mxREAL);
  prhsi[1]   = mxCreateNumericArray(dL, sL, mxUINT8_CLASS , mxREAL);
  src        = (double*)mxGetPr(prhsi[0]);
  label      = (unsigned char*)mxGetPr(prhsi[1]);
  for(long i = 0; i < nL; i++) src[i]   = src0[i]; 
  for(long i = 0; i < nL; i++) label[i] = label0[i]; 
          
  nvox = dims[0]*dims[1]*dims[2];
  
  for(int i = 0; i < nvox; i++) {
    max_vol = MAX(src[i], max_vol);
  }

  offset = 0.2*max_vol;
  
  /* add offset to ensure that CSF values are much larger than background noise */
  for (int i=0; i<nvox; i++) {
    if (label[i] > 0) src[i] += offset;
  }

  /* initial labeling using Kmeans */
  if (init>0) {
    mask = (unsigned char *)mxMalloc(sizeof(unsigned char)*nvox);
    if(mask == NULL) {
      mexErrMsgTxt("Memory allocation error\n");
      exit(EXIT_FAILURE);
    }
    for (int i=0; i<nvox; i++)
      mask[i] = (src[i]>0) ? 255 : 0;

    thresh            = 0;
    thresh_kmeans_int = 128;
    iters_nu          = 0;    /* bias correction works better inside Amap */

    /* 
     * kmeans.c: 
     *   double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters,
     *     double *voxelsize, int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve, double bias_fwhm) 
     */
  
    /* initial Kmeans estimation with 6 classes */
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, thresh, thresh_kmeans_int, iters_nu, KMEANS, bias_fwhm);
    /* final Kmeans estimation with 3 classes */
    max_vol = Kmeans( src, label, mask, 25, n_classes, voxelsize, dims2, thresh, thresh_kmeans_int, iters_nu, NOPVE, bias_fwhm);

    mxFree(mask);
  }
  

  /* 
   * Amap.c: 
   *   void Amap(double *src, unsigned char *label, unsigned char *prob, double *mean, int n_classes, int niters, 
   *     int sub, int *dims, int pve, double weight_MRF, double *voxelsize, int niters_ICM, double offset, double bias_fwhm, 
   *     double *fmeans, double *fstd)
   */
  Amap(src, label, prob, mean, n_classes, niters, sub, dims2, pve, weight_MRF, voxelsize, iters_icm, offset, bias_fwhm, verb, fmeans, fstds);
  

  /* Pve.c:
   *   void Pve5(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
   *   void Pve6(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
   */
  if(pve==6) Pve6(src, prob, label, mean, dims2);
  if(pve==5) Pve5(src, prob, label, mean, dims2);
  
  
  /* new dynamic output for segmentation mean and std values */
  if ( nlhs>1 ) {
    for (int i=0; i<n_classes; i++) fmeanso[i] = fmeans[i];
  }
  if ( nlhs>2 ) {
    for (int i=0; i<n_classes; i++) fstdso[i] = fstds[i];
  }
  if ( nlhs>3 ) {
    for (int i=0; i<nvox; i++) srco[i] = src[i] - offset;
  }
  
  /* clear internal variables */ 
  /*
  mxDestroyArray(prhsi[0]);
  mxDestroyArray(prhsi[1]);
   */
}

