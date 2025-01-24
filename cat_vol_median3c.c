/* Median Filter
 * ________________________________________________________________________
 * Median Filter for 3d single image D. Bi is used to mask voxels for the 
 * filter process, whereas Bn is used to mask voxels that are used as 
 * neighbors in the filter process. 
 *
 *  This is a special subversion to filter label maps!
 *
 *  M = median3c(D[,Bi,Bn,iter,nb])
 *
 *  D      (single)   3d matrix for filter process 
 *  Bi     (logical)  3d matrix that mark voxel that should be filtered
 *  Bn     (logical)  3d matrix that mark voxel that are used as neighbors 
 *  iter   (double)   number of interations (<=10, default=1)
 *  nb     (double)   number of neighbors (<=10, default=1); 
 *
 *  Examples: 
 *   1)
 *     A = round(smooth3(rand(50,50,3,'single')*3));
 *     B = false(size(A)); B(5:end-4,5:end-4,:)=true; 
 *     C = cat_vol_median3c(A,B,B,2,2); 
 *     ds('d2smns','',1,A+B,C,2);
 *
 * ______________________________________________________________________
 *
 * Christian Gaser, Robert Dahnke
 * Structural Brain Mapping Group (https://neuro-jena.github.io)
 * Departments of Neurology and Psychiatry
 * Jena University Hospital
 * ______________________________________________________________________
 * $Id$ 
 */

#include "mex.h"   
#include "math.h"
#include "float.h"

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

float abs2(float n) {	if (n<0) return -n; else return n; }        

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<1) mexErrMsgTxt("ERROR:median3c: not enough input elements\n");
  if (nrhs>5) mexErrMsgTxt("ERROR:median3c: too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:median3c: not enough output elements\n");
  if (nlhs>1) mexErrMsgTxt("ERROR:median3c: too many output elements\n");

  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int) mxGetNumberOfElements(prhs[0]);

  int it,iter,nb; 
  if ( dL != 3 || mxIsSingle(prhs[0])==0)        mexErrMsgTxt("ERROR:median3c: first input must be a single 3d matrix\n");
  if ( nrhs>1) {
    if ( mxGetNumberOfDimensions(prhs[1]) != 3 ) mexErrMsgTxt("ERROR:median3c: second input must be 3d - to use a later parameter use ''true(size( input1 ))''\n");
    if (mxIsLogical(prhs[1])==0)                 mexErrMsgTxt("ERROR:median3c: second input must be a logical 3d matrix\n");
  } else if ( nrhs>2) {
    if ( mxGetNumberOfDimensions(prhs[2]) != 3 ) mexErrMsgTxt("ERROR:median3c: third input must be 3d - to use a later parameter use ''true(size( input1 ))'\n");
    if ( mxIsLogical(prhs[2])==0)                mexErrMsgTxt("ERROR:median3c: third input must be a logical 3d matrix\n"); 
  }
  if (nrhs<4)  iter    = 1; else {double *diter    = (double *) mxGetPr(prhs[3]); iter = (int) fmin(10,fmax(1,round( diter[0] )));}
  if (nrhs<5)  nb      = 1; else {double *dnb      = (double *) mxGetPr(prhs[4]); nb   = (int) fmin(10,fmax(1,round( dnb[0]   )));}
  

  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  int NV[256],ind,ni,n;
  bool *Bi, *Bn;
        
  /* in- and output */
  float *D = (float *) mxGetPr(prhs[0]);
  if (nrhs>1)  Bi = (bool *) mxGetPr(prhs[1]); 
  if (nrhs>2)  Bn = (bool *) mxGetPr(prhs[2]); 

  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *M  = (float *) mxGetPr(plhs[0]);
  float *M2 = (float *) mxGetPr(plhs[1]);

  /* for interation */
  for (long ix=0;ix<nL;ix++) M2[ix]=D[ix];

  /* filter process */
  for (int it=1;it<=iter;it++) {
    
    for (int z=0;z<sL[2];z++) for (int y=0;y<sL[1];y++) for (int x=0;x<sL[0];x++) {
      ind = index(x,y,z,sL);
      if ((nrhs==1 || (nrhs>1 && Bi[ind])) && M2[ind]>=0 && !mxIsNaN(M2[ind]) && !mxIsInf(M2[ind])) {
        n = 0;
        /* go through all elements in a 3x3x3 box */
        for (int i=0;i<255;i++) NV[i]=0;  
        for (int i=-nb;i<=nb;i++) for (int j=-nb;j<=nb;j++) for (int k=-nb;k<=nb;k++) {
          /* check borders */ 
          if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2])) {
            ni = index(x+i,y+j,z+k,sL);
            /* check masks and NaN or Infinities */
            if ((nrhs>2 && Bn[ni]==0) || M2[ni]<0 || mxIsNaN(M2[ni]) || mxIsInf(M2[ni]) ) ni = ind;
            NV[(int) M2[ni]]++;
          }
        }
        /* find maximum occurent value */
        M[ind] = 0.0; n=0; for (int i=0;i<255;i++) {if ((NV[i]>=0) && (NV[i]>n)) {n=NV[i]; M[ind]=(float) i;}}
       }
      else {
        M[ind] = M2[ind];
      }
    }

    /* final update for iteration */
    if (it<iter)
      for (long ix=0;ix<nL;ix++) M2[ix]=M[ix];
  }
}


