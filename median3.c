/*
 * Robert Dahnke
 * $Id$ 
 *
 */

/* Median Filter
 * ________________________________________________________________________
 * Median Filter for 3d single image D. Bi is used to mask voxels for the 
 * filter process, whereas Bn is used to mask voxels that are used as 
 * neighbors in the filter process. Both mask can changed by intensity 
 * threshold (Bi_low,Bi_high,Bn_low,Bn_high) for D.
 *
 *  M = median3(D[,Bi,Bn,sf,Bi_low,Bi_high,Bn_low,Bn_high])
 *
 *  D      (single)   3d matrix for filter process 
 *  Bi     (logical)  3d matrix that mark voxel that should be filtered
 *  Bn     (logical)  3d matrix that mark voxel that are used as neighbors 
 *  sf     (double)   threshold that is used to filter the result
 *                      sf=0: no filter
 *                      sf<0: only smaller changes
 *                      sf>0: only bigger changes
 *  Bi_low  (double)  low  threshold in D for filtering (add to Bi)
 *  Bi_high (double)  high threshold in D for filtering (add to Bi)
 *  Bn_low  (double)  low  threshold in D for neighbors (add to Bn)
 *  Bn_high (double)  high threshold in D for neighbors (add to Bn)
 *
 * Used slower quicksort for median calculation, because the faster median 
 * of the median application implemenation leads to wrong results. 
 *
 * TODO: check all input elements... 
 * ________________________________________________________________________
 * Robert Dahnke 2011_01
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) + 1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) + 1;        
  *x = i % sy + 1;
}

/* qicksort */
void swap(float *a, float *b)
{
  float t=*a; *a=*b; *b=t;
}

void sort(float arr[], int beg, int end)
{
  if (end > beg + 1)
  {
    float piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    sort(arr, beg, l);
    sort(arr, r, end);
  }
}

float abs2(float n) {	if (n<0) return -n; else return n; }        

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<1) mexErrMsgTxt("ERROR:median3: not enought input elements\n");
  if (nrhs>8) mexErrMsgTxt("ERROR:median3: to many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:median3: to less output elements\n");
  if (nlhs>1) mexErrMsgTxt("ERROR:median3: to many output elements\n");
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int) sL[0];
  const int     y  = (int) sL[1];
  const int     xy = x*y;
  
  if ( dL != 3 || mxIsSingle(prhs[0])==0)                mexErrMsgTxt("ERROR:median3: first input must be a single 3d matrix\n");
  if ( nrhs>1 && mxGetNumberOfDimensions(prhs[1]) != 3 ) mexErrMsgTxt("ERROR:median3: second input must be 3d - to use a later parameter use ''true(size( input1 ))''\n");
  if ( nrhs>2 && mxGetNumberOfDimensions(prhs[1]) != 3 ) mexErrMsgTxt("ERROR:median3: third input must be 3d - to use a later parameter use ''true(size( input1 ))'\n");
  if ( nrhs>1 && mxIsLogical(prhs[1])==0)                mexErrMsgTxt("ERROR:median3: second input must be a logical 3d matrix\n");
  if ( nrhs>2 && mxIsLogical(prhs[2])==0)                mexErrMsgTxt("ERROR:median3: third input must be a logical 3d matrix\n"); 

  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = {0,  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  float NV[27], bi, bn, sf, bil, bih, bnl, bnh;
  int u,v,w,nu,nv,nw,ni,i,n; 
  bool *Bi, *Bn;
  
  /* in- and output */
  float*D = (float *)mxGetPr(prhs[0]);
  if (nrhs>=2) Bi  = (bool *) mxGetPr(prhs[1]); 
  if (nrhs>=3) Bn =  (bool *) mxGetPr(prhs[2]); 
  if (nrhs<4) sf  = 0; 
  else        sf  = (float) *mxGetPr(prhs[3]);

  if (nrhs<5) bil = -FLT_MAX;   
  else        bil = (float) *mxGetPr(prhs[4]);
  if (nrhs<6) bih =  FLT_MAX;   
  else        bih = (float) *mxGetPr(prhs[5]);
  if (nrhs<7) bnl = -FLT_MAX;   
  else        bnl = (float) *mxGetPr(prhs[6]);  
  if (nrhs<8) bnh =  FLT_MAX;   
  else        bnh = (float) *mxGetPr(prhs[7]);

  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float *M = (float *)mxGetPr(plhs[0]);
  
  /* filter process */
  for (i=0;i<nL;i++) {
    if ((nrhs==1 || (nrhs>=2 && Bi[i])) && D[i]>=bil && D[i]<=bih) {
      ind2sub(i,&u,&v,&w,xy,x);
      for (n=0;n<=26;n++) {
        ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
#if defined(_WIN32)
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || 
              (nrhs>=3 && Bn[ni]==0) || D[ni]<bnl ||  D[ni]>bnh || _isnan(D[ni]) || D[ni]==FLT_MAX || D[i]==-FLT_MAX ) ni=i;
#else
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || 
              (nrhs>=3 && Bn[ni]==0) || D[ni]<bnl ||  D[ni]>bnh || isnan(D[ni]) || D[ni]==FLT_MAX || D[i]==-FLT_MAX ) ni=i;
#endif
        NV[n]=D[ni];  
      }
      sort(NV,0,26);
      M[i]=NV[13];
    }
    else {
      M[i]=D[i];
    }
  }

  /* selective filter settings - only big changes (only change extremly noisy data) */
  if (sf>0) {
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (abs2(D[i]-M[i])<sf) ) M[i]=D[i];
    }
  }
  /* selective filter settings - only small changes */
  if (sf<0) { 
    for (i=0;i<nL;i++) {
      if ( (nrhs>=2 && Bi[i]) && D[i]>bil && D[i]<bih && (abs2(D[i]-M[i])>-sf) ) M[i]=D[i];
    }
  }
 
}


