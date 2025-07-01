/* function [mn,std,min,max,sum,n,median?] = cat_vol_ROIval(Ya,Yv)
 * _____________________________________________________________________
 * Estimation of mean, standard deviation, minimum, maximum, sum, number
 * and median (not yet implemented) of Yv in a ROI described by Ya.
 *
 * Example: 
 *  1) 
 *    A = rand(50,50,3,'single');  
 *    for i=1:size(A,2), A(:,i,:) = A(:,i,:) + (( i/size(A,2) ) - 0.5); end
 *    L = zeros(size(A),'single'); L(round(numel(L)*0.367))=1; 
 *    L(round(numel(L)*0.533))=2; L(round(numel(L)*0.616))=3; 
 *    [D,I] = cat_vbdist(L); L=L(I); L = uint8(round(L));
 *    [MN,STD,MIN,MAX,SM,N,MD] = cat_vol_ROIval(A,L); 
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
#ifdef _MSC_VER
  #define FINFINITY (FLT_MAX+FLT_MAX);
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FINFINITY 1.0f/0.0f;
  #define FNAN 0.0f/0.0f
#endif

  
/* qicksort for median */
/*
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
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* main input handling */
  if (nrhs!=2)                  
    mexErrMsgTxt("ERROR:cat_vol_ROIval: requires 2 maps. A 3d uint8 atlas map and a single value map.\n");
  if (mxIsUint8(prhs[0])==0 || mxGetNumberOfDimensions(prhs[0])!=3 )    
    mexErrMsgTxt("ERROR:cat_vol_ROIval: 1st input must be an 3d uint8 matrix.\n");
  if (mxIsSingle(prhs[1])==0  || mxGetNumberOfDimensions(prhs[0])!=3 )   
    mexErrMsgTxt("ERROR:cat_vol_ROIval: 2nd input must be an 3d single matrix.\n");
  
  /* main information about input data (size, dimensions, ...) */
  const unsigned int n  = mxGetNumberOfElements(prhs[0]);
  const unsigned int n1 = mxGetNumberOfElements(prhs[1]);
  const mwSize *sL0 = mxGetDimensions(prhs[0]); 
  const mwSize *sL1 = mxGetDimensions(prhs[0]);
  if (n!=n1 || sL0[0]!=sL1[0] || sL0[1]!=sL1[1] || sL0[2]!=sL1[2])
    mexErrMsgTxt("ERROR:cat_vol_ROIval: 1nd and 2nd input must have the same size.\n");

  /* input data */
  const unsigned char *Ya = (unsigned char *) mxGetPr(prhs[0]);
  const float         *Yv = (float *) mxGetPr(prhs[1]); 
  
  unsigned int ti, tii, maxYa=0, maxYa2=0;


  for (int i=0;i<n;i++) {if (maxYa<Ya[i]) maxYa=Ya[i];}
  const int sV[2] = {1,256}; 
  const mwSize sVm[2] = {1,2}, sdsv[2] = {1,256}; 
  
  /* output data */
  plhs[0]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  plhs[1]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  plhs[2]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  plhs[3]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  plhs[4]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  plhs[5]    = mxCreateNumericArray(2,sVm,mxSINGLE_CLASS,mxREAL);
  /* plhs[6]    = mxCreateNumericArray(2,sV,mxSINGLE_CLASS,mxREAL); */
  float*rmn  = (float *)mxGetPr(plhs[0]);
  float*rstd = (float *)mxGetPr(plhs[1]);
  float*rmin = (float *)mxGetPr(plhs[2]);
  float*rmax = (float *)mxGetPr(plhs[3]);
  float*rsum = (float *)mxGetPr(plhs[4]); 
  float*rn   = (float *)mxGetPr(plhs[5]);
  /* f
  float*rmed = (float *)mxGetPr(plhs[6]);
  float Yt[16777216]; for (i=0;i<16777216;i++) Yt[i]=0;
  */ 
  
  
  /* initialization of output data with 0 */
  for (int id=0;id<maxYa;id++) {
    rn[id]   = 0.0;
    rmn[id]  = 0.0;
    rstd[id] = 0.0;
    rmin[id] = FLT_MAX;
    rmax[id] = -FLT_MAX;
    rsum[id] = 0.0;
    /*
    rmed[id]=0; 
    */
  }

  /* estimation */
  for (int i=0;i<n;i++) {
    int id= (int) (Ya[i]) - 1; 
    if (id>=0) {
      rn[id]     = rn[id]   + 1; 
      rsum[id]   = rsum[id] + Yv[i];
      if (rmax[id]<Yv[i]) 
        rmax[id] = Yv[i];
      if (rmin[id]>Yv[i]) 
        rmin[id] = Yv[i]; 
    }
  }
  
  return;
  /* correct minimum and maximum for empty ROIs */
  for (int id=0;id<maxYa;id++) {
    if (rn[id]==0) {
      rmax[id]=FNAN;
      rmin[id]=FNAN;
    }
  }

  /* mean */
  for (int id=0;id<maxYa;id++) {
    if (rn[id]>0)
      rmn[id] = rsum[id] / rn[id];
    else
      rmn[id] = FNAN;
  }
  
  /* standard deviation */
  for (int i=0;i<n;i++) {
    int id= (int) (Ya[i]) - 1; 
    if (id>=0 & rn[Ya[i]]>0) {
      rstd[id] += powf(Yv[i]-rmn[Ya[i]-1],2);
    }
  }
  for (int id=0;id<maxYa;id++) {
   if (rn[id]>1)
     rstd[id] = sqrtf( 1.0 / (rn[id] - 1.0) * rstd[id]);
   else {
     if (rn[id]==1)
       rstd[id] = 1.0;
     else
       rstd[id] = FNAN;
    }
  }

 /* median */
 /*
 for (id=0;id<maxYa;id++) {
    if (rn[id]>2) {
      
      ti=0; 
      for (i=0;i<n;i++) {
        if (Ya[i]>0 && id==(Ya[i]-1) && ti<16777216) {
          Yt[ti] = 1; //Yv[i]; 
          ti++;
        }
      }
      //printf("%d-%d-%f\n",ti,(int) rn[id],Yt[ti]);   
      //for (tii=0;tii<ti;tii++) 
      // rmed[id]=rmed[id] + Yt[tii];
        
      //sort(Yt,0,ti);
  
      //rmed[id]=Yt[(unsigned char) floor( ( (double) rn[id] )/2 )]/2 + 
      //         Yt[(unsigned char)  ceil( ( (double) rn[id] )/2 )]/2;
    }
    else
      rmed[id]=rmn[id];
  }
  */
}