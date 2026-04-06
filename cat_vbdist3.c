/* voxel-based euclidean distance calculation
 * ________________________________________________________________________
 * Calculates the euclidean distance to Ysegment > 0.5 moderated by a second
 * map Yspeed with values between 0 and 1, where 0 voxels are not visited at 
 * all. The speedmap allows to quantify asymetric conditions, e.g., for 
 * sulcal gaps.
 * The closest voxel is determined by a quasi eikonal time distance.
 * Unvisited points were set to Infinity. 
 * 
 *  [Ydistance,Yindex,Ytime] = vbdist3(Ysegment[,Yspeed[,vx_vol]])
 *  
 *  Ysegment (numeric/logical) input image with zero for non elements 
 *              values are converted to single and clamped to [0,1]
 *  Yspeed   (numeric/logical) speed map or mask to limit the distance calculation roughly.
 *              Values are used for quasi-Eikonal time propagation.
 *  Ydistance (single)  euclidean distance image to the closest object point
 *  Yindex    (uint32)  index of nearest point
 *  Ytime     (single)  time-distance image using the speed map
 *
 * Examples: 
 * (1) distance from two points with a simple mask 
 *   A = zeros(50,50,3,'single'); A(15,25,2)=1; A(35,25,2)=1; 
 *   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
 *   D = cat_vbdist(A,B); 
 *   ds('d2sm','',1,A/3+B/3+1/3,D/20,2);
 *
 * (2) not working mask definition
 *   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
 *   B(10:end-10,10:end-10,:) = false; 
 *   D = cat_vbdist(A,B); 
 *   ds('d2sm','',1,A/3+B/3+1/3,D/20,2);
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
/* #include "matrix.h" */
#include "float.h"
#include <stdio.h>
#include <stdlib.h>


#ifdef _MSC_VER
  #define FINFINITY (FLT_MAX+FLT_MAX);
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FINFINITY 1.0f/0.0f;
  #define FNAN 0.0f/0.0f
#endif

/* estimate minimum of A and its index in A */
void pmin(float A[], int sA, float *minimum, int *index)
{
  *minimum = FLT_MAX; 
  *index = 0; 
  for(int i=0; i<sA; i++) {
    if ((A[i]>0.0f) && (*minimum>A[i]))
    { 
      *minimum = A[i]; 
      *index   = i;
    }
  }
}

/* estimate x,y,z position of index i in an array size sx and sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy) {
/* integer division rounds torwards zero */
  *z = i / sxy + 1;
   i = i % sxy;
  *y = i / sy + 1;        
  *x = i % sy + 1;
}

/* copy numeric or logical input to single values clamped between 0 and 1 */
void copyToSingle(const mxArray *A, float *B, int n) {
  if (mxIsComplex(A)) mexErrMsgTxt("ERROR:cat_vbdist: input must be real.\n");
  switch (mxGetClassID(A)) {
    case mxLOGICAL_CLASS: {
      mxLogical *src = mxGetLogicals(A);
      for (int i=0; i<n; i++) B[i] = src[i] ? 1.0f : 0.0f;
      return;
    }
    case mxSINGLE_CLASS: {
      float *src = (float*)mxGetData(A);
      for (int i=0; i<n; i++) {
        float v = src[i];
        B[i] = (v < 0.0f) ? 0.0f : ((v > 1.0f) ? 1.0f : v);
      }
      return;
    }
    case mxDOUBLE_CLASS: {
      double *src = mxGetPr(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxINT8_CLASS: {
      signed char *src = (signed char*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxUINT8_CLASS: {
      unsigned char *src = (unsigned char*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxINT16_CLASS: {
      short *src = (short*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxUINT16_CLASS: {
      unsigned short *src = (unsigned short*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxINT32_CLASS: {
      int *src = (int*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxUINT32_CLASS: {
      unsigned int *src = (unsigned int*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxINT64_CLASS: {
      long long *src = (long long*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    case mxUINT64_CLASS: {
      unsigned long long *src = (unsigned long long*)mxGetData(A);
      for (int i=0; i<n; i++) {
        double v = src[i];
        B[i] = (float)((v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v));
      }
      return;
    }
    default:
      mexErrMsgTxt("ERROR:cat_vbdist: input must be a real numeric or logical matrix.\n");
  }
}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1) mexErrMsgTxt("ERROR:cat_vbdist: not enough input elements\n");
  if (nlhs>3) mexErrMsgTxt("ERROR:cat_vbdist: too many output elements\n");
  if (!(mxIsNumeric(prhs[0]) || mxIsLogical(prhs[0]))) mexErrMsgTxt("ERROR:cat_vbdist: first  input (object) must be a real numeric or logical matrix\n");
  if (nrhs==2 && !(mxIsNumeric(prhs[1]) || mxIsLogical(prhs[1]))) mexErrMsgTxt("ERROR:cat_vbdist: second input (mask) must be a real numeric or logical matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0) mexErrMsgTxt("ERROR:cat_vbdist: third input (vx_vol) must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) mexErrMsgTxt("ERROR:cat_vbdist: third input (vx_vol) must have 3 Elements"); 
  if ( mxGetNumberOfElements(prhs[0])>=2147483647 ) mexErrMsgTxt("ERROR:cat_vbdist: Matrix size is too big for currently used integer datatype.\n");
  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sizeYsegment = mxGetDimensions(prhs[0]);
  const int     dimY = mxGetNumberOfDimensions(prhs[0]);
  const int     numY = mxGetNumberOfElements(prhs[0]);
  if (nrhs==2 && mxGetNumberOfElements(prhs[1]) != (size_t)numY) mexErrMsgTxt("ERROR:cat_vbdist: second input (mask) must have same size as first input\n");
  const int     sizeX  = (int)sizeYsegment[0];
  const int     sizeY  = (int)sizeYsegment[1];
  const int     sizeXY = sizeX*sizeY;
  
  /* get parameters or set defaults */
  static const double default_vx_vol[3] = {1.0, 1.0, 1.0};
  const double *vx_vol;
  if (nrhs<3) {
    vx_vol = default_vx_vol;
  } 
  else {
    vx_vol = mxGetPr(prhs[2]);
  }

  /* define the different neibhor distances*/
  float vx1 = (float)fabs(vx_vol[0]), vx2 = (float)fabs(vx_vol[1]), vx3 = (float)fabs(vx_vol[2]);
  const float   vx12  = sqrtf(  vx1*vx1  + vx2*vx2); /* xy - voxel size */
  const float   vx13  = sqrtf(  vx1*vx1  + vx3*vx3); /* xz - voxel size */
  const float   vx23  = sqrtf(  vx2*vx2  + vx3*vx3); /* yz - voxel size */
  const float   vx123 = sqrtf( vx12*vx12 + vx3*vx3); /* xyz - voxel size */
  
  /* indices of the neighbor relativeNeighborIndexes (index distance) and euclidean distance relativeNeighborDistances */
  const int   numNeighbors = 14;
  const int   relativeNeighborIndexes[14] = {  0, -1,-sizeX+1, -sizeX,-sizeX-1,  -sizeXY+1,-sizeXY,-sizeXY-1,  -sizeXY+sizeX+1,-sizeXY+sizeX,-sizeXY+sizeX-1,  -sizeXY-sizeX+1,-sizeXY-sizeX,-sizeXY-sizeX-1};  
  const float relativeNeighborDistances[14] = {0.0, vx1, vx12, vx2, vx12,    vx13, vx3,  vx13,     vx123,  vx23,   vx123,     vx123,  vx23,   vx123};
  float       neighborDistances[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float       neighborTimes[14] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  float       maxNeighborDistance = FLT_MAX;
  float       maxNeighborTime = FLT_MAX;
  int         maxNeighborIndices = 0;
  int         maxNeighborTimeIndex = 0;
  
  /* get data */
  float *Ysegment = (float *)mxCalloc((size_t)numY, sizeof(float));
  float *Yspeed   = (float *)mxCalloc((size_t)numY, sizeof(float));
  copyToSingle(prhs[0], Ysegment, numY);
  if (nrhs>1) {
    copyToSingle(prhs[1], Yspeed, numY);
  } else {
    for (int i=0; i<numY; i++) {
      Yspeed[i] = 1.0f;
    }
  }

  /* data output depending on nlhs */
  plhs[0] = mxCreateNumericArray(dimY,sizeYsegment,mxSINGLE_CLASS,mxREAL);
  if (nlhs>1) plhs[1] = mxCreateNumericArray(dimY,sizeYsegment,mxUINT32_CLASS,mxREAL);
  if (nlhs>2) plhs[2] = mxCreateNumericArray(dimY,sizeYsegment,mxSINGLE_CLASS,mxREAL);

  float *Ydistance = (float *)mxGetPr(plhs[0]);
  unsigned int *Yindex;
  float *Ytime;

  /* prepare additional output parameter with the index used for the assignment  
   * of the closes voxel defined by the speedmap and the time distance to this 
   * voxel. 
   */
  if (nlhs>1) {
    Yindex = (unsigned int  *)mxGetPr(plhs[1]);
  } 
  else {
    Yindex = (unsigned int *)mxCalloc((size_t)numY, sizeof(unsigned int));
  }
  if (nlhs>2) {
    Ytime = (float *)mxGetPr(plhs[2]);
  } 
  else {
    Ytime = (float *)mxCalloc((size_t)numY, sizeof(float));
  }
  
  /* initialize the index, distance and time maps */
  for (int i=0;i<numY;i++) {
    Yindex[i] = (unsigned int)i;
    if (Ysegment[i] >= 0.5f) {
      Ydistance[i] = 0.0f;
      Ytime[i] = 0.0f;
    } 
    else {
      Ydistance[i] = FLT_MAX;
      Ytime[i] = FLT_MAX;
    }
  }
  
  /* local coordinates of a voxel i (u,v,w = x,y,z) and the temporar neigbhor
   * displacement nu,nv,nw.
   */
  int u,v,w,nu,nv,nw,ni = 0; 
  /* forward direction that consider all points smaller than i */
  for (int i=0;i<numY;i++) {
    if ( (Ydistance[i]>0) && (Yspeed[i] > 0.0f) ) {
      ind2sub(i, &u, &v, &w, sizeXY, sizeX);
      
      /* read neighbor values */
      for (int n=0;n<numNeighbors;n++) {
        ni = i + relativeNeighborIndexes[n];
        ind2sub(ni, &nu, &nv, &nw, sizeXY, sizeX);
        if ( (ni<0) || (ni>=numY) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        neighborDistances[n] = Ydistance[ni] + relativeNeighborDistances[n];
        float speed = 0.5f * (Yspeed[i] + Yspeed[ni]);
        neighborTimes[n] = (speed > 0.0f) ? Ytime[ni] + relativeNeighborDistances[n] / speed : FLT_MAX;
      }

      /* find minimum distance/time within the neighborhood */
      pmin(neighborDistances, numNeighbors, &maxNeighborDistance, &maxNeighborIndices);
      pmin(neighborTimes, numNeighbors, &maxNeighborTime, &maxNeighborTimeIndex);

      /* update values by using the time-index */
      if (maxNeighborIndices > 0) {
        Yindex[i] = Yindex[i + relativeNeighborIndexes[ maxNeighborTimeIndex ] ];
        ind2sub((int)Yindex[i], &nu, &nv, &nw, sizeXY, sizeX); 
        Ydistance[i] = sqrtf(powf((u-nu)*vx1,2) + powf((v-nv)*vx2,2) + powf((w-nw)*vx3,2));
      }
      if (maxNeighborTimeIndex>0 && maxNeighborTime < Ytime[i]) {
        Ytime[i] = maxNeighborTime;
      }
    }
  }
    
  /* backward direction that consider all points larger than i */
  for (int i=numY-1;i>=0;i--) {
    if ( (Ydistance[i]>0) && (Yspeed[i] > 0.0f) ) {
              
      ind2sub(i, &u, &v, &w, sizeXY, sizeX);

      /* read neighbor values */
      for (int n=0;n<numNeighbors;n++) {
        ni = i - relativeNeighborIndexes[n];
        ind2sub((int)ni, &nu, &nv, &nw, sizeXY, sizeX);
        if ( (ni<0) || (ni>=numY) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        neighborDistances[n] = Ydistance[ni] + relativeNeighborDistances[n];
        float speed = 0.5f * (Yspeed[i] + Yspeed[ni]);
        neighborTimes[n] = (speed > 0.0f) ? Ytime[ni] + relativeNeighborDistances[n] / speed : FLT_MAX;
      }

      /* find minimum distance/time within the neighborhood */
      pmin(neighborDistances, numNeighbors, &maxNeighborDistance, &maxNeighborIndices);
      pmin(neighborTimes, numNeighbors, &maxNeighborTime, &maxNeighborTimeIndex);

      /* update values by using the time-index */
      if (maxNeighborIndices > 0) {
        Yindex[i] = Yindex[i - relativeNeighborIndexes[ maxNeighborTimeIndex ]];
        ind2sub((int)Yindex[i], &nu, &nv, &nw, sizeXY, sizeX); 
        Ydistance[i] = sqrtf(powf((u-nu)*vx1,2.0f) + powf((v-nv)*vx2,2.0f) + powf((w-nw)*vx3,2.0f));
      }
      if (maxNeighborTimeIndex>0 && maxNeighborTime < Ytime[i]) {
        Ytime[i] = maxNeighborTime;
      }
    }
  }

  /* final settings */
  for (int i=0;i<numY;i++) {
    if ( Ydistance[i] == FLT_MAX ) {
      Ydistance[i] = INFINITY; /* set unvisited points to infinity */
      Yindex[i] = (unsigned int)i; /* set index to itself for unvisited points */
      Ytime[i] = INFINITY; /* set unvisited points to infinity */
    }
  }

  /* free internal input variables */  
  mxFree(Ysegment);
  mxFree(Yspeed);

  /* set outputs and correct index values for matlab */
  if (nlhs<2) mxFree(Yindex); else { for (int i=0;i<numY;i++) Yindex[i]=Yindex[i]+1; } 
  if (nlhs<3) mxFree(Ytime);

}