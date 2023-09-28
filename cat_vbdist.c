/* voxel-based euclidean distance calculation
 * ________________________________________________________________________
 * Calculates the euclidean distance without PVE to an object in P with a 
 * boundary of 0.5 for all voxels within a given mask M that should define
 * a convex hull with direct connection between object and estimation 
 * voxels. Unvisited points were set to FLT_MAX. 
 * 
 *  [D,I,L] = vbdist(P[,M])
 *  
 *  P (single)  input image with zero for non elements 
 *  M (logical) mask to limit the distance calculation roughly, e.g., to 
 *              the brain defined by a convex hull
 *              WARNING: Voxels have to see objects within the mask!
 *  D (single)  distance image
 *  L (uint8)   label map
 *  I (uint32)  index of nearest point
 *
 * Examples: 
 *  % (1) distance from two points with a simple mask 
 *   A = zeros(50,50,3,'single'); A(15,25,2)=1; A(35,25,2)=1; 
 *   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
 *   D = cat_vbdist(A,B); 
 *   ds('d2sm','',1,A/3+B/3+1/3,D/20,2);
 *
 *  % (2) not working mask definition
 *   B = false(size(A)); B(5:end-5,5:end-5,:) = true; 
 *    B(10:end-10,10:end-10,:) = false; 
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
  *minimum = FLT_MAX; *index = 0; /* printf("%d ",sizeof(A)/8); */
  for(int i=0; i<sA; i++) {
    if ((A[i]>0.0) && (*minimum>A[i]))
    { 
      *minimum = A[i]; 
      *index   = i;
    }
  }
}


/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1)                                       mexErrMsgTxt("ERROR:cat_vbdist: not enough input elements\n");
  if (nlhs>4)                                       mexErrMsgTxt("ERROR:cat_vbdist: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:cat_vbdist: first  input (object) must be an 3d single matrix\n");
  if (nrhs==2 && mxIsLogical(prhs[1])==0)           mexErrMsgTxt("ERROR:cat_vbdist: second input (mask) must be an 3d logical matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0)            mexErrMsgTxt("ERROR:cat_vbdist: third input (vx_vol) must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) mexErrMsgTxt("ERROR:cat_vbdist: third input (vx_vol) must have 3 Elements"); 
  if (nrhs==4 && mxIsDouble(prhs[3])==0)            mexErrMsgTxt("ERROR:cat_vbdist: fourth input (debug) must be one double value\n");
  if ( mxGetNumberOfElements(prhs[0])>=2147483647 ) mexErrMsgTxt("ERROR:cat_vbdist: Matrix size is too big for currently used integer datatype.\n");
  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;
  
  /* get parameters or set defaults */
  bool debug; 
  const mwSize sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double  *S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1.0; S[1]=1.0; S[2]=1.0;} else {S = mxGetPr(prhs[2]);}
  if (nrhs<4) {debug=false;} else {double*debugd = (double *) mxGetPr(prhs[3]); debug= (*debugd) > 0;}

  float s1 = fabs((float)S[0]),s2 = fabs((float)S[1]),s3 = fabs((float)S[2]);
  const float   s12  = (float) sqrt( (double)  s1*s1  + s2*s2); /* xy - voxel size */
  const float   s13  = (float) sqrt( (double)  s1*s1  + s3*s3); /* xz - voxel size */
  const float   s23  = (float) sqrt( (double)  s2*s2  + s3*s3); /* yz - voxel size */
  const float   s123 = (float) sqrt( (double) s12*s12 + s3*s3); /* xyz - voxel size */
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[] = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[] = {0.0, s1, s12, s2, s12,    s13, s3,  s13,     s123,  s23,   s123,     s123,  s23,   s123};
  const int   sN = sizeof(NI)/4; /* division by 4 to get from the number of bytes to the number of elements */ 
  float       DN[sN];
  float       DNm = FLT_MAX;
  int         DNi;
  
  /* data intern */
  mxArray *Y[2];
  Y[0] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);
  Y[1] = mxCreateNumericArray(dL,sL,mxUINT8_CLASS,mxREAL);
  unsigned int  *I = (unsigned int  *)mxGetPr(Y[0]);
  unsigned char *L = (unsigned char *)mxGetPr(Y[1]);
  
  /* data output depending on nlhs */
              plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  if (nlhs>1) plhs[1] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);
  if (nlhs>2) plhs[2] = mxCreateNumericArray(dL,sL,mxUINT8_CLASS,mxREAL);
  
  float *D; unsigned int * IO; unsigned char * LO; 
  D = (float *)mxGetPr(plhs[0]);
  if (nlhs>1) IO = (unsigned int  *)mxGetPr(plhs[1]);
  if (nlhs>2) LO = (unsigned char *)mxGetPr(plhs[2]);
  
  /* get data */
  float *V = (float*)mxGetPr(prhs[0]);
  bool  *R; if (nrhs>1) R=(bool *)mxGetPr(prhs[1]); 
  bool e255 = false; 
 
  
  /* intitialisation */
  for (int i=0;i<nL;i++) 
  {
    L[i]=(unsigned char) round(V[i]);
    I[i]=(unsigned int)i;
  
    /* FLT_MAX should be better FINFINITY but this will cause errors in different CAT functions */
    /* such as cat_run_job.APRGs:196 where the old call "Ymx  = Ybgd./(Ybd+Ybgd);" was working but 
     * the new version caused a treshold in an image with INF that lead to an empty brainmask 
     * Although, this can be solved by "Ymx  = min(1,Ybgd./(Ybd+Ybgd));" I do not know other 
     * problematic call now and this should maybe changed later. 
     */
    if (L[i]>=0.5) D[i]=0.0; else D[i]=FLT_MAX; 
    if (L[i]>255.0)  
    {
      if (e255==false) 
      {
        printf("Warning: First parameter of vbdist > 255!\n"); 
        e255 = true;
      }
      L[i] = 255;
    }
  }

  if (debug) {
    printf(" Debug:  debug = %d\n",debug); 
    printf(" Size: nL=%d, sN=%d\n",nL,sN); 
    printf(" Neighbour distances: x=%1.2f,y=%1.2f,z=%1.2f - xy=%1.2f,yz=%1.2f,xz=%1.2f - xyz=%1.2f\n",s1,s2,s3,s12,s23,s13,s123); 
  }
  
  int u,v,w,nu,nv,nw; 
  /* forward direction that consider all points smaller than i */
  for (int i=0;i<nL;i++) 
  {
    if ( (D[i]>0) && (nrhs==1 || (nrhs>1 && R[i]==true) ) )
    {
      ind2sub(i,&u,&v,&w,xy,x);
      
      /* read neighbor values */
      for (int n=0;n<sN;n++)
      {
        int ni = i + NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i+NI[DNi]];
        I[i] = (unsigned int) I[i+NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = (float)sqrt(pow((double)(u-nu)*s1,2) + pow((double)(v-nv)*s2,2) + pow((double)(w-nw)*s3,2));
      }
    }
  }
    
  /* backward direction that consider all points larger than i */
  for (int i=nL-1;i>=0;i--) 
  {
    if ( (D[i]>0) && (nrhs==1 || (nrhs>1 && R[i]==true) ) )
    {
              
      ind2sub(i,&u,&v,&w,xy,x);

      /* read neighbor values */
      for (int n=0;n<sN;n++)
      {
        int ni = i - NI[n];
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) ) ni=i;
        DN[n] = D[ni] + ND[n];
      }

      /* find minimum distance within the neighborhood */
      pmin(DN,sN,&DNm,&DNi);

      /* update values */
      if (DNi>0) {
        L[i] = L[i-NI[DNi]];
        I[i] = (unsigned int)  I[i-NI[DNi]];
        D[i] = DNm; 
        ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
        D[i] = (float)sqrt(pow((double)(u-nu)*s1,2) + pow((double)(v-nv)*s2,2) + pow((double)(w-nw)*s3,2));
      }
    }
  }

  /* clear internal variables */
  /*
  mxDestroyArray(plhs[1]);
  mxDestroyArray(plhs[2]);
   */
  
  /* set outputs and correct index values */
  if (nlhs>1) for (int i=0;i<nL;i++) IO[i]=I[i]+1;
  if (nlhs>2) for (int i=0;i<nL;i++) LO[i]=L[i];

}


