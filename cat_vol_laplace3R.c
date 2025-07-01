/* laplace calculation
 * _____________________________________________________________________
 * Filter SEG within the intensity range of low and high until the changes
 * are below TH.   
 * 
 *  L = cat_vol_laplace3R(SEG,R,TH)
 * 
 *  SEG  .. 3D single input matrix
 *  R    .. 3D boolean volume to describe the filter area
 *  TH   .. threshold to control the number of iterations
 *          maximum change of an element after iteration
 *
 * Example: 
 *   A = zeros(50,50,3,'single'); A(10:end-9,10:end-9,2)=0.5; 
 *   A(20:end-19,20:end-19,2)=1;
 *   B = A==0.5; 
 *   C = cat_vol_laplace3R(A,0,1,0.001); ds('d2smns','',1,A,C,2); 
 * ______________________________________________________________________
 *
 * Christian Gaser, Robert Dahnke
 * Structural Brain Mapping Group (https://neuro-jena.github.io)
 * Departments of Neurology and Psychiatry
 * Jena University Hospital
 * ______________________________________________________________________
 * $Id$ 
 */

/*
 * TODO: change of L1 and L2 by pointer
 */

#include "mex.h"   
#include "math.h"
#include "float.h"
/* #include "matrix.h" */

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

float abs2(float n) {	if (n<0) return -n; else return n; }

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<2)                                       mexErrMsgTxt("ERROR:laplace3R: not enough input elements\n");
  if (nrhs>4)                                       mexErrMsgTxt("ERROR:laplace3R: too many input elements\n");
  if (nlhs>1)                                       mexErrMsgTxt("ERROR:laplace3R: too many output elements\n");
  if (nlhs<1)                                       mexErrMsgTxt("ERROR:laplace3R: not enough output elements\n");
  if (mxIsSingle(prhs[0])==0)                       mexErrMsgTxt("ERROR:laplace3R: 1st input must be an 3d single matrix\n");
  if (mxIsLogical(prhs[1])==0)                      mexErrMsgTxt("ERROR:laplace3R: 2nd input must be an 3d logical matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0)            mexErrMsgTxt("ERROR:laplace3R: 3rd input must be an double matrix\n");
  if (nrhs==4 && mxIsDouble(prhs[3])==0)            mexErrMsgTxt("ERROR:laplace3R: 4th input (voxelsize) must be a double matrix\n");
  if (nrhs==4 && mxGetNumberOfElements(prhs[3])!=3) mexErrMsgTxt("ERROR:laplace3R: 4th input (voxelsize) must have 3 Elements");
  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;
  
  /* input data */
  float*SEG = (float *)mxGetPr(prhs[0]);
  bool *M   = (bool  *)mxGetPr(prhs[1]); 
  float TH  = (float) mxGetScalar(prhs[2]); if ( TH>=0.5 || TH<0.000001 ) mexErrMsgTxt("ERROR:laplace3R: threshhold must be >0.000001 and smaller than 0.5\n");
  const mwSize sS[2] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1.0; S[1]=1.0; S[2]=1.0;} else {S = mxGetPr(prhs[2]);}
 
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   sN = 6;
  const int   NI[6]  = { -1, 1, -x, x, -xy, xy};  
  
  /* output data */
  mxArray *hlps[2];
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  hlps[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  hlps[1] = mxCreateLogicalArray(dL,sL);

  float  *L1 = (float *)mxGetPr(plhs[0]);
  float  *L2 = (float *)mxGetPr(hlps[0]);
  bool   *LN = (bool  *)mxGetPr(hlps[1]);
  
  /* intitialisiation */
  for (int i=0;i<nL;i++) {
    if ( mxIsNaN(SEG[i]) ) { L1[i] = FLT_MAX; } else { L1[i] = SEG[i]; }
    L2[i] = L1[i];
    LN[i] = M[i];
  }

  int u,v,w,nu,nv,nw,ni,iter=0,maxiter=2000;
  float Nn, diff, maxdiffi, maxdiff=1.0;
  while ( maxdiff > TH && iter < maxiter) {
    maxdiffi=0; iter++;
    for (int i=0;i<nL;i++) {
      if ( M[i] && LN[i] ) {  
        ind2sub(i,&u,&v,&w,xy,x);

        /* read neighbor values */
        L2[i]=0.0; Nn=0.0;
        for (int n=0;n<sN;n++) {
          ni = i + NI[n]; 
          ind2sub(ni,&nu,&nv,&nw,xy,x);
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX) )==false) 
            {L2[i] = L2[i] + L1[ni]; Nn++;}
        }
        if (Nn>0) {L2[i]/=Nn;} else {L2[i]=L1[i];}
        
        diff  = abs2( L1[i] - L2[i] ); /*printf("%f %f %f\n",L1[i],L2[i],diff); */
        if ( diff>(TH/10.0) ) { 
          for (int n=0;n<sN;n++) {
            ni = i + NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
            if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX) )==false) 
              LN[ni] = true; /* if i change his neigbors has to be recalculated */
          }
        }
      
        LN[i]=false;
        if ( maxdiffi<diff ) maxdiffi=diff; 
      }
    }
    maxdiff = maxdiffi;

    /* update of L1 */
    for (int i=0;i<nL;i++) { L1[i]=L2[i]; }
  }
  /* printf("%d\n",iter); */
  
  /* clear internal variables */
  /*
  mxDestroyArray(hlps[0]);
  mxDestroyArray(hlps[1]);
   */
}


