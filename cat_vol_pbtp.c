/* quasi-euclidean distance calculation
 * _____________________________________________________________________________
 * [GMT,RPM] = cat_vol_pbtp(SEG,WMD,CSFD)
 *
 * SEG  = (single) segment image with low and high boundary bd
 * GMT  = (single) thickness image
 * RPM  = (single) radial position map
 * WMD  = (single) CSF distance map
 * CSFD = (single) CSF distance map
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
 * TODO:
 *  - eikonal distance for subsegmentation (region growing)
 *  - own labeling (
 */

#include "mex.h"   
#include "math.h"
#include <stdlib.h>
/* #include "matrix.h" */

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

/* comparison function for qsort */
int compare_float(const void *a, const void *b) {
  float fa = *(const float*)a;
  float fb = *(const float*)b;
  return (fa > fb) - (fa < fb);
}

/* 
Get all values of the voxels which are in WMD-range (children of this voxel).
The largest values will describe the local thickness of the band connected to the current voxel.
We use the maximum then to estimate the median of the highest neibors. 
However, we also use the siblings of the voxel (with similar WMD) to estimate a mean value.
This mean is then used to have a smoother thickness estimation.
*/
float pmax(const float GMT[], const float RPM[], const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA) {
  float n=0.0f, maximum=WMD, high_avg=0.0f;

  /* the pure maximum */
  /* ... the maximum causes overestimations that helps for gyri but is 
   *     problematic in sulci and in case of blood vessels and meninges */
  if ( 0==1 ) {
    for (int i=0;i<sA;i++) {
      if (  ( GMT[i] < 1e15f ) && ( maximum < GMT[i] ) &&                  /* thickness/WMD of neighbors should be larger */
            ( SEG[i] > 1.0f ) && ( SEGI>1.0f && SEGI<3.0f ) &&           /* projection range */
            ( ( ( RPM[i] - ND[i] * 1.2f ) <= WMD ) ) &&                    /* upper boundary - maximum distance */
            ( ( ( RPM[i] - ND[i] * 0.5f ) >  WMD ) || ( SEG[i]<1.5f ) ) &&  /* lower boundary - minimum distance - corrected values outside */
            ( ( ( (SEGI * max(1.0f,min(1.2f,SEGI-1.5f)) ) >= SEG[i] ) ) || ( SEG[i]<1.5f ) ) )  /* for high values will project data over sulcal gaps */
        { maximum = GMT[i]; n++; } 
    }
  }
  
  /* the median of the highest values*/
  float high_vals[14] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  int count = 0;
  for (int i=0;i<sA;i++) {
    if ( ( GMT[i] < 1e15f ) && ( (maximum + .5 ) < GMT[i] ) && /* the maximum is here the WMD of the voxel itself + .5 */
         ( SEG[i] >= 1.0f ) && ( SEGI>1.0f && SEGI<3.0f ) && 
         ( ( (RPM[i] - ND[i] * 1.2f ) <= WMD ) ) && 
         ( ( (RPM[i] - ND[i] * 0.5f ) >  WMD ) || ( SEG[i]<1.5f ) ) &&
         ( ( ( (SEGI * max(1.0f,min(1.2f,SEGI-1.5f)) ) >= SEG[i] ) ) || ( SEG[i]<1.5f ) ) ) 
      { high_vals[count++] = GMT[i]; } 
  }
  if ( count > 0 ) {
    qsort(high_vals, count, sizeof(float), compare_float);
    /* Use central 50% of sorted values for robust averaging */
    n=0.0f; 
    for (int i=(int)round( (float)count * 0.25f); i<(int)round( (float)count * 0.75f); i++) {
      high_avg = high_avg + high_vals[i];
      n++; 
    }
    if ( n > 0 )
      high_avg /= n;  
    else
      high_avg = maximum; 
  }
  else
    high_avg = maximum; 

  return high_avg;
}



/* estimate x,y,z position of index i in an array size sx and sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy) {
/* integer division rounds torwards zero */
  *z = i / sxy + 1;
   i = i % sxy;
  *y = i / sy + 1;        
  *x = i % sy + 1;
}



/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enough input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: too many input elements.\n");
  if (nlhs<1) mexErrMsgTxt("ERROR: not enough output elements.\n");
  if (nlhs>2) mexErrMsgTxt("ERROR: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first  input must be an 3d single matrix\n");
  
  /* main information about input data (size, dimensions, ...) */
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  if (dL != 3) mexErrMsgTxt("ERROR: first input must be a 3d single matrix\n");
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
  const float   s2 = sqrtf(2.0f);
  const float   s3 = sqrtf(3.0f);
  const int     nr = nrhs;
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   sN = 14;
  const int   NI[14]   = {   0,  -1,-x+1,  -x,-x-1,  -xy+1, -xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[14]   = {0.0f,1.0f,  s2,1.0f,  s2,     s2,1.0f,   s2,       s3,   s2,     s3,       s3,   s2,     s3};
  float       GMTN[14] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  float       WMDN[14] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  float       SEGN[14] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
  float       DNm = 0.0f;
  float       du  = 0.0f, dv  = 0.0f, dw  = 0.0f, dnu = 0.0f, dnv = 0.0f, dnw = 0.0f, d = 0.0f, dcf = 0.0f; 
  float       WMu = 0.0f, WMv = 0.0f, WMw = 0.0f, GMu = 0.0f, GMv = 0.0f, GMw = 0.0f, SEGl = 0.0f, SEGu = 0.0f, tmpfloat = 0.0f;
  int         ni = 0, u = 0, v = 0, w = 0, nu = 0, nv = 0, nw = 0, WMC = 0, CSFC = 0;
    
  /* main volumes */
  mxArray *hlps[1];  
  hlps[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);

  /* input variables */
  float*SEG  = (float *)mxGetPr(prhs[0]);
  if (mxIsSingle(prhs[1])==0 || mxGetNumberOfElements(prhs[1]) != nL) mexErrMsgTxt("ERROR: second input must be a single matrix with same size as first input\n");
  if (mxIsSingle(prhs[2])==0 || mxGetNumberOfElements(prhs[2]) != nL) mexErrMsgTxt("ERROR: third input must be a single matrix with same size as first input\n");
  const float*WMD_in  = (float *)mxGetPr(prhs[1]);
  const float*CSFD_in = (float *)mxGetPr(prhs[2]);
  
  /* create internal copies to preserve input arrays */
  float *WMD  = (float *)mxMalloc(nL * sizeof(float));
  float *CSFD = (float *)mxMalloc(nL * sizeof(float));
  for (int i=0; i<nL; i++) {
    WMD[i]  = WMD_in[i];
    CSFD[i] = CSFD_in[i];
  }

  /* output variables */
  float *GMT  = (float *)mxGetPr(plhs[0]);
  float *RPM  = (float *)mxGetPr(hlps[0]);
           
  
  /* intitialisiation */
  for (int i=0;i<nL;i++) {
    GMT[i] = WMD[i];
    RPM[i] = WMD[i];
    /* proof distance input */
    if ( SEG[i]>2.5f ) WMC++;
    if ( SEG[i]<1.5f ) CSFC++;
  }
  if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
 
  
  
/* thickness calculation
 * ===================================================================== */
  for (int i=0;i<nL;i++) {
    if (SEG[i]>1.0f && SEG[i]<3.0f) {
      ind2sub(i,&u,&v,&w,xy,x);
      
      /* read neighbour values */
      for (int n=0;n<sN;n++) {
        ni = i + NI[n];
        if (ni < 0 || ni >= nL) ni = i;
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        GMTN[n] = GMT[ni]; 
        WMDN[n] = RPM[ni]; 
        SEGN[n] = SEG[ni];
      }

      /* find minimum distance within the neighborhood */
      DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN);
      GMT[i] = DNm;
    }
  }
  
  for (int i=nL-1;i>=0;i--) {
    if (SEG[i]>1.0f && SEG[i]<3.0f) {
      ind2sub(i,&u,&v,&w,xy,x);
      
      /* read neighbour values */
      for (int n=0;n<sN;n++) {
        ni = i - NI[n];
        if (ni < 0 || ni >= nL) ni = i;
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        GMTN[n] = GMT[ni]; 
        WMDN[n] = RPM[ni]; 
        SEGN[n] = SEG[ni];
      }

      /* find minimum distance within the neighborhood */
      DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN);
      if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm;
    }
  }
 

  for (int i=0;i<nL;i++) if (SEG[i]<1.5f || SEG[i]>2.5f ) GMT[i] = 0.0f; 



/* final settings...
 * ===================================================================== */
  float CSFDc = 0.0f, GMTi = 0.0f, CSFDi = 0.0f;
  for (int i=0;i<nL;i++) { 
    GMT[i] = max(0.0f, min(CSFD[i] + WMD[i], GMT[i]));
    if isinf(GMT[i]) GMT[i] = 0.0f; 
    if (SEG[i]>=1.5f && SEG[i]<=2.5f) {
      GMTi   = CSFD[i] + WMD[i];  
      CSFDi  = GMT[i]  - WMD[i];
    
      if ( CSFD[i]>CSFDi )  CSFD[i] = CSFDi;          
      else                  GMT[i]  = GMTi;
    }
  }

 
/* estimate RPM
 * ===================================================================== */
  for (int i=0;i<nL;i++) {
    if ( SEG[i]>=2.5f)   
      RPM[i] = 1.0f; 
    else {
      if ( SEG[i]<=1.5f || GMT[i]==0.0f ) 
        RPM[i] = 0.0f;
      else {
        RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
        if (RPM[i]>1.0f) RPM[i] = 1.0f;
        if (RPM[i]<0.0f) RPM[i] = 0.0f; 
      }
    } 
  }
  
  /* create second output */
  if ( nlhs > 1 ) {
    plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
    float *RPMO = (float *)mxGetPr(plhs[1]);
    for (int i=0;i<nL;i++) RPMO[i] = RPM[i]; 
  }
  
  /* clear internal variables */
  mxFree(WMD);
  mxFree(CSFD);
  //mxDestroyArray(hlps[0]);
}


