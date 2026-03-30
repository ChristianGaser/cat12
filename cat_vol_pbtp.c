/* quasi-euclidean distance calculation
 * _____________________________________________________________________________
 * [GMT,RPM] = cat_vol_pbtp(SEG,WMD,CSFD[,opt])
 *
 * SEG  = (single) segment image with low and high boundary bd
 * GMT  = (single) thickness image
 * RPM  = (single) radial position map
 * WMD  = (single) CSF distance map
 * CSFD = (single) CSF distance map
 *
 * opt.bd   = (single) [low,high] boundary values (default 1.5 and 2.5)
 * opt.CSFD = calculate CSFD
 * opt.PVE  = use PVE information (0=none,1=fast,2=exact)
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

struct opt_type {
  int   CSFD;                         /* use CSFD */
  int   PVE;                          /* 0, 1=fast, 2=exact */
  float LB, HB, LLB, HLB, LHB, HHB;   /* boundary */
  int   sL[3];
  } opt;

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
  float n=0.0, maximum=WMD, high_avg=0.0;

  /* the pure maximum */
  for (int i=0;i<sA;i++) {
    if (  ( GMT[i] < 1e15 ) && ( maximum < GMT[i] ) &&                  /* thickness/WMD of neighbors should be larger */
          ( SEG[i] >= 1.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) &&           /* projection range */
          ( ( ( RPM[i] - ND[i] * 1.2 ) <= WMD ) ) &&                    /* upper boundary - maximum distance */
          ( ( ( RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&  /* lower boundary - minimum distance - corrected values outside */
          ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1.5)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) )  /* for high values will project data over sulcal gaps */
      { maximum = GMT[i]; n++; } 
  }
  
  /* the median of the highest values*/
  float high_vals[14]; int count = 0;
  for (int i=0;i<sA;i++) {
    if ( ( GMT[i] < 1e15 ) && ( (maximum - 1) < GMT[i] ) && 
         ( SEG[i] >= 1.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) && 
         ( ( (RPM[i] - ND[i] * 1.2 ) <= WMD ) ) && 
         ( ( (RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&
         ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1.5)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) ) 
      { high_vals[count++] = GMT[i]; } 
  }
  if ( count > 0 ) {
    qsort(high_vals, count, sizeof(float), compare_float);
    /* Use central 50% of sorted values for robust averaging */
    n=0.0; 
    for (int i=(int)round(count*.25); i<(int)round(count*.75); i++) {
      high_avg = high_avg + high_vals[i];
      n++; 
    }
    if ( n>0 )
      high_avg /= n;  
    else
      high_avg = maximum; 
  }
  else
    high_avg = maximum; 

  return high_avg;
}




/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  if (i < 0) i = 0;
  else if (i >= snL) i = snL - 1;

  *z = (int)floor( (double)i / (double)sxy ) ; 
  i = i % sxy;
  *y = (int)floor( (double)i / (double)sy ) ;        
  *x = i % sy ;
}



/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enough input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: too many input elements.\n");
  if (nlhs<1) mexErrMsgTxt("ERROR: not enough output elements.\n");
  if (nlhs>2) mexErrMsgTxt("ERROR: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0) mexErrMsgTxt("ERROR: second input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[2])==0) mexErrMsgTxt("ERROR: third input must be an 3d single matrix\n");
 
  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  const mwSize *sW = mxGetDimensions(prhs[1]); 
  const mwSize *sC = mxGetDimensions(prhs[2]); 
  mwSize sSEG[3] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
  if (mxGetNumberOfDimensions(prhs[1]) != dL || mxGetNumberOfDimensions(prhs[2]) != dL ||
      sW[0] != sL[0] || sW[1] != sL[1] || (dL>2 && sW[2] != sL[2]) ||
      sC[0] != sL[0] || sC[1] != sL[1] || (dL>2 && sC[2] != sL[2])) {
    mexErrMsgTxt("ERROR: input matrices must have identical dimensions\n");
  }
  const float   s2 = sqrt(2.0);
  const float   s3 = sqrt(3.0);
  const int     nr = nrhs;
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   sN = 14;
  const int   NI[14]  = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[14]  = {0.0,1.0,  s2,1.0,  s2,     s2,1.0,   s2,       s3,   s2,     s3,       s3,   s2,     s3};
  float       DN[14],DI[14],GMTN[14],WMDN[14],SEGN[14],DNm;
  float       du, dv, dw, dnu, dnv, dnw, d, dcf, WMu, WMv, WMw, GMu, GMv, GMw, SEGl, SEGu, tmpfloat;
  int         ni,u,v,w,nu,nv,nw, tmpint, WMC=0, CSFC=0;
    
  /* main volumes - actual without memory optimization ... */
  mxArray *hlps[1];  
  hlps[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
/* not yet defined  
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[3] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[4] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);  
*/  

  /* input variables */
  float*SEG  = (float *)mxGetPr(prhs[0]);
  const float*WMD_in  = (float *)mxGetPr(prhs[1]);
  const float*CSFD_in = (float *)mxGetPr(prhs[2]);
  
  /* create internal copies to preserve input arrays */
  float *WMD  = (float *)mxMalloc(nL * sizeof(float));
  float *CSFD = (float *)mxMalloc(nL * sizeof(float));
  for (int i=0; i<nL; i++) {
    WMD[i]  = WMD_in[i];
    CSFD[i] = CSFD_in[i];
  }

  if (nrhs > 3) {
    mxArray *field;
    field = mxGetField(prhs[3], 0, "CSFD");
    if (field) opt.CSFD = (int)mxGetScalar(field);
    field = mxGetField(prhs[3], 0, "PVE");
    if (field) opt.PVE = (int)mxGetScalar(field);
    field = mxGetField(prhs[3], 0, "LB");
    if (field) opt.LB = (float)mxGetScalar(field);
    field = mxGetField(prhs[3], 0, "HB");
    if (field) opt.HB = (float)mxGetScalar(field);
  } else {
    opt.CSFD = 1; opt.PVE = 2; opt.LB = 1.5; opt.HB = 2.5;
  }
  opt.LLB=floor(opt.LB), opt.HLB=ceil(opt.LB), opt.LHB=floor(opt.HB), opt.HHB=ceil(opt.HB);
  
  /* output variables */
  float        *GMT  = (float *)mxGetPr(plhs[0]);
  float        *RPM  = (float *)mxGetPr(hlps[0]);
           
  
  /* intitialisiation */
  for (int i=0;i<nL;i++) {
    GMT[i] = WMD[i];
    RPM[i] = WMD[i];
    /* proof distance input */
    if ( SEG[i]>=opt.HB ) WMC++;
    if ( SEG[i]<=opt.LB ) CSFC++;
  }
  if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
  if (CSFC==0) opt.CSFD = 0;

  
  
/* thickness calculation
 ======================================================================= */
  for (int i=0;i<nL;i++) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);
      
      /* read neighbour values */
      for (int n=0;n<sN;n++) {
        ni = i + NI[n];
        if (ni < 0 || ni >= nL) ni = i;
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      /* find minimum distance within the neighborhood */
      DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN);
      GMT[i] = DNm;
    }
  }
  
  for (int i=nL-1;i>=0;i--) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,&u,&v,&w,nL,xy,x);
      
      /* read neighbour values */
      for (int n=0;n<sN;n++) {
        ni = i - NI[n];
        if (ni < 0 || ni >= nL) ni = i;
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      /* find minimum distance within the neighborhood */
      DNm = pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN);
      if ( GMT[i] < DNm && DNm>0 ) GMT[i] = DNm;
    }
  }
  
  for (int i=0;i<nL;i++) if (SEG[i]<opt.LB || SEG[i]>opt.HB) GMT[i] = 0.0; 

 
  


/* final settings...
 ======================================================================= */
  float CSFDc = 0.0, GMTi, CSFDi; /* 0.125 */
  for (int i=0;i<nL;i++) { 
    /* GMT[i] = min(CSFD[i] + WMD[i],GMT[i]); */
    if (SEG[i]>=opt.LB && SEG[i]<=opt.HB) {
      GMTi   = CSFD[i] + WMD[i];  
      CSFDi  = GMT[i]  - WMD[i];
    
      if ( CSFD[i]>CSFDi )  CSFD[i] = CSFDi;          
      else                  GMT[i]  = GMTi;
    }
  }

 
/* estimate RPM
 ======================================================================= */
  for (int i=0;i<nL;i++) {
    if ( SEG[i]>=opt.HB )   
      RPM[i] = 1.0; 
    else {
      if ( SEG[i]<=opt.LB || GMT[i]==0.0 ) 
        RPM[i] = 0.0;
      else {
        RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
        if (RPM[i]>1.0) RPM[i] = 1.0;
        if (RPM[i]<0.0) RPM[i] = 0.0; 
      }
    } 
  }
  
  /* create second output */
  if ( nlhs > 1 ) {
    plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
    float        *RPMO  = (float *)mxGetPr(plhs[1]);
    for (int i=0;i<nL;i++) RPMO[i] = RPM[i]; 
  }
  
  /* clear internal variables */
  mxFree(WMD);
  mxFree(CSFD);
  mxDestroyArray(hlps[0]);
}


