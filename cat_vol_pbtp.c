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

/* get all values of the voxels which are in WMD-range (children of this voxel) */
float pmax(const float GMT[], const float RPM[], const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA) {
  float T[27];  
  float n=0.0, maximum=WMD; 
  for (int i=0;i<27;i++) T[i]=-1;

  /* the pure maximum */
  /* (GMT[i]<1e15) && (maximum < GMT[i]) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.5)>WMD) && (SEGI)>=SEG[i] && SEG[i]>1 && SEGI>1.66) */
  for (int i=0;i<=sA;i++) {
    if (  ( GMT[i] < 1e15 ) && ( maximum < GMT[i] ) &&                  /* thickness/WMD of neighbors should be larger */
          ( SEG[i] >= 1.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) &&           /* projection range */
          ( ( ( RPM[i] - ND[i] * 1.2 ) <= WMD ) ) &&                    /* upper boundary - maximum distance */
          ( ( ( RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&  /* lower boundary - minimum distance - corrected values outside */
          ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1.5)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) )  /* for high values will project data over sulcal gaps */
      { maximum = GMT[i]; }
  }

  
  /* the mean of the highest values*/
  float maximum2=maximum; float m2n=0.0; 
  for (int i=0;i<=sA;i++) {
    if ( ( GMT[i] < 1e15 ) && ( (maximum - 1) < GMT[i] ) && 
         ( SEG[i] >= 1.0 ) && ( SEGI>1.2 && SEGI<=2.75 ) && 
         ( ( (RPM[i] - ND[i] * 1.2 ) <= WMD ) ) && 
         ( ( (RPM[i] - ND[i] * 0.5 ) >  WMD ) || ( SEG[i]<1.5 ) ) &&
         ( ( ( (SEGI * max(1.0,min(1.2,SEGI-1.5)) ) >= SEG[i] ) ) || ( SEG[i]<1.5 ) ) ) 
      { maximum2 = maximum2 + GMT[i]; m2n++; } 
  }
  if ( m2n > 0.0 )  maximum = (maximum2 - maximum) / m2n;

  return maximum;
}




/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  /* not here ... 
   *  if (i<0) i=0; 
   *  if (i>=snL) i=snL-1;
  */
  
  *z = (int)floor( (double)i / (double)sxy ) ; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) ;        
  *x = i % sy ;
}



/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enough input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: too many input elements.\n");
  if (nlhs<1) mexErrMsgTxt("ERROR: not enough output elements.\n");
  if (nlhs>2) mexErrMsgTxt("ERROR: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first  input must be an 3d single matrix\n");
 
  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  mwSize sSEG[3] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
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
  float*WMD  = (float *)mxGetPr(prhs[1]);
  float*CSFD = (float *)mxGetPr(prhs[2]);

  /*if ( nrhs>1) {
    tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"CSFD"));  printf("X=%d", tmpint); if ( tmpint!=NULL && (tmpint>=0 && tmpint<=1) ) opt.CSFD = tmpint;   else opt.CSFD  = 1;
    tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"PVE"));   printf("X=%d", tmpint); if ( tmpint!=NULL && (tmpint>=0 && tmpint<=2) ) opt.PVE  = tmpint;   else opt.PVE   = 2;
    tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"LB"));  printf("X=%d", tmpfloat); if ( tmpfloat!=NULL )                           opt.LB   = tmpfloat; else opt.LB    = 1.5;
    tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"HB"));  printf("X=%d", tmpfloat); if ( tmpfloat!=NULL )                           opt.HB   = tmpfloat; else opt.HB    = 2.5;
  } 
  else */{ opt.CSFD = 1;opt.PVE = 2;opt.LB = 1.5;opt.HB = 2.5; }
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
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
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
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
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
    if (SEG[i]>=opt.LB && SEG[i]<=opt.LB) {
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
  /*
  mxDestroyArray(hlps[0]);
   */
}


