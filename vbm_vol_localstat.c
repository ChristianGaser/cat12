/* Local Mean, Minimum, Maximum, SD, and Peak estimation
 * ________________________________________________________________________
 * S = vbm_vol_localstat(V,M,nb,stat)
 * 
 * V    (single)    input volume
 * M    (logical)   mask volume
 * nb   (double)    neigbhour distance (1 .. 10)
 * stat (double)    1-mean, 2-min, 3-max, 4-std, 5-peak
 *
 * ________________________________________________________________________
 * Robert Dahnke 2012_01
 * Center of Neuroimaging 
 * University Jena
 *
 * $Id$ 
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))


float abs2(float n) {	if (n<0) return -n; else return n; }        
float pow2(float n) { n * n;}

/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs<2) mexErrMsgTxt("ERROR:vbm_vol_localstat: not enough input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR:vbm_vol_localstat: too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:vbm_vol_localstat: not enough output elements\n");
  if (nlhs>5) mexErrMsgTxt("ERROR:vbm_vol_localstat: too many output elements\n");

  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = (int) mxGetNumberOfElements(prhs[0]);
  const mwSize *sB = mxGetDimensions(prhs[0]);
  const int     dB = mxGetNumberOfDimensions(prhs[0]);
  const int     nB = (int) mxGetNumberOfElements(prhs[0]);
  int nh, st; 
  
  if ( dL != 3 || mxIsSingle(prhs[0])==0)                mexErrMsgTxt("ERROR:vbm_vol_localstat: first input must be a single 3d matrix\n");
  if ( dB != 3 || mxIsLogical(prhs[1])==0 || nL != nB )  mexErrMsgTxt("ERROR:vbm_vol_localstat: second input must be a logical 3d matrix with equal size than input 1\n");
  if (nrhs<3)  nh = 1; else  nh = (int) *mxGetPr(prhs[2]);
  if ( nh > 10 )                                         mexErrMsgTxt("ERROR:vbm_vol_localstat: number of neighbors is limited to 10. (Use reduce resolution instead.) \n");
  if (nrhs<4)  st = 1; else st = (int) *mxGetPr(prhs[3]);
  if ( st<1 || st>6 )                                    mexErrMsgTxt("ERROR:vbm_vol_localstat: fourth input has to be 1=mean, 2=min, 3=max, 4=std. \n");
  if (nrhs==5 && mxIsDouble(prhs[4])==0)                 mexErrMsgTxt("ERROR:vbm_vol_localstat: fifth input (vx_vol) must be an double matrix\n");
  if (nrhs==5 && mxGetNumberOfElements(prhs[4])!=3)      mexErrMsgTxt("ERROR:vbm_vol_localstat: fifth input (vx_vol) must have 3 Elements"); 
  
  /* vx_vol */
  const int sS[] = {1,3}; mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL); double*S = mxGetPr(SS);
  if (nrhs<4) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[3]);}
  float s1 = abs2((float)S[0]),s2 = abs2((float)S[1]),s3 = abs2((float)S[2]);
   
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
const int NVs=(int) (2*nh+1)*(2*nh+1)*(2*nh+1);
//printf("%d",st);

  float NVstdth,nx,NVnmax,stdth=0.90; //1-1/nh; 1 - 1/(2*nh+1); //*(2*nh+1));
  float NV[9261],NVn[9261],GV[9261],DN[9261], NVmn, NVstd; // nmax ==10
  
  int i,j,k,ind,ni,x,y,z,n,nn,HIST[1000]; //,HIST1[1000],HIST2[1000];
        
  /* in- and output */
  float *D = (float *) mxGetPr(prhs[0]);
  bool  *B = (bool  *) mxGetPr(prhs[1]);
  
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M   = (float *) mxGetPr(plhs[0]);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M2  = (float *) mxGetPr(plhs[1]);
  
  int n1i,n2i; 
  int HISTmax=40; //HISTmax1=40; HISTmax2=40; //(int) (((float) (NVs))/20); if (HISTmax>1000) HISTmax=1000; 
  const int     sx  = (int) sL[0];
  const int     sy  = (int) sL[1];
  const int     sxy = x*y;
 
  /* filter process */
  for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
    ind = index(x,y,z,sL);
    if (B[ind]) {
      n = 0;
      /* go through all elements in a nh x nh x nh box */
      for (i=-nh;i<=nh;i++) for (j=-nh;j<=nh;j++) for (k=-nh;k<=nh;k++) {
        /* check borders, masks, NaN or Infinities */
        if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2]) ) {
          ni     = index(x+i,y+j,z+k,sL);
          DN[ni] = (float) sqrt( (double) ((i * i) + 
                                 (j * j) + 
                                 (k * k)) );
          // ( B[ni]==1 && isnan(D[ni])==0 && D[ni]!=FLT_MAX && D[ind]!=-FLT_MAX && (DN[ni]<=(float)nh)
          if ( D[ni]>0 && (DN[ni]<=(float)nh) ) { //&& (DN[ni]<=(float)nh) ) {
            NV[n] = D[ni];
            n++;
          }
        }
      }
     
      NVmn = 0;nx=0; 
      if (st==1) { M[ind]=   0.0; for (nn=0;nn<n;nn++) { M[ind]+=NV[nn]; nx++;}; M[ind]/=nx;};
      if (st==2) { M[ind]=D[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]<M[ind]) M[ind]=NV[nn];};}; 
      if (st==3) { M[ind]=D[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]>M[ind]) M[ind]=NV[nn];};}; 
      if (st==4) {
        M[ind] = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind];
        NVstd  = 0.0; for (nn=0;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
        M[ind] = (float)sqrt((double)(NVstd/(nx-1)));
      };
      if (st==5) {
        /* max in histogram */
        for (nn=0;nn<200;nn++) {HIST[nn]=0;}
        for (nn=0;nn<n;nn++) {HIST[(int) round( NV[nn]* (float) HISTmax) ]++;}
        M[ind]=0; for (nn=0;nn<200;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;}
        M[ind]=M[ind]/(float) HISTmax;
      };  
      if (st==6) {
        /* max in histogram */
        NVmn = D[ni]; //= 0.0; for (nn=0;nn<n;nn++) { NVmn  +=  NV[nn];}; NVmn/=n;
        
        for (nn=0;nn<200;nn++) {HIST[nn]=0;}
        for (nn=0;nn<n;nn++) {HIST[(int) round( NV[nn]* (float) HISTmax) ]++;}
        //for (nn=0;nn<n;nn++) if (NV[nn]<MVmn) {HIST1[(int) round( NV[nn]* (float) HISTmax) ]++;}
        M[ind]=0; for (nn=(int) round( NVmn * (float) HISTmax);nn<200;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;}
        //M[ind]=0; for (nn=0;nn<200;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;}
        M[ind]=M[ind]/(float) HISTmax;
        M2[ind]=0; for (nn=0;nn<(int) round( NVmn * (float) HISTmax);nn++) { if (HIST[nn]>M2[ind]) M2[ind]=(float) nn;}
        M2[ind]=M2[ind]/(float) HISTmax;
      };  
      
      if ((M[ind]==-FLT_MAX) || (D[ind]==FLT_MAX) || (isnan(D[ind])) ) M[ind]=0;
     /*
      H[ind]=-FLT_MAX; L[ind]=FLT_MAX; NVmn = 0; NVstd = 0; nx=0; 
      for (nn=0;nn<n;nn++) {
        if (NV[nn]>H[ind]) H[ind]=NV[nn];                                  // max  
        if (NV[nn]<L[ind]) L[ind]=NV[nn];                                  // min  
        NVmn+=NV[nn]; nx++;}; NVmn/=nx;                                    // mean 
        NVstd+=pow(NV[nn]-NVmn,2);}; NVstd=sqrt(NVstd/(nx-1));             // std   
      }
    */
    /*  
      NVmn = 0; NVstd = 0; 
      if (n>=0) {
        for (nn=0;nn<n;nn++) {NVmn+=NV[nn];}; NVmn/=n;                                         
        for (nn=0;nn<n;nn++) {NVstd+=(NV[nn]-NVmn)*(NV[nn]-NVmn);}; NVstd=sqrt(NVstd/(n-1));  
      
        if ( stdth>0 ) {
          NVnmax=0; for (i=0;i<n;i++) {
            NVn[i]=abs2(NVn[i]-NVmn); if (NVnmax<NVn[i]) NVnmax=NVn[i];
            // NVn[i]=abs2(NVn[i]-NVmn) + GV[i]; if (NVnmax<NVn[i]) NVnmax=NVn[i];
          }
          NVstdth=NVnmax*stdth;
          
          H[ind]=-FLT_MAX; L[ind]=FLT_MAX; NVmn = 0; NVstd = 0; nx=0; 
          for (nn=0;nn<n;nn++) {if (NVn[nn]<NVstdth) {if (NV[nn]>H[ind]) H[ind]=NV[nn];};}                               
          for (nn=0;nn<n;nn++) {if (NVn[nn]<NVstdth) {if (NV[nn]<L[ind]) L[ind]=NV[nn];};}                              
          for (nn=0;nn<n;nn++) {if (NVn[nn]<NVstdth) {NVmn+=NV[nn]; nx++;};}; NVmn/=nx;                                 
          for (nn=0;nn<n;nn++) {if (NVn[nn]<NVstdth) {NVstd+=(NV[nn]-NVmn)*(NV[nn]-NVmn);};}; NVstd=sqrt(NVstd/(nx-1)); 
        
          if (H[ind]==-FLT_MAX) H[ind]=0;
          if (L[ind]== FLT_MAX) L[ind]=0;
          if (isnan(NVmn))      NVmn=0;
          if (isnan(NVstd))     NVstd=0;
          
        }

      
      if (H[ind]==-FLT_MAX) H[ind]=0;
      if (L[ind]== FLT_MAX) L[ind]=0;
      if (isnan(NVmn))      NVmn=0;
      if (isnan(NVstd))     NVstd=0;

      M[ind] = NVmn;
      SD[ind] = NVstd;
           */    
     

    }
    else {
      M[ind] = 0; /*D[ind];*/
     //SD[ind] = 0;      
    }
    
    
  } 
  
  
  
}


