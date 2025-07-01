/* Local Mean, Minimum, Maximum, SD, and Peak estimation
 * ________________________________________________________________________
 * Estimates specific functions in a volume V within a mask region M. 
 * For each voxel v of M, the values of the neigbors of v that belong to M 
 * and are within a distance smaller than nb where used (i.e. masked voxels
 * within a sphere of radius nb). 
 *
 * If V contains NaNs, -INFs or INFs these values are ignored and added to 
 * the mask M.  Masked voxels in the output volume were defined depending 
 * on the mskval  variable, i.e. as zeros (0, default), input (1), NANs (2),
 * -INF (3), or INF(4).  
 *
 * The function was designed to extract the inhomogeneity in noisy data in 
 * a well-known area, such a tissue class in structural images. In general  
 * the mean estimated in a local neighborhood nb or after some interations
 * iter. However, tissues are often affected by partial volume effects near
 * the tissue boundaries and minimum/maximum can be quite helpful to reduce
 * such effects. For noise estimation the local variance/standard deviation 
 * is also quite useful. 
 * Besides the mean value the local peak of the histogram can also work 
 *
 *  S = cat_vol_localstat(V,M[,nb,stat,iter,filter0,verb])
 * 
 *  V    (single)    input volume
 *  M    (logical)   mask volume
 *  nb   (double)    neigbhour distance (1 .. 10)
 *  stat (double)    1-mean, 2-min, 3-max, 4-std  
 *                   5-peak1, 6-peak2, 7-peak3    (experimental)
 *                   8-median                     
 *                   9-hist                       (experimental)
 * iter             number of iterations          (default=1)
 * filter0 (double) originally values <=0 were ignored (default=0)
 * mskval  (double) setting of masked voxels
 *                   (0-zeros,1-input,2-NAN,3--INF,4-INF)
 * verb (double)    verbose output for debugging
 *
 *
 * Examples: 
 *  Here are some simple samples to outline the subfunctions. The mask area
 *  is defined by NaN. The simulated data of A is between -1 and 1 and B is 
 *  a locial mask.
 *
 *  == input variables ==
 *  A  = rand(20,20,3,'single') - 1;           
 *    for i=1:size(A,2), A(:,i,:) = A(:,i,:) + (( i/size(A,2) ) - 0.5); end
 *  B  = smooth3(smooth3(rand(size(A))))>0.5;
 *
 *   
 *  == function calls ==
 *  (1) MEAN:     values around 0
 *      C = cat_vol_localstat(A,B,2,1,2,2); ds('d2smns','',1,A,C,2); 
 *
 *  (2) MINIMUM:  values trending torwards -1
 *      C = cat_vol_localstat(A,B,2,2,2,2); ds('d2smns','',1,A,C,2); 
 *
 *  (3) MAXIMUM:  values trending torwards 1
 *      C = cat_vol_localstat(A,B,2,3,2,2); ds('d2smns','',1,A,C,2); 
 *
 *  (4) STANDARD DEVIATION:  values about 0.5 
 *      C = cat_vol_localstat(A,B,2,4,1,2); ds('d2smns','',1,A,C,2); 
 *
 *  (8) MEDIAN:   values around 0
 *      C = cat_vol_localstat(AL,B,2,8,1,2); ds('d2smns','',1,A,C,2); 
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
#include <stdio.h>
#include <stdlib.h>

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifdef _MSC_VER
  #define FINFINITY (FLT_MAX+FLT_MAX);
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FINFINITY 1.0f/0.0f;
  #define FNAN 0.0f/0.0f
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))


/* qicksort for median */
void swap_float(float *a, float *b)
{
  float t=*a; *a=*b; *b=t;
}

void sort_float(float arr[], int beg, int end)
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
        swap_float(&arr[l], &arr[--r]);
    }
    swap_float(&arr[--l], &arr[beg]);
    sort_float(arr, beg, l);
    sort_float(arr, r, end);
  }
}



/* floating point versions */
float min(float a, float b) { if (a<b) return a; else return b; }
float max(float a, float b) { if (a>b) return a; else return b; }
float abs2(float n) { if (n<0) return -n; else return n; }        
float pow2(float n) { return n*n; }



/* main function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check size of in-/output */
  if (nrhs<2) mexErrMsgTxt("ERROR:cat_vol_localstat: Not enough input elements\n");
  if (nrhs>8) mexErrMsgTxt("ERROR:cat_vol_localstat: Too many input elements\n");
  if (nlhs<1) mexErrMsgTxt("ERROR:cat_vol_localstat: Not enough output elements\n");
  if (nlhs>1) mexErrMsgTxt("ERROR:cat_vol_localstat: Too many output elements\n");

  
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     dB = mxGetNumberOfDimensions(prhs[0]);
  const int     nB = mxGetNumberOfElements(prhs[0]);
  
  
  /* test properties of input */ 
  if ( dL != 3 || mxIsSingle(prhs[0])==0)                mexErrMsgTxt("ERROR:cat_vol_localstat: First input must be a single 3d matrix\n");
  if ( dB != 3 || mxIsLogical(prhs[1])==0 || nL != nB )  mexErrMsgTxt("ERROR:cat_vol_localstat: Second input must be a logical 3d matrix with equal size than input 1\n");
  
  
  /* set control parameter and defaults */   
  int nh, st, iter, mskval, filter0, verb;
  if (nrhs<3)  nh      = 1; else {double *dnh      = (double *) mxGetPr(prhs[2]); nh      = (int) round( dnh[0]     );}   
  if (nrhs<4)  st      = 1; else {double *dst      = (double *) mxGetPr(prhs[3]); st      = (int) round( dst[0]     );}
  if (nrhs<5)  iter    = 1; else {double *diter    = (double *) mxGetPr(prhs[4]); iter    = (int) round( diter[0]   );}  
  if (nrhs<6)  mskval  = 0; else {double *dmskval  = (double *) mxGetPr(prhs[5]); mskval  = (int) round( dmskval[0] );}  
  if (nrhs<7)  filter0 = 0; else {double *dfilter0 = (double *) mxGetPr(prhs[6]); filter0 = (int) ( dfilter0[0]>0   );}  
  if (nrhs<8)  verb    = 0; else {double *dverb    = (double *) mxGetPr(prhs[7]); verb    = (int) ( dverb[0]>0      );} 
  /* check control parameter */
  if ( nh > 10 )                mexErrMsgTxt("ERROR:cat_vol_localstat: Number of neighbors is limited to 10. (Use reduce resolution instead.) \n");
  if ( st<1 || st>9 )           mexErrMsgTxt("ERROR:cat_vol_localstat: Fourth input has to be 1=mean, 2=min, 3=max, 4=std, 5=peak1, 6=peak2, 7=peak3, 8=median, 9=hist. \n");
  if ( (st>4 && st<8) || st>8 ) mexErrMsgTxt("ERROR:cat_vol_localstat: Experimental function! \n");
  if ( iter<1 )                 printf("WARNING:cat_vol_localstat: Number of iteration is 0 and nothing is done except masking! \n");
  if ( mskval<0 || mskval>4 )   mexErrMsgTxt("ERROR:cat_vol_localstat: mskval has to be 0=zeros,1=input,2=NAN,3=-INF,4=INF. \n");

  
  /* input varialbes */
  const float *D = (float *) mxGetPr(prhs[0]);
  const bool  *B = (bool  *) mxGetPr(prhs[1]);
  
  
  /* create and initialise output variables and internal variables */
  mxArray *hlps[4], *msks[1];
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M   = (float *) mxGetPr(plhs[0]); /* only real output */
  hlps[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M2  = (float *) mxGetPr(hlps[0]); /* helping variable */
  hlps[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *M3  = (float *) mxGetPr(hlps[1]); /* helping variable */
  hlps[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); float *MD  = (float *) mxGetPr(hlps[2]); /* helping variable */
  msks[0] = mxCreateLogicalArray(dL,sL);                       bool  *MB  = (bool  *) mxGetPr(msks[0]); /* helping variable */
  /* initialise */
  for (int i=0;i<nL;i++) M[i]  = D[i]; 
  for (int i=0;i<nL;i++) M2[i] = (float)0.0; 
  for (int i=0;i<nL;i++) M3[i] = (float)0.0;  
  for (int i=0;i<nL;i++) MD[i] = D[i]; /* define for first iteration */
  for (int i=0;i<nL;i++) MB[i] = B[i]; /* define for first iteration */
  /* update input variables -	correction of non-visited or other incorrect voxels */
  for (int i=0;i<nL;i++) { 
    if ( mxIsInf(D[i]) || mxIsNaN(D[i] )) MB[i] = false; 
    if ( filter0==0 && D[i]==0 )          MB[i] = false; 
  }
  
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  float nx; 
  float NV[9261],DN[9261], NVmd, NVmn, NVstd; /* nmax == 10, i.e. 21^3 */
  float stdd[3],stdp[3],stdn[3];
  int   stddc[3],stdpc[3],stdnc[3];
  int   di,i,j,k,ind,ni,x,y,z,n,nn,md; /*,HIST1[1000],HIST2[1000]; */
  int   HISTmax=256; /* HISTmax1=40; HISTmax2=40; (int) (((float) (NVs))/20); if (HISTmax>1000) HISTmax=1000; */ 
  int   HISTmin=0;

  
  /*
   * Display initial parameter
   */
  if ( verb ) {
    printf("\ncat_vol_localstat.c debuging mode:\n  Initialize Parameter: \n");
    printf("    size(B) = %d %d %d\n",(int)sL[0],(int)sL[1],(int)sL[2]); 
    printf("    nb      = %d\n",nh); 
    printf("    stat    = %d\n",st); 
    printf("    iter    = %d\n",iter); 
    printf("    mskval  = %d\n",mskval); 
    printf("    fitler0 = %d\n",filter0); 
    printf("    verb    = %d\n",verb); 
  }
 
  
  /* prepare special variables for histogram based subfunctions */
  if (st==7) {
    for (i=0;i<nL;i++) { if (HISTmin>D[i]) HISTmin=(int)D[i]; };  HISTmin--;
    for (i=0;i<nL;i++) { if (HISTmax<D[i]) HISTmax=(int)D[i]; };  HISTmax++;
    for (i=0;i<nL;i++) {
      MD[i] = (float)ROUND(max(min(MD[i],(float)HISTmax),(float)HISTmin));
    }
  }
  int HISTn  = HISTmax - HISTmin + 1;
  float *HIST; HIST = (float *) malloc(sizeof(float)*HISTn);
  for (nn=0;nn<HISTn;nn++) HIST[nn]=0.0;

  

  
  /* main filter process */
  if ( verb ) printf("  Filtering ");   
  for (int iteri=0;iteri<iter;iteri++) {
    if ( verb ) printf(".");   
      
    for (z=0;z<sL[2];z++) for (y=0;y<sL[1];y++) for (x=0;x<sL[0];x++) {
      ind = index(x,y,z,sL);
      if ( MB[ind] && mxIsNaN(MD[ind])==false && mxIsInf(MD[ind])==false  ) {
        n  = 0; /* number of accessed voxels */

        /* go through all elements in a nh x nh x nh box and get their values */
        for (i=-nh;i<=nh;i++) for (j=-nh;j<=nh;j++) for (k=-nh;k<=nh;k++) {
          /* check borders, masks, NaN or Infinities */
          if ( ((x+i)>=0) && ((x+i)<sL[0]) && ((y+j)>=0) && ((y+j)<sL[1]) && ((z+k)>=0) && ((z+k)<sL[2])) {
            ni = index(x+i,y+j,z+k,sL);
            if ( ( MB[ni] || st==7 ) && mxIsNaN(MD[ni])==false && mxIsInf(MD[ni])==false ) {
              DN[n] = sqrtf( ( (float) ((i * i) + (j * j) + (k * k)) ) );
              if ( DN[n]<=(float) nh )  { /* && (DN[ni]<=(float)nh) ) { */
                NV[n] = MD[ni];
                n++;
              }
            }
          }
        }

        NVmn = 0.0; nx = 0.0; 

        /* mean */
        if (st==1) { M[ind]=0.0; for (nn=0;nn<n;nn++) { M[ind]+=NV[nn]; nx++;}; M[ind]/=nx;};

        /* mean, min, max */
        if (st==10) {
          M[ind]=0.0; 
          M2[ind]=MD[ind]; 
          M3[ind]=MD[ind]; 
          for (nn=0;nn<n;nn++) { 
            M[ind]+=NV[nn]; nx++;
            if (NV[nn]<M2[ind]) M2[ind]=NV[nn];
            if (NV[nn]>M3[ind]) M3[ind]=NV[nn];
          }
          M[ind]/=nx;
        }

        /* minimum */
        if (st==2) { M[ind]=MD[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]<M[ind]) M[ind]=NV[nn];};}; 

        /* maximum */
        if (st==3) { M[ind]=MD[ind]; for (nn=0;nn<n;nn++) { if (NV[nn]>M[ind]) M[ind]=NV[nn];};}; 

        /* standard deviation */
        if (st==4) {
          M[ind]  = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind];
          M2[ind] = M[ind];
          NVstd   = 0.0; for (nn=0;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
          M[ind]  = sqrtf((float)(NVstd/(nx-1.0)));
        };


        /* ===============================================================
         * experimental histogram functions 
         * ===============================================================
         */ 
        /* max in histogram 1 */
        if (st==5) {
          for (nn=0;nn<HISTn;nn++) {HIST[nn]=0;}
          for (nn=0;nn<n;nn++) {HIST[(int) ROUND( NV[nn]* (float) HISTmax) ]++;}
          M[ind]=0.0; for (nn=0;nn<HISTn;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) HIST[nn]; }
          M[ind]=M[ind]/(float) HISTmax;
        };  

        /* max in histogram 2 */
        if (st==6) {
          NVmn = MD[ni]; /*= 0.0; for (nn=0;nn<n;nn++) { NVmn  +=  NV[nn];}; NVmn/=n; */

          for (nn=0;nn<HISTn;nn++) {HIST[nn]=0;}
          for (nn=0;nn<n;nn++) {HIST[(int) ROUND( NV[nn]* (float) HISTmax) ]++;}
          /*for (nn=0;nn<n;nn++) if (NV[nn]<MVmn) {HIST1[(int) ROUND( NV[nn]* (float) HISTmax) ]++;} */
          M[ind]=0; for (nn=(int) ROUND( NVmn * (float) HISTmax);nn<HISTn;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;}
          /*M[ind]=0; for (nn=0;nn<200;nn++) { if (HIST[nn]>M[ind]) M[ind]=(float) nn;} */
          M[ind]=M[ind]/(float) HISTmax;
          M2[ind]=0.0; for (nn=0;nn<(int) ROUND( NVmn * (float) HISTmax);nn++) { if (HIST[nn]>M2[ind]) M2[ind]=(float) nn;}
          M2[ind]=M2[ind]/(float) HISTmax;
        };  

        /* max in histogram 3 */
        if (st==7) {
          for (nn=0;nn<HISTn;nn++) {HIST[nn]=0.0;} /* init histogram */
          for (nn=0;nn<n;nn++) {HIST[(int) (NV[nn] - HISTmin)]++;}  /* estimate histogram */
          NVmn=0.0; M[ind]=0.0; for (nn=0;nn<HISTn;nn++) { if (HIST[nn]>NVmn) {NVmn=HIST[nn]; M[ind]=(float) nn;}}
          M[ind]=M[ind] + HISTmin;
        };   

        
        

        /* median */
        if (st==8) {
          if (n>nh*nh*nh) { 
            sort_float(NV,0,n); 
            md=(int)ROUND(n*0.5);
            NVmd   = NV[md];
            M[ind] = NV[md];
          }
        }



        /* ===============================================================
         * experimental noise/signal functions 
         * ===============================================================
         */ 
        if (st==9) {
          /* estimation of noise and strukture a subregions (normally the WM) */

          /* mean (or better median) intensity in the area */
          /*  M[ind] = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind]; */
          if (n>1) { if (n==2) {
              NVmd = (NV[0] + NV[1]) / 2.0f;  
            }
            else {
              sort_float(NV,0,n); 
             NVmd = NV[(int)(n/2.0)];
            }
          }
          M[ind] = 0.0; for (nn=0;nn<n;nn++) { M[ind]+= NV[nn]; nx++;}; M[ind]/=nx; NVmn=M[ind];
          NVstd  = 0.0; for (nn=0;nn<n;nn++) { NVstd += (NV[nn]-NVmn)*(NV[nn]-NVmn);};
          M[ind] = (float)sqrt((double)(NVstd/(nx-1.0)));

          /* estimation of noise for each dimention */
          stdd[0]=0.0; stdd[1]=0.0;stdd[2]=0.0;
          stdp[0]=0.0; stdp[1]=0.0;stdp[2]=0.0;
          stdn[0]=0.0; stdn[1]=0.0;stdn[2]=0.0;

          stddc[0]=0; stddc[1]=0;stddc[2]=0;
          stdpc[0]=0; stdpc[1]=0;stdpc[2]=0;
          stdnc[0]=0; stdnc[1]=0;stdnc[2]=0;

          for (di=0;di<=2;di++) { 
            for (i=-nh*2;i<=nh*2;i++) {
              if (di==0)
                ni = index(x+i,y,z,sL);
              else {
                if (di==1)
                  ni = index(x,y+i,z,sL);
                else
                  ni = index(x,y,z+1,sL);
              }

              if ( (x+i)>=0 && (x+i)<sL[0] ) {
                if (MB[ni] && mxIsNaN(MD[ni])==false && mxIsInf(MD[ni])==false ) {
                  stdd[di] += (MD[ni]-NVmn)*(MD[ni]-NVmn); stddc[di]++;
                  if (MD[ind]>NVmd) {
                    stdp[di] += (MD[ni]-NVmn)*(MD[ni]-NVmn); stdpc[di]++;
                  }
                  else {
                    stdn[di] += (MD[ni]-NVmn)*(MD[ni]-NVmn); stdnc[di]++;
                  }
                }
              }
            }

            if ( stddc[di]>1 ) {stdd[di]=sqrtf((float)(stdd[di]/(stddc[di]-1)));} else {stdd[di] = 0.0;}
            if ( stdpc[di]>1 ) {stdp[di]=sqrtf((float)(stdp[di]/(stdpc[di]-1)));} else {stdp[di] = 0.0;}
            if ( stdnc[di]>1 ) {stdn[di]=sqrtf((float)(stdn[di]/(stdnc[di]-1)));} else {stdn[di] = 0.0;}
          }
         /* sort_float(stdd,0,2); */

          M[ind]=stdd[0];
          M2[ind]=stdd[1];
          /* M[ind]  = (stdd[0] + stdd[1] + stdd[2])/3 * 2; */
          /* M[ind]  = (stdp[0] + stdp[1] + stdp[2])/3 * 2; */
          /* M2[ind] = (stdn[0] + stdn[1] + stdn[2])/3 * 2 -  M[ind]; */ 

  /* if ((M[ind]==-FLT_MAX) || (MD[ind]==FLT_MAX) || (mxIsNaN(MD[ind])) ) M[ind]=0.0; */
        }
      } 
    }
  
    /* update for next iteration */
    for (i=0;i<nL;i++) MD[i]=M[i];  
 
  }
  
  /* final setting of output values */
  /* https://www.gnu.org/software/libc/manual/html_node/Infinity-and-NaN.html */
  switch ( mskval ) { 
    case 0: for (i=0;i<nL;i++) if ( MB[i]==false ) M[i]=0.0; break;
    case 1: for (i=0;i<nL;i++) if ( MB[i]==false ) M[i]=D[i]; break;
    case 2: for (i=0;i<nL;i++) if ( MB[i]==false ) M[i]=FNAN; break;
    case 3: for (i=0;i<nL;i++) if ( MB[i]==false ) M[i]=-FINFINITY; break;
    case 4: for (i=0;i<nL;i++) if ( MB[i]==false ) M[i]=FINFINITY; break;
  }
  
  if ( verb ) printf(" done. \n");   
}
