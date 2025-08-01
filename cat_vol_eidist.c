/* 
 * This function estimates the Euclidean distance D to the closest boundary
 * voxel I, given by the 0.5 isolevel in B with values between 0 and 1, 
 * defining subvoxel position by partial volume effects. Larger object 
 * voxel overrule smaller voxels.  
 *
 * To align the closest voxel a modified Eikonal distances T is estimated 
 * on the field L, ie. L is a speed map and T is the travel time.
 * For a correct side alignment a harder option "csf" is used that avoids
 * the diffusion to voxels with greater values of L (L[i]>L[ni] for
 * diffusion).
 *
 * Voxels with NaN and -INF are ignored and produce NaN in D, whereas 
 * in I the default index is used, because I has to be a integer matrix 
 * where NaN is not available. With setnan=0, the result of NAN voxel in 
 * D is changed to INF.
 * 
 *  [D,I,T] = cat_vol_eidist(B,L,[vx_vol,euclid,csf,setnan,verb])
 * 
 *  D         Euclidean distance map to the nearest Boundary point in the
 *            Eikonal field (3d-single-matrix)
 *            (air distance)
 *  T         Eikonal distance map in the Eikonal field (3d-single-matrix)
 *            (way length or time) 
 *  I         index map      (3d-uint32-matrix)
 *  B         boundary map   (3d-single-matrix)
 *  L         speed map      (3d-single-matrix)
 *  vx_vol    voxel size     (1x3-double-matrix): default=[1 1 1]
 *            ... not tested yet
 *  csf       option         (1x1-double-value): 0-no; 1-yes; default=1
 *  euclid    option         (1x1-double-value): 0-no; 1-yes; default=1
 *            output euclidean or speed map as first value
 *  setnan    option         (1x1-double-value): 0-no; 1-yes; default=1
 *  verb      option         (1x1-double-value): 0-no, 1-yes; default=0
 * 
 *
 * Small test examples:
 * 1) "Face" with distance from the eyes and an odd nose
 *   A=zeros(50,50,3,'single'); A(20:30,5:15,2)=10; A=smooth3(A); 
 *   A(20:30,35:45,2)=1; A(1:5,1:25,:)=nan; A(1:5,26:50,:)=-inf; 
 *   F=ones(size(A),'single'); F(10:40,20,:)=0.5; F(40,10:40,:)=0;
 *   [D,I,T]=cat_vol_eidist(A,F,[1 1 1],1,1); 
 *   ds('d2smns','',1,A - F,D/10,2); title('Euclidean distance')
 *   ds('d2smns','',1,A - F,T/10,2); title('Eikonal distance')
 *
 * 2) 1D-examples to count distance
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,6,:)=0.2; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *
 *   [D,I]=cat_vol_eidist(A,F,[1 1 1],1,1,0,1); ds('x2','',1,A,D,D,I,5)
 *
 * See also compile.m
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
 * ToDo:
 * - check anisotropic distance estimation 
 */
 

#include "mex.h"   
#include "math.h"
#include "float.h"
#include "limits.h"
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


/* 
 * Estimate x,y,z position of index i in an array size sx,sxy.
 * The index is given in c-notation for a matrix with n elements i is 
 * element of 0,1,...,n-1. Also the x,y,z axis are described in the 
 * c-notation and a matrix A with 20x30x40 elements will have coordinates 
 * with x=0,1,...,19; y=0,1,...,30; and z=0,1,...40. The index i=0 is 
 * identical with x=0,y=0,z=0; whereas the last index in the matrix A 
 * is i=n-1=(20*30*40)-1 and has the coordinates x=19,y=29, and z=39.
 *
 */


// estimation of the xyz values based on the index value 
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  // handling of boundaries 
  if (i<0) i=0; 
  if (i>=snL) i=snL-1;
  
  *z = (int)floor( (double)i / (double)sxy ); 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy );        
  *x = i % sy ;
}


/* 
 * Estimate index i of a voxel x,y,z in an array size s.
 * See also for ind2sub.
 */
int sub2ind(int x, int y, int z, int s[]) {
  // handling on boundaries 
  if (x<0) x=0; if (x>s[0]-1) x=s[0]-1; 
  if (y<0) y=0; if (y>s[1]-1) y=s[1]-1; 
  if (z<0) z=0; if (z>s[2]-1) z=s[2]-1; 
  
  //   z * (number of voxels within a slice) 
  // + y * (number of voxels in a column)  
  // + x ( what is the position within the column )      
  return (z)*s[0]*s[1] + (y)*s[0] + (x);
}


float pow_float(float x, float y) {
  return (float) pow((double) x,(double) y); 
}

float sqr_float(float x) {
  return x*x; 
}

float sqrt_float(float x) {
  return (float) sqrt((double) x); 
}

float floor_float(float x) {
  return (float) floor((double) x); 
}

/* 
 * Read out the linear interpolated value of a volume SEG with the size
 * s on the position x,y,z (c-notation). See also ind2sub for details of
 * the c-notation.
 */
float isoval(float SEG[], float x, float y, float z, int s[]){

  int i;
  float seg=0.0, n=0.0;
  float fx = floor_float(x),   fy = floor_float(y),   fz = floor_float(z);
  float cx = floor_float(x+1), cy = floor_float(y+1), cz = floor_float(z+1);
    
  float wfx = cx-x, wfy = cy-y, wfz = cz-z;
  float wcx = x-fx, wcy = y-fy, wcz = z-fz;

  /* value of the 8 neighbors and there distance weight */
  float N[8], W[8];  
  N[0]=SEG[sub2ind((int)fx,(int)fy,(int)fz,s)];  W[0]=wfx * wfy * wfz; 
  N[1]=SEG[sub2ind((int)cx,(int)fy,(int)fz,s)];  W[1]=wcx * wfy * wfz;
  N[2]=SEG[sub2ind((int)fx,(int)cy,(int)fz,s)];  W[2]=wfx * wcy * wfz;
  N[3]=SEG[sub2ind((int)cx,(int)cy,(int)fz,s)];  W[3]=wcx * wcy * wfz;
  N[4]=SEG[sub2ind((int)fx,(int)fy,(int)cz,s)];  W[4]=wfx * wfy * wcz;
  N[5]=SEG[sub2ind((int)cx,(int)fy,(int)cz,s)];  W[5]=wcx * wfy * wcz;
  N[6]=SEG[sub2ind((int)fx,(int)cy,(int)cz,s)];  W[6]=wfx * wcy * wcz; 
  N[7]=SEG[sub2ind((int)cx,(int)cy,(int)cz,s)];  W[7]=wcx * wcy * wcz;
    
  for (int i=0; i<8; i++) {
    if ( mxIsNaN(N[i])==false || mxIsInf(N[i])==false )
      seg = seg + N[i] * W[i]; n+= W[i];
  }
  if ( n>0.0 )
    return seg/n; 
  else 
    return FNAN;
}


/* 
 * MAINFUNCTION [D,I] = cat_vol_eidist(B,L,vx_vol,euclid,csf,setnan,verb])
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /* 
   * Check input 
   */
  if (nrhs<2) mexErrMsgTxt("ERROR:cat_vol_eidist: not enough input elements\n");
  if (nrhs>7) mexErrMsgTxt("ERROR:cat_vol_eidist: too many input elements.\n");
  if (nlhs>3) mexErrMsgTxt("ERROR:cat_vol_eidist: too many output elements.\n");
  if (mxIsSingle(prhs[0])==false)
    mexErrMsgTxt("ERROR:cat_vol_eidist: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==false)
    mexErrMsgTxt("ERROR:cat_vol_eidist: second input must be an 3d single matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==false)
    mexErrMsgTxt("ERROR:cat_vol_eidist: third  input must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3)
    printf("ERROR:cat_vol_eidist: third input must have 3 Elements");
  if (nrhs==4 && mxIsDouble(prhs[3])==false &&  mxGetNumberOfElements(prhs[3])!=1) 
    printf("ERROR:cat_vol_eidist: fourth input must be one double value"); 
  if (nrhs==5 && mxIsDouble(prhs[4])==false &&  mxGetNumberOfElements(prhs[4])!=1)
    printf("ERROR:cat_vol_eidist: fifth input must be one double value");
  if (nrhs==6 && mxIsDouble(prhs[5])==false &&  mxGetNumberOfElements(prhs[5])!=1)
    printf("ERROR:cat_vol_eidist: sixth input must be one double value");
  if (nrhs==7 && mxIsDouble(prhs[6])==false &&  mxGetNumberOfElements(prhs[6])!=1)
    printf("ERROR:cat_vol_eidist: seventh input must be one double value");
    
  /* 
   * set default variables 
   */
  int csf, euclid, setnan, verb;    
  float nanres;

  if (nrhs>=4) csf    = (int)*mxGetPr(prhs[3]); else csf    = 1;
  if (nrhs>=5) euclid = (int)*mxGetPr(prhs[4]); else euclid = 1;
  if (nrhs>=6) setnan = (int)*mxGetPr(prhs[5]); else setnan = 1;
  if (nrhs>=7) verb   = (int)*mxGetPr(prhs[6]); else verb   = 0;
  
  if (setnan>=1) nanres = FNAN; else nanres = FINFINITY; 
  
  /* 
   * main input variables B (boundary map) and L (speed map) 
   */
  float*B  = (float *)mxGetPr(prhs[0]); 
  float*L  = (float *)mxGetPr(prhs[1]);             
  
  /* 
   * main information about input data (size, dimensions, ...) 
   */
  const mwSize *sL = mxGetDimensions(prhs[0]); 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     nD = mxGetNumberOfElements(prhs[1]);
  const int     x  = (int)sL[0];
  const int     y  = (int)sL[1];
  const int     xy = x*y;
  int           sizeL[3] = {(int)sL[0],(int)sL[1],(int)sL[2]}; 

  if (nL!=nD) mexErrMsgTxt("ERROR:cat_vol_eidist: images must have the same number of elements.\n");
  
  /* 
   * Voxel dimension and standard distance in the N26 neighborhood.
   * Indices of the neighbor NI (index distance) and euclidean distance ND. 
   * Set default voxel size=[1 1 1] or use input.
   */
  const mwSize sS[2] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double   *S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1; S[1]=1; S[2]=1;} else S=mxGetPr(prhs[2]);
  
  float s1 = (float)fabs(S[0]);
  float s2 = (float)fabs(S[1]);
  float s3 = (float)fabs(S[2]); /* x,y,z - voxel size */
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy  - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz  - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz  - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* nL - voxel size */
  
  const int   NI[26]  = { 1, -1, x, -x, xy, -xy, -x-1, -x+1, x-1, x+1,      
                       -xy-1, -xy+1, xy-1, xy+1, -xy-x, -xy+x, xy-x, xy+x,
                       -xy-x-1, -xy-x+1, -xy+x-1, -xy+x+1, xy-x-1, xy-x+1,
                       xy+x-1,xy+x+1};  
  const float ND[26]  = { s1, s1, s2, s2, s3, s3, s12, s12,s12,s12,
                        s13, s13, s13, s13, s23, s23, s23, s23,  
                        s123, s123, s123, s123, s123, s123, s123, s123};
  const int kllv = (sL[0]+sL[1]+sL[2]);

  
  /* 
   * other variables 
   */
  float dinu, dinv, dinw, dcf, WMu, WMv, WMw, WM, DIN, DINE;  /* distance and intensity variables */
  int   ni,u,v,w,nu,nv,nw,iu,iv,iw;                       /* nL and index-values of a voxel, one neighbor and the nearest boundary voxel */
  int   nC=sL[0]*sL[1]*sL[2], nCo=INT_MAX;                    /* runtime variables */
  int   kll=0;                                                /* runtime variables */
  int   fast=1;                                               /* stop if number of unvisited points stay constant */

  
  /* 
   * Create main output volumes and variables D (distance map) and I (index map)
   */ 
  mxArray *hlps[2]; /* helping matrixes */
  plhs[0]         = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); 
  hlps[0]         = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL); 
  hlps[1]         = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  float        *D = (float *)mxGetPr(plhs[0]);        
  unsigned int *I = (unsigned int *)mxGetPr(hlps[0]);	/* index map */
  float        *T = (float *)mxGetPr(hlps[1]); 
  
  /*
   * Display Initial Parameter
   */
  if ( verb ) {
    printf("\ncat_vol_eidist.c debuging mode:\n  Initialize Parameter: \n");
    printf("    size(B) = %d %d %d\n",(int)sL[0],(int)sL[1],(int)sL[2]); 
    printf("    vx_vol  = %0.0f %0.0f %0.0f\n",s1,s2,s3); 
    printf("    euclid  = %d\n",(int) euclid); 
    printf("    setnan  = %d\n",(int) setnan); 
  }
      
 
  /* 
   * Check input values.
   */
  int vx=0;
  for (int i=0;i<nL;i++) { 
    if ( B[i]>=0.5 ) vx++;                                              /* count object voxel */
    
    if ( (L[i]<=0.0 || mxIsNaN(L[i])) && B[i]<0.5) B[i] = FNAN;         /* ignore voxel that cannot be visited */
    
    if ( mxIsInf(B[i]) && B[i]<0.0 ) B[i] = FNAN;                       /* voxel to ignore */
    if ( B[i]>1.0 ) B[i] = 1.0;                                         /* normalize object definition limitation */
    if ( B[i]<0.0 ) B[i] = 0.0;                                         /* normalize object definition limitation */
    I[i] = (unsigned int) i;                                            /* initialize index map */
    if ( nlhs>2 ) T[i] = 0.0;
    
    if ( L[i]<=0.0 || mxIsNaN(B[i]) || mxIsNaN(L[i]) ) 
      D[i] = FNAN;
    else 
      D[i] = FINFINITY;
  }
  
    
  /* 
   * Check if there is an object and return FINFINITY and i=i+1 (for matlab) 
   * for all points if there is no object.
   */
  if ( vx==0 ) { 
    for (int i=0;i<nL;i++) { D[i]=nanres; I[i]=(unsigned int)i+1; } 
    printf("WARNING:cat_vol_eidist: Found no object for distance estimation!\n");
    return;
  }
    

  /* 
   * Eikonal distance initialization:
   * In 3D this can be really complex, especially for anisotropic volumes. 
   * So only a simple approximation by the PVE of a voxel is used. 
   * RD202310: However, these was not really working in all (interpolated) 
   * images and was simplyfied.
   */ 
  if ( verb ) printf("  Initialize Eikonal distance estimation and index alignment \n");
  for (int i=0;i<nL;i++) { 
    if ( D[i]>0.0 && mxIsNaN(D[i])==false) { 
  
      if ( B[i]>=0.5 ) {
        D[i] = 0.0; 
        ind2sub(i,&u,&v,&w,nL,xy,x);
        
        for (int n=0;n<26;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 

          if ( ( (ni<0) || (ni>=nL) || 
               (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || 
               (ni==i) )==false && B[ni]<0.5) {
              /* PVE idea not robust for interpolation and GM variance */
              DIN = ND[n]; /* * (B[ni] - 0.5) / ( B[ni] - B[i] ); */ /* RD202310: simplification */
              
              if ( fabs(D[ni])>DIN ) {
                D[ni] = -DIN; 
                I[ni] = I[i];
                if ( nlhs>2 ) T[ni] = DIN;
              }
           
          }
        }
      }
    }
  }
 
  
  /* 
   * iterative Eikonal distance estimation
   */ 
  if ( verb ) printf("  Eikonal distance estimation and index alignment \n");
  while ( nC>0 && kll<kllv && (nC!=nCo || fast==false)) {

    kll++; 
    nCo=nC;
    nC=0;

    /*
     * Forward direction:
     * - push based, i.e., for each voxel (index i) we estimate the  
     *   weighed distance of its (suited) neighbors
     */
    for (int i=0;i<nL;i++) { 
      if ( D[i]<=0.0 && mxIsNaN(D[i])==false) { 
        if ( D[i]<0.0 && mxIsInf(D[i])==false ) D[i]*=-1.0; /* mark voxel as visited */

        /* get xyz coordinates of the current voxel */
        ind2sub(i,&u,&v,&w,nL,xy,x);

        /* run mapping for all neighbours  */
        for (int n=0;n<26;n++) {

          /* get index and xyz coordinates of the current voxel */
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 


          /* Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || 
               (abs(nw-w)>1) || (ni==i) || mxIsNaN(D[ni]) || B[ni]>=0.5 )==false &&
               ( csf==0 || L[ni]<=fabs(L[i]+0.5) ) ) { 

            /* new distance for the neighbor */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);

            /* use DIN, if the new value is smaller than the actual value of the neighbor */
            if ( fabs(D[ni]) > DIN ) {
              if (D[ni]>0.0) nC++; /* count changes */
              D[ni] = -DIN; /* set new thickness - negative value so that the neigbor will also push */ 
              I[ni] = I[i]; /* index value of the current voxel (closest boundary index) */
              if ( nlhs>2 ) T[ni] = fabs(D[i]) + ND[n]; /* DIN; */ /* second output map with simple voxel distance */
            }
          }
        }
      }
    }


    /*
     * Backward direction:
     * Same as forward, with the small difference of demarking the start 
     * voxel at the end.
     */
    for (int i=nL-1;i>=0;i--) { 
      if ( D[i]<=0.0 && mxIsNaN(D[i])==false) { 
        if ( D[i]<0.0 && mxIsInf(D[i])==false ) D[i]*=-1.0; /* mark voxel as visited */

        ind2sub(i,&u,&v,&w,nL,xy,x);

        for (int n=0;n<26;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 

          /* Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || 
               (abs(nw-w)>1) || (ni==i) || mxIsNaN(D[ni]) || B[ni]>=0.5 )==false && 
               ( csf==0 || L[ni]<=fabs(L[i]+0.5) ) ) { 

            /* new distance */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);

            /* use DIN, if the actual value is larger */
            if ( fabs(D[ni])>DIN ) {
              if (D[ni]>0.0) nC++;
              D[ni] = -DIN; 
              I[ni] = I[i];
              if ( nlhs>2 ) T[ni] = fabs(D[i]) + ND[n]; /* DIN; */
            }
          }
        }

        /* 
         * Demarking the start voxels
         */
        if (D[i]==0.0) D[i]=-FINFINITY; 

      }
    }

    if ( verb && kll<30 ) 
      printf("    nC=%10d, kll=%4d, kllv=%4d\n",nC,kll,kllv);
    if ( nC==nCo ) { csf=0; nCo++; } /* further growing??? */
  }

 

  /* 
   * Correction of non-visited points due to miss-estimations of the 
   * exact boundary.
   */
  if ( verb ) printf("  Correction of unvisited points \n");  
  for (int i=0;i<nL;i++) {
    if ( mxIsInf(D[i]) && D[i]<0.0 ) D[i]=0.0;
    if ( D[i]<0.0 ) D[i]=-FINFINITY; 
  }


  /*
   * Euclidean distance estimation based on the nearest voxel in I.
   */
  for (int i=0;i<nL;i++) if (D[i]<0.0 && mxIsNaN(D[i])==false && mxIsInf(D[i])==false) D[i]*=-1.0;

  if ( verb ) printf("  Euclidean distance estimation \n"); 
  if ( euclid ) {
    for (int i=0;i<nL;i++) { 
      if ( mxIsNaN(B[i])==false && mxIsInf(B[i])==false && D[i]>0.0 && I[i]!=(unsigned int)i ) { 

        ni = (int) I[i];

        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        ind2sub(i ,&iu,&iv,&iw,nL,xy,x);

        /* standard euclidean distance between i and closest object point I[i] */
        dinu = (float)iu - (float)nu; dinu *= s1;
        dinv = (float)iv - (float)nv; dinv *= s2;
        dinw = (float)iw - (float)nw; dinw *= s3;
        DIN  = sqrt_float(sqr_float(dinu) + sqr_float(dinv) + sqr_float(dinw) - (3*0.5) ); /* 0.5 is the boundary vs. grid-distance */  

        /* For voxels that are not too close to the object the exact 
         * Euclidean distance should be estimated. For closer points
         * the previous distance value is good enough.
         */
        if ( 1 ) {
          /* Estimation of the exact PVE boundary by estimation of a point
           * next to the boundary voxel.
           */
          dinu /= DIN; dinv /= DIN; dinw /= DIN; /* normal vector in normal space */
          dinu /= s1;  dinv /= s2;  dinw /= s3;  /* normal vector for anisotropic space */
          WM = isoval(B,(float)nu + dinu,(float)nv + dinv,(float)nw + dinw,sizeL);

          if ( B[ni]!=WM ) {
            /* estimate new point before border to get the gradient based on this the exact HB */
            dcf = (B[ni] - 0.5) / ( B[ni] - WM );
            WMu = (float)nu + dinu*dcf; 
            WMv = (float)nv + dinv*dcf; 
            WMw = (float)nw + dinw*dcf; 
            // WM  = isoval(B,WMu,WMv,WMw,sizeL);

            /* new exact distance to interpolated boundary */ 
            dinu = (float)iu - WMu; dinu *= s1;
            dinv = (float)iv - WMv; dinv *= s2;
            dinw = (float)iw - WMw; dinw *= s3;
            DINE = sqrt_float(sqr_float(dinu) + sqr_float(dinv) + sqr_float(dinw)); 

            if ( false ) { // WM<0.45 || WM>0.55 ) { // 0.4 0.6
              WMu = (float)nu + 0.5*dinu*dcf; 
              WMv = (float)nv + 0.5*dinv*dcf; 
              WMw = (float)nw + 0.5*dinw*dcf; 
              WM  = isoval(B,WMu,WMv,WMw,sizeL);
            }
            
            /* use the new estimated distance only if it is a meanful */
            if ( DINE>0.0 && mxIsNaN(DINE)==false && mxIsInf(DINE)==false ) {
              D[i] = fmin( DIN , DINE );
            }
            else {
              /* use the voxelboundary corrected euclidean distance DIN
               * in case of larger values, or otherwise use the initial
               * value used before ... for higher values the speed map
               * lead to non euclidean distances! 
               */
              D[i] = DIN; 
            }
          }
        }
        else {
        /* simple euclidean distance without PVE for tests */
          D[i] = DIN; 
        }
      }
    }
  }

          
  
  /*
   * Final corrections
   */
  if ( verb ) printf("  Final corrections \n");  
  for (int i=0;i<nL;i++) { 
    /* correct for matlab index */
    if ( mxIsNaN(I[i])==false && I[i]>0 ) I[i]++; else I[i]=1;                        

    /* correction of non-visited or other incorrect voxels */ 
    if ( D[i]<0.0 || mxIsNaN(D[i]) || mxIsInf(D[i]) ) D[i]=0; // nanres; 

  } 
  
  
  /* final alignments to dynamic output variables */
  if (nlhs>1) {
    plhs[1]           = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL); 
    unsigned int *IO  = (unsigned int *) mxGetPr(plhs[1]);
    for (int i=0;i<nL;i++) IO[i] = I[i]; 
  }
  if (nlhs>2) {
    plhs[2]   = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); 
    float *TO = (float *) mxGetPr(plhs[2]);
    for (int i=0;i<nL;i++) TO[i] = T[i]; 
  }
  
  /* clear internal variables */
  /*
  mxDestroyArray(hlps[0]);
  mxDestroyArray(hlps[1]);
  mxDestroyArray(hlps[2]);
  */
  
  if ( verb ) printf("done. \n");   
}
