/* 
 * This function estimates the Euklidean distance D to the closest  
 * boundary voxel I given by the 0.5 isolevel in B that should contain
 * values between 0 and 1 to describe the boundary with PVE. 
 *
 * To align the closest voxel a modified Eikonal distances is estimated 
 * on the field L. 
 * For a correct side alignment a harder option "csf" is used that avoid
 * the diffusion to voxels with greater value of L (L[i]>L[ni] for
 * diffusion).
 *
 * Voxels with NaN and -INF were ignored, and produce NaN in D, whereas 
 * in I the default index is used, becuase I has to be a integer matrix 
 * where NaN is not available. With setnan=0, the result of NAN voxel in 
 * D is changed to INF.
 * 
 *  [D,I] = vbm_vol_eidist(B,L,[vx_vol,euklid,csf,setnan,verb])
 * 
 *  D         distance map   (3d-single-matrix)
 *  I         index map      (3d-uint32-matrix)
 *  B         boundary map   (3d-single-matrix)
 *  L         speed map      (3d-single-matrix)
 *  vx_vol    voxel size     (1x3-double-matrix): default=[1 1 1]
 *  csf       option         (1x1-double-value): 0-no; 1-yes; default=1
 *  euklid    option         (1x1-double-value): 0-no; 1-yes; default=1
 *  setnan    option         (1x1-double-value): 0-no; 1-yes; default=1
 *  verb      option         (1x1-double-value): 0-no, 1-yes; default=0
 * 
 *
 * Small test examples:
 * 1) 
 *   A=zeros(50,50,50,'single'); A(20:30,5:15,20:30)=1; A(20:30,35:45,20:30)=1; 
 *   A=smooth3(A); A(1:5,1:25,:)=nan; A(1:5,26:50,:)=-inf; A(45:50,26:50,:)=inf;
 *   F=ones(size(A),'single'); F(10:40,20,:)=0.5; F(40,10:40,:)=0;
 *
 * 2) 1D-examples
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,6,:)=0.2; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *   A=zeros(10,20,10,'single'); A(:,1:5,:)=1; A(:,15:end,:)=nan; F=ones(size(A),'single');
 *
 *   [D,I]=vbm_vol_eidist(A,F,[1 1 1],1,1,0,1);
 *   ds('x2','',1,A,D,D,I,25)
 * _____________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id$ 
 */


/*
 * ToDo:
 * - the main goal is the replacement of the vbm_vol_vbdist and vbm_vol_eikonal3 functions
 * - variable class input
 * - 
 */
 

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"
#include "limits.h"
#include "algorithm.h" 

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
 * There is a small test output in the main function after the parameter 
 * initialization ( search for "ind2subtest" ).
 */
void ind2sub(int i, int *x, int *y, int *z, int snL, int sxy, int sy) {
  if (i<0) i=0; 
  if (i>=snL) i=snL-1;
  
  *z = (int)floor( i / (double)sxy ) ; 
   i = i % (sxy);
  *y = (int)floor( i / (double)sy ) ;        
  *x = i % sy ;
}
/* 
 * Estimate index i of a voxel x,y,z in an array size s.
 * See also for ind2sub.
 */
int sub2ind(int x,int y, int z, int s[]) {
  if (x<0) x=0; if (x>s[0]-1) x=s[0]-1; 
  if (y<0) y=0; if (y>s[1]-1) y=s[1]-1; 
  if (z<0) z=0; if (z>s[2]-1) z=s[2]-1; 
  
  /* int i=(z)*s[0]*s[1] + (y)*s[0] + (x) - 1; */
  return (z)*s[0]*s[1] + (y)*s[0] + (x);
}


float min(float a, float b) {
  if (a<b) return a; else return b; 
}


float fpow(float x, float y) {
  return (float) pow((double) x,(double) y); 
}


float fsqrt(float x) {
  return (float) sqrt((double) x); 
}


/* 
 * Read out the linear interpolatet value of a volume SEG with the size
 * s on the position x,y,z (c-notation). See also ind2sub for details of
 * the c-noation.
 */
float isoval(float SEG[], float x, float y, float z, int s[]){

  int i;
  float fx = floor(x),   fy = floor(y),	  fz = floor(z);
  float cx = floor(x+1), cy = floor(y+1), cz = floor(z+1);
    
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
    
  float seg=0, n=0;
  for (i=0; i<8; i++) {
    if ( mxIsNaN(N[i])==0 || mxIsInf(N[i])==0 )
      seg = seg + N[i] * W[i]; n+= W[i];
  }
  if ( n>0 )
    return seg/n; 
  else 
    return FNAN;
}


/* 
 * MAINFUNCTION [D,I] = vbm_vol_eidist(B,L,vx_vol,euklid,csf,setnan,verb])
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /* 
   * Check input 
   */
  if (nrhs<2) mexErrMsgTxt("ERROR:vbm_vol_eidist:  not enought input elements\n");
  if (nrhs>7) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many input elements.\n");
  if (nlhs>3) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: second input must be an 3d single matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: third  input must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) {printf("ERROR:vbm_vol_eidist: third input must have 3 Elements"); };
  if (nrhs==4 && mxIsDouble(prhs[3])==0 &&  mxGetNumberOfElements(prhs[3])!=1) {printf("ERROR:vbm_vol_eidist: fourth input must be one double value"); }; 
  if (nrhs==5 && mxIsDouble(prhs[4])==0 &&  mxGetNumberOfElements(prhs[4])!=1) {printf("ERROR:vbm_vol_eidist: fifth input must be one double value"); }; 
  if (nrhs==6 && mxIsDouble(prhs[5])==0 &&  mxGetNumberOfElements(prhs[5])!=1) {printf("ERROR:vbm_vol_eidist: sixth input must be one double value"); }; 
  if (nrhs==7 && mxIsDouble(prhs[6])==0 &&  mxGetNumberOfElements(prhs[6])!=1) {printf("ERROR:vbm_vol_eidist: seventh input must be one double value"); }; 
  
    
  /* 
   * set default variables 
   */
  int csf=1;    double*dcsf;    if (nrhs>=4) {dcsf=mxGetPr(prhs[3]);    csf    = (int) dcsf[0]>0.5;};    
  int euklid=1; double*deuklid; if (nrhs>=5) {deuklid=mxGetPr(prhs[4]); euklid = (int) deuklid[0]>0.5;};    
  int setnan=1; double*dsetnan; if (nrhs>=6) {dsetnan=mxGetPr(prhs[5]); setnan = (int) dsetnan[0]>0.5;};    
  int verb=0;   double*dverb;   if (nrhs>=7) {dverb=mxGetPr(prhs[6]);   verb   = (int) dverb[0]>0.5;};    
  float nanres;  if ( setnan>=1 ) nanres = FNAN; else nanres = FINFINITY; 

  
  /* 
   * main input variables B (boundary map) and L (speed map) 
   */
  float*B	 = (float *)mxGetPr(prhs[0]);	
  float*L  = (float *)mxGetPr(prhs[1]);	
                 
  
  /* 
   * main informations about input data (size, dimensions, ...) 
   */
  const mwSize *sL = mxGetDimensions(prhs[0]); int sizeL[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     nD = mxGetNumberOfElements(prhs[1]);
  if (nL!=nD) mexErrMsgTxt("ERROR:vbm_vol_eidist: images must have the same number of elements.\n");
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;

  
  /* 
   * Voxel dimension and standard distance in the N26 neighborhood.
   * Indices of the neighbor NI (index distance) and euclidean distance ND. 
   * Set default voxel size=[1 1 1] or use input.
   */
  const int sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double   *S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[2]);}
  
  float s1 = fabs((float)S[0]),s2 = fabs((float)S[1]),s3 = fabs((float)S[2]); /* x,y,z - voxel size */
  const float   s12  = sqrt( s1*s1  + s2*s2); /* xy  - voxel size */
  const float   s13  = sqrt( s1*s1  + s3*s3); /* xz  - voxel size */
  const float   s23  = sqrt( s2*s2  + s3*s3); /* yz  - voxel size */
  const float   s123 = sqrt(s12*s12 + s3*s3); /* nL - voxel size */
  
  const int   NI[]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3,  s12, s12,s12,s12,   s13,  s13, s13, s13,   s23,  s23, s23, s23,     s123,   s123,   s123,   s123,   s123,  s123,  s123,  s123};
  const int kllv = (sL[0]+sL[1]+sL[2]);

  
  /* 
   * other variables 
   */
  float 	dinu, dinv, dinw, dcf, WMu, WMv, WMw, WM, DIN;	/* distance and intensity variables */
  int     i,n,ii,ni,u,v,w,nu,nv,nw,iu,iv,iw;          	/* nL and index-values of a voxel, one neighbor and the nearest boundary voxel */
  int	  	nC=sL[0]*sL[1]*sL[2], nCo=FINFINITY;                /* runtime variables */
  int     kll=0;                                               /* runtime variables */
  int     fast=1;                                             /* stop if number of unvisited points stay constant */

  
  /* 
   * Create main output volumes and variables D (distance map) and I (index map)
   */ 
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); 
  plhs[1] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL); 
  float *D  				= (float *)mxGetPr(plhs[0]);				
  unsigned int *I  = (unsigned int *)mxGetPr(plhs[1]);	
  
 
  /*
   * Display Initial Parameter
   */
  if ( verb ) printf("\nvbm_vol_eidist.c debuging mode:\n  Initialize Parameter: \n");
  if ( verb ) printf("    size(B) = %d %d %d\n",sL[0],sL[1],sL[2]); 
  if ( verb ) printf("    vx_vol  = %0.0f %0.0f %0.0f\n",s1,s2,s3); 
  if ( verb ) printf("    euklid  = %d\n",(int) euklid); 
  if ( verb ) printf("    setnan  = %d\n",(int) setnan); 
  
  
  /* 
   * Small test of the ind2sub and sub2ind functions
   * ind2subtest
   */
  if ( 0 ) {
    i=0; ind2sub(i,&u,&v,&w,nL,xy,x); printf("i=%10d > nL=%3d,%3d,%3d\n",i,u,v,w);
    i=1; ind2sub(i,&u,&v,&w,nL,xy,x); printf("i=%10d > nL=%3d,%3d,%3d\n",i,u,v,w);
    i=sL[0]*sL[1]*sL[2] - 1 ; ind2sub(i,&u,&v,&w,nL,xy,x); printf("i=%10d > nL=%3d,%3d,%3d\n",i,u,v,w);
    i=sL[0]*sL[1]*sL[2]  ; ind2sub(i,&u,&v,&w,nL,xy,x); printf("i=%10d > nL=%3d,%3d,%3d\n",i,u,v,w);

    u=0; v=0; w=0; i=sub2ind(u,v,w,sizeL); printf("i=%10d < nL=%3d,%3d,%3d\n",i,u,v,w);
    u=1; v=0; w=0; i=sub2ind(u,v,w,sizeL); printf("i=%10d < nL=%3d,%3d,%3d\n",i,u,v,w);
    u=sL[0]-1; v=sL[1]-1; w=sL[2]-1; i=sub2ind(u,v,w,sizeL); printf("i=%10d < nL=%3d,%3d,%3d\n",i,u,v,w);
    u=sL[0]; v=sL[1]-1; w=sL[2]-1; i=sub2ind(u,v,w,sizeL); printf("i=%10d < nL=%3d,%3d,%3d\n",i,u,v,w);
  }
    
 
  /* 
   * Check input values.
   */
  int vx=0;
 	for (i=0;i<nL;i++) { 
    if ( B[i]>=0.5 ) vx++;                                              /* count object voxel */
    
    if ( (L[i]<=0 || mxIsNaN(L[i])) && B[i]<0.5) B[i] = FNAN;           /* ignore voxel that canot be visited */
    
    if ( mxIsInf(B[i]) && B[i]<0 ) B[i] = FNAN;                         /* voxel to ignore */
    if ( B[i]>1 ) B[i] = 1;                                             /* normalize object */
    if ( B[i]<0 ) B[i] = 0;                                             /* normalize object */
    I[i] = (unsigned int) i;                                            /* initialize index map */
    
    if ( L[i]<=0 || mxIsNaN(B[i]) || mxIsNaN(L[i]) ) 
      D[i] = FNAN;
    else 
      D[i] = FINFINITY;
  }
  
    
  /* 
   * Check if there is a object and return FINFINITY and i=i+1 (for matlab) 
   * for all points if there is no object.
   */
 	if ( vx==0 ) { 
    for (i=0;i<nL;i++) { D[i]=nanres; I[i]=(unsigned int)i+1; } 
    printf("WARNING:vbm_vol_eidist: Found no object for distance estimation!\n");
    return;
  }
    

  /* 
   * Eikonal distance initialisation:
   * In 3D this can be really complex, especialy for anisotripic volumes. 
   * So only a simple approximation by the PVE of a voxle is used. 
   */ 
 if ( verb ) printf("  Initialize Eikonal distance estimation and index alignment \n");
 for (i=0;i<nL;i++) { 
    if ( D[i]>0 && mxIsNaN(D[i])==0) { 
  
		 	if ( B[i]>=0.5 ) {
        D[i] = 0; 
      
        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(i,&u,&v,&w,nL,xy,x);
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 

          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==i) )==0 && B[ni]<0.5) {
            if ( B[ni]==0 ) {
              DIN = ND[n] * (1.5-B[i]);
              
              if ( fabs(D[ni])>DIN ) {
                D[ni] = -DIN; 
                I[ni] = I[i];
              }
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
	while ( nC>0 && kll<kllv && (nC!=nCo || fast==0)) {
    
    kll++; 
    nCo=nC;
    nC=0;
    
    
    /*
     * Forward direction:
     */
    for (i=0;i<nL;i++) { 
      if ( D[i]<=0 && mxIsNaN(D[i])==0) { 
        if ( D[i]<0 && mxIsInf(D[i])==0 ) D[i]*=-1; /* mark voxel as visited */
        
        ii=(int)i;
        ind2sub(ii,&u,&v,&w,nL,xy,x);

        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 
          
          
          /*
           * Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==ii) || 
            mxIsNaN(D[ni]) || B[ni]>=0.5 )==0 && ( csf==0 || L[ni]<=(L[i]+0.1) )  ) { 

            /* new distance */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);
            
            /* use DIN, if the actual value is larger	*/
            if ( fabs(D[ni])>DIN ) {
              if (D[ni]>0) nC++;
              D[ni] = -DIN; 
              I[ni] = I[i];
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
    for (i=nL-1;i>=0;i--) { 
      if ( D[i]<=0 && mxIsNaN(D[i])==0) { 
        if ( D[i]<0 && mxIsInf(D[i])==0 ) D[i]*=-1; /* mark voxel as visited */
      
        ii=(int)i;
        ind2sub(ii,&u,&v,&w,nL,xy,x);

        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,nL,xy,x); 
          
          /*
           * Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==ii) || 
            mxIsNaN(D[ni]) || B[ni]>=0.5 )==0 && ( csf==0 || L[ni]<=(L[i]+0.1) )  ) { 

            /* new distance */
            DIN = fabs(D[i]) + ND[n] / (FLT_MIN + L[ni]);
            
            /* use DIN, if the actual value is larger	*/
            if ( fabs(D[ni])>DIN ) {
              if (D[ni]>0) nC++;
              D[ni] = -DIN; 
              I[ni] = I[i];
            }
          }
        }
        
        
        
        /* 
         * Demarking the start voxels
         */
        if (D[i]==0) D[i]=-FINFINITY; 
        
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
  if ( verb ) printf("  Correction of non-visited points \n");  
	for (i=0;i<nL;i++) {
    if ( mxIsInf(D[i]) && D[i]<0 ) D[i]=0;
		if ( D[i]<0 ) D[i]=-FINFINITY; 
	}
  
     
  /*
   * Euclidean distance estimation based on the nearest voxel in I.
   */
  if ( verb ) printf("  Euclidean distance estimation \n"); 
  if ( 1 && euklid ) {
    for (i=0;i<nL;i++) { 
      if ( mxIsNaN(B[i])==0 && mxIsInf(B[i])==0 && D[i]>0 && I[i]!=(unsigned int)i ) { 

        ni = (int) I[i];
        
        ind2sub(ni,&nu,&nv,&nw,nL,xy,x);
        ind2sub(i ,&iu,&iv,&iw,nL,xy,x);

        /* standard euclidean distance between i and closes object point I[i] */
        dinu = (float)iu-(float)nu; dinu *= s1;
        dinv = (float)iv-(float)nv; dinv *= s2;
        dinw = (float)iw-(float)nw; dinw *= s3;
        DIN  = fsqrt(fpow(dinu,2) + fpow(dinv,2) + fpow(dinw,2)); 

        
        /* For voxels that are not to close to the object the exact 
         * Euklidean distance should be estimated. For closer points
         * the previous distance value is good enought.
         */
        if ( 1 || DIN>s123 ) {
          /* Estimation of the exact PVE boundary by estimation of a point
           * next to the boundary voxel.
           */
          dinu /= DIN; dinv /= DIN; dinw /= DIN; /* normal vector in normal space */
          dinu /= s1;  dinv /= s2;  dinw /= s3;  /* normal vector for anisotropic space */
          WM = isoval(B,((float)nu) + dinu,((float)nv) + dinv,((float)nw) + dinw,sizeL);
          
          if ( B[ni]!=WM ) {
            /* estimate new point bevor border to get the gradient based on this the exact HB */
            dcf = (B[ni] - 0.5) / ( B[ni] - WM );
            WMu = ((float)nu) + dinu*dcf; WMv = ((float)nv) + dinv*dcf; WMw = ((float)nw) + dinw*dcf; 
            WM  = isoval(B,WMu,WMv,WMw,sizeL);
            if ( WM<0.4 || WM>0.6 ) {
              WM = isoval(B,((float)nu) + 2*dinu,((float)nv) + 2*dinv,((float)nw) + 2*dinw,sizeL);
              dcf = (B[ni] - 0.5) / ( B[ni] - WM );
              WMu = ((float)nu) + dinu*dcf; WMv = ((float)nv) + dinv*dcf; WMw = ((float)nw) + dinw*dcf; 
              WM  = isoval(B,WMu,WMv,WMw,sizeL);
            }
            
            /* new exact distance to interpolated boundary */ 
            if ( WM>0.4 && WM<0.6 ) {
              dinu = iu-WMu; dinv = iv-WMv; dinw = iw-WMw;
              DIN  = fsqrt(fpow(dinu*s1,2) + fpow(dinv*s2,2) + fpow(dinw*s3,2));   

              if ( DIN>0 && mxIsNaN(DIN)==0 && mxIsInf(DIN)==0 )
                D[i] = DIN;
            }
            /* some test output for error handling
            else {
              D[i]= -1; 
              printf("-1: %0.4f  -  %3d,%3d,%3d   -  %3.2f,%3.2f,%3.2f  -  %3.2f,%3.2f,%3.2f   \n",WM,nu,nv,nw,dinu,dinv,dinw,WMu,WMv,WMw);
            }
            */
          }
          /* some test output for error handling
          else {
            D[i] = -2;
            printf("-2: %0.4f  -  %3d,%3d,%3d   -  %3.2f,%3.2f,%3.2f  -  %3.2f,%3.2f,%3.2f   \n",WM,nu,nv,nw,dinu,dinv,dinw,((float)nu) + dinu,((float)nv) + dinv,((float)nw) + dinw);
          }
          */
        }
      }
    }
  }
       
          
	/*
   * Final corrections
   */
  if ( verb ) printf("  Final corrections \n");  
	for (i=0;i<nL;i++) { 
    /* correct for matlab index */
		if ( mxIsNaN(I[i])==0 && I[i]>0 ) I[i]++; else I[i]=1;                        

    /* correction of non-visited or other incorrect voxels */	
    if ( D[i]<0 || mxIsNaN(D[i]) || mxIsInf(D[i]) ) D[i]=nanres; 
    
	} 
  if ( verb ) printf("done. \n");   
}
