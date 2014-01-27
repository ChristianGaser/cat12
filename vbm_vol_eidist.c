/* 
 * This function estimates the Euklidean distance D to the closest  
 * boundary voxel I given by the 0.5 isolevel in B that should contain
 * values between 0 and 1 to describe the boundary with PVE. 
 *
 * To align the closest voxel a modified Eikonal distances is estimated 
 * on the field L. 
 *
 * Voxels with NaN and -INF were ignored, and produce NaN in D, whereas 
 * in I the default index is used, becuase I has to be a integer matrix 
 * where NaN is not available. With setnan=0, the result of NAN voxel in 
 * D is changed to INF.
 * 
 *  [D,I] = vbm_vol_eidist(B,L,[vx_vol,euklid,setnan,verb])
 * 
 *  D         distance map   (3d-single-matrix)
 *  I         index map      (3d-uint32-matrix)
 *  B         boundary map   (3d-single-matrix)
 *  L         speed map      (3d-single-matrix)
 *  vx_vol    voxel size     (1x3-double-matrix): default=[1 1 1]
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
 *   [D,I]=vbm_vol_eidist(A,F,[1 1 1],1,0,1);
 *   ds('x2','',1,A,D,D,I,25)
 * _____________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id$ 
 */
 

#include "mex.h"   
#include "matrix.h"
#include "math.h"
#include "float.h"
#include "limits.h"

#ifdef _MSC_VER
  #define FINFINITY (FLT_MAX+FLT_MAX);
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FINFINITY 1.0f/0.0f;
  #define FNAN 0.0f/0.0f
#endif

float min(float a, float b) {	if (a<b) return a; else return b; }

/* 
 * Estimate x,y,z position of index i in an array size sx,sxy.
 */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy) {
  *z = (int)floor( i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( i / (double)sy ) +1;        
  *x = i % sy + 1;
}



/* 
 * Estimate index i of voxel x,y,z in an array size s.
 */
int sub2ind(int x,int y, int z, int s[]) {
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}


/* 
 * Read out the linear interpolatet value of a volume SEG with the size
 * s on the position x,y,z.
 */
float isoval(float SEG[], float x, float y, float z, int s[]){
	if (x<0) x=0; if (x>s[0]) x=(float)s[0];
	if (y<0) y=0; if (y>s[1]) y=(float)s[1];
	if (z<0) z=0; if (z>s[2]) z=(float)s[2];

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
    
  float seg=0;
  for (i=0; i<8; i++) seg = seg + N[i]*W[i];
  return seg;  
}


/* 
 * MAINFUNCTION [D,I] = vbm_vol_eidist(B,L,vx_vol,euklid,verb)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* 
   * Check input 
   */
  if (nrhs<2) mexErrMsgTxt("ERROR:vbm_vol_eidist:  not enought input elements\n");
  if (nrhs>6) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many input elements.\n");
  if (nlhs>3) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: second input must be an 3d single matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: third  input must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) {printf("ERROR:vbm_vol_eidist: third input must have 3 Elements"); };
  if (nrhs==4 && mxIsDouble(prhs[3])==0 &&  mxGetNumberOfElements(prhs[3])!=1) {printf("ERROR:vbm_vol_eidist: fourth input must be one double value"); }; 
  if (nrhs==5 && mxIsDouble(prhs[4])==0 &&  mxGetNumberOfElements(prhs[4])!=1) {printf("ERROR:vbm_vol_eidist: fifth input must be one double value"); }; 
  if (nrhs==6 && mxIsDouble(prhs[5])==0 &&  mxGetNumberOfElements(prhs[5])!=1) {printf("ERROR:vbm_vol_eidist: sixth input must be one double value"); }; 
  
    
  /* set default variables */
  int euklid=1; double*deuklid; if (nrhs>=4) {deuklid=mxGetPr(prhs[3]); euklid = (int) deuklid[0]>0.5;};    
  int setnan=1; double*dsetnan; if (nrhs>=5) {dsetnan=mxGetPr(prhs[4]); setnan = (int) dsetnan[0]>0.5;};    
  int verb=0;   double*dverb;   if (nrhs>=6) {dverb=mxGetPr(prhs[5]);   verb   = (int) dverb[0]>0.5;};    

  float nanres;
  if ( setnan>=1 ) nanres = FNAN; else nanres = FINFINITY; 

  /* main input variables B (boundary map) and L (speed map) */
  float*B	 = (float *)mxGetPr(prhs[0]);	
  float*L  = (float *)mxGetPr(prhs[1]);	
                 
  
  /* main informations about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]); int sizeL[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;

  /* 
   * Voxel dimension:
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
  const float   s123 = sqrt(s12*s12 + s3*s3); /* xyz - voxel size */
  const int     nr = nrhs;
  
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3,  s12, s12,s12,s12,   s13,  s13, s13, s13,   s23,  s23, s23, s23,     s123,   s123,   s123,   s123,   s123,  s123,  s123,  s123};
  const int kllv = (sL[0]+sL[1]+sL[2]);
  
  /* other variables */
  float 	dinu, dinv, dinw, dcf, WMu, WMv, WMw, WM, DN, DIN;	/* distance and intensity variables */
  int     i,j,n,ii,ni,DNi,u,v,w,nu,nv,nw,iu,iv,iw;          	/* xyz and index-values of a voxel, one neighbor and the nearest boundary voxel */
  int	  	nC=1, nCVL=0, kll=0, nCVo, kllo;                    /* runtime variables */
  
  
  /* 
   * Create main output volumes and variables D (distance map) and I (index map)
   */ 
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); 
  plhs[1] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL); 
  float*D  				= (float *)mxGetPr(plhs[0]);				
  unsigned int*I  = (unsigned int *)mxGetPr(plhs[1]);	
  
 
  /*
   * Display Initial Parameter
   */
  if ( verb ) printf("vbm_vol_eidist.c verbing mode:\n  Initialize Parameter: \n");
  if ( verb ) printf("    size(B) = %d %d %d\n",sL[0],sL[1],sL[2]); 
  if ( verb ) printf("    vx_vol  = %0.0f %0.0f %0.0f\n",s1,s2,s3); 
  if ( verb ) printf("    euklid  = %d\n",(int) euklid); 
  if ( verb ) printf("    setnan  = %d\n",(int) setnan); 
  
    
  /* 
   * Check input values.
   */
  int vx=0;
 	for (i=0;i<nL;i++) { 
    if ( B[i]>=0.5 ) vx++;                                              /* count object voxel */
    
    if ( (L[i]<=0 || mxIsNaN(L[i])) && B[i]<0.5) B[i] = FNAN;          /* ignore voxel that canot be visited */
    
    if ( mxIsInf(B[i]) && B[i]<0 ) B[i] = FNAN;                        /* voxel to ignore */
    if ( B[i]>1 ) B[i] = 1;                                             /* normalize object */
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
    for (i=0;i<nL;i++) { D[i]=nanres; I[i]=i+1; } 
    if ( verb ) printf("done.\n");
    printf("WARNING:vbm_vol_eidist: Found no object for distance estimation!\n");
    return;
  }
    
  
  /* 
   *  
   * Eikonal distance initialisation:
   * In 3D this can be really complex, especialy for anisotripic volumes. 
   * So only a simple approximation by the PVE of a voxle is used. 
   */ 
 if ( verb ) printf("\n  Initialize Eikonal distance estimation and index alignment: ");
 for (i=0;i<nL;i++) { 
    if ( D[i]>0 && mxIsNaN(D[i])==0) { 

      
		 	if ( B[i]>=0.5 ) {
        D[i] = 0; 
      
        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(i,&u,&v,&w,xy,x);
          ind2sub(ni,&nu,&nv,&nw,xy,x); 

          if ( ( (ni<0) || (ni>nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==i) )==0 && B[ni]<0.5) {
            if ( B[i]==1 )
              DIN = ND[n] * (0.5 + min(0.5-B[ni],0.5));
            else
              DIN = ND[n] * (0.5 + min(1-B[i],0.5));

            if ( fabs(D[ni])>DIN ) {
              D[ni] = -DIN; 
              I[ni] = I[i];
            }
          }
        }
      }
    }
	}

 
  /* 
   * iterative Eikonal distance estimation
   */ 
  if ( verb ) printf("done.\n  Eikonal distance estimation and index alignment: ");
	while ( nC>0 && kll<kllv ) {
    
    kll++; 
    nC=0;
    
    /*
     * Forward direction:
     */
    for (i=0;i<nL;i++) { 
      if ( D[i]<=0 && mxIsNaN(D[i])==0) { 
        if ( D[i]<0 && mxIsInf(D[i])==0 ) D[i]*=-1; /* mark voxel as visited */
        
        ii=(int)i;
        ind2sub(ii,&u,&v,&w,xy,x);

        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,xy,x); 
          
          /*
           * Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==ii) || 
            mxIsNaN(D[ni]) || B[ni]>=0.5 )==0  && L[ni]<=(L[i]+0.01)  ) { 

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
     * voxel at the end
     */
    for (i=nL;i>=0;i--) { 
      if ( D[i]<=0 && mxIsNaN(D[i])==0) { 
        if ( D[i]<0 && mxIsInf(D[i])==0 ) D[i]*=-1; /* mark voxel as visited */
        
        ii=(int)i;
        ind2sub(ii,&u,&v,&w,xy,x);

        for (n=0;n<27;n++) {
          ni = i + NI[n];
          ind2sub(ni,&nu,&nv,&nw,xy,x); 
          
          /*
           * Only process this part for real neighbor indices and
           * only if L is lower (to avoid region-growing over sulci)! 
           */
          if ( ( (ni<0) || (ni>nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (ni==ii) || 
            mxIsNaN(D[ni]) || B[ni]>=0.5 )==0  && L[ni]<=(L[i]+0.01)  ) { 

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

    if ( verb ) printf("\n    nC=%10d, kll=%4d, kllv=%4d",nC,kll,kllv);
	}

  
  /* 
   * Correction of non-visited points due to miss-estimations of the 
   * exact boundary.
  */
  if ( verb ) printf("\n  Correction of non-visited points: ");  
	for (i=0;i<nL;i++) {
    if ( mxIsInf(D[i]) && D[i]<0 ) D[i]=0;
		if ( D[i]<0 ) D[i]=-FINFINITY; 
	}
  
     
  /*
   * Euclidean distance estimation based on the nearest voxel in I.
   */
  if ( verb ) printf("done.\n  Euclidean distance estimation: "); 
  if ( euklid ) {
    for (i=0;i<nL;i++) { 
      if ( mxIsNaN(B[i])==0 && D[i]>0 && I[i]!=i) { 

        ni = (int) I[i];
        
        ind2sub(ni,&nu,&nv,&nw,xy,x);
        ind2sub(i ,&iu,&iv,&iw,xy,x);

        /* standard euclidean distance */
        dinu = (float)nu-(float)iu; dinu *= s1;
        dinv = (float)nv-(float)iv; dinv *= s2;
        dinw = (float)nw-(float)iw; dinw *= s3;
        DIN  = sqrt(pow(dinu,2) + pow(dinv,2) + pow(dinw,2)); 

        if ( 1 ) {
          /* difference distance with normaliced vectors */
          dinu /= DIN; dinv /= DIN; dinw /= DIN; /* normal vector in normal space */
          dinu /= s1;  dinv /= s2;  dinw /= s3;  /* normal vector for anisotropic space */
          WM = isoval(B,(float)nu + dinu,(float)nv + dinv,(float)nw + dinw,sizeL);

          if ( B[ni]!=WM && B[i]!=WM ) {
            /* estimate new point bevor border to get the gradient based on this the exact HB */
            dcf = (B[ni] - 0.5) / ( B[ni] - WM );
            WMu = (float)nu + dinu*dcf; WMv = (float)nv + dinv*dcf; WMw = (float)nw + dinw*dcf; 
            WM  = isoval(B,WMu,WMv,WMw,sizeL);

            /* new exact distance to interpolatet boundary */
            if ( WM>0.45 && WM<0.65 ) {
              dinu = (float)iu-WMu; dinv = (float)iv-WMv; dinw = (float)iw-WMw;
              D[i]  = sqrt(pow(dinu*s1,2) + pow(dinv*s2,2) + pow(dinw*s3,2));   
            }
          }
        }
      }
    }
  }
       
          
	/*
   * Final corrections
   */
  if ( verb ) printf("done.\n  Final corrections: ");  
	for (i=0;i<nL;i++) { 
    /* correct for matlab index */
		if ( mxIsNaN(I[i])==0 && I[i]>0 ) I[i]++; else I[i]=1;                        

    /* correction of non-visited or other incorrect voxels */	
    if ( D[i]<0 || mxIsNaN(D[i]) || mxIsInf(D[i]) ) D[i]=nanres; 
    
	} 
  if ( verb ) printf("done.\n");  
}