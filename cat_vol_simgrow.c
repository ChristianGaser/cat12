/* ______________________________________________________________________
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
/* #include "matrix.h" */

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

/* HELPFUNCTIONS */

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
  *z = (int)floor( (double)i / (double)sxy ) +1; 
   i = i % (sxy);
  *y = (int)floor( (double)i / (double)sy ) +1;        
  *x = i % sy + 1;
}

int sub2ind(int x,int y, int z, const int s[]) {
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}

float abs2(float n) {	if (n<0.0) return -n; else return n; }
float sign(float n) {	if (n<0.0) return 1; else return 0.0; }
float max2(float a, float b) { if (a>b) return a; else return b; }

/* MAINFUNCTION */
/* ROI, F, 1 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3)                                                      mexErrMsgTxt("ERROR:cat_vol_simgrow: not enough input elements\n");
  if (nrhs>5)                                                      mexErrMsgTxt("ERROR:cat_vol_simgrow: too many input elements.\n");
  if (nlhs>2)                                                      mexErrMsgTxt("ERROR:cat_vol_simgrow: too many output elements.\n");
  if (mxIsSingle(prhs[0])==0)                                      mexErrMsgTxt("ERROR:cat_vol_simgrow: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0)                                      mexErrMsgTxt("ERROR:cat_vol_simgrow: second input must be an 3d single matrix\n");
  if (mxIsDouble(prhs[2])==0 || mxGetNumberOfElements(prhs[2])!=1) mexErrMsgTxt("ERROR:cat_vol_simgrow: third input must one double value\n");
  
  if (nrhs==4 && mxIsDouble(prhs[3])==0)                           mexErrMsgTxt("ERROR:cat_vol_simgrow: fourth  input must be an double matrix\n");
  if (nrhs==4 && mxGetNumberOfElements(prhs[3])!=3)                mexErrMsgTxt("ERROR:cat_vol_simgrow: fourth input must have 3 Elements");
  if (nrhs==5 && mxIsDouble(prhs[4])==0)                           mexErrMsgTxt("ERROR:cat_vol_simgrow: fifth input must be an double matrix\n");
  if (nrhs==5 && mxGetNumberOfElements(prhs[4])!=2)                mexErrMsgTxt("ERROR:cat_vol_simgrow: fifth input must have 2 Elements");
    
  /* main information about input data (size, dimensions, ...) */
  const mwSize *sL = mxGetDimensions(prhs[0]);
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;

  const mwSize sSS[2] = {1,3}, sdsv[2] = {1,2}; 
  mxArray *SS  = mxCreateNumericArray(2,sSS, mxDOUBLE_CLASS,mxREAL); double*S  = mxGetPr(SS);
  mxArray *dsv = mxCreateNumericArray(2,sdsv,mxDOUBLE_CLASS,mxREAL); double*dd = mxGetPr(dsv);
  float dI=0.0; double*SEGd; if (nrhs>=3) {SEGd=mxGetPr(prhs[2]); dI=(float) SEGd[0];}; 
  if (nrhs<4) {S[0]=1.0; S[1]=1.0; S[2]=1.0;} else {S=mxGetPr(prhs[3]);}
  if (nrhs<5) {dd[0]=0.1; dd[1]=10.0;}    else {dd=mxGetPr(prhs[4]);}
    
  /* indices of the neighbor Ni (index distance) and euclidean distance NW */
  const int   NI[26]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x, 
                         -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
 
  int         ni,u,v,w,nu,nv,nw;
    
  /* input variables */
  float*ALAB = (float *)mxGetPr(prhs[0]);	/* label map */
  float*SEG  = (float *)mxGetPr(prhs[1]);	/* tissue map */
  
  /* main volumes - actual without memory optimization ... */
  mxArray *hlps[1]; 
  hlps[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);  
  float*DIST  = (float *)mxGetPr(hlps[0]);	/* distance map */
  
  /* output variables */
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* label map */ 
  if (nlhs>0) plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); /* tissue map (speed) */
  /* output variables */
  float*DISTO;  
  float*SLAB  = (float *)mxGetPr(plhs[0]);	/* label map */
  if (nlhs>0) DISTO = (float *)mxGetPr(plhs[1]);	/* distance map for ouput */
  
  
  int	  nCV   = 0; 		/* # voxel of interest (negative voxel that have to processed) */ 
	int 	kll   = 0;
	int   nC    = 1;
	int 	kllv  = 2000;
	float Dmax  = 0.2;
	float DISTN = 0.0;
	
  
	/* initialisation of parameter volumes */
	for (int i=0;i<nL;i++) { 
		SLAB[i] = ALAB[i]; 
    
    if ( mxIsNaN(SLAB[i]) ) SLAB[i] = 0.0;
   
		if ( SLAB[i] == 0.0 ) {  
      DIST[i] = FLT_MAX;
    }
		else {
			if ( SLAB[i] == -FLT_MAX || mxIsInf(SLAB[i]) ) {
        DIST[i] = -FLT_MAX;
      } 
			else { 									 
        DIST[i] = 0.0; 
        nCV++;
      }
    }
  }
  
  
  /* diffusion */ 
	while ( nCV>0 && kll<kllv && nC>0 ) {
		kll++; nC=0;
	
		for (int i=0;i<nL;i++) { 
			if ( (DIST[i]<=0.0) && ( DIST[i]!=-FLT_MAX ) ) { 
				if ( (DIST[i]<0.0) && ( DIST[i]!=-FLT_MAX ) )  DIST[i]=-DIST[i]; 
        nCV--; /* demark points - also the with zero distance */
      
				ind2sub(i,&u,&v,&w,xy,x); 
				for (int n=0;n<26;n++) {
					ni = i+NI[n]; ind2sub(ni,&nu,&nv,&nw,xy,x);
					
					if ( ((ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) )==false &&
               (abs2(SEG[i]-SEG[ni])<(dI*SEG[i]*SEG[i])) && ALAB[ni]==0.0 ) {
            DISTN = DIST[i] + abs2(SEG[i]-SEG[ni]) + 0.001;
            
						if ( (ALAB[ni]==0.0) && (DIST[ni]>0) && (DIST[ni]!=-FLT_MAX) && (abs2(DIST[ni])>abs2(DISTN)) && DISTN<Dmax)  {	
							if (DIST[ni]>0.0) nCV++; nC++;
							DIST[ni] = -DISTN;
							SLAB[ni] = SLAB[i];
						}
					}				
				}
				if (DIST[i]==0.0) DIST[i]=-FLT_MAX; /* demark start points */
			}
		}
	}


	for (int i=0;i<nL;i++) { if (DIST[i]==-FLT_MAX) DIST[i]=0.0; }
  if (nlhs>0) for (int i=0;i<nL;i++) { DISTO[i]=DIST[i]; } 

  
  /* clear internal variables */
  /*
  mxDestroyArray(hlps[0]);
  */
}
