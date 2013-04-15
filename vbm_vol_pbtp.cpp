/* quasi-euclidean distance calculation
 * _____________________________________________________________________________
 * [GMT,RPM,WMD,CSFD,II] = eqdist(SEG,WMD,CSFD[,opt])
 *
 * SEG  = (single) segment image with low and high boundary bd
 * GMT  = (single) thickness image
 * RPM  = (single) radial position map
 * WMD  = (single) CSF distance map
 * CSFD = (single) CSF distance map
 * II  = (uint32) index of the inner (WM)  boundary voxel
 *
 * opt.bd   = (single) [low,high] boundary values (default 1.5 and 2.5)
 * opt.CSFD = calculate CSFD
 * opt.PVE  = use PVE informations (0=none,1=fast,2=exact)
 *
 * TODO:
 *  - eikonal distance for subsegmentation (region growing)
 *  - own labeling (
 * _____________________________________________________________________________
 * Robert Dahnke 2009_11
 * Center of Neuroimaging 
 * University Jena
 */

#include "mex.h"   
#include "matrix.h"
#include "math.h"


// estimate minimum of A and its index in A    
void pmin(const float A[], int sA, float & minimum, int & index) {
  minimum=INFINITY; index=0; 
  for(int i=0;i<sA;i++) {
    if ((A[i]>0) && (minimum>A[i])) { 
      minimum = A[i];
      index   = i;
    }
  }
}
float mean3(float a, float b, float c) {
  return (a+b+c)/3;
}


// simple maximum and absolute value
float max(float n1, float n2) { if (n1>n2) return n1;	else return n2; }
float abs2(float n) {	if (n<0) return -n; else return n; }

/*
void pmax(const float GMT[],const float RPM[],const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA, float & maximum, int & index) {
  float T[27]; for (int i=0;i<27;i++) T[i]=-1; float n=0.0; maximum=WMD; index=0; float maximum2=WMD; 
  //maximum=WMD; index=0; //printf("%d ",sizeof(A)/8);
  for(int i=0;i<=sA;i++) {
    
    if (  (GMT[i]<INFINITY) && (maximum < GMT[i]) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.65)>WMD) && SEG[i]>1.5 && SEGI>=2)     //work: +-0.4 to +-0.5
      {
        n++; 
        maximum = GMT[i]; 
        index = i;
      }
//    if (  (GMT[i]<INFINITY) && (GMT[i]>WMD) && ((RPM[i]-ND[i]-0.3)<WMD) && ((RPM[i]-ND[i]+0.9)>WMD) && SEG[i]>1.5 && SEGI>=2) 
//      {
//        maximum2 = (maximum2 + GMT[i])/2; 
//      } 
    //if (maximum<=maximum2+0.5) maximum = 0.9*maximum + 0.1*maximum2;
  }
}

*/


// get all values ot the voxels witch are in WMD-range (children of this voxel)  
void pmax(const float GMT[],const float RPM[],const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA, float & maximum, int & index) {
  float T[27]; for (int i=0;i<27;i++) T[i]=-1; float n=0.0; maximum=WMD; index=0; 
  
  for(int i=0;i<=sA;i++) {
    if (  (GMT[i]<INFINITY) && (maximum < GMT[i]) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.5)>WMD) && (SEGI)>=SEG[i] && SEG[i]>1 && SEGI>1.66)  // vorher ohne // 1.15 und 0.6
    //if (  (GMT[i]<INFINITY) && (maximum < GMT[i]) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.65)>WMD) )     
      { maximum = GMT[i]; index = i; }
  }
  float maximum2=maximum; float m2n=0; 
  for(int i=0;i<=sA;i++) {
    if (  (GMT[i]<INFINITY) && (GMT[i]>WMD) && ((RPM[i]-ND[i]*1.25)<=WMD) && ((RPM[i]-ND[i]*0.5)>WMD) && (SEGI)>=SEG[i] && SEG[i]>1.00 && SEGI>1.66)   // 1.00 0.5 // 1.50 0.75 && SEG[i]>=1.5 && SEGI>1.5
      { maximum2 = maximum2 + GMT[i]; m2n++; } 
  }
  if ( m2n > 0 )  maximum = (maximum2 - maximum)/m2n; 
  }




// estimate x,y,z position of index i in an array size sx,sxy=sx*sy...
void ind2sub(int i,int &x,int &y, int &z, int sxy, int sy) {
  z = floor( i / (sxy) ) +1; 
  i = i % (sxy);
  y = floor( i / sy ) +1;        
  x = i % sy + 1;
}

int sub2ind(int x,int y, int z, int s[]) {
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}

// read out linear interpolatet value of the volume
float isoval(float*SEG,float x, float y, float z, int sSEG[]){
  float fx = floor(x),   fy = floor(y),   fz = floor(z);
  float cx = floor(x+1), cy = floor(y+1), cz = floor(z+1);
  
  float wfx = cx-x, wfy = cy-y, wfz = cz-z;
  float wcx = x-fx, wcy = y-fy, wcz = z-fz;

  // value of the 8 neighbors and there distance weight
  float N[8], W[8];  
  N[0]=SEG[sub2ind(int(fx),int(fy),int(fz),sSEG)];  W[0]=wfx * wfy * wfz; 
  N[1]=SEG[sub2ind(int(cx),int(fy),int(fz),sSEG)];  W[1]=wcx * wfy * wfz;
  N[2]=SEG[sub2ind(int(fx),int(cy),int(fz),sSEG)];  W[2]=wfx * wcy * wfz;
  N[3]=SEG[sub2ind(int(cx),int(cy),int(fz),sSEG)];  W[3]=wcx * wcy * wfz;
  N[4]=SEG[sub2ind(int(fx),int(fy),int(cz),sSEG)];  W[4]=wfx * wfy * wcz;
  N[5]=SEG[sub2ind(int(cx),int(fy),int(cz),sSEG)];  W[5]=wcx * wfy * wcz;
  N[6]=SEG[sub2ind(int(fx),int(cy),int(cz),sSEG)];  W[6]=wfx * wcy * wcz; 
  N[7]=SEG[sub2ind(int(cx),int(cy),int(cz),sSEG)];  W[7]=wcx * wcy * wcz;
    
  float seg=0;
  for (int i=0; i<8; i++) seg = seg + N[i]*W[i];
  return seg;  
}


struct opt_type {
	int   CSFD;													// use CSFD
	int   PVE;													// 0, 1=fast, 2=exact
	float LB, HB, LLB, HLB, LHB, HHB;  	// boundary
	int   sL[3];
	// ...
	} opt;



// main function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<3) mexErrMsgTxt("ERROR: not enought input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR: to many input elements.\n");
  if (nlhs>2) mexErrMsgTxt("ERROR: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR: first  input must be an 3d single matrix\n");
 
  
  // main informations about input data (size, dimensions, ...)
  const mwSize *sL = mxGetDimensions(prhs[0]); int sSEG[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;
  const float   s2 = sqrt(2.0);
  const float   s3 = sqrt(3.0);
  const int     nr = nrhs;
  
  // indices of the neighbor Ni (index distance) and euclidean distance NW
  const int   NI[]  = {  0, -1,-x+1, -x,-x-1,  -xy+1,-xy,-xy-1,  -xy+x+1,-xy+x,-xy+x-1,  -xy-x+1,-xy-x,-xy-x-1};  
  const float ND[]  = {0.0,1.0,  s2,1.0,  s2,     s2,1.0,   s2,       s3,   s2,     s3,       s3,   s2,     s3};
  const int   sN  = sizeof(NI)/4;  
  float       DN[sN],DI[sN],GMTN[sN],WMDN[sN],SEGN[sN],DNm;
  
  float 			du, dv, dw, dnu, dnv, dnw, d, dcf, WMu, WMv, WMw, GMu, GMv, GMw, SEGl, SEGu, tmpfloat;
  int         ni,DNi,u,v,w,nu,nv,nw, tmpint, WMC=0, CSFC=0;
    
  // main volumes - actual without memory optimation ...
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[2] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[3] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL);
  plhs[4] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL);  
  
  // input variables
  float*SEG  = (float *)mxGetPr(prhs[0]);
  float*WMD  = (float *)mxGetPr(prhs[1]);
  float*CSFD = (float *)mxGetPr(prhs[2]);
  
  /*if ( nrhs>1) {
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"CSFD"));  printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=1) ) opt.CSFD = tmpint;   else opt.CSFD 	= 1;
		tmpint   = (int)mxGetScalar(mxGetField(prhs[1],1,"PVE"));   printf("X=%d", tmpint);	if ( tmpint!=NULL && (tmpint>=0 && tmpint<=2) ) opt.PVE  = tmpint; 	 else opt.PVE		= 2;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"LB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.LB   = tmpfloat; else opt.LB  	= 1.5;
		tmpfloat = (float)mxGetScalar(mxGetField(prhs[1],1,"HB"));  printf("X=%d", tmpfloat);	if ( tmpfloat!=NULL ) 													opt.HB   = tmpfloat; else opt.HB 		= 2.5;
	} 
	else */{ opt.CSFD = 1;opt.PVE = 2;opt.LB = 1.5;opt.HB	= 2.5; }
	opt.LLB=floor(opt.LB), opt.HLB=ceil(opt.LB), opt.LHB=floor(opt.HB), opt.HHB=ceil(opt.HB);
  
  // output variables
  float        *GMT  = (float *)mxGetPr(plhs[0]);
  float        *RPM  = (float *)mxGetPr(plhs[1]);
  
  // intitialisiation
  for (unsigned int i=0;i<nL;i++) {
  	GMT[i] = WMD[i];
    RPM[i] = WMD[i];
		// proof distance input
    if ( SEG[i]>=opt.HB ) WMC++;
    if ( SEG[i]<=opt.LB ) CSFC++;
  }
	if (WMC==0)  mexErrMsgTxt("ERROR: no WM voxel\n");
	if (CSFC==0) opt.CSFD = 0;

  
  
// thickness calcuation
// =============================================================================
  for (int i=0;i<nL;i++) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,u,v,w,xy,x);
      
      // read neighbor values
      for (int n=0;n<sN;n++) {
        ni = i + NI[n];
        ind2sub(ni,nu,nv,nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      // find minimum distance within the neighborhood
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN,DNm,DNi);
      GMT[i] = DNm;
    }
  }
  for (int i=nL-1;i>=0;i--) {
    if (SEG[i]>opt.LLB && SEG[i]<opt.HHB) {
      ind2sub(i,u,v,w,xy,x);
      
      // read neighbor values
      for (int n=0;n<sN;n++) {
        ni = i - NI[n];
        ind2sub(ni,nu,nv,nw,xy,x);
        if ( (ni<0) || (ni>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
        GMTN[n] = GMT[ni]; WMDN[n] = RPM[ni]; SEGN[n] = SEG[ni];
      }

      // find minimum distance within the neighborhood
      pmax(GMTN,WMDN,SEGN,ND,WMD[i],SEG[i],sN,DNm,DNi);
      if (GMT[i] < DNm) GMT[i] = DNm;
    }
  }
 
 
 	for (int i=0;i<nL;i++) if (SEG[i]<opt.LB || SEG[i]>opt.HB) GMT[i]=0; //WMD[i]

 
  
  

// final setings...
// =============================================================================
	float CSFDc = 0, GMTi, CSFDi; // 0.125
	for (int i=0;i<nL;i++) { 
		GMT[i] = GMT[i];
		if (SEG[i]>=opt.LB & SEG[i]<=opt.LB) {
			GMTi   = CSFD[i] + WMD[i];	
			CSFDi  = GMT[i] - WMD[i];
		
			if ( CSFD[i]>CSFDi )	CSFD[i] = CSFDi; 					
			else   		 			  		GMT[i]  = GMTi;
		}
	}
 
 
// median filter
// =============================================================================
	if (1) {
		int n1i,n2i; 
		for (int i=0;i<nL;i++) 
		{
			ind2sub(i,u,v,w,xy,x);
			n1i=i-1; ind2sub(n1i,nu,nv,nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2 ) n1i=i; 
			n2i=i+1; ind2sub(n2i,nu,nv,nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2 ) n2i=i;
			//if (GMT[n1i]==0), n1i=n2i; if (GMT[n2i]==0), n2i=n2i; 
			if 			( GMT[n1i]>GMT[i] && GMT[i]<GMT[n2i] && GMT[n1i]<GMT[n2i] ) RPM[i] = GMT[n1i];
			else if ( GMT[n1i]>GMT[i] && GMT[i]<GMT[n2i] && GMT[n1i]>GMT[n2i] ) RPM[i] = GMT[n2i];
			else																									 							RPM[i] = GMT[i];
		}
		
		for (int i=0;i<nL;i++) 
		{
			ind2sub(i,u,v,w,xy,x);
			n1i=i-x; ind2sub(n1i,nu,nv,nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2) n1i=i;
			n2i=i+x; ind2sub(n2i,nu,nv,nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2) n2i=i;
			if 			( RPM[n1i]>RPM[i] && RPM[i]<RPM[n2i] && RPM[n1i]<RPM[n2i] ) GMT[i] = RPM[n1i];
			else if ( RPM[n1i]>RPM[i] && RPM[i]<RPM[n2i] && RPM[n1i]>RPM[n2i] ) GMT[i] = RPM[n2i];
			else																																GMT[i] = RPM[i];
		}
		
		for (int i=0;i<nL;i++) 
		{
			ind2sub(i,u,v,w,xy,x);
			n1i=i-xy; ind2sub(n1i,nu,nv,nw,xy,x); if ( (n1i<0) || (n1i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2) n1i=i;
			n2i=i+xy; ind2sub(n2i,nu,nv,nw,xy,x); if ( (n2i<0) || (n2i>=nL) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || SEG[n1i]<2) n2i=i;
			if 			( GMT[n1i]>GMT[i] && GMT[i]<GMT[n2i] && GMT[n1i]<GMT[n2i] ) RPM[i] = GMT[n1i];
			else if ( GMT[n1i]>GMT[i] && GMT[i]<GMT[n2i] && GMT[n1i]>GMT[n2i] ) RPM[i] = GMT[n2i];
			else																																RPM[i] = GMT[i];
		}
		
		for (int i=0;i<nL;i++) {
			//if ( (GMT[i]-RPM[i])>-0.25 && (GMT[i]-RPM[i])<0.25 ) 
			if (SEG[i]>=2) GMT[i]=RPM[i]; 
		}
  }
 
 
// estimate RPM
// =============================================================================
	for (int i=0;i<nL;i++) {
		if ( SEG[i]>=opt.HB ) 	RPM[i]=1.0; 
		else {
			if ( SEG[i]<=opt.LB || GMT[i]==0.0 ) RPM[i]=0.0;
			else {
				RPM[i] = (GMT[i] - WMD[i]) / GMT[i];
				if (RPM[i]>1) RPM[i]=1;
				if (RPM[i]<0) RPM[i]=0;	
			}
		} 
	}
 
}


