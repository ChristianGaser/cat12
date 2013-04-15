/* Eikonal based nearest neighbor relation for a euklidean distance calcultion
 * _____________________________________________________________________________
 * [D,I] = vbm_vol_eidist(L,B,range,side)
 * 
 * L      (single)  input image
 * B      (logical) mask
 * range  
 * side
 *
 * Examples:
 *
 *
 *
 * _____________________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id$ 
 */
 
 #include "mex.h"   
#include "matrix.h"
#include "math.h"

// HELPFUNCTIONS:
// estimate x,y,z position of index i in an array size sx,sxy=sx*sy...
void ind2sub(int i,int &x,int &y, int &z, int sxy, int sy) {
	if (i<0) {
	  x=0;y=0;z=0; 
	  }
	else {
 	 	z = floor( i / (sxy) ) +1; 
  	i = i % (sxy);
  	y = floor( i / sy ) +1;        
  	x = i % sy + 1;
  }
}
int sub2ind(int x,int y, int z, const int s[]) {
	if (x<0) x=0; if (x>s[0]) x=s[0];
	if (y<0) y=0; if (y>s[1]) y=s[1];
	if (z<0) z=0; if (z>s[2]) z=s[2];
	
  int i=(z-1)*s[0]*s[1] + (y-1)*s[0] + (x-1);
  if (i<0 || i>s[0]*s[1]*s[2]) i=1;
  return i;
}

float max(float a, float b) { if (a<b) return b; else return a; }
float min(float a, float b) { if (a<b) return a; else return b; }
float abs2(float n) {	if (n<0) return -n; else return n; }
float sign(float n) {	if (n<0) return 1; else return 0; }



// read out linear interpolatet value of the volume
float isoval(float SEG[], float x, float y, float z, int s[]){
	if (x<0) x=0; if (x>s[0]) x=s[0];
	if (y<0) y=0; if (y>s[1]) y=s[1];
	if (z<0) z=0; if (z>s[2]) z=s[2];

  int i;
  float fx = floor(x),   fy = floor(y),		fz = floor(z);
  float cx = floor(x+1), cy = floor(y+1), cz = floor(z+1);
    
  float wfx = cx-x, wfy = cy-y, wfz = cz-z;
  float wcx = x-fx, wcy = y-fy, wcz = z-fz;

  // value of the 8 neighbors and there distance weight
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



// MAINFUNCTION [D,I] = ed(L,bd,r,s) = eidist(L,boundary,range,side)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs<1) mexErrMsgTxt("ERROR:vbm_vol_eidist:  not enought input elements\n");
  if (nrhs>4) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many input elements.\n");
  if (nlhs>3) mexErrMsgTxt("ERROR:vbm_vol_eidist: to many output elements.\n");
  if (mxIsSingle(prhs[0])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: first  input must be an 3d single matrix\n");
  if (mxIsSingle(prhs[1])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: second input must be an 3d single matrix\n");
  if (nrhs==3 && mxIsDouble(prhs[2])==0) mexErrMsgTxt("ERROR:vbm_vol_eidist: third  input must be an double matrix\n");
  if (nrhs==3 && mxGetNumberOfElements(prhs[2])!=3) {printf("ERROR:vbm_vol_eidist: third input must have 3 Elements"); throw 1; };
  if (nrhs==4 && mxIsDouble(prhs[3])==0 &&  mxGetNumberOfElements(prhs[3])!=1) {printf("ERROR:vbm_vol_eidist: fourth input must be one double value"); throw 1; }; 
  
  // input variables
  float*B	 = (float *)mxGetPr(prhs[0]);	// label map
  float*L  = (float *)mxGetPr(prhs[1]);	// label map
  
  //if (nrhs<2) float*bd = (float *)mxGetPr(prhs[1]); else float*bd  = 0; // tissue map
  //if (nrhs<3) float*r  = (float *)mxGetPr(prhs[2]); else float*r[] = {-INFINITY,0};
  //if (nrhs<4) float*s  = (float *)mxGetPr(prhs[3]); else float*s   = 1; 
  float bd  = 0.5;
  float r[] = {-INFINITY,0.5};
  float s   = 1;
  
  
  // main informations about input data (size, dimensions, ...)
  const mwSize *sL = mxGetDimensions(prhs[0]); int sSEG[] = {sL[0],sL[1],sL[2]}; 
  const int     dL = mxGetNumberOfDimensions(prhs[0]);
  const int     nL = mxGetNumberOfElements(prhs[0]);
  const int     x  = sL[0];
  const int     y  = sL[1];
  const int     xy = x*y;

  const int sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<3) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[2]);}
  bool euklid=true; double*deuklid; if (nrhs==4) {deuklid=mxGetPr(prhs[3]); euklid = (bool) deuklid[0]>0.5;};    
  
  float s1 = abs2(float(S[0])),s2 = abs2(float(S[1])),s3 = abs2(float(S[2]));
  const float   s12  = sqrt( s1*s1  + s2*s2); // xy - voxel size
  const float   s13  = sqrt( s1*s1  + s3*s3); // xz - voxel size
  const float   s23  = sqrt( s2*s2  + s3*s3); // yz - voxel size
  const float   s123 = sqrt(s12*s12 + s3*s3); // xyz - voxel size
  const int     nr = nrhs;
  // printf("%1.2f,%1.2f,%1.2f - %1.2f,%1.2f,%1.2f - %1.2f",s1,s2,s3,s12,s23,s13,s123);
  
    
  // indices of the neighbor Ni (index distance) and euclidean distance NW
  const int   NI[]  = {  1, -1,  x, -x, xy,-xy, -x-1,-x+1,x-1,x+1, -xy-1,-xy+1,xy-1,xy+1, -xy-x,-xy+x,xy-x,xy+x,  -xy-x-1,-xy-x+1,-xy+x-1,-xy+x+1, xy-x-1,xy-x+1,xy+x-1,xy+x+1};  
  const float ND[]  = { s1, s1, s2, s2, s3, s3,  s12, s12,s12,s12,   s13,  s13, s13, s13,   s23,  s23, s23, s23,     s123,   s123,   s123,   s123,   s123,  s123,  s123,  s123};
  //const float ND[]  = {1.0,1.0,1.0,1.0,1.0,1.0,   s2,  s2, s2, s2,    s2,   s2,  s2,  s2,    s2,   s2,  s2,  s2,       s3,     s3,     s3,     s3,     s3,    s3,    s3,    s3};
 
  float 	dinu, dinv, dinw, dcf, WMu, WMv, WMw, WM, DN, DIN;	// distance and intensity variables 
  int     i,j,n,ii,ni,DNi,u,v,w,nu,nv,nw,iu,iv,iw; 						// xyz and index-values of a voxel, one neighbor and the nearest boundary voxel
  int	  	nCV=0, nC=1, nCVL=0, kll=0, nCVo, kllo, kllv = (sL[0]+sL[1]+sL[2]);	// runtime variables   
 
    
  
  
  // main volumes - actual without memory optimation ...
  plhs[0] = mxCreateNumericArray(dL,sL,mxSINGLE_CLASS,mxREAL); // label map (single/float) 
  plhs[1] = mxCreateNumericArray(dL,sL,mxUINT32_CLASS,mxREAL); // index map (uint32)
  
  // output variables
  float*D  				= (float *)mxGetPr(plhs[0]);				// distance map
  unsigned int*I  = (unsigned int *)mxGetPr(plhs[1]);	// index map
  
	
	// check if there is a object and return INFINITY and i=i+1 (for matlab) for all points if there is no object
	int vx=0;	for (i=0;i<nL;i++) { if ( B[i]>=bd ) vx++; }
	if (vx==0) { for (i=0;i<nL;i++) { D[i]=INFINITY; I[i]=i+1; } return;}
		

	// initialisation of parameter volumes
	for (i=0;i<nL;i++) { 
	// only temporary to avoid unregular voxel
		ind2sub(i,u,v,w,xy,x); 
		if ( (i<0) || (i>=nL) || u<1 || v<1 || w<1 || u>(sL[0]) || v>(sL[1]) || w>(sL[2])  ) D[i] = -INFINITY; 
//		if ( (i<0) || (i>=nL) || u<1 || v<1 || w<1 || u>=(sL[0]-1) || v>=(sL[1]-1) || w>=(sL[2]-1)  ) D[i] = -INFINITY; 

	// regular part		
		I[i]=(unsigned int)i; 
		if ( L[i]<=0 || B[i]==-INFINITY) {
			D[i] = -INFINITY;
		} // neg object
		if ( D[i]!=-INFINITY) {
		 	if ( B[i]>bd )	{D[i]	 = 0; nCV++;}
			else       		 	{D[i]	 = INFINITY;}	
		}
	}

	// diffusion 
	while ( nCV>0 && nC>0 && kll<kllv) {
		kll++; nC=0;
		for (i=0;i<nL;i++) { 
			// only new values (-n), that are not -INF and that are not to far from the actual iteration (kll)
		try {
	 		if ( D[i]<=0 && D[i]!=-INFINITY && abs2(D[i])<=(s3*float(kll)+0.5) && B[i]>=0 && B[i]<=1 && L[i]>0 ) { // && L[i]<=1 
				if (D[i]<0) D[i]=-D[i]; 
				nCV--;  // demark points - also the with zero distance
				ii=(int)I[i];
				ind2sub(i,u,v,w,xy,x); 
				ind2sub(ii,iu,iv,iw,xy,x);
				for (n=0;n<26;n++) {
					ni=i+NI[n]; ind2sub(ni,nu,nv,nw,xy,x);
					if ( ( ( (ni<0) || (ni>=nL) || (abs2(nu-u)>1) || (abs2(nv-v)>1) || (abs2(nw-w)>1) || (ni==ii) )==0 ) && ( L[ni]<=(L[i]+0.00) ) ) { 
            if (euklid) {
              // standard euclidean distance ... but we want the real boundary
              dinu = float(nu)-float(iu); dinu *= s1;
              dinv = float(nv)-float(iv); dinv *= s2;
              dinw = float(nw)-float(iw); dinw *= s3;
              DIN  = sqrt(pow(dinu,2) + pow(dinv,2) + pow(dinw,2)); 

              if ( DIN>0 ) {
                // difference distance with normaliced vectors
                dinu /= DIN; dinv /= DIN; dinw /= DIN; // normal vector in normal space
                dinu /= s1;  dinv /= s2;  dinw /= s3;  // normal vector for anisotropic space
                WM = isoval(B,float(iu) + dinu,float(iv) + dinv,float(iw) + dinw,sSEG);
                if (B[ii]!=WM && B[ni]!=WM) {
                  // estimate new point bevor border to get the gradient based on this the exact HB
                  dcf = (B[ii] - bd) / ( B[ii] - WM );
                  WMu = float(iu) + dinu*dcf; WMv = float(iv) + dinv*dcf; WMw = float(iw) + dinw*dcf; 
                  WM = isoval(B,WMu,WMv,WMw,sSEG);
                  // new exact distance to interpolatet boundary
                  if ( WM>0.40 && WM<0.60 ) {
                    dinu = float(nu)-WMu; dinv = float(nv)-WMv; dinw = float(nw)-WMw;
                    DIN  = sqrt(pow(dinu*s1,2) + pow(dinv*s2,2) + pow(dinw*s3,2));   
                  }
                }
                else {
                  DIN=DIN-0.5;
                }
              }
            }
            else {
              DIN=abs2(D[ni]+ND[n]);
            }
						// use DN?	
						if ( abs2(D[ni])>DIN && D[ni]!=-INFINITY ) {
							if (D[ni]>0) nCV++; nC++;
							D[ni] = -DIN; 
							I[ni] = I[i];
						}
					}
				}
				if (D[i]==0) D[i]=-INFINITY; // demark start points
			}
		}
		catch (...) {
		  printf(":");
		}
		
		}
	}

	// correction of non-visited points due to miss-estimations of the exact boundary... use values of our neighbors 
	kll=0; nCV=0;
	for (i=0;i<nL;i++) {
		if (D[i]==INFINITY) nCV++; 
		if (D[i]<0 && D[i]!=-INFINITY) D[i]=-D[i]; 
	}
	nC=nCV; 
	while ( nCV>0 && nC>0 && kll<kllv) {
		kll++; nC=0;
		for (i=0;i<nL;i++) { 
			if ( (D[i]==INFINITY ) || (D[i]<=0 && D[i]!=-INFINITY) && abs2(D[i])<=(s3*float(kll)+0.5) && B[i]>=0 && B[i]<=1 && L[i]>0 && L[i]<=1) { 
				if (D[i]<0) D[i]=-D[i]; 
				nCV--;  // demark points - also the with zero distance
				ind2sub(i,u,v,w,xy,x); 
				for (n=0;n<26;n++) {
					ni=i+NI[n]; ind2sub(ni,nu,nv,nw,xy,x);
					if ( ( (ni<0) || (ni>=nL) || (abs2(nu-u)>1) || (abs2(nv-v)>1) || (abs2(nw-w)>1) && D[ni]!=INFINITY )==0 ) { 
						DN = abs2(D[ni]) + ND[n]; 
						if ( abs2(D[i])>DN ) { // use DN?
							if (D[i]>0) nCV++; nC++;
							D[i] = -DN; 
							I[i] = I[ni];
						}
					}
				}
			}
		}
	}

	// last correction
	for (i=0;i<nL;i++) { 
		if (D[i]<0 && D[i]!=-INFINITY) D[i]=-D[i]; 
		if (I[i]>0) I[i]++; else I[i]=1;						// correct for matlab index
		if (D[i]==-INFINITY || D[i]==NAN) D[i]=0; 	// correction of non-visited or other incorrect voxels
	} 
}