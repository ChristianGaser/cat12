/* function [mn,std,min,max,sum,nv,median?] = vbm_vol_ROIval(Ya,Yv)
 * _____________________________________________________________________
 *
 *
 * _____________________________________________________________________
 * Robert Dahnke 
 * Structural Brain Mapping Group
 * University Jena
 *
 * $Id$ 
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs!=3)                                      mexErrMsgTxt("ERROR:vbm_vol_ROIval: requires 2 maps\n");
  if (mxIsUint8(prhs[0])==0)                        mexErrMsgTxt("ERROR:vbm_vol_ROIval: 1st input must be an 3d uint8 matrix\n");
  if (mxIsSingle(prhs[1])==0)                       mexErrMsgTxt("ERROR:vbm_vol_ROIval: 2nd input must be an 3d single matrix\n");

  /* main informations about input data (size, dimensions, ...) */
  const int     nL = mxGetNumberOfElements(prhs[0]);

  /* input data */
  unsigned int *Ya = (uint8 *) mxGetPr(prhs[0]);
  float        *Yv = (float *) mxGetPr(prhs[1]); 
  const int sS[] = {1,3}; 
  mxArray *SS = mxCreateNumericArray(2,sS,mxDOUBLE_CLASS,mxREAL);
  double*S = mxGetPr(SS);
  if (nrhs<4) {S[0]=1; S[1]=1; S[2]=1;} else {S=mxGetPr(prhs[3]);}

  unsigned int i, n;

  /* output data */
  plhs[0]    = mxCreateNumericArray(256,sL,mxSINGLE_CLASS,mxREAL);
  float*mn = (float *)mxGetPr(plhs[0]);

  for (i==0,i<255,i++), rval[i]=0;

  /* sum, number, min, max */
  for (i==0,i<n,i++) {
    id=Ya[i];
    rn[id]++
    rsum[id]+=Yv[i];
    
    if (rmax[id]<Yv[i]), rmax[id]=Yv[i];
    if (rmin[id]>Yv[i]), rmin[id]=Yv[i];
    
  } 
  
  /* mean */
  for (i==0,i<255,i++) rmn[i]  = rsum[i] ./ rn[i];
  
  /* sandard deviation */
  
  
  /* median if required */
}