/*
 * Christian Gaser
 * $Id: upfirdn2dMex.c 226 2009-12-03 07:22:53Z gaser $ 
 *
 *
 * The core algorithm is based on upfirdn.cc in octave
 * Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France)
 *
 */

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"

void
upfirdn (double* x, double* y, double* h, int p, int q, const int* dims_x, int length_h)
{
  int rx, cx, Ly, c, m, n, lm, k, ix, ih;
  double r, accum;
  
  rx = dims_x[0];
  cx = dims_x[1];
  
  r = p/(double)q;
  Ly = ceil ((double) ((rx-1)*p + length_h) / (double) (q));  

  for (c = 0; c < cx; c++)
  {
    m = 0;
    while (m < Ly)
    {
      n = floor(m/r);
      lm = (m * q) % p;
      k = 0;
      accum = 0.0;
      do
      {
        ix = n - k;
        if (ix >= rx)
        {
          k ++;
          continue;
        }
              
        ih = k * p + lm;
        if ((ih >= length_h) | (ix < 0))
          break;

        accum += h[ih] * x[ix + c*dims_x[0]];
        k++;
      }
      while (1);
          
      y[m + c*Ly] = accum;
      m ++;
    }
  }
  
}
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations */
double *x, *y, *h;
int p, q, ndim_x;
int dims_y[2], length_h;
const int *dims_x, *dims_h;

/* check inputs */
if (nrhs!=4)
  mexErrMsgTxt("4 inputs required.");
else if (nlhs>1)
  mexErrMsgTxt("Too many output arguments.");
  
if (!mxIsDouble(prhs[0]))
	mexErrMsgTxt("First argument must be double.");

/* get input image */
x = (double*)mxGetPr(prhs[0]);

ndim_x = mxGetNumberOfDimensions(prhs[0]);
if (ndim_x!=2)
  mxErrMsgTxt("Image does not have 2 dimensions.");
  
dims_x = mxGetDimensions(prhs[0]);
if (dims_x[0] == 1)
  mxErrMsgTxt("Vector must be transposed.");

/* filter */
h = (double*)mxGetPr(prhs[1]);

dims_h = mxGetDimensions(prhs[1]);
if ((dims_h[0]>1) && (dims_h[1]>1))
  mxErrMsgTxt("Filter is not a vector.");

if (dims_h[0] == 1)
  length_h = dims_h[1];
else
  length_h = dims_h[0];

/* get parameters */
p = (int)(mxGetScalar(prhs[2]));
q = (int)(mxGetScalar(prhs[3]));

/* dimensions for output */
dims_y[0] = ceil ((double) ((dims_x[0]-1)*p + length_h) / (double) (q));
dims_y[1] = dims_x[1];

/*Allocate memory and assign output pointer*/
plhs[0] = mxCreateNumericArray(2,dims_y,mxDOUBLE_CLASS, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
y = mxGetPr(plhs[0]);

upfirdn (x, y, h, p, q, dims_x, length_h);

return;
}
