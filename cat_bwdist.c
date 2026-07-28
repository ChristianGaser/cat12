/*
 * cat_bwdist.c
 *
 * Euclidean distance transform of a binary object, distance only.
 * Toolbox free replacement for bwdist of the Image Processing Toolbox.
 *
 *   D = cat_bwdist(P [,vx_vol])
 *
 *   P      (single) object map, the object is given by all voxels P>0.5
 *   vx_vol (double) 1x3 voxel size in mm, default [1 1 1]
 *   D      (single) Euclidean distance in mm
 *
 * Separable exact algorithm of Felzenszwalb & Huttenlocher (2012), which
 * computes the squared distance with three one dimensional passes, one per
 * dimension, each in linear time.  In contrast to the vector propagation
 * used by cat_vbdist it does not have to keep track of the nearest object
 * voxel, which is what makes it about ten times faster.  It therefore
 * cannot return the index or label maps, and is meant for the many callers
 * that only threshold the distance (e.g. the morphological operations).
 *
 * Differences to bwdist of the Image Processing Toolbox:
 *   - the object is P>0.5 and not every nonzero voxel, which follows the
 *     convention of cat_vbdist and makes this a drop in replacement for it
 *   - only the distance is returned, there is no second output IDX
 *   - only the euclidean metric, no cityblock/chessboard/quasi-euclidean
 *   - the voxel size can be given, so the distance is in mm and not in
 *     voxels, and anisotropic voxels are handled correctly
 *
 * The one dimensional transform computes
 *     d[q] = min_p ( f[p] + h^2 (q-p)^2 )
 * as the lower envelope of the parabolas rooted at every p, where h is the
 * voxel size along that dimension.
 *
 * ______________________________________________________________________
 * Christian Gaser, Robert Dahnke
 * Structural Brain Mapping Group (https://neuro-jena.github.io)
 * Departments of Neurology and Psychiatry
 * Jena University Hospital
 * ______________________________________________________________________
 * $Id$
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#ifndef INFINITY
#define INFINITY (1.0f/0.0f)
#endif

/* Large finite sentinel for "no object seen yet".  Deliberately not
   INFINITY: the parabola intersection below subtracts two such values, and
   INF-INF is NaN, while a finite object compared against an infinite one
   gives an intersection of -INF, which satisfies s <= z[0] = -INF and makes
   the envelope index k walk off the front of the array. */
#define FAR 1.0e18f

/* one dimensional squared distance transform of a single line */
static void dt1d(float *f, float *d, int *v, float *z, int n, float h2)
{
  int q, k = 0;
  float s;

  v[0] = 0;
  z[0] = -INFINITY;
  z[1] =  INFINITY;

  for (q = 1; q < n; q++) {
    /* intersection of the parabolas rooted at q and at v[k] */
    s = ((f[q] + h2*(float)q*(float)q) - (f[v[k]] + h2*(float)v[k]*(float)v[k]))
        / (2.0f*h2*((float)q - (float)v[k]));
    while (k > 0 && s <= z[k]) {   /* k>0 guard: never index below v[0]/z[0] */
      k--;
      s = ((f[q] + h2*(float)q*(float)q) - (f[v[k]] + h2*(float)v[k]*(float)v[k]))
          / (2.0f*h2*((float)q - (float)v[k]));
    }
    k++;
    v[k]   = q;
    z[k]   = s;
    z[k+1] = INFINITY;
  }

  k = 0;
  for (q = 0; q < n; q++) {
    while (z[k+1] < (float)q) k++;
    d[q] = h2*((float)q - (float)v[k])*((float)q - (float)v[k]) + f[v[k]];
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mwSize *sz;
  mwSize nel;
  int nx, ny, nz, x, y, z_, i, nmax;
  double *vx;
  float h2x, h2y, h2z;
  float *P, *D, *f, *d, *zz;
  int *v;

  if (nrhs < 1 || nrhs > 2) mexErrMsgTxt("ERROR:cat_bwdist: one or two inputs required.");
  if (nlhs > 1)             mexErrMsgTxt("ERROR:cat_bwdist: only one output, use cat_vbdist if the index map is needed.");
  if (!mxIsSingle(prhs[0])) mexErrMsgTxt("ERROR:cat_bwdist: first input must be an 3d single matrix.");
  if (mxGetNumberOfDimensions(prhs[0]) != 3)
                            mexErrMsgTxt("ERROR:cat_bwdist: first input must be an 3d single matrix.");

  sz = mxGetDimensions(prhs[0]);
  nx = (int)sz[0]; ny = (int)sz[1]; nz = (int)sz[2];
  nel = mxGetNumberOfElements(prhs[0]);

  h2x = h2y = h2z = 1.0f;
  if (nrhs == 2) {
    if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 3)
      mexErrMsgTxt("ERROR:cat_bwdist: vx_vol must be a double matrix with 3 elements.");
    vx = (double *)mxGetData(prhs[1]);
    h2x = (float)(vx[0]*vx[0]);
    h2y = (float)(vx[1]*vx[1]);
    h2z = (float)(vx[2]*vx[2]);
  }

  /* mxGetData and not mxGetPr: with the interleaved complex API mxGetPr is
     only valid for double arrays and returns garbage for single ones */
  plhs[0] = mxCreateNumericArray(3, sz, mxSINGLE_CLASS, mxREAL);
  D = (float *)mxGetData(plhs[0]);
  P = (float *)mxGetData(prhs[0]);

  /* squared distance, initialised to 0 on the object and to FAR elsewhere */
  for (i = 0; i < (int)nel; i++) D[i] = (P[i] > 0.5f) ? 0.0f : FAR;

  nmax = nx; if (ny > nmax) nmax = ny; if (nz > nmax) nmax = nz;
  f  = (float *)mxMalloc(sizeof(float) * nmax);
  d  = (float *)mxMalloc(sizeof(float) * nmax);
  zz = (float *)mxMalloc(sizeof(float) * (nmax + 1));
  v  = (int   *)mxMalloc(sizeof(int)   * nmax);

  /* dimension 1 (contiguous) */
  for (z_ = 0; z_ < nz; z_++)
    for (y = 0; y < ny; y++) {
      float *col = D + (mwSize)z_*nx*ny + (mwSize)y*nx;
      for (x = 0; x < nx; x++) f[x] = col[x];
      dt1d(f, d, v, zz, nx, h2x);
      for (x = 0; x < nx; x++) col[x] = d[x];
    }

  /* dimension 2 */
  for (z_ = 0; z_ < nz; z_++)
    for (x = 0; x < nx; x++) {
      float *base = D + (mwSize)z_*nx*ny + x;
      for (y = 0; y < ny; y++) f[y] = base[(mwSize)y*nx];
      dt1d(f, d, v, zz, ny, h2y);
      for (y = 0; y < ny; y++) base[(mwSize)y*nx] = d[y];
    }

  /* dimension 3 */
  for (y = 0; y < ny; y++)
    for (x = 0; x < nx; x++) {
      float *base = D + (mwSize)y*nx + x;
      for (z_ = 0; z_ < nz; z_++) f[z_] = base[(mwSize)z_*nx*ny];
      dt1d(f, d, v, zz, nz, h2z);
      for (z_ = 0; z_ < nz; z_++) base[(mwSize)z_*nx*ny] = d[z_];
    }

  mxFree(f); mxFree(d); mxFree(zz); mxFree(v);

  /* squared distance -> distance */
  for (i = 0; i < (int)nel; i++) D[i] = sqrtf(D[i]);
}
