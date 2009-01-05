/*
 * Christian Gaser
 * $Id$ 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "optimizer3d.h"
#include "diffeo3d.h"
#include "Amap.h"

#define TINY 5e-2

struct dartel_prm {
  int rform;         /* regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy */
  double rparam[6];  /* regularization parameters */
  double lmreg;      /* LM regularization */
  int cycles;        /* number of cycles for full multi grid (FMG) */
  int its;           /* Number of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial */
};

/* First order hold resampling - trilinear interpolation */
/* modified version from spm_vol_utils.c from John Asburner */
void resample_trilinear(int m, double *vol, double *out, double *x, double *y, double *z,
  int dim[3], int offset)
{
	int i;
	double k111,k112,k121,k122,k211,k212,k221,k222;
	double dx1, dx2, dy1, dy2, dz1, dz2;
	int off1, off2, xcoord, ycoord, zcoord;

	for (i=0; i<m; i++)
	{
		if ((z[i]>=1) && (z[i]<dim[0]) &&
			  (y[i]>=1) && (y[i]<dim[1]) &&
			  (x[i]>=1) && (x[i]<dim[2]))
		{
		

			xcoord = (int)floor(x[i]); dx1=x[i]-xcoord; dx2=1.0-dx1;
			ycoord = (int)floor(y[i]); dy1=y[i]-ycoord; dy2=1.0-dy1;
			zcoord = (int)floor(z[i]); dz1=z[i]-zcoord; dz2=1.0-dz1;

      off1 = xcoord-1 + dim[0]*(xcoord-1 + dim[1]*(xcoord-1)) + offset;
      k222 = vol[off1]; k122 = vol[off1+1]; off2 = off1+dim[0];
      k212 = vol[off2]; k112 = vol[off2+1]; off1+= dim[0]*dim[1];
      k221 = vol[off1]; k121 = vol[off1+1]; off2 = off1+dim[0];
      k211 = vol[off2]; k111 = vol[off2+1];

			out[i] =  (((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
				+ (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1;

		}
		else out[i] = 0.0;

	}
}

void WarpPriors(unsigned char *prob, unsigned char *priors, unsigned char *mask, float *flow, int *dims, int loop)
{
  int vol_samp, vol, vol2, vol3, vol_samp2, vol_samp3, i, j, k, l, m;
  int size_samp[4], samp;
  double buf[3], ll[3], samp_1; 
  float *f, *g, *v, *flow1, *scratch, *mask_tmp, max;
  double *flow2x, *flow2y, *flow2z, *xs, *ys, *zs;
  int it, it0, it1, it_scratch, ndims4, dims_samp[3]; 
  int x, y, z, index_samp, xsamp, ysamp, ysamp2, zsamp, zsamp2;
  int yoffset, area_samp, zoffset, area, index;
  int z_area, y_dims, ix, iy, iz;
   
  int code = 2;    /* multinomial */
  int rform = 0;   /* linear energy */
  double lmreg = 0.01;
  static double param[3] = {1.0, 1.0, 1.0};

  samp = 4;
  
  /* only use gm/wm */
  ndims4 = 2;

  /* define grid dimensions */
  for(j=0; j<3; j++) dims_samp[j] = (int) ceil((dims[j]-1)/((double) samp))+1;
  
  /* find grid point conversion factor */
  samp_1 = 1.0/((double) samp);

  area = dims[0]*dims[1];
  vol  = dims[0]*dims[1]*dims[2];
  vol2 = 2*vol;
  vol3 = 3*vol;
  area_samp = dims_samp[0]*dims_samp[1];
  vol_samp  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol_samp2 = 2*vol_samp;
  vol_samp3 = 3*vol_samp;
  
  f      = (float *)malloc(sizeof(float)*vol_samp3);
  g      = (float *)malloc(sizeof(float)*vol_samp3);
  v      = (float *)malloc(sizeof(float)*vol_samp3);
  flow1  = (float *)malloc(sizeof(float)*vol_samp3);
  flow2x = (double *)malloc(sizeof(float)*vol);
  flow2y = (double *)malloc(sizeof(float)*vol);
  flow2z = (double *)malloc(sizeof(float)*vol);
  xs = (double *)malloc(sizeof(float)*dims[0]);
  ys = (double *)malloc(sizeof(float)*dims[1]);
  zs = (double *)malloc(sizeof(float)*dims[2]);
  mask_tmp = (float *)malloc(sizeof(float)*vol);

  for (i=0; i<dims[0]; i++) xs[i] = (double)i/(double)samp;
  for (i=0; i<dims[1]; i++) ys[i] = (double)i/(double)samp;
  for (i=0; i<dims[2]; i++) zs[i] = (double)i/(double)samp;

  /* for warping we need floating values of mask */
  for (i = 0; i < vol; i++) mask_tmp[i] = (float)mask[i];

  /* initialize flow field */
  for (i = 0; i < vol_samp3; i++) v[i] = flow[i];

  struct dartel_prm* prm = (struct dartel_prm*)malloc(sizeof(struct dartel_prm)*10);
  
  /* first three entries of param are equal */
  for (j = 0; j < loop; j++)
    for (i = 0; i < 3; i++)  prm[j].rparam[i] = param[i];

  /* some entries are equal */
  for (j = 0; j < loop; j++) {
    prm[j].rform = rform;
    prm[j].cycles = 3;
    prm[j].its = 3;
    prm[j].code = code;
    prm[j].lmreg = lmreg;
  }

  prm[0].rparam[3] = 4.0;   prm[0].rparam[4] = 2.0;    prm[0].rparam[5] = 1e-6; prm[0].k = 0; 
  prm[1].rparam[3] = 2.0;   prm[1].rparam[4] = 1.0;    prm[1].rparam[5] = 1e-6; prm[1].k = 0; 
  prm[2].rparam[3] = 1.0;   prm[2].rparam[4] = 0.5;    prm[2].rparam[5] = 1e-6; prm[2].k = 1; 
  prm[3].rparam[3] = 0.5;   prm[3].rparam[4] = 0.25;   prm[3].rparam[5] = 1e-6; prm[3].k = 2; 
  prm[4].rparam[3] = 0.25;  prm[4].rparam[4] = 0.125;  prm[4].rparam[5] = 1e-6; prm[4].k = 4; 
  prm[5].rparam[3] = 0.125; prm[5].rparam[4] = 0.0625; prm[5].rparam[5] = 1e-6; prm[5].k = 6; 

  for (i = 0; i < vol_samp3; i++) f[i] = 0.0;

  /* loop over neighborhoods of the grid points */
  for (k=-samp; k<=samp; k++) {
    for (l=-samp; l<=samp; l++) {
      for (m=-samp; m<=samp; m++) {
        for (z = 0; z < dims_samp[2]; z++) {
          zsamp = z*samp + k;
          if ((zsamp>=0) & (zsamp<dims[2])) {
            zsamp2 = zsamp*area;
            zoffset = z*area_samp;
            for (y=0; y<dims_samp[1]; y++) {
              ysamp = y*samp + l;
              if ((ysamp>=0) & (ysamp<dims[1])) {
                ysamp2 = ysamp*dims[0];
                yoffset = zoffset + y*dims_samp[0];
                for (x=0; x<dims_samp[0]; x++) {
                  xsamp = x*samp + m;
                  if ((xsamp>=0) & (xsamp<dims[0])) {
                    index_samp = yoffset + x;
                    index = zsamp2 + ysamp2 + xsamp;
                    f[index_samp          ] += ((float)priors[index       ]);
                    f[index_samp+vol_samp ] += ((float)priors[index + vol ]);
                    f[index_samp+vol_samp2] += ((float)priors[index + vol2]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }



  max = -1e15;
  for (i=0; i<vol_samp3; i++) max = MAX(f[i], max);
  for (i=0; i<vol_samp3; i++) f[i] /= max;

  for(i=0; i < 3; i++) size_samp[i] = dims_samp[i];
    
  size_samp[3] = ndims4;
  it = 0;
  for (it0 = 0; it0 < loop; it0++) {
    it_scratch = iteration_scratchsize((int *)size_samp, prm[it0].code, prm[it0].k);
    scratch = (float *)malloc(sizeof(float)*it_scratch);

    for (i = 0; i < vol_samp3; i++) g[i] = 0.0;

    /* load probabilities and use updated mask */
    /* loop over neighborhoods of the grid points */
    for (k=-samp; k<=samp; k++) {
      for (l=-samp; l<=samp; l++) {
        for (m=-samp; m<=samp; m++) {
          for (z = 0; z < dims_samp[2]; z++) {
            zsamp = z*samp + k;
            if ((zsamp>=0) & (zsamp<dims[2])) {
              zsamp2 = zsamp*area;
              zoffset = z*area_samp;
              for (y=0; y<dims_samp[1]; y++) {
                ysamp = y*samp + l;
                if ((ysamp>=0) & (ysamp<dims[1])) {
                  ysamp2 = ysamp*dims[0];
                  yoffset = zoffset + y*dims_samp[0];
                  for (x=0; x<dims_samp[0]; x++) {
                    xsamp = x*samp + m;
                    if ((xsamp>=0) & (xsamp<dims[0])) {
                      index_samp = yoffset + x;
                      index = zsamp2 + ysamp2 + xsamp;
                      g[index_samp          ] += ((float)prob[index + vol ])*((float)mask[index]);
                      g[index_samp+vol_samp ] += ((float)prob[index + vol2])*((float)mask[index]);
                      g[index_samp+vol_samp2] += ((float)prob[index       ])*((float)mask[index]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    max = -1e15;
    for (i=0; i<vol_samp3; i++) max = MAX(g[i], max);
    for (i=0; i<vol_samp3; i++) g[i] /= max;

    for (it1 = 0; it1 < prm[it0].its; it1++) {
      it++;
      iteration(size_samp, prm[it0].k, v, f, g, (float *)0, prm[it0].rform, prm[it0].rparam, prm[it0].lmreg, 
        prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);              
      printf("%02d:\t%.2f\t%6.2f\t%6.2f\t%6.2f\n", it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
      fflush(stdout);
      for (i = 0; i < vol_samp3; i++) v[i] = flow[i];
    }
    free(scratch);

    expdef(size_samp, 6, -1, v, flow, flow1, (float *)0, (float *)0); 

    resample_trilinear(vol, (double *)flow, flow2x, xs, ys, zs, dims, 0);    
    resample_trilinear(vol, (double *)flow, flow2y, xs, ys, zs, dims, vol_samp);
    resample_trilinear(vol, (double *)flow, flow2z, xs, ys, zs, dims, vol_samp2);

    /* warp mask */
    for (i = 0; i < vol; i++) {
      sampn(dims, mask_tmp, 1, vol, flow2x[i]*(double)samp-1.0, flow2y[i]*(double)samp-1.0, flow2z[i]*(double)samp-1.0, buf);
      mask[i] = (unsigned char) ROUND(buf[0]);
    }

  }

  /* rescue flow field */
  for (i = 0; i < vol_samp3; i++) flow[i] = v[i];

  free(flow1);
  free(flow2x);
  free(flow2y);
  free(flow2z);
  free(mask_tmp);
  free(xs);
  free(ys);
  free(zs);
  free(v);
  free(f);
  free(g);

}