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

struct dartel_prm {
  int rform;         /* regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy */
  double rparam[6];  /* regularization parameters */
  double lmreg;      /* LM regularization */
  int cycles;        /* number of cycles for full multi grid (FMG) */
  int its;           /* Number of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial */
};


void WarpPriors(unsigned char *prob, unsigned char *priors, unsigned char *mask, float *flow, int *dims, int loop)
{
  int vol, vol0, vol2, vol3, i, j, k, l, m, is, size[4], samp;
  double buf[3], ll[3], sum_buf, samp_1; 
  float *f, *g, *v, *flow1, *scratch, max;
  int it, it0, it1, it_scratch, ndims4, dims_samp[3]; 
  int x, y, z, index_samp, xsamp, ysamp, ysamp2, zsamp, zsamp2;
  int yoffset, area_samp, zoffset, area, index;
  int z_area, y_dims, ix, iy, iz;
  unsigned char *mask_samp;
   
  int code = 2;    /* multinomial */
  int rform = 0;   /* linear energy */
  double lmreg = 0.01;
  static double param[3] = {1.0, 1.0, 1.0};

  samp = 3;
  
  /* only use gm/wm */
  ndims4 = 3;

  /* define grid dimensions */
  for(j=0; j<3; j++) dims_samp[j] = (int) ceil((dims[j]-1)/((double) samp))+1;
  
  /* find grid point conversion factor */
  samp_1 = 1.0/((double) samp);

  area = dims[0]*dims[1];
  vol0 = dims[0]*dims[1]*dims[2];
  area_samp = dims_samp[0]*dims_samp[1];
  vol  = dims_samp[0]*dims_samp[1]*dims_samp[2];
  vol2 = 2*vol;
  vol3 = 3*vol;
  
  f     = (float *)malloc(sizeof(float)*vol3);
  g     = (float *)malloc(sizeof(float)*vol3);
  v     = (float *)malloc(sizeof(float)*vol3);
  flow1 = (float *)malloc(sizeof(float)*vol3);
  mask_samp = (unsigned char*)malloc(sizeof(unsigned char)*vol);

  /* initialize flow field */
  for (i = 0; i < vol3; i++) v[i] = flow[i];

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

  for (i = 0; i < vol3; i++) f[i] = 0.0;

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
                    f[index_samp     ] += ((float)priors[index         ]);
                    f[index_samp+vol ] += ((float)priors[index + vol0  ]);
                    f[index_samp+vol2] += ((float)priors[index + vol0*2]);
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
  for (i=0; i<vol3; i++) max = MAX(f[i], max);
  for (i=0; i<vol3; i++) f[i] /= max;

  for(i=0; i < 3; i++) size[i] = dims_samp[i];
    
  size[3] = ndims4;
  it = 0;
  for (it0 = 0; it0 < loop; it0++) {
    it_scratch = iteration_scratchsize((int *)size, prm[it0].code, prm[it0].k);
    scratch = (float *)malloc(sizeof(float)*it_scratch);

    for (i = 0; i < vol3; i++) g[i] = 0.0;

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
                      g[index_samp     ] += ((float)prob[index + vol0  ])*((float)mask[index]);
                      g[index_samp+vol ] += ((float)prob[index + vol0*2])*((float)mask[index]);
                      g[index_samp+vol2] += ((float)prob[index         ])*((float)mask[index]);
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
    for (i=0; i<vol3; i++) max = MAX(g[i], max);
    for (i=0; i<vol3; i++) g[i] /= max;

    for (it1 = 0; it1 < prm[it0].its; it1++) {
      it++;
      iteration(size, prm[it0].k, v, f, g, (float *)0, prm[it0].rform, prm[it0].rparam, prm[it0].lmreg, 
        prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);              
      printf("%02d:\t%.2f\t%6.2f\t%6.2f\t%6.2f\n", it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
      fflush(stdout);
      for (i = 0; i < vol3; i++) v[i] = flow[i];
    }
    free(scratch);

    /* warp mask */
    expdef(size, 6, -1, v, flow, flow1, (float *)0, (float *)0); 
    for (i = 0; i < vol; i++) {
      sampn(size, f, 3, vol, (double)flow[i]-1.0, (double)flow[vol+i]-1.0, (double)flow[vol2+i]-1.0, buf);
      sum_buf = 0.0;
      for(j=0; j<3; j++) sum_buf += buf[j];
      mask_samp[i] = (unsigned char) MIN(ROUND(255.0*sum_buf), 255);
    }

    /* loop over image points */
    for (z = 0; z < dims[2]; z++) {
      z_area=z*area;
      for (y = 0; y < dims[1]; y++) {
        y_dims=y*dims[0];
        for (x = 0; x < dims[0]; x++)  {
	  
          i = x + y_dims + z_area;
          
          /* find the interpolation factors */
          ix = (int)(samp_1*x);
          iy = (int)(samp_1*y);
          iz = (int)(samp_1*z);
          index_samp = iz*area_samp + iy*dims_samp[0] + ix;
          
          mask[i] = mask_samp[index_samp];
        }
      }
    }
  }

  for (i = 0; i < vol3; i++) flow[i] = v[i];

  free(flow1);
  free(mask_samp);
  free(v);
  free(f);
  free(g);

}