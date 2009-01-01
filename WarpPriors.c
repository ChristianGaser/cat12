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

struct dartel_prm {
  int rform;         // regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy
  double rparam[6];  // regularization parameters
  double lmreg;      // LM regularization
  int cycles;        // number of cycles for full multi grid (FMG)
  int its;           // Number of relaxation iterations in each multigrid cycle
  int k;             // time steps for solving the PDE
  int code;          // objective function: 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial
};


void WarpPriors(unsigned char *prob, unsigned char *priors, unsigned char *mask, float *flow, int *dims, int loop)
{
  int vol, vol2, vol3, i, j;
  
  vol  = dims[0]*dims[1]*dims[2];
  vol2 = 2*vol;
  vol3 = 3*vol;

  float *f     = (float *)malloc(sizeof(float)*vol3);
  float *g     = (float *)malloc(sizeof(float)*vol3);
  float *v     = (float *)malloc(sizeof(float)*vol3);
  float *flow1 = (float *)malloc(sizeof(float)*vol3);

  for (i = 0; i < vol3; i++) v[i] = 0.0;

  int code = 2;    // multinomial
  int rform = 0;   // linear energy
  double lmreg = 0.01;
  static double param[3] = {1.0, 1.0, 1.0};

  struct dartel_prm* prm = (struct dartel_prm*)malloc(sizeof(struct dartel_prm)*10);
  // first three entrys of param are equal
  for (j = 0; j < loop; j++)
    for (i = 0; i < 3; i++)  prm[j].rparam[i] = param[i];

  // some entry are equal
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

  // change order to gm/wm/csf
  for (i = 0; i < vol; i++) {
    f[i     ] = ((float)priors[i+vol ])/255.0;
    f[i+vol ] = ((float)priors[i+vol2])/255.0;
    f[i+vol2] = ((float)priors[i     ])/255.0;
    g[i     ] = ((float)  prob[i+vol ])/255.0;
    g[i+vol ] = ((float)  prob[i+vol2])/255.0;
    g[i+vol2] = ((float)  prob[i     ])/255.0;
  }

  int size[4];
  for(i=0; i < 3; i++) size[i] = dims[i];
    
  // only use gm/wm
  size[3] = 2;
  int it = 0, it0, it1; 
  double ll[3]; 
  for (it0 = 0; it0 < loop; it0++) {
    int it_scratch = iteration_scratchsize((int *)size, prm[it0].code, prm[it0].k);
    float *scratch = (float *)malloc(sizeof(float)*it_scratch);
    for (it1 = 0; it1 < prm[it0].its; it1++) {
      it++;
      iteration(size, prm[it0].k, v, f, g, (float *)0, prm[it0].rform, prm[it0].rparam, prm[it0].lmreg, 
        prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);              
      printf("%02d:\t%6g\t%6g\t%6g\t%6g\n", it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
      fflush(stdout);
      for (i = 0; i < vol3; i++) v[i] = flow[i];
    }
    free(scratch);
  }

  expdef(size, 6, -1, v, flow, flow1, (float *)0, (float *)0); 
  for (i = 0; i < vol; i++) {
    double buf[3], sum_buf = 0.0; 
    int j;
    sampn(size, f, 3, vol,
            (double)flow[i]-1.0, (double)flow[vol+i]-1.0, (double)flow[vol2+i]-1.0, buf);
    for(j=0; j<3; j++) sum_buf += buf[j];
    // set result to zero if sum of probabilities is too low
    if (sum_buf < 0.25) {
      prob[i] = 0;
      prob[i+vol] = 0;
      prob[i+vol2] = 0;
    }
  }

  free(flow1);
  free(v);
  free(f);
  free(g);

}