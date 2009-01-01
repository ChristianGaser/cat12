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
#include "PveAmap.h"

void PveAmap(double *src, unsigned char *priors, unsigned char *mask, unsigned char *prob, double *mean, double *separations, int *dims)
{

  int thresh, thresh_kmeans_int, vol, i, j;
  int n_loops, update_label;
  unsigned char *label, *mask_init;
  double max_src;
  float *flow;
  
  vol = dims[0]*dims[1]*dims[2];
  label     = (unsigned char*)malloc(sizeof(unsigned char)*vol);
  mask_init = (unsigned char*)malloc(sizeof(unsigned char)*vol);
  flow  = (float *)malloc(sizeof(float)*vol*3);
  
  /* compute mask based on sum of tissue priors for GM/WM/CSF */
  for (j=0; j<3; j++)
    for (i=0; i<vol; i++) 
      mask_init[i] = priors[i] + priors[i+vol] + priors[i+2*vol];
    
  Niters = 3;
  thresh_brainmask = 0.01;

  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  max_src = Kmeans( src, label, mask_init, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, pve);
      
  if (pve) {
    update_label = 0;
    Pve5(src, prob, label, mean, dims, update_label);
  } else {
    for(i=0; i<vol; i++) {
      switch(label[i]) {
      case 0:
        prob[i] = 0; prob[i+vol] = 0; prob[i+2*vol] = 0;
        break;
      case 1:
        prob[i] = 255; prob[i+vol] = 0; prob[i+2*vol] = 0;
        break;
      case 2:
        prob[i] = 0; prob[i+vol] = 255; prob[i+2*vol] = 0;
        break;
      case 3:
        prob[i] = 0; prob[i+vol] = 0; prob[i+2*vol] = 255;
        break;
      }
    }
  }

  n_loops = 1;
//  WarpPriors(prob, priors, mask_init, flow, dims, n_loops);
  
  Amap( src, label, prob, mean, n_pure_classes, Niters, Nflips, subsample, dims, pve);

  if (pve) {
    printf("Calculate Partial Volume Estimate.\n");
    update_label = 1;
    Pve5(src, prob, label, mean, dims, update_label);
  }

  free(mask_init);
  free(label);
  free(flow);
}