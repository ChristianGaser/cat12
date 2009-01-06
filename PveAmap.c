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

  int thresh, thresh_kmeans_int, vol, i;
  int n_loops, update_label, sum_priors;
  unsigned char *label;
  double max_src, max_mask;
  float *flow;
  
  vol = dims[0]*dims[1]*dims[2];
  label     = (unsigned char*)malloc(sizeof(unsigned char)*vol);
  flow  = (float *)malloc(sizeof(float)*vol*3);
  
  /* initialize flow field with zeros */
  for (i = 0; i < (vol*3); i++) flow[i] = 0.0;
  
  /* check maximum of mask to indicate whether it's defined or not */
  max_mask = -1e15;
  for (i=0; i<vol; i++) max_mask = MAX(mask[i], max_mask);

  /* compute mask based on sum of tissue priors for GM/WM/CSF if not given */
  if(max_mask == 0) {
    for (i=0; i<vol; i++) {
      sum_priors = (int)priors[i] + (int)priors[i+vol] + (int)priors[i+2*vol];
      if(sum_priors > 255) mask[i] = 255;
      else mask[i] = (unsigned char) sum_priors;
    }
  }
    
  Niters = 5;
  thresh_brainmask = 0.01;
  pve = 1;

  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  max_src = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, pve);
      
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

  n_loops = 3;
//  WarpPriors(prob, priors, mask, flow, dims, n_loops);
  
  for(i=0; i<vol; i++)
    if(mask[i] < 32) src[i] = 0.0;
  
  Amap( src, label, prob, mean, n_pure_classes, Niters, Nflips, subsample, dims, pve);

  if (pve) {
    printf("Calculate Partial Volume Estimate.\n");
    update_label = 1;
    Pve5(src, prob, label, mean, dims, update_label);
  }

  n_loops = 6;
  WarpPriors(prob, priors, mask, flow, dims, n_loops);

  for(i=0; i<vol; i++) {
    if(mask[i] < 32) {
      prob[i      ] = 0;
      prob[i+vol  ] = 0;
      prob[i+vol*2] = 0;
    }
  }
  
  free(label);
  free(flow);
}