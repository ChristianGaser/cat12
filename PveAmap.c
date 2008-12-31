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

void PveAmap(double *src, unsigned char *label, unsigned char *mask, unsigned char *prob, double *mean, double *separations, int *dims)
{

  int thresh, thresh_kmeans_int;
  double max_src;
    
  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  max_src = Kmeans( src, label, mask, 25, n_pure_classes, separations, dims, thresh, thresh_kmeans_int, iters_nu, pve);
    
  Amap( src, label, prob, mean, n_pure_classes, Niters, Nflips, subsample, dims, pve);

  if (pve) {
    printf("Calculate Partial Volume Estimate.\n");
    Pve5(src, prob, label, mean, dims);
  }

}