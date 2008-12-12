/*
 * Tskmeans.c
 *
 * Tree structure k-means algorithm
 *
 * Jagath C. Rajapakse (raja@cns.mpg.de) 23-07-97
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#ifndef MAX
#define MAX( x, y )  ( ((x) >= (y)) ? (x) : (y) )
#endif

double EstimateKmeans(double *src, unsigned char *label, unsigned char *mask, int nc, double *mean, int ni, int BG, int *dims, int thresh_mask, int thresh_kmeans, double max_src)
/* perform k-means algorithm give initial mean estimates */    
{
  int i, j, j0, x, y, z, v;
  int count;
  long histo[256], lut[256], vol, area;
  long z_area, y_dims;
  double diff, dmin, dx, xnorm, sum;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  /* build intensity histogram */
  for (i = 0; i < 256; i++) histo[i] = 0;
  for (z=0;z<dims[2];z++) {
    z_area = z*area;
    for (y=0;y<dims[1];y++) {
      y_dims = y*dims[0];
      for (x=0;x<dims[0];x++) {
         v = (int)round(255.0*src[z_area + y_dims + x]/max_src);
         if (v < BG) continue;
         if ((thresh_mask > 0) && ((int)mask[z_area + y_dims + x] < thresh_kmeans))
           continue;
         if (v < 0) v = 0;
         if (v > 255) v = 255;	
         histo[v]++;
      }
    }
  }

  /* use only value in histogram where cumsum is between 1..99% */
  long cumsum[256];
  cumsum[0] = histo[0];
  for (i = 1; i < 256; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 256; i++) cumsum[i] = (long) round(1000.0*(double)cumsum[i]/(double)cumsum[255]);
  for (i = 0; i < 256; i++) if ((cumsum[i] <= 10) || (cumsum[i] >= 990)) histo[i] = 0;

  /* loop through */
  diff = 1e15;  count = 0;
  while (diff > 1.0 && count < ni) {

    /* assign class labels */
    for (i = BG; i < 256; i++) {
      dmin = 256.0 * 256.0;
      for (j = 0; j < nc; j++) {
	    dx = (double) i - mean[j];
	    dx *= dx;
	    if (dx < dmin) {
	      lut[i] = j;
	      dmin = dx;
	    }
      }
    }

    /* find the new cluster centers */
    diff = 0;
    for (i = 0; i < nc; i++) {
      xnorm = 0.0; sum = 0.0;
      for (j = BG; j < 256; j++)
	    if (lut[j] == i) {
	      xnorm += histo[j];
	      sum +=  j * histo[j];
	    }
      sum = xnorm > 0 ? sum /= xnorm: 0.0;
      dx = sum - mean[i];
      mean[i] = sum;
      dx *= dx;
      diff += dx;
    }
    count++;
  }

  /* assign final labels to voxels */
  for (i=BG; i<256; i++) {
    dmin = 1e15;
    j0 = 0;
    for (j = 0; j < nc; j++) {
      if (fabs((double) i - mean[j]) < dmin) {
	    dmin = fabs((double)i - mean[j]);
	    j0 = j;
      }
    }
    lut[i] = j0;
  }
  
  if (BG == 1) lut[0] = 0;

  /* adjust for the background label */
  diff = 0;
  
  for (z=0;z<dims[2];z++) {
    z_area = z*area;
    for (y=0;y<dims[1];y++) {
      y_dims = y*dims[0];
      for (x=0;x<dims[0];x++) {
         v = (int)round(255.0*src[z_area + y_dims + x]/max_src);
         if (v >= BG) {
           if (v < 0) v = 0;
           if (v > 255) v = 255;
           label[z_area + y_dims + x] = (unsigned char)(lut[v] + BG);	
           diff += SQR((double)v - mean[lut[v]]);
           if ((thresh_mask > 0) && ((int)mask[z_area + y_dims + x] < thresh_mask))
               label[z_area + y_dims + x] = 0;	
         }
         else label[z_area + y_dims + x] = 0;	
      }
    }
  }

  /* return square error */
  return(diff);
}

double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int nclusters, int BG, double *separations, int *dims, int labelto3, int thresh_mask, int thresh_kmeans, int iters_nu)
{
  int i, j, l, k, x, y, z;
  double e, emin, eps, *nu, *src_bak, th_src, val_nu;
  double last_err = 1e10;
  double max_src = -1e10;
  double mean_nu = 0.0;
  long n[nclusters];
  double mean[nclusters];
  double var[nclusters];
  double mu[nclusters];
  double Mu[nclusters];
  int val, nc;
  long vol, count, area, z_area, y_dims;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  src_bak = (double *)malloc(sizeof(double)*vol);
  
  if (iters_nu > 0) 
    nu = (double *)malloc(sizeof(double)*vol);
  
  /* find maximum and mean inside mask */
  count = 0;
  for (i = 0; i < vol; i++) {
    if (mask[i] > 0) {
      max_src = MAX(src[i], max_src);
    }
  }
  
  /* for reducing 5 labels to 3 restrict initial segmentation to 3 classes */
  int nc_initial = nclusters;
  if (labelto3) nc_initial = 3;
  
  /* go through all sizes of cluster beginning with two clusters */
  for (nc=2; nc<=nc_initial; nc++) {

    if (nc == 2) {
      /* initialize for the two cluster case; */
      n[0]=0; mean[0] = 0.0; var[0] = 0.0;

      for (z=0;z<dims[2];z++) {
        z_area = z*area;
        for (y=0;y<dims[1];y++) {
          y_dims = y*dims[0];
          for (x=0;x<dims[0];x++) {
            val = (int)round(255.0*src[z_area + y_dims + x]/max_src);
            if (val < BG) continue;
            n[0]++;
            mean[0] += (double) val;
            var[0] += (double) val*(double) val;
          }
        }
      } 

      Mu[0] = n[0] != 0 ? mean[0]/n[0]: 0.0;
      var[0] = n[0] > 1 ? (var[0] - n[0]*Mu[0]*Mu[0])/(n[0] - 1.0) : 1.0;
      eps = 0.5*sqrt(var[0]);
    }
    else {
      /* find the deviant (epsilon) for the node being divided */
      eps = Mu[0];
      for (i=0; i<nc-2; i++)
        if (Mu[i+1] - Mu[i] < eps)
          eps = Mu[i+1] - Mu[i];
      if (255 - Mu[nc-2] < eps)
        eps = 255 - Mu[nc-2];
      eps = eps*0.5;
    }

    /* go through low order clustering */
    emin = 1e15;
    for (k=0; k<nc-1; k++) {
      for (i=nc-1; i>k+1; i--) mean[i] = Mu[i-1];
      mean[k+1] = Mu[k] + eps;  mean[k] = Mu[k] - eps;
      for (i=BG; i<k; i++) mean[i] = Mu[i];
      e = EstimateKmeans(src, label, mask, nc, mean, NI, BG, dims, thresh_mask, thresh_kmeans, max_src);
      if (e < emin) {
        emin = e;
        for (i=0; i<nc; i++) 
          mu[i] = mean[i];
      }
    }
    for (i=0; i<nc; i++) Mu[i] = mu[i];     
  }

  /* only use values above the mean of the lower two cluster for nu-estimate */
  th_src = max_src*(double)((mu[0]+mu[1])/2.0)/255.0;

  /* extend initial 3 clusters to 5 clusters by averaging clusters */
  if (labelto3 == 1) {
    mu[4] = mu[2];
    mu[2] = mu[1];
    mu[3] = (mu[4]+mu[2])/2;
    mu[1] = (mu[0]+mu[2])/2;
  }
  
  /* extend initial 3 clusters to 6 clusters by averaging clusters */
  if (labelto3 == 2) {
    mu[5] = mu[2];
    mu[3] = mu[1];
    mu[1] = mu[0];
    mu[4] = (mu[3]+mu[5])/2;
    mu[2] = (mu[1]+mu[3])/2;
    mu[0] = mu[1]/2;
  }

  /* find the final clustering and correct for nu */
  if (iters_nu > 0) {
    int count_err = 0;
    for (j = 0; j <= iters_nu; j++) {  
      count = 0;
      mean_nu = 0.0;
      for (i = 0; i < vol; i++) {
        nu[i] = 0.0;
        /* only use values above threshold where mask is defined for nu-estimate */
        if ((src[i] > th_src) && (mask[i] > thresh_kmeans)) {
          val_nu = src[i]/mu[label[i]-BG];
          if ((isfinite(val_nu))) {
            nu[i] = val_nu;
            mean_nu += val_nu;
            count++;
          }
        }

      }
      
      /* correct nu input to a mean of 1 to remain original intensity range */
      mean_nu /= (double)count;
      for (i=0; i<vol; i++)
        nu[i] /= mean_nu;
      
      /* spline estimate: start with distance of 1500 end end up with 500 */
      splineSmooth(nu, 0.01, MAX(500,1500.0/(j+1)), 4, separations, dims);
      
      /* apply nu correction to source image */
      for (i=0; i<vol; i++) {
        if (nu[i] > 0)
          src[i] /= nu[i];
      }
      
      /* update k-means estimate */
      e = EstimateKmeans(src, label, mask, nclusters, mu, NI, BG, dims, thresh_mask, thresh_kmeans, max_src);

      if (e > last_err)  count_err++;
      else count_err = 0;

      /* interrupt if last error was for the last 2 iterations larger 
       * or change is < 0.25% */
      if ((count_err > 1) || (((last_err-e)/e < 0.0025) && ((last_err-e)/e > 0))) {
        /* rescue old values from previous iteration */
        for (i=0; i<vol; i++) src[i] = src_bak[i];
        break;
      }

      /* save old values */
      for (i=0; i<vol; i++) src_bak[i] = src[i];

      last_err = e;
    
      fprintf(stderr,"iters:%2d error: %6.2f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",j+1, e/(nc*dims[0]*dims[1]*dims[2]));
    
    }
  } else {
    e = EstimateKmeans(src, label, mask, nclusters, mu, NI, BG, dims, thresh_mask, thresh_kmeans, max_src);
  }
  
  max_src = -1e10;
  for (i = 0; i < vol; i++)
    max_src = MAX(src[i], max_src);

  fprintf(stderr,"\nK-Means: ");
  for (i=0; i<nclusters; i++) fprintf(stderr,"%3.3f ",max_src*mu[i]/255.0); 
  fprintf(stderr,"\terror: %3.3f\n",e/(nc*dims[0]*dims[1]*dims[2]));    

  free(src_bak);
  
  if (iters_nu > 0)
    free(nu);
      
  return(max_src);    
}













