/*
 * Christian Gaser
 * $Id$ 
 *
 */

/* This code is a substantially modified version of Amap.C 
 * from Jagath C. Rajapakse
 * 
 * Original author : Jagath C. Rajapakse
 *
 * See:
 * Statistical approach to single-channel MR brain scans
 * J. C. Rajapakse, J. N. Giedd, and J. L. Rapoport
 * IEEE Transactions on Medical Imaging, Vol 16, No 2, 1997
 *
 * Comments to raja@cns.mpg.de, 15.10.96
 */

/* The likelihood and PVE calculations are a modified version from
 * the PVE software bundle:
 * Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
 * Technology, 2002 - 2004.
 * P.O. Box 553, FIN-33101, Finland
 * E-mail: jussi.tohka@tut.fi
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Amap.h"

/* calculate the mean and variance for every class on a grid size SUBxSUBxSUB */
static void GetMeansVariances(double *src, unsigned char *label, int nc, struct point *r, int sub, int *dims, double *thresh)
{
  int i, j, ind;
  int area, narea, nvol, zsub, ysub, xsub, yoffset, zoffset;
  int zsub2, ysub2;
  int nix, niy, niz, k, l, m, z, y, x, label_value;
  double val;
  struct ipoint *ir;
  int label_value_BG;

  area = dims[0]*dims[1];

  /* define grid dimensions */
  nix = (int) ceil((dims[0]-1)/((double) sub))+1;
  niy = (int) ceil((dims[1]-1)/((double) sub))+1;
  niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
  narea = nix*niy;
  nvol  = nix*niy*niz;
  
  ir = (struct ipoint*)malloc(sizeof(struct ipoint)*nc*nvol);
  if(ir == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  for (i=0; i<nc; i++) {
    for (j=0; j<nvol; j++) {
      ind = (i*nvol)+j; 
      ir[ind].n = 0;
      ir[ind].s = 0.0;
      ir[ind].ss = 0.0;
    }
  }
  

  /* loop over neighborhoods of the grid points */
  for (k=-sub; k<=sub; k++) {
    for (l=-sub; l<=sub; l++) {
      for (m=-sub; m<=sub; m++) {
        for (z = 0; z < niz; z++) {
          zsub = z*sub + k;
          if (zsub>=0 & zsub<dims[2]) {
            zsub2 = zsub*area;
            zoffset = z*narea;
            for (y=0; y<niy; y++) {
              ysub = y*sub + l;
              if (ysub>=0 & ysub<dims[1]) {
                ysub2 = ysub*dims[0];
                yoffset = zoffset + y*nix;
                for (x=0; x<nix; x++) {
                  xsub = x*sub + m;
                  if (xsub>=0 & xsub<dims[0]) {
                    label_value = (int)label[zsub2 + ysub2 + xsub];
                    label_value_BG = label_value - 1;
                    if (label_value_BG < 0) continue;
                    val = src[zsub2 + ysub2 + xsub];
                    
                    /* exclude values out of quartile 1-99% */
                    if ((val<thresh[0]) || (val>thresh[1])) continue;
                    ind = ((label_value_BG)*nvol)+yoffset+x;
                    ir[ind].n++;
                    ir[ind].s += val; ir[ind].ss += val*val;
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* find means and standard deviations */
  for (i=0; i<nc; i++) {
    for (j=0; j<nvol; j++) {
      ind = (i*nvol)+j;
      if (ir[ind].n > G) {
        r[ind].mean = ir[ind].s/ir[ind].n;
        if (ir[ind].n == 1)
          r[ind].var = 0.0;
        else    r[ind].var = (ir[ind].ss -ir[ind].n*SQR(r[ind].mean))/(ir[ind].n-1);
      } else r[ind].mean = 0.0;
    }
  }

  free(ir);

  return;
}

/* Computes likelihood of value given parameters mean and variance */ 
double ComputeGaussianLikelihood(double value, double mean , double var)

{ 
  return(exp(-(SQR(value - mean))/(2 * var))/SQRT2PI/sqrt(var));
}

/* -------------------------------------------------------------------
Computes the likelihoods for the mixed classes. Returns the likelihood.
var1,var2 are the variances of pdfs representing pure classes. measurement_var 
is the measurement noise. So the model for the variable y (representing the 
intensity value) that is composed of t * tissue1 and (1 - t)* tissue2 becomes :

y = t*x1 + (1 - t)*x2 + xm,
x1 ~ N(mean1,var1) , x2 ~ N(mean2,var2) , xm ~ N(0,measurement_var).

Note: The numerical integration routine used by the 
function is primitive , but so is the mankind...

*/

double ComputeMarginalizedLikelihood(double value, double mean1 , double mean2, 
                                       double var1, double var2, 
                                       unsigned int nof_intervals)

{ 
  double lh, tmean, tvar, delta, step;
  
  step = 1.0 / (double) nof_intervals;
  lh = 0;
  
  for(delta = 0.0; delta <= 1.0; delta += step) {
    tmean = delta * mean1 + ( 1 - delta ) * mean2;
    tvar = SQR(delta) * var1 + SQR(1 - delta) * var2;
    lh += ComputeGaussianLikelihood(value, tmean, tvar)*step;
  }
  
  return(lh);
 }


/* Find maximum argument out of the n possibilities */

unsigned char MaxArg(double *pval, unsigned char n)
{
  double maximum;
  unsigned char i, index;
  
  maximum = pval[0];
  index = 1;
  for(i = 1; i < n; i++) {
    if(pval[i] > maximum) {
      index = i + 1;
      maximum = pval[i];
    }
  }
  return(index);
}

/* Compute initial PVE labeling based on marginalized likelihood */
void ComputeInitialPveLabel(double *src, unsigned char *label, struct point *r, int nc, int sub, int *dims)
{
  
  int x, y, z, z_area, y_dims, index, label_value;
  int i, ix, iy, iz, ind, ind2, nix, niy, niz, narea, nvol;
  long area, vol;
  double val, sub_1, mean[nc], var[nc], d_pve[nc];
  
  area = dims[0]*dims[1];
  vol = area*dims[2];

  /* find grid point conversion factor */
  sub_1 = 1.0/((double) sub);

  /* define grid dimensions */
  nix = (int) ceil((dims[0]-1)/((double) sub))+1;
  niy = (int) ceil((dims[1]-1)/((double) sub))+1;
  niz = (int) ceil((dims[2]-1)/((double) sub))+1; 

  narea = nix*niy;
  nvol = nix*niy*niz;
  
  /* loop over image points */
  for (z = 1; z < dims[2]-1; z++) {
    z_area=z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims=y*dims[0];
      for (x = 1; x < dims[0]-1; x++)  {
	  
        index = x + y_dims + z_area;
        label_value = (int)label[index];
        if (label_value < 1) continue;
        val = src[index];
          
        /* find the interpolation factors */
        ix = (int)(sub_1*x), iy = (int)(sub_1*y), iz = (int)(sub_1*z);
        ind = iz*narea + iy*nix + ix;
          
        for (i=0; i<nc; i++) {
          ind2 = (i*nvol) + ind;            
          if (r[ind2].mean > 0.0) {
            mean[1+i*2] = r[ind2].mean;
            var[1+i*2]  = r[ind2].var;
          }
        }

        if (fabs(mean[CSFLABEL]) > TINY) {
          d_pve[CSFLABEL] = ComputeGaussianLikelihood(val, mean[CSFLABEL], var[CSFLABEL]);
        } else d_pve[CSFLABEL] = HUGE;

        if (fabs(mean[GMLABEL]) > TINY) {
          d_pve[GMLABEL] = ComputeGaussianLikelihood(val, mean[GMLABEL], var[GMLABEL]);
        } else d_pve[GMLABEL] = HUGE;

        if (fabs(mean[WMLABEL]) > TINY) {
          d_pve[WMLABEL] = ComputeGaussianLikelihood(val, mean[WMLABEL], var[WMLABEL]);
        } else d_pve[WMLABEL] = HUGE;

        if (fabs(mean[CSFLABEL]) > TINY) {
          d_pve[BKGCSFLABEL] = ComputeMarginalizedLikelihood(val, 0.0, mean[CSFLABEL],
                                        0.1*MIN3(var[CSFLABEL],var[GMLABEL],var[WMLABEL]), var[CSFLABEL], 100 );
        } else d_pve[BKGCSFLABEL] = HUGE;

        if ((fabs(mean[WMLABEL]) > TINY) && (fabs(mean[GMLABEL]) > TINY)) {
          d_pve[WMGMLABEL] = ComputeMarginalizedLikelihood(val, mean[WMLABEL], mean[GMLABEL],
                                        var[WMLABEL], var[GMLABEL], 100 );
        } else d_pve[WMGMLABEL] = HUGE;
            
        if ((fabs(mean[CSFLABEL]) > TINY) && (fabs(mean[GMLABEL]) > TINY)) {
          d_pve[GMCSFLABEL] = ComputeMarginalizedLikelihood(val, mean[GMLABEL], mean[CSFLABEL],
                                        var[GMLABEL], var[CSFLABEL], 100 );
        } else d_pve[GMCSFLABEL] = HUGE;

        label[index] = (unsigned char) MaxArg(d_pve, MAX_NC);
      }
    }
  }   
} 

void Normalize(double* pval, char n)
{
  double sca = 0.0;
  int i;

  for(i = 0;i < n;i++) 
    sca += pval[i];
 
  if(fabs(sca) >  TINY) {       /* To avoid divisions by zero */
    for(i = 0;i < n;i++) {
      pval[i] /= sca;
    }
  }
}

void ComputeMrfProbability(double *mrf_probability, double *exponent, unsigned char *label, int x, int y , int z, int *dims,
                             int nc, double beta)
{
  int i,j,k;
  unsigned char label1, label2;  
  double distance;
  double slice_width_sq[3] = {1.0, 1.0, 1.0};
  int similarity_value;
  int same = -5;
  int similar = -1;
  int different = 5; 
  
  /* To determine if it's possible to get out of image limits. 
     If not true (as it usually is) this saves the trouble calculating this 27 times */
  for(label1 = 0;label1 < nc;label1++)
    exponent[label1] = 0;
  
  for(i = -1; i < 2; i++) {
    for(j = -1; j < 2; j++) {
      for(k = -1; k < 2; k++) {
        if( i != 0 || j != 0 || k != 0 ) {
           
          label2 = label[(x+i)+dims[0]*(y+j)+dims[0]*dims[1]*(z+k)];
               
          for(label1 = 1;label1 < nc+1;label1++) { 
            if(label1 == label2) similarity_value = same;
            else if(abs(label1 - label2)<2) similarity_value = similar;
            else similarity_value = different;

            distance = sqrt(slice_width_sq[0] * abs(i) +
                            slice_width_sq[1] * abs(j) +
                            slice_width_sq[2] * abs(k));

            exponent[label1-1] += similarity_value/distance;             
          }
        }   
      }
    }
  }

  for(label1 = 0;label1 < nc;label1++)
    mrf_probability[label1] = exp(-(beta*exponent[label1])); 
  
} 


/* Iterative conditional mode */
void ICM(double *src, unsigned char *label, double *mean, double *var, int nc, int *dims, double beta, int iterations)
{
  
  int i, iter, x, y, z, z_area, y_dims, index, label_value, sum_changed;
  long area, vol;
  double val, d_pve[nc], mrf_probability[nc];
  double exponent[nc];
  unsigned char *prob, new_label;
    
  area = dims[0]*dims[1];
  vol = area*dims[2];
  
  prob = (unsigned char*)malloc(sizeof(unsigned char)*vol*nc);
  if(prob == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  /* loop over image points */
  for (z = 1; z < dims[2]-1; z++) {
    z_area=z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims=y*dims[0];
      for (x = 1; x < dims[0]-1; x++)  {
	  
        index = x + y_dims + z_area;
        label_value = (int)label[index];
        if (label_value < 1) continue;
        val = src[index];
          
        for (i=0; i<nc; i++) 
          d_pve[i] = ComputeGaussianLikelihood(val, mean[i], var[i]);        
        Normalize(d_pve, nc);

        for (i=0; i<nc; i++)
          prob[index+i*vol] = (unsigned char)round(255.0*d_pve[i]);
          
        label[index] = (unsigned char) MaxArg(d_pve, nc);
      }
    }
  }   

  for(iter=0; iter < iterations; iter++) {
    sum_changed = 0;
    /* loop over image points */
    for (z = 1; z < dims[2]-1; z++) {
      z_area=z*area;
      for (y = 1; y < dims[1]-1; y++) {
        y_dims=y*dims[0];
        for (x = 1; x < dims[0]-1; x++)  {
	  
          index = x + y_dims + z_area;
          if(label[index] > 0) {
            ComputeMrfProbability(mrf_probability, exponent, label, x, y, z, dims, nc, beta);
          
            for (i=0; i<nc; i++)
              mrf_probability[i] *= (double)prob[index+i*vol];

            new_label = (unsigned char) MaxArg(mrf_probability, nc);
            if (new_label != label[index]) {
              sum_changed++;
              label[index] = new_label;
            }
          }
        }
      }
    }
    printf("ICM: %2d voxels changed: %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iter, sum_changed);
    fflush(stdout);
    if(sum_changed < 10) break;
  }   
  printf("\n");
  free(prob);
} 

void EstimateSegmentation(double *src, unsigned char *label, unsigned char *prob, struct point *r, double *mean, double *var, int nc, int niters, int sub, int *dims, double *thresh)
{
  int i;
  int area, narea, nvol, vol, z_area, y_dims, index, ind;
  double sub_1, dmin, val, beta[1];
  double d[nc], alpha[nc], log_alpha[nc], log_var[nc];
  double pvalue[nc], psum;
  int nix, niy, niz, iters, count_change;
  int x, y, z, label_value, xBG;
  int ix, iy, iz, ind2;
  double ll, ll_old, change_ll;

  MrfPrior(label, nc, alpha, beta, 0, dims);    

  area = dims[0]*dims[1];
  vol = area*dims[2];

  /* find grid point conversion factor */
  sub_1 = 1.0/((double) sub);

  /* define grid dimensions */
  nix = (int) ceil((dims[0]-1)/((double) sub))+1;
  niy = (int) ceil((dims[1]-1)/((double) sub))+1;
  niz = (int) ceil((dims[2]-1)/((double) sub))+1; 

  narea = nix*niy;
  nvol = nix*niy*niz;

  for (i=0; i<nc; i++) log_alpha[i] = log(alpha[i]);
    
  ll_old = HUGE;
  count_change = 0;
  
  for (iters = 0; iters<niters; iters++)  {
      
    ll = 0.0;
    
    /* get means for grid points */
    GetMeansVariances(src, label, nc, r, sub, dims, thresh);    

    /* loop over image points */
    for (z = 1; z < dims[2]-1; z++) {
      z_area=z*area;
      for (y = 1; y < dims[1]-1; y++) {
        y_dims=y*dims[0];
        for (x = 1; x < dims[0]-1; x++)  {
	  
          index = x + y_dims + z_area;
          label_value = (int) label[index];
          if (label_value < 1) continue;
          val = src[index];
          
          /* find the interpolation factors */
          ix = (int)(sub_1*x);
          iy = (int)(sub_1*y);
          iz = (int)(sub_1*z);
          ind = iz*narea + iy*nix + ix;
          
          for (i=0; i<nc; i++) {
            ind2 = (i*nvol) + ind;  
            if (r[ind2].mean > TINY) {
              mean[i] = r[ind2].mean;
              var[i]  = r[ind2].var;
              log_var[i] = log(var[i]);
            }
          }
          
          /* compute energy at each point */
          dmin = HUGE; xBG = 1; 
          psum = 0.0;
          for (i=0; i<nc; i++) {
            if (fabs(mean[i]) > TINY) {
              
              d[i] = 0.5*(SQR(val-mean[i])/var[i]+log_var[i])-log_alpha[i];
              pvalue[i] = exp(-d[i])/SQRT2PI;
              psum += pvalue[i];
            } else d[i] = HUGE;
            if ( d[i] < dmin) {
              dmin = d[i];
              xBG = i;
            }
          }
          	  
          /* scale p-values to a sum of 1 */
          if (psum > TINY) {
            for (i=0; i<nc; i++) pvalue[i] /= psum;
            ll -= log(psum);
          } else  for (i=0; i<nc; i++) pvalue[i] = 0.0;
         
          for (i=0; i<nc; i++)
            prob[(vol*i) + index] = (unsigned char)ROUND(255*pvalue[i]);
         
          /* if the class has changed modify the label */
          if (xBG + 1 != label_value) label[index] = (unsigned char) (xBG + 1); 
         
        }
      }
    }

    ll /= (double)vol;
    change_ll = (ll_old - ll)/fabs(ll);
    printf("iters:%3d log-likelihood: %7.5f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters, ll);
    fflush(stdout);
    ll_old = ll;
    
    /* break if log-likelihood has not changed significantly two iterations */
    if (change_ll < TH_CHANGE) count_change++;
    if (count_change > 1) break;    
  }

  printf("\nFinal Mean*Std: "); 
  for (i=0; i<nc; i++) printf("%.3f*%.3f  ",mean[i],sqrt(var[i])); 
  printf("\n"); 
  fflush(stdout);

}


/* perform adaptive MAP on given src and initial segmentation label */
void Amap(double *src, unsigned char *label, unsigned char *prob, double *mean, int nc, int niters, int sub, int *dims, int pve, double weight_MRF)
{
  int i, nix, niy, niz;
  int area, nvol, vol;
  int histo[65536];
  double var[MAX_NC];
  double thresh[2];
  double min_src = HUGE, max_src = -HUGE;
  int cumsum[65536];
  struct point *r;
      
  area = dims[0]*dims[1];
  vol = area*dims[2];
 
  for (i=0; i<vol; i++) {
    min_src = MIN(src[i], min_src);
    max_src = MAX(src[i], max_src);
  }

  /* build histogram */
  for (i = 0; i < 65536; i++) histo[i] = 0;
  for (i=0; i<vol; i++) {
    if ((int)label[i] < 1) continue;
    histo[(int)ROUND(65535.0*(src[i]-min_src)/(max_src-min_src))]++;
  }

  /* find values between 1% and 99% quartile */
  cumsum[0] = histo[0];
  for (i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 65536; i++) cumsum[i] = (int) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
  for (i = 0; i < 65536; i++) if (cumsum[i] >= 10) break;
  thresh[0] = (double)i/65535.0*(max_src-min_src);
  for (i = 65535; i > 0; i--) if (cumsum[i] <= 990) break;
  thresh[1] = (double)i/65535.0*(max_src-min_src);
 
  /* define grid dimensions */
  nix = (int) ceil((dims[0]-1)/((double) sub))+1;
  niy = (int) ceil((dims[1]-1)/((double) sub))+1;
  niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
  nvol  = nix*niy*niz;

  r = (struct point*)malloc(sizeof(struct point)*MAX_NC*nvol);
  if(r == NULL) {
    fprintf(stderr,"Memory allocation error\n");
    exit(EXIT_FAILURE);
  }
    
  EstimateSegmentation(src, label, prob, r, mean, var, nc, niters, sub, dims, thresh);

  /* Use marginalized likelihood to estimate initial 6 classes */
  if (pve) {
    GetMeansVariances(src, label, nc, r, sub, dims, thresh);    
    ComputeInitialPveLabel(src, label, r, nc, sub, dims);
    nc += 3;
    EstimateSegmentation(src, label, prob, r, mean, var, nc, niters, sub, dims, thresh);
  }

  /* use much smaller beta for if no pve is selected */
  if(!pve) weight_MRF /= 20.0;
  
  ICM(src, label, mean, var, nc, dims, weight_MRF, 50);

  free(r);

  return;    
}

