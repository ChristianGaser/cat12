#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Amap.h"

/* use always 6 classes */
#define Kb 6
#define MAXK 30

void Bayes(double *src, unsigned char *label, unsigned char *priors, unsigned char *prob, double *separations, int *dims, int correct_nu)
{
  int i, j, k, l, k1, subit, n_loops, subsample_warp;
  int z_area, y_dims, count, subsample, masked_smoothing;
  double ll = -HUGE, llr=0.0, *nu, fwhm[3];
  double mn[MAXK], vr[MAXK], mg[MAXK], mn2[MAXK], vr2[MAXK], mg2[MAXK];
  double mom0[MAXK], mom1[MAXK], mom2[MAXK], mgm[MAXK];
  double q[MAXK], bt[MAXK], b[MAXK], qt[MAXK];
  double tol1 = 1e-3, bias_fwhm, sum;
  float *flow;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  int K = 0, K2 = 0;
  
  /* multiple gaussians are not yet working */
  int ngauss[6] = {1,1,1,1,1,1};
  int iters_EM[5] = {2, 10, 10, 10, 10};

  int histo[65536], lkp[MAXK];
  double mn_thresh, mx_thresh;
  double min_src = HUGE, max_src = -HUGE;
  int cumsum[65536], order_priors[6] = {2,0,1,3,4,5};

  subsample_warp = ROUND(9.0/(separations[0]+separations[1]+separations[2]));
  
  flow  = (float *)malloc(sizeof(float)*vol*3);
  /* initialize flow field with zeros */
  for (i = 0; i < (vol*3); i++) flow[i] = 0.0;

  if (correct_nu) nu = (double *)malloc(sizeof(double)*vol);

  bias_fwhm = 60.0;
  if (correct_nu) {
    /* use larger filter */
    for(i=0; i<3; i++) fwhm[i] = 2*bias_fwhm;
       
    /* estimate mean */
    count = 0; sum = 0.0;
    for (i = 0; i < vol; i++) {
      if(src[i]>0) {
        sum += src[i];
        count++;
      }
    }

    if (count==0) return;

    sum /= (double)count;


    for (i = 0; i < vol; i++) {
      if(src[i]>0)
        nu[i] = src[i] - sum;
      else nu[i] = 0;
    }
        
    /* use subsampling for faster processing */
    subsample = 2;
    masked_smoothing = 1;
    smooth_subsample_double(nu, dims, separations, fwhm, masked_smoothing, subsample);
        
    /* and correct bias */
    for (i = 0; i < vol; i++)
      if(src[i]>0)
        src[i] -= nu[i];

  }

  for (i=0; i<vol; i++) {
    min_src = MIN(src[i], min_src);
    max_src = MAX(src[i], max_src);
  }

  /* build histogram */
  for (i = 0; i < 65536; i++) histo[i] = 0;
  for (i=0; i<vol; i++) {
    if (src[i] == 0) continue;
    histo[(int)ROUND(65535.0*(src[i]-min_src)/(max_src-min_src))]++;
  }

  /* find values between 0.1% and 99.9% quartile */
  cumsum[0] = histo[0];
  for (i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 65536; i++) cumsum[i] = (int) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
  for (i = 0; i < 65536; i++) if (cumsum[i] >= 1) break;
  mn_thresh = (double)i/65535.0*(max_src-min_src);
  for (i = 65535; i > 0; i--) if (cumsum[i] <= 999) break;
  mx_thresh = (double)i/65535.0*(max_src-min_src);

  /* K = sum(ngauss) */
  for (k1=0; k1<Kb; k1++)   K  += ngauss[k1]; 
  for (k1=0; k1<Kb-1; k1++) K2 += ngauss[k1]; 

  /* build lkp */
  l = 0;
  for (k1=0; k1<Kb; k1++) {
    for (j=0; j<ngauss[k1]; j++) {
      lkp[l] = k1; 
      l++;
    }
  }

  /* initialize mean, var, mg */
  for (k=0; k<K; k++) {
    mn[k] = mx_thresh * 1.0/(double)(k+1);
    mg[k] = 1.0/(double)K;
    vr[k] = mx_thresh*mx_thresh + TINY;
  }
          
  /* start with a few EM iterations and after nu-corection use more iterations */
  for (j=0; j < 3; j++) {
    for (subit=0; subit<iters_EM[j]; subit++) {
      double oll = ll;
      ll = llr;
      for (k=0; k<K; k++) {
        mom0[k] = 0.0;
        mom1[k] = 0.0;
        mom2[k] = 0.0;
      }
      for (k1=0; k1<Kb; k1++) mgm[k1] = 0.0;
    
      for (k=0; k<K; k++) {
        mg2[k] = mg[k];
        mn2[k] = mn[k];
        vr2[k] = vr[k];
      }
      
      for (i=0; i<vol; i++) {
        if((src[i]>mn_thresh) && (src[i]<mx_thresh)) {
          double s = TINY;
          int k2 = 0;
          for (k1=0; k1<Kb; k1++) {
            bt[k1] = (double) priors[i+(vol*k1[order_priors])]; 
            for (l=0; l<ngauss[k1]; l++) {
              s += bt[k1]*mg[k2];
              k2++;
            }
          }

          for (k1=0; k1<Kb; k1++) {
            bt[k1] /= s;
            mgm[k1] += bt[k1];
          }
          double sq = TINY;
          double qmax = -HUGE;
          int kmax;
          
          for (k=0; k<K; k++) {
            q[k] = mg[k]*bt[lkp[k]]*exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
            sq += q[k];
            if (q[k] > qmax) {
              qmax = q[k];
              if (k>K2-1) kmax = 0; else kmax = k + 1;
            } 
          }
          
          /* prepare prob for warping */
          sq = TINY;
          for (k1=0; k1<Kb; k1++) qt[k1] = 0.0;
          for (k=0; k<K; k++) { 
            double p1 = mg[k]*bt[lkp[k]]*exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
            qt[lkp[k]] += p1;
            sq += qt[lkp[k]];
          }
          for (k1=0; k1<Kb; k1++)
            prob[i+(vol*k1[order_priors])] = (unsigned char)ROUND(255.0*qt[k1]/sq);

          label[i] = kmax;
          ll += log10(sq);
          for (k=0; k<K; k++) {
            double p1 = q[k]/sq; mom0[k] += p1;
            p1 *= src[i];        mom1[k] += p1;
            p1 *= src[i];        mom2[k] += p1;
          }      
        }
      } 

      if ((subit==1) && (j>0)) WarpPriors(prob, priors, flow, dims, 6, 0, subsample_warp);
          
      for (k=0; k<K; k++) {
        mg[k] = (mom0[k]+TINY)/(mgm[lkp[k]]+TINY);
        mn[k] = mom1[k]/(mom0[k]+TINY);
        vr[k] = (mom2[k]-SQR(mom1[k])/mom0[k]+1e6*TINY)/(mom0[k]+TINY);
        vr[k] = MAX(vr[k],TINY);
        /* rescue previous values in case of nan */
        if ((isnan(mg[k])) || (isnan(mn[k])) || (isnan(vr[k]))) {
          mg[k] = mg2[k];
          mn[k] = mn2[k];
          vr[k] = vr2[k];
          break;
        }
      }
      
      printf("%7.4f\b\b\b\b\b\b\b",ll/vol);    
      printf("\n");
      fflush(stdout);
      if((ll-oll)<tol1*vol) break;
    }
    	  
    if(correct_nu) {
      for (i = 0; i < vol; i++) {
        nu[i] = 0.0;
        /* only use values above threshold for nu-estimate */
        if ((src[i] > mn_thresh) && (label[i] < 4)) {
          double val_nu = src[i]-mn[label[i]-1];
          if ((finite(val_nu)))
            nu[i] = val_nu;
        }
      }

      /* smoothing of residuals */
      for(i=0; i<3; i++) fwhm[i] = bias_fwhm;
        
      /* use subsampling for faster processing */
      subsample = 2;
      masked_smoothing = 0;

      smooth_subsample_double(nu, dims, separations, fwhm, masked_smoothing, subsample);

      /* apply nu correction to source image */
      for (i=0; i<vol; i++) {
        if (nu[i] > 0.0)
          src[i] -= nu[i];
      }
    }    
  }

  for (i=0; i<vol; i++) {
    if (src[i] != 0) {
      for (k1=0; k1<Kb; k1++) {
        bt[k1] = (double) priors[i+(vol*k1[order_priors])]; 
        qt[k1] = 0.0;
      }
      double s = TINY;
      for (k=0; k<K; k++) { 
        b[k] = bt[lkp[k]]*mg[k];
        s += b[k];
      }
      double sq = TINY;
      for (k=0; k<K; k++) { 
        double p1 = exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));
        qt[lkp[k]] += p1*b[k]/s;
        sq += qt[lkp[k]];
      }
      double qmax = -HUGE;
      int kmax;
      double psum = 0.0;
      for (k1=0; k1<Kb; k1++) {
        double p1 = qt[k1]/sq;
        /* save probabilities */
        prob[i+(vol*k1[order_priors])] = (unsigned char)ROUND(255.0*p1);
        psum += p1;
        if (p1 > qmax) { 
          qmax = p1;     
          kmax = k1 + 1;
        }   
      }
      /* label only if sum of all probabilities is > 0 */
      if (psum > 0)
        label[i] = kmax;
      else
        label[i] = 0;
    } else label[i] = 0;
  }

/*  for (k=0; k<Kb; k++) 
    printf("%g %g\n",mn[k],sqrt(vr[k]));
*/    

  free(flow);
  if (correct_nu) free(nu);
    
  return;
}  
