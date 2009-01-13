#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "Amap.h"

void Bayes(double *src, unsigned char *label, unsigned char *priors, unsigned char *mask, double *separations, int *dims, int correct_nu)
{
  int i, j, k, k1, subit;
  int z_area, y_dims;
  double ll = -FLT_MAX, llr=0.0, *nu;
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  int K = 0, K2 = 0;
  int Kb = 6; /* use always 6 classes */
  
  /* multiple gaussians are not yet working */
//  int ngauss[6] = {2,2,2,4,2,2};
  int ngauss[6] = {1,1,1,1,1,1};

  int histo[65536];
  double mn_thresh, mx_thresh;
  double min_src = FLT_MAX, max_src = -FLT_MAX;
  int cumsum[65536], order_priors[6] = {2,0,1,3,4,5};

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

  if (correct_nu) nu = (double *)malloc(sizeof(double)*vol);

  /* K = sum(ngauss) */
  for (k1=0; k1<Kb; k1++)   K  += ngauss[k1]; 
  for (k1=0; k1<Kb-1; k1++) K2 += ngauss[k1]; 

  /* build lkp */
  int lkp[K];
  int l = 0;
  for (k1=0; k1<Kb; k1++) {
    for (j=0; j<ngauss[k1]; j++) {
      lkp[l] = k1; 
      l++;
    }
  }

  /* initialize mean, var, mg */
  double mn[K], vr[K], mg[K], mn2[K], vr2[K], mg2[K];
  for (k=0; k<K; k++) {
    mn[k] = mx_thresh * drand48();
    mg[k] = 1.0/(double)K;
    vr[k] = mx_thresh*mx_thresh + TINY;
  }
        
  double mom0[K], mom1[K], mom2[K], mgm[Kb];
  double q[K], bt[Kb], b[K], qt[Kb];
  double tol1 = 1e-4;
  
  /* start with a few EM iterations and after nu-corection use more iterations */
  int iters_EM[5] = {2, 10, 10, 10, 10};
  for (j=0; j < 5; j++) {
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
          double qmax = -FLT_MAX;
          int kmax;
          
          for (k=0; k<K; k++) {
            q[k] = mg[k]*bt[lkp[k]]*exp(SQR(src[i]-mn[k])/(-2*vr[k]))/(SQRT2PI*sqrt(vr[k]));            
            sq += q[k];
            if (q[k] > qmax) {
              qmax = q[k];
              if (k>K2-1) kmax = 0; else kmax = k + 1;
            } 
          }
          
          if (kmax < 4) label[i] = kmax; else label[i] = 0;
          ll += log10(sq);
          for (k=0; k<K; k++) {
            double p1 = q[k]/sq; mom0[k] += p1;
            p1 *= src[i];        mom1[k] += p1;
            p1 *= src[i];        mom2[k] += p1;
          }      
        } 
      }  
          
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
      
      if((ll-oll)<tol1*vol) break;
      printf("%5.4f\n",ll/vol);    
    }
    	  
    if(correct_nu) {
      for (i = 0; i < vol; i++) {
        nu[i] = 0.0;
        /* only use values above threshold where mask is defined for nu-estimate */
        if ((src[i] > mn_thresh) && (mask[i] > 64)) {
          double val_nu = src[i]/mn[label[i]-1];
          if ((isfinite(val_nu))) {
            nu[i] = val_nu;
          }
        }
      }
      printf("Fitting splines ...\n");
      /* spline estimate with increasing spatial resolution */
      splineSmooth(nu, 0.01, MAX(500,1500.0/(j+1)), 4, separations, dims);
      
      /* apply nu correction to source image */
      for (i=0; i<vol; i++) {
        if (nu[i] > 0.0)
          src[i] /= nu[i];
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
      double qmax = -FLT_MAX;
      int kmax;
      double psum = 0.0;
      for (k1=0; k1<Kb; k1++) {
        double p1 = qt[k1]/sq;
        psum += p1;
        if (p1 > qmax) { 
          qmax = p1;     
          kmax = k1 + 1;
        }   
      }
      /* label only if sum of all probabilities is > 0.1 */
      if ((psum > 0) && (kmax < 4))
        label[i] = kmax;
      else
        label[i] = 0;
    } else label[i] = 0;
  }

  for (k=0; k<Kb; k++) 
    printf("%g %g\n",mn[k],sqrt(vr[k]));
    
  if (correct_nu) free(nu);
    
  return;
}  
