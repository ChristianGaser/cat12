/*
 * Amap.C
 * 
 * author : Jagath C. Rajapakse
 *
 * Statistical approach to single-channel MR brain scans
 * J. C. Rajapakse, J. N. Giedd, and J. L. Rapoport
 * IEEE Transactions on Medical Imaging, Vol 16, No 2, 1997
 *
 * Comments to raja@cns.mpg.de, 15.10.96
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#define SQRT2PI 2.506628
#define G 30

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#define MAX_NC 100

#ifndef ROUND
#define ROUND( x ) ((int) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif
void MrfPrior(unsigned char *label, int nc, double *alpha, double *beta, int BG, int init, int *dims);

struct point {
        double mean;
        double var;
};

struct ipoint {
        int n;
        double s;
        double ss;
};

/* calculate the mean and variance for every class on a grid size SUBxSUBxSUB */
static void get_means(double *src, unsigned char *label, int nc, struct point *r, int sub, int BG, int *dims, double mn_thresh, double mx_thresh)
{
        int i, j, ind;
        int area, narea, nvol, zsub, ysub, xsub, yoffset, zoffset;
        int zsub2, ysub2;
        int nix, niy, niz, k, l, m, z, y, x, labval;
        double val;
        struct ipoint *ir;
        int labval_BG;

        area = dims[0]*dims[1];

        /* define grid dimensions */
        nix = (int) ceil((dims[0]-1)/((double) sub))+1;
        niy = (int) ceil((dims[1]-1)/((double) sub))+1;
        niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
        narea = nix*niy;
        nvol = nix*niy*niz;
        
        ir = (struct ipoint*)malloc(sizeof(struct ipoint)*nc*nvol);

        for (i=0; i<nc; i++) {
                for (j=0; j<nvol; j++) {
                        ind = (i*nvol)+j; 
                        ir[ind].n = 0;
                        ir[ind].s = 0.0;
                        ir[ind].ss = 0.0;
                }
        }
        

        /* loop over neighborhoods of the grid points */
#ifdef _OPENMP
        # pragma omp parallel default(none) \
        private(k,l,m,z, y, x, zsub, zsub2, ysub, ysub2, xsub, zoffset, yoffset, labval, val, ind, labval_BG) \
        shared(r,nc,sub, niz, dims, area, narea, nix, label, BG, src, mn_thresh, mx_thresh, nvol, ir, niy)
#endif
        for (k=-sub; k<=sub; k++) {
                for (l=-sub; l<=sub; l++) {
                        for (m=-sub; m<=sub; m++) {
#ifdef _OPENMP
        # pragma omp for
#endif
                                for (z = 0; z < niz; z++)        {
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
                                                                                labval = (int)label[zsub2 + ysub2 + xsub];
                                                                                labval_BG = labval - BG;
                                                                                if (labval_BG < 0) continue;
                                                                                val = src[zsub2 + ysub2 + xsub];
                                                                                
                                                                                /* exclude values out of quartile 1/99% */
                                                                                if ((val<mn_thresh) || (val>mx_thresh)) continue;
                                                                                ind = ((labval_BG)*nvol)+yoffset+x;
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
#ifdef _OPENMP
        # pragma omp parallel private(i, j, ind) 
#endif
        for (i=0; i<nc; i++) {
#ifdef _OPENMP
        # pragma omp for
#endif
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

/* perform adaptive MAP on given src and initial segmentation label */
void Amap(double *src, unsigned char *label, unsigned char *prob, double *mean, int nc, int BG, int niters, int nflips, int sub, int *dims, double weight_MRF)
{
        int i, index;
        int area, narea, nvol, vol, z_area, y_dims;
        int histo[65536];
        double sub_1, beta[1], dmin, val;
        double var[MAX_NC], d[MAX_NC], alpha[MAX_NC], log_alpha[MAX_NC], log_var[MAX_NC];
        double pvalue[MAX_NC], psum, error;
        int nix, niy, niz, iters;
        int x, y, z, labval, ind, xBG, flag;
        int ix, iy, iz, iBG, ind2;
        double first, mn_thresh, mx_thresh;
        double min_src = FLT_MAX, max_src = -FLT_MAX;
        int cumsum[65536];
        struct point *r;
                        
#ifdef _OPENMP
        #pragma omp parallel
        #pragma omp master
        printf("Number of threads: %d\n",omp_get_num_threads());
#endif

        MrfPrior(label, nc, alpha, beta, BG, 0, dims);
        
        /* use pre-defined MRF prior */
        if (weight_MRF < 1.0) {
	        beta[0] = weight_MRF;
        	printf("weighted MRF prior beta: %g\n",beta[0]);
        }
        for (i=0; i<nc; i++)
                log_alpha[i] = log(alpha[i]);

        area = dims[0]*dims[1];
        vol = area*dims[2];
 
        for (i=0; i<vol; i++) {
                min_src = MIN(src[i], min_src);
                max_src = MAX(src[i], max_src);
        }

        /* build histogram */
        for (i = 0; i < 65536; i++) histo[i] = 0;
        for (i=0; i<vol; i++) {
                if ((int)label[i] < BG) continue;
                histo[(int)ROUND(65535.0*(src[i]-min_src)/(max_src-min_src))]++;
        }

        /* find values between 1% and 99% quartile */
        cumsum[0] = histo[0];
        for (i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
        for (i = 0; i < 65536; i++) cumsum[i] = (int) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
        for (i = 0; i < 65536; i++) if (cumsum[i] >= 10) break;
        mn_thresh = (double)i/65535.0*(max_src-min_src);
        for (i = 65535; i > 0; i--) if (cumsum[i] <= 990) break;
        mx_thresh = (double)i/65535.0*(max_src-min_src);
 
        /* find grid point conversion factor */
        sub_1 = 1.0/((double) sub);

        /* define grid dimensions */
        nix = (int) ceil((dims[0]-1)/((double) sub))+1;
        niy = (int) ceil((dims[1]-1)/((double) sub))+1;
        niz = (int) ceil((dims[2]-1)/((double) sub))+1; 

        narea = nix*niy;
        nvol = nix*niy*niz;

        r = (struct point*)malloc(sizeof(struct point)*nc*nvol);

        error = 0.0;
        for (iters = 0; iters<=niters; iters++)        {
                int flips = 0;
                
                /* get means for grid points */
                get_means(src, label, nc, r, sub, BG, dims, mn_thresh, mx_thresh);                

                /* loop over image points */
                for (z = 1; z < dims[2]-1; z++) {
                        z_area=z*area;
                        for (y = 1; y < dims[1]-1; y++) {
                                y_dims=y*dims[0];
                                for (x = 1; x < dims[0]-1; x++)        {
	        
                                        index = x + y_dims + z_area;
                                        labval = (int)label[index];
                                        if (labval < BG) continue;
                                        val = src[index];
                                        
                                        /* find the interpolation factors */
                                        ix = (int)(sub_1*x), iy = (int)(sub_1*y), iz = (int)(sub_1*z);
                                        ind = iz*narea + iy*nix + ix;
                                        
                                        for (i=0; i<nc; i++) {
                                                ind2 = (i*nvol) + ind;                                                
                                                if (r[ind2].mean > 0.0) {
                                                        mean[i] = r[ind2].mean;
                                                        var[i]  = r[ind2].var;
                                                        log_var[i] = log(var[i]);
                                                }
                                        }

                                        /* compute energy at each point */
                                        flag = 0;
                                        dmin = 1e15; xBG = BG; 
                                        psum = 0.0;
                                        for (i=0; i<nc; i++) {
                                                first=0.0;
                                                iBG = i+BG;
                                                if ((int)label[index-1] == iBG)       first++;
                                                if ((int)label[index+1] == iBG)       first++;
                                                if ((int)label[index-dims[0]] == iBG) first++;
                                                if ((int)label[index+dims[0]] == iBG) first++;
                                                if ((int)label[index-area] == iBG)    first++;
                                                if ((int)label[index+area] == iBG)    first++;

                                                if (mean[i] > 0.0) {
                                                        if (var[i] > 0.0) {
                                                                d[i] = 0.5*(SQR(val-mean[i])/var[i]+log_var[i])-log_alpha[i]-beta[0]*first;
                                                                pvalue[i] = exp(-d[i])/SQRT2PI;
                                                                psum += pvalue[i];
                                                        } else flag = 1;
                                                } else d[i] = 1e15;
                                                if ( d[i] < dmin) {
                                                        dmin = d[i]; xBG = i;
                                                }
                                        }
                                        
                                        if (flag == 1) {
                                                dmin = 1e15; xBG = BG;
                                                for (i=0; i<nc; i++) {
                                                        if (mean[i] > 0.0) d[i] = SQR(val-mean[i])-log_alpha[i]-beta[0]*first;
                                                        else d[i] = 1e15;
                                                        if ( d[i] < dmin) {
                                                                dmin = d[i];
                                                                xBG = i;
                                                        }
                                                }
                                        }
	        
                                 /* scale p-values to a sum of 1 */
                                 if (psum > 0.0) {
                                         for (i=0; i<nc; i++) pvalue[i] /= psum;
                                 } else  for (i=0; i<nc; i++) pvalue[i] = 0.0;;
                                 
                                 for (i=0; i<nc; i++)
                                         prob[(vol*i) + index] = (unsigned char)ROUND(255*pvalue[i]);
                                 
                                 /* if the class has changed increment flips and change the label */
                                 if (xBG + BG != labval) {flips++; label[index] = (unsigned char) (xBG + BG); }
                                 
                         }
                 }
         }
                
                printf("iters:%3d flips:%6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters, flips);
                
                if (flips <= nflips) break;                
        }

        printf("\nFinal Mean*Std: "); 
        for (i=0; i<nc; i++) printf("%5.3f*%5.3f        ",mean[i],sqrt(var[i])); 
        printf("\n"); 

        free(r);

        return;                
}
