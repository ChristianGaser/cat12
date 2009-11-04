/*
 * Christian Gaser
 * $Id$ 
 *
 */

#define SPLINESMOOTH 0
#define SQRT2PI 2.506628
#define G 6

#define MAX_NC 6
#define TH_COLOR 1
#define TH_CHANGE 0.00001
#define TINY 1e-15 
#ifndef HUGE
#define HUGE 1e15 
#endif

#define NOPVE 0
#define MARGINALIZED 1
#define KMEANS 2

#define BKGCSFLABEL 0
#define CSFLABEL    1
#define GMCSFLABEL  2
#define GMLABEL     3
#define WMGMLABEL   4
#define WMLABEL     5

#define NOLABEL 0
#define LABEL 1
#define PVELABEL 2

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifndef MIN3
#define MIN3(a,b,c) (MIN(a,MIN(b,c)))
#endif

extern double Kmeans(double *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *separations, int *dims, int thresh_mask, int thresh_kmeans, int iters_nu, int pve);

extern void Bayes(double *src, unsigned char *label, unsigned char *priors, unsigned char *mask, double *separations, int *dims, int correct_nu);

extern void WarpPriors(unsigned char *prob, unsigned char *priors, unsigned char *mask, float *flow, int *dims, int loop, int samp);

extern void Amap(double *src, unsigned char *label, unsigned char *prob, double *mean, int nc, int niters, int sub, int *dims, int pve, double weight_MRF);

extern void Pve6(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims, int update_label);

extern void MrfPrior(unsigned char *label, int nc, double *alpha, double *beta, int init, int *dims);

extern void sampn(int dm[], float f[], int n, int mm, double x, double y, double z, double v[]);

extern int splineSmooth( double *src, double lambda, double distance, int subsample, double *separations, int *dims);

struct point {
  double mean;
  double var;
};

struct ipoint {
  int n;
  double s;
  double ss;
};
