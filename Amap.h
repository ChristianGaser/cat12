/*
 * Christian Gaser
 * $Id$ 
 *
 */

#define SQRT2PI 2.506628
#define G 6

#define MAX_NC 5
#define TH_COLOR 1
#define TH_CHANGE 0.00001
#define TINY 1e-15 
#ifndef HUGE
#define HUGE 1e15 
#endif

#define NOPVE 0
#define MARGINALIZED 1
#define KMEANS 2
#define BAYES 3

#define CSFLABEL   0
#define GMCSFLABEL 1
#define GMLABEL    2
#define WMGMLABEL  3
#define WMLABEL    4

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

struct point {
  double mean;
  double var;
};

struct ipoint {
  int n;
  double s;
  double ss;
};
