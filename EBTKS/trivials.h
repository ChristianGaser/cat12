/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: trivials.h,v $
$Revision: 1.2 $
$Author: jason $
$Date: 2002/03/20 21:42:44 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef MIN
#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))
#endif

#ifndef ROUND
#define ROUND(x) ((int) ((x) + 0.5))
#endif

#ifndef SIGN
#define SIGN(x)   (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0))
#endif

#ifndef SQR
#define SQR(x)    ((x)*(x))
#endif

#ifndef TRIVIALS_H
#define TRIVIALS_H

#include <math.h>
#include <cstdlib>
#include "miscTemplateFunc.h"

#ifdef FALSE
#undef FALSE
#undef TRUE
#endif

const char FALSE = 0;
const char TRUE  = 1;

const int SUCCESS = 1;
const int FAILURE = 0;

//#define NAN -MAXDOUBLE

#define toggleOn(value, bit)        ((value) & (1 << (bit)))
#define setToggle(value, bit)       ((value) = (value) & (1 << (bit)))
#define clearToggle(value, bit)     ((value) = (value) & ~(1 << (bit)))

inline int intCompareAsc(const void *x, const void *y) {
  return *(int *) x - *(int *) y; }
inline int floatCompareAsc(const void *x, const void *y) {
  float diff = *(float *) x - *(float *) y; return SIGN(diff); }
inline int doubleCompareAsc(const void *x, const void *y) {
  double diff = *(double *) x - *(double *) y; return SIGN(diff); }

inline double gauss(double mean, double std) {
  double v1 = 2.0*rand() - 1.0;
  double v2 = 2.0*rand() - 1.0;
  double s = v1*v1 + v2*v2;

  while (s >= 1.0) {
    v1 = 2.0*rand() - 1.0;
    v2 = 2.0*rand() - 1.0;
    s  = v1*v1 + v2*v2;
  }

  return mean + std*v1*sqrt(-2.0*log(s)/s);
}

#endif
