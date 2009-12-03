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
$RCSfile: Complex.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _COMPLEX_H
#define _COMPLEX_H

// A few functions that define (bogus) math ops for complex

#include <complex.h>
#ifdef __GNUC__
  typedef complex<double> dcomplex;
#else
  typedef complex dcomplex;
#endif

static _errorCount = 100;

int operator < (const dcomplex&, const dcomplex&) {
  if (_errorCount) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount--;
  }
  return 0;
}

int operator <= (const dcomplex&, const dcomplex&) {
  if (_errorCount) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount--;
  }
  return 0;
}

int operator > (const dcomplex&, const dcomplex&) {
  if (_errorCount) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount--;
  }
  return 0;
}

int operator >= (const dcomplex&, const dcomplex&) {
  if (_errorCount) {
    cerr << "Comparison of dcomplex numbers undefined" << endl;
    _errorCount--;
  }
  return 0;
}

#endif
