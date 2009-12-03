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
$RCSfile: miscTemplateFunc.h,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2004/12/08 16:42:20 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MISC_TEMPLATE_FUNC_H
#define _MISC_TEMPLATE_FUNC_H

#include <stdlib.h>
#include "dcomplex.h"
#include "fcomplex.h"

inline double asDouble(dcomplex value) { return sqrt(std::norm(value)); }
inline double asDouble(fcomplex value) { return sqrt(std::norm(value)); }
inline double asDouble(signed char value) { return double(value); }
inline double asDouble(unsigned char value) { return double(value); }
inline double asDouble(short value) { return double(value); }
inline double asDouble(int value) { return double(value); }
inline double asDouble(unsigned int value) { return double(value); }
inline double asDouble(float value) { return double(value); }
inline double asDouble(double value) { return double(value); }

// Value swapping
#if defined(__GNUC__) && (__GNUC__ <= 2)
template <class Type>
inline void swap(Type& x, Type& y) { Type tmp(x); x = y; y = tmp; }
#endif /* (bert) end of newly conditionalized code. */

// Power functions
inline double intPower(double x, int y) {
  if (!y) return(1.0);
  if (!x) return(0.0);
  if (x == 1.0) return(1.0);

  register double result = x;
  for (register unsigned i = abs(y) - 1; i; i--) 
    result *= x;

  return (y < 0) ? 1.0/result : (double) result;
}

#if defined(__GNUC__) && (__GNUC__ <= 2) /* (bert) removed unneeded min/max */
// Two argument min/max
template <class Type>
inline Type min(const Type& x, const Type& y) { return (x < y) ? x : y; }
template <class Type>
inline Type max(const Type& x, const Type& y) { return (x > y) ? x : y; }
#endif /* (bert) end of newly conditionalized code */

// Three argument min/max
template <class Type>
inline Type min(const Type& x, const Type& y, const Type& z) { 
  Type value = (x < y) ? x : y; return (value < z) ? value : z; }
template <class Type>
inline Type max(const Type& x, const Type& y, const Type& z) { 
  Type value = (x > y) ? x : y; return (value > z) ? value : z; }

// Clamping
template <class Type>
inline Type clamp(const Type& value, const Type& a, const Type& b) {
  if (b < a) {
    if (value > a)
      return a;
    if (value < b)
      return b;
  }
  else {
    if (value < a)
      return a;
    if (value > b)
      return b;
  }

  return value;
}

template <class Type>
Type nextPowerOf2(Type value)
{
  double x;

  if (value >= 0) {
    x = 1;
    while (x < value)
      x *= 2;
  }
  else {
    x = -1;
    while (x > value)
      x *= 2;
  }

  return Type(x);
}

// Explicit template instantiation of the above function, to
// avoid unpleasant compiler warnings.
//
template <> 
inline unsigned nextPowerOf2(unsigned value)
{
  double x;

  x = 1;
  while (x < value)
    x *= 2;

  return unsigned(x);
}

template <class Type> 
int isPowerOf2(Type value)
{
  if (!value)
    return 0;

  if (value < 0)
    value = -value;

  Type x = 1;
  while (x < value)
    x *= 2;

  return x == value;
}

// Explicit template instantiation of the above function, to
// avoid unpleasant compiler warnings.
//
template <> 
inline int isPowerOf2(unsigned value)
{
  if (!value)
    return 0;

  unsigned x = 1;
  while (x < value)
    x *= 2;

  return x == value;
}
#endif
