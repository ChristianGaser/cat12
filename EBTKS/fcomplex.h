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
$RCSfile: fcomplex.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _FCOMPLEX_H
#define _FCOMPLEX_H

#include <complex>		/* (bert) changed from complex.h */
typedef std::complex<float> fcomplex;
#if 0				/* (bert) removed code here */

//Code written by Georges Aboutanos
//Based on Rose and Stroustrup and Sun complex.h.

#include <errno.h>
#include <math.h>
#include <dcomplex.h>


/****************************class fcomplex********************/
class fcomplex {
  float re, im;
public:
  fcomplex() {re=0.0; im=0.0; }
  fcomplex(float r, float i=0.0) { re=r; im=i; }
  fcomplex(const dcomplex& a) { re = real(a); im = imag(a); }

  operator dcomplex () { return dcomplex(re, im); }
  
  friend fcomplex conj(fcomplex);
  friend float    arg(fcomplex);
  friend float    imag(fcomplex);
  friend float    real(fcomplex);
  friend float    norm(fcomplex);
  friend float    abs(fcomplex);
  friend float    magnitude(fcomplex);
  friend float    magln(fcomplex);
  friend float    mag10(fcomplex);
  friend float    power(fcomplex);
  friend fcomplex operator+(fcomplex, fcomplex);
  friend fcomplex operator-(fcomplex);
  friend fcomplex operator-(fcomplex, fcomplex);
  friend fcomplex operator*(fcomplex, fcomplex);
  friend fcomplex operator/(fcomplex, fcomplex);
  friend int operator==(fcomplex, fcomplex);
  friend int operator!=(fcomplex, fcomplex);
  
  void operator+=(fcomplex);
  void operator-=(fcomplex);
  void operator*=(fcomplex);
  void operator/=(fcomplex);

  friend ostream& operator << (ostream& os, fcomplex value);
  friend istream& operator >> (istream& is, fcomplex& value);
};

inline fcomplex operator+(fcomplex a1, fcomplex a2)
{
   return fcomplex(a1.re+a2.re,a1.im+a2.im);
}

inline fcomplex operator-(fcomplex a1, fcomplex a2)
{
   return fcomplex(a1.re-a2.re,a1.im-a2.im);
}

inline fcomplex operator-(fcomplex a)
{
   return fcomplex(-a.re, a.im);
}

inline fcomplex operator*(fcomplex a1, fcomplex a2)
{
   return fcomplex((a1.re*a2.re)-(a1.im*a2.im),(a1.re*a2.im)+(a1.im*a2.re));
}

inline fcomplex operator/(fcomplex a1, fcomplex a2)
{
  float denominator = a2.re*a2.re + a2.im*a2.im;
  return fcomplex(((a1.re*a2.re) + (a1.im*a2.im))/denominator,
		  ((a1.im*a2.re) - (a1.re*a2.im))/denominator);
}

inline fcomplex conj(fcomplex a)
{
   return fcomplex(a.re, -a.im);
}

inline float arg(fcomplex a)
{
   return atan2(a.im, a.re);
}

inline float real(fcomplex a)
{
   return a.re;
}

inline float imag(fcomplex a)
{
   return a.im;
}

inline float norm(fcomplex a)
{
  return a.re*a.re + a.im*a.im;
}

inline float abs(fcomplex a)
{
  return float(sqrt(a.re*a.re + a.im*a.im));
}

inline float magnitude(fcomplex a) 
{ 
  return abs(a); 
}

inline float magln(fcomplex a)
{
  return float(log(sqrt(a.re*a.re + a.im*a.im)));
}

inline float mag10(fcomplex a)
{
  return float(log10(sqrt(a.re*a.re + a.im*a.im)));
}

inline float power(fcomplex a)
{
  return a.re*a.re + a.im*a.im;
}

inline int operator==(fcomplex a, fcomplex b)
{
   return (a.re==b.re && a.im==b.im);
}

inline int operator!=(fcomplex a, fcomplex b)
{
   return (a.re!=b.re || a.im!=b.im);
}

inline void fcomplex::operator+=(fcomplex a)
{
   re += a.re;   im += a.im;
}

inline void fcomplex::operator-=(fcomplex a)
{
   re -= a.re;   im -= a.im;
}

inline void fcomplex::operator*=(fcomplex a)
{
   re = (re*a.re)-(im*a.im); im = (re*a.im)+(im*a.re);
}

inline void fcomplex::operator/=(fcomplex a)
{
  float denominator = a.re*a.re + a.im*a.im;
  re = ((re*a.re) + (im*a.im))/denominator;
  im = ((im*a.re) - (re*a.im))/denominator;
}
#endif /* (bert) end of removed code */

int operator < (const fcomplex&, const fcomplex&);
int operator <= (const fcomplex&, const fcomplex&);
int operator > (const fcomplex&, const fcomplex&);
int operator >= (const fcomplex&, const fcomplex&);
std::ostream& operator << (std::ostream& os, fcomplex value);
std::istream& operator >> (std::istream& is, fcomplex&);
#endif

