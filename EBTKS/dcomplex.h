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
$RCSfile: dcomplex.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _DCOMPLEX_H
#define _DCOMPLEX_H

#include <complex>		/* (bert) changed from complex.h */
typedef std::complex<double> dcomplex;

int operator < (const dcomplex&, const dcomplex&);
int operator <= (const dcomplex&, const dcomplex&);
int operator > (const dcomplex&, const dcomplex&);
int operator >= (const dcomplex&, const dcomplex&);
int operator && (const dcomplex&, const dcomplex&);
int operator || (const dcomplex&, const dcomplex&);

#endif
