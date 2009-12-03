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
$RCSfile: Polynomial.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "SimpleArray.h"

#ifndef USE_FLMAT
#define USE_FLMAT
#endif

#include "Matrix.h"

#ifndef FloatArray
typedef SimpleArray<float> FloatArray;
#endif
#ifndef DblArray
typedef SimpleArray<double> DblArray;
#endif
#ifndef FlMat
typedef Mat<int> IntMat;
#endif
#ifndef FlMat
typedef Mat<float> FlMat;
#endif
#ifndef DblMat
typedef Mat<double> DblMat;
#endif

/************************************************************************
 *
 * Polynomial class.
 *
 * This class stores the coefficients of an n-th order, m-dimensional
 * polynomial fitted in a least squares sense to the points specified
 * in the arguments of the constructors. The private members <_coef> 
 * and  <_expComb> keep all coefficients and all required combinations 
 * of exponents for each dimension, respectively. 
 *
 * I.e., a two-dimensional polynomial
 *
 *                                  2        2      3
 *        f(x, y) = a0 + a1*x + a2*x y + a3*y + a4*y
 *
 * is represented as
 *
 *        _coef:    [a0 a1 a2 a3 a4]
 *                  
 *        _expComb: [0 1 2 0 0]  (exponents of x-coord)
 *                  [0 0 1 2 3]  (exponents of y-coord)
 *
 ************************************************************************/

class Polynomial {
  IntMat   _expComb;
  DblArray _coef;
  unsigned _nDimensions; // Provided for the sake of speed
  unsigned _nCoef;       // Provided for the sake of speed

public:
// Constructors/destructor
  // General n-D polynomials, dimension and coefficients specified.
  //  
  // The coefficients coef = {a_i} should follow this example:
  //
  //                     0 0 0      1 0 0     
  //    f(x, y, z) = a0 x y z + a1 x y z 
  //
  //                     0 1 0      1 1 0
  //    	   + a2 x y z + a3 x y z 
  //
  //                     0 0 1      1 0 1
  //		   + a4 x y z + a5 x y z 
  //
  //                     0 1 1      1 1 1 
  //		   + a6 x y z + a7 x y z 
  //
  // To create a polynomial with mvo = 1, a3, a5, a6, and a7 should be 0,
  // but must still be specified.
  Polynomial(unsigned dim, const DblArray& coef);

  // General n-D polynomials, exponent combinations and coefficients specified.
  //  
  // To create an MVO=2 2D polynomial, create an integer  matrix will all exponent
  // combinations, like
  //
  // expComb = [ 0 0 0 1 1 2 ]
  //           [ 0 1 2 0 1 0 ]
  //
  // and the corresponding coefficients array:
  //
  // coef    = [ . . . . . . ]
  //
  Polynomial(const IntMat& expComb, const DblArray& coef);

  // 1D fitted polynomials
  Polynomial(unsigned mvo, const FloatArray& x, const DblArray& f);
  Polynomial(unsigned mvo, const IntArray& x, const DblArray& f);

  // 2D fitted polynomials
  Polynomial(unsigned mvo, const FloatArray& x, const FloatArray& y, const DblArray& f);
  Polynomial(unsigned mvo, const IntArray& x, const IntArray& y, const DblArray& f);

  // n-D fitted polynomials
  Polynomial(unsigned mvo, const FlMat& x, const DblArray& f);
  Polynomial(unsigned mvo, const IntMat& x, const DblArray& f);

// Evaluation functions
  double operator () (float x) const;
  double operator () (float x, float y) const;
  double operator () (float x, float y, float z) const;
  double operator () (const FloatArray& x) const;

private:
  void _allExpComb(unsigned dimension, unsigned maxOrder);
  void _pruneExpComb(unsigned mvo);
  void _fit(const FlMat&, const DblArray&);
};

#endif
