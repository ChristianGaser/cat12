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
$RCSfile: Spline.h,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2003/04/16 17:59:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef SPLINE_3D_H
#define SPLINE_3D_H

#include "SimpleArray.h"
#include "Matrix.h"
#include "OrderedCltn.h"
#ifdef HAVE_INTERFACE
#include "MRegion.h"
#endif

#ifndef DblArray
typedef SimpleArray<double> DblArray;
#endif
#ifndef FloatArray
typedef SimpleArray<float> FloatArray;
#endif
#ifndef DblMat
typedef Mat<double> DblMat;
#endif
#ifndef FlMat
typedef Mat<float> FlMat;
#endif

typedef double (*RadialFunction)(double r2);

/*******************
 * Spline base class
 *******************/

class Spline {
protected:
  unsigned _nDimensions;
  unsigned _nKnots;
  FlMat    _knots;
  unsigned _nCoef;
  DblArray _coef;
  Boolean  _fitted;
  DblMat   _AtA, _AtF, _Aj, _Ajt;
  DblMat   _J;
  double   _lambda;
  float   *_tempPoint;

public:
  static Boolean verbose;

  static double rToTheThird(double r2) { return r2 * sqrt(r2); }
  static double rSquareLogR(double r2) { if (r2) r2 *= log(r2)/2; return r2; }
  static double r(double r2)           { return sqrt(r2); }
  
// Constructors/destructor
  Spline(unsigned nDimensions);
  Spline(const FloatArray& xKnots);
  Spline(const FloatArray& xKnots, const FloatArray& yKnots);
  Spline(const FloatArray& xKnots, const FloatArray& yKnots, const FloatArray& zKnots) ;
  Spline(const FlMat& knots);
#ifdef HAVE_INTERFACE
  Spline(const MRegion& knots);
#endif
  Spline(const OrderedCltn& knots); // OC of MPoint3D
  virtual ~Spline() { delete [] _tempPoint; }

  // Returns a pointer to an internal array where the application can
  // write the coordinates of data points, both for use in the fitting
  // (addDataPoint) as well as for evaluation (operator double()).
  float *point() { return _tempPoint; }

  void    lambda(double lambda) { _lambda = lambda; }
  Boolean J(const DblMat& J);

  // Fitting functions which allow data points to be specified individually
  virtual Boolean clearDataPoints(); 
  Boolean         addDataPoint(double value) {
    return addDataPoint(_tempPoint, value); }
  virtual Boolean addDataPoint(const float *point, double value) = 0;
  virtual Boolean fit();

  virtual Boolean fit(const FlMat& coord, const DblArray& values);
  virtual Boolean fit(const OrderedCltn& coord, const DblArray& values) {
    return fit(coord2FlMat(coord), values); }
#ifdef HAVE_INTERFACE
  virtual Boolean fit(const MRegion& coord, const DblArray& values) {
    return fit(coord2FlMat(coord), values); }
#endif

  int operator ! () const      { return !_fitted; }
  operator double () const { return operator () (_tempPoint); }

  virtual double operator () (const float *point = 0) const = 0;

  double operator () (const int *point) const {
    if (point)
      for (unsigned d = 0; d < _nDimensions; d++)
	_tempPoint[d] = float(point[d]);
    return operator () (_tempPoint);
  }
  double operator () (float x) const { return operator () (&x); }
  double operator () (float x, float y) const {
    _tempPoint[0] = x; _tempPoint[1] = y; return operator () (_tempPoint); }
  double operator () (float x, float y, float z) const {
    _tempPoint[0] = x; _tempPoint[1] = y; _tempPoint[2] = z; 
    return operator () (_tempPoint); }
#ifdef HAVE_INTERFACE
  double operator () (const MPoint& point) const {
    _tempPoint[0] = point.x; _tempPoint[1] = point.y; return operator () (_tempPoint); }
  double operator () (const MPoint3D& point) const {
    _tempPoint[0] = point.x; _tempPoint[1] = point.y; _tempPoint[2] = point.z;
    return operator () (_tempPoint); }
#endif

  virtual void saveState(char *filename = "matlab.mat") const;
  virtual DblArray getCoefficients(void) const { return _coef; }
  virtual Boolean  putCoefficients(const DblArray& coef) { 
    if (coef.size() != _nCoef) return FALSE;
    else { _coef = coef; _fitted = TRUE; return TRUE; } }

#ifdef HAVE_INTERFACE
  static FlMat coord2FlMat(const MRegion& coord);
#endif
  static FlMat coord2FlMat(const OrderedCltn& coord); // OC of MPoint3D!

protected:
  virtual void _initialize();

private:
  void _prune(); // Eliminate duplicate knots
};

/****************
 * B-spline class
 ****************/

/*********************
 * Radial spline class
 *********************/

class RSpline : public Spline {
protected:
  RadialFunction _radialFunction;
  double         _coordScale, _coordScale2; // Scale factor for coords (, squared)

public:
// Constructors/destructor
  RSpline(const FloatArray& xKnots, double scale = 1.0);
  RSpline(const FloatArray& xKnots, const FloatArray& yKnots, double scale = 1.0);
  RSpline(const FloatArray& xKnots, const FloatArray& yKnots, const FloatArray& zKnots,
	  double scale = 1.0);
  RSpline(const FlMat& knots, double scale = 1.0);
#ifdef HAVE_INTERFACE
  RSpline(const MRegion& knots, double scale = 1.0);
#endif
  RSpline(const OrderedCltn& knots, double scale = 1.0); // OC of MPoint3D
  virtual ~RSpline() {}

  void radialFunction(RadialFunction F) { _radialFunction = F; }

  virtual Boolean addDataPoint(const float *point, double value);

  virtual double operator () (const float *point = 0) const;

protected:
  double _r2(const float *p1, const float *p2) const {
    double r2 = 0.0;
    for (unsigned i = 0; i < _nDimensions; i++) { /* (bert) changed d to i */
      double d = double(*p1++) - *p2++;
      r2 += d*d;
    }
    return _coordScale2*r2;
  }
  virtual void _initialize();
};

/*************************
 * Thin-plate spline class
 *************************/

class TPSpline : public RSpline {
protected:
  DblMat _Psi;
  DblMat _Psi_r;
  DblMat _Ft; 

public:
// Constructors/destructor
  TPSpline(const FloatArray& xKnots, double scale = 1.0);
  TPSpline(const FloatArray& xKnots, const FloatArray& yKnots, double scale = 1.0);
  TPSpline(const FloatArray& xKnots, const FloatArray& yKnots, const FloatArray& zKnots,
	   double scale = 1.0);
  TPSpline(const FlMat& knots, double scale = 1.0);
#ifdef HAVE_INTERFACE
  TPSpline(const MRegion& knots, double scale = 1.0);
#endif
  TPSpline(const OrderedCltn& knots, double scale = 1.0); // OC of MPoint3D
  virtual ~TPSpline();

  virtual Boolean clearDataPoints();
  virtual Boolean addDataPoint(const float *point, double value);
  virtual Boolean fit();

  virtual double operator () (const float *point = 0) const;

protected:
  virtual void _initialize();
};

#endif
