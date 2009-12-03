/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, John G. Sled, 
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
$RCSfile: TBSpline.h,v $
$Revision: 1.3 $
$Author: bert $
$Date: 2003/04/16 18:04:03 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : TBSpline.h,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tensor cubic B-splines in N dimensions with limited bending energy
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : $Log: TBSpline.h,v $
@MODIFIED   : Revision 1.3  2003/04/16 18:04:03  bert
@MODIFIED   : Removed meaningless const qualifiers on operator[] to avoid compiler warnings
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/20 21:42:44  jason
@MODIFIED   : Now compiles with gcc 3
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  2001/11/09 16:37:25  jason
@MODIFIED   : First, non-compiling import
@MODIFIED   :
@MODIFIED   : Revision 1.4  1997/02/06 15:56:29  jgsled
@MODIFIED   : Optimized using lookup tables.
@MODIFIED   :
 * Revision 1.3  1996/12/10  04:48:46  jgsled
 * intermediate version -- doesn't work
 *
 * Revision 1.2  1996/12/09  15:50:54  jgsled
 * modified to save memory during spline evalutation
 *
 * Revision 1.1  1996/08/23  20:02:40  jgsled
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/29  16:04:07  jgsled
 * last version before modification for inclusion in the AI source tree
 *
 * Revision 1.3  1996/07/09  15:14:00  jgsled
 * Working version
 *
 * Revision 1.2  1996/04/23  13:39:13  jgsled
 * Optimized addDataPoint function, modified to be compatible with changes to
 * Spline.h
 *
 * Revision 1.1  1996/04/21  17:21:13  jgsled
 * Tensor B Spline class with limited bending energy
 * Initial version: results have been verified with matlab for one and two
 * dimensions
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid_tbspline_h[] = "$Header: /software/source/libraries/EBTKS/include/TBSpline.h,v 1.3 2003/04/16 18:04:03 bert Exp $";
#endif

#ifndef TBSPLINE_H
#define TBSPLINE_H

#include "SimpleArray.h"
#include "Matrix.h"
#include "MatrixSupport.h"
#include "Spline.h"

#define DEBUG_TBSPLINE

#ifndef DblArray
typedef SimpleArray<double> DblArray;
#endif
#ifndef FloatArray
typedef SimpleArray<float> FloatArray;
#endif
#ifndef IntArray
typedef SimpleArray<int> IntArray;
#endif
#ifndef DblMat
typedef Mat<double> DblMat;
#endif
#ifndef FlMat
typedef Mat<float> FlMat;
#endif


/*************************
 * Tensor cubic B-spline class
 *************************/

class TBSpline : public Spline {
protected:
  TBSpline *_thisTBSpline;  // used to circumvent constant restrictions
  Boolean  _haveData;
  DblMat   _domain;
  double   _distance;  // knot spacing
  double   _scale;     // normalizing factor
  DblMat   _knots;
  DblMat   _zero;      // location of knot zero
  IntArray _n;         // number of basis functions in each dimension
  int      _nProduct;  // product of elements of n
  int      _four;      // 4^nDimensions
  IntArray _smallN;    // { 4, 4, 4, 4 ...}  size of small tensor

protected:
  // working matrices for intermediate results
  DblMat   _terms;          // 4 by nDimensions
  IntArray _blockIndex;     // nDimensions
  DblArray _values;         // 4^nDimensions
  IntArray _locations;      // 4^nDimensions
  int *_dloc_i;
  int *_dloc_j;

public:
  static double _default_lambda;

  // Constructors/destructor
  TBSpline(const DblMat &domain, 
	                  // an n by 2 matrix specifying a region in R^n
	                  // upon which the spline is defined
	                  // eg. [ x1 x2; y1 y2 ] 
	   double distance,  // the distance between adjacent knots 
	   double lambda = _default_lambda, // 0 < lambda < inf 
                                             // lambda is a limit on the
                                             // bending energy of the spline
           int allocate_flag = TRUE);  // TRUE -> normal usage
                                       // FALSE -> cannot call add_data
  ~TBSpline(void)   { delete _tempPoint; };

  // functions to allow incremental addition of new data points
  virtual Boolean clearDataPoints();
  virtual Boolean addDataPoint(const float *point, double value); 
  virtual Boolean fit();

  // evaluation functions
  operator Boolean() const { return _fitted; }
  virtual double operator () (const float *point) const;
  double operator () (float x) const { return operator () (&x); }
  double operator () (float x, float y) const {
    _tempPoint[0] = x; _tempPoint[1] = y; return operator () (_tempPoint); }
  double operator () (float x, float y, float z) const {
    _tempPoint[0] = x; _tempPoint[1] = y; _tempPoint[2] = z; 
    return operator () (_tempPoint); }  

  virtual void saveState(char *filename) const;

#ifdef DEBUG_TBSPLINE
  void test_J(void);
  void test_solver(void);
#endif

protected:
  void bendingEnergyTensor(const IntArray &n, DblMat &J);
  DblMat bendingEnergy(int size, int order);
  inline double cube(double x) const { return x*x*x; };
  DblMat solveSymmetricSystem(DblMat &A, DblMat b, int *info);
};

// tensor index
class TIndex   
{
protected:
  int _flat;
  IntArray _index;
  IntArray _n;
  IntArray _step;
  int _nDimensions;
  void init(const IntArray &n);
public:
  TIndex(const IntArray &n);  // initial value is zero
  TIndex(const IntArray &index, const IntArray &n);  // initial value is index
  virtual int operator[] (unsigned i) { return _index[i]; }; // (bert) lose const
  int flat(void) { return _flat; };
  virtual void operator ++ (void);
  virtual void reset(void) { _index.clear(0); };
};


// index of small tensor within large tensor
class TSubIndex : public TIndex  
{
protected:
  IntArray _smallIndex;
  IntArray _smallN;
  IntArray _smallStep;
  void init(const IntArray &smallN);
public:
  TSubIndex(const TIndex &index, const IntArray &smallN);
  TSubIndex(const TIndex &index, const IntArray &smallIndex,
	    const IntArray &smallN);
  virtual int operator[] (unsigned i) { return _smallIndex[i]; }; // (bert) lose const
  virtual void operator ++ ();
  void reset(void) { _index.clear(0); _smallIndex.clear(0); };
};

#ifdef DEBUG_TBSPLINE
void test_index(void);
#endif


/*************************************
 * Tensor cubic B-spline volume class
 *************************************/
// A derived class of TBSpline optimized for voxelated volumes

#define VDIM  3

class TBSplineVolume : public TBSpline {
protected:
  double *_1dSpline[VDIM]; // 1D splines evaluated at voxel centers
  unsigned short *_offset[VDIM];   // offsets into _coef
  IntArray _sizes;  // dimensions of volume

public:
  TBSplineVolume(const double start[VDIM],  // center of (0,0,0)th voxel 
                 const double step[VDIM],   // voxel separations
                 const int count[VDIM],     // number of voxels along each axis
                 double distance,  // the distance between adjacent knots 
                 double lambda = _default_lambda, // 0 < lambda < inf 
                                             // lambda is a limit on the
                                             // bending energy of the spline
                 int allocate_flag = TRUE);  // TRUE -> normal usage
                                             // FALSE -> cannot call add_data
  TBSplineVolume(const DblMat &domain, 
	                  // an n by 2 matrix specifying a region in R^n
	                  // upon which the spline is defined
	                  // eg. [ x1 x2; y1 y2 ] 
                 const double start[VDIM],  // center of (0,0,0)th voxel 
                 const double step[VDIM],   // voxel separations
                 const int count[VDIM],     // number of voxels along each axis
                 double distance,  // the distance between adjacent knots 
                 double lambda = _default_lambda, // 0 < lambda < inf 
                                             // lambda is a limit on the
                                             // bending energy of the spline
                 int allocate_flag = TRUE);  // TRUE -> normal usage
                                             // FALSE -> cannot call add_data

  virtual Boolean addDataPoint(int x, int y, int z, double value);
  virtual void saveState(char *filename) const;
  // evaluation functions
  operator Boolean() const { return _fitted; }
  double operator () (int x, int y, int z) const;

protected:
  DblMat computeDomain(const double start[VDIM], const double step[VDIM],
                       const int count[VDIM]); 
  void createLookup(const double start[VDIM],  
                               // center of (0,0,0)th voxel 
                             const double step[VDIM],   // voxel separations
                             const int count[VDIM]);
};

#endif











