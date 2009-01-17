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
$RCSfile: MPoint.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef MPOINT_H
#define MPOINT_H

#include <math.h>
#include <stdlib.h>
#include <iostream>		/* (bert) changed from iostream.h */
#include "trivials.h"
#include "MTypes.h"
#include "Pool.h"

//class MRectangle;

/**************
 * MPoint class
 **************/

class MPoint {
  static Pool<MPoint> _pool;

public:
  short x;
  short y;

  MPoint(short x1 = 0, short y1 = 0)     { x = x1; y = y1; }
  MPoint(const MPoint& point)        { x = point.x; y = point.y; }

  void *operator new(size_t)         { return _pool.alloc(); }
  void  operator delete(void *point) { _pool.free(point); }

  MPoint& operator = (const MPoint& point) {
    x = point.x; y = point.y; return *this; }

  MPoint& operator () (short x1, short y1) {
    x = x1; y = y1; return *this; }

  double  distanceTo(const MPoint& point) const {
    return sqrt((double) SQR(x - point.x) + SQR(y - point.y)); }
  double  distanceTo(short x1, short y1) const {
    return sqrt((double) SQR(x - x1) + SQR(y - y1)); }
  double  sqrDistanceTo(const MPoint& point) const {
    return (double) SQR(x - point.x) + SQR(y - point.y); }
  double  sqrDistanceTo(short x1, short y1) const {
    return (double) SQR(x - x1) + SQR(y - y1); }

  MPoint  halfwayTo(const MPoint& point) const {
    return MPoint((x + point.x)/2, (y + point.y)/2); }

  Boolean connected4with(const MPoint& point) const {
    return (((x == point.x) && (abs(y - point.y) <= 1)) ||
	    ((y == point.y) && (abs(x - point.x) <= 1))); }

  Boolean connected8with(const MPoint& point) const {
    return ((abs(x - point.x) <= 1) && (abs(y - point.y) <= 1)); }

  MPoint  above() const { return MPoint(x, y - 1); }
  MPoint  below() const { return MPoint(x, y + 1); }
  MPoint  left() const  { return MPoint(x - 1, y); }
  MPoint  right() const { return MPoint(x + 1, y); }
  Boolean in(int x0, int y0, int x1, int y1) const {
    return ((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1)); }
  Boolean in(const MPoint& p0, const MPoint& p1) const {
    return in(p0.x, p0.y, p1.x, p1.y); }
  MPoint& magnify(int mag);
  MPoint  magnify(int mag) const { MPoint result(*this); return result.magnify(mag); }

  MPoint& operator += (const MPoint& point) { x += point.x; y += point.y; return *this; }
  MPoint& operator -= (const MPoint& point) { x -= point.x; y -= point.y; return *this; }
  MPoint  operator +  (const MPoint& point) const { 
    return MPoint(x + point.x, y + point.y); }
  MPoint  operator -  (const MPoint& point) const { 
    return MPoint(x - point.x, y - point.y); }
  MPoint  operator *  (const MPoint& point) const { 
    return MPoint(x * point.x, y * point.y); }
  Boolean operator < (const MPoint& point) const { 
    return((x < point.x) && (y < point.y)); }
  Boolean operator > (const MPoint& point) const { 
    return((x > point.x) && (y > point.y)); }
  Boolean operator == (const MPoint& point) const { 
    return((x == point.x) && (y == point.y)); }
  Boolean operator != (const MPoint& point) const { 
    return((x != point.x) || (y != point.y)); }
};

// Don't know why these can't be implemented as friends to the template class
// (See also Array.h)
std::istream& operator >> (std::istream& is, MPoint *&point);
std::istream& operator >> (std::istream& is, MPoint& point);
std::ostream& operator << (std::ostream& os, const MPoint& point);


/****************
 * MPoint3D class
 ****************/

class MPoint3D {
//  static Pool<MPoint3D> _pool; *** Crashes on MPoint3D. Couldn't figure out why ***

public:
  short x;
  short y;
  short z;

  MPoint3D(short x1 = 0, short y1 = 0, short z1 = 0) { x = x1; y = y1; z = z1; }
  MPoint3D(const MPoint3D& point)              { x = point.x; y = point.y; z = point.z; }

  //void *operator new(size_t)         { return _pool.alloc(); }
  //void  operator delete(void *point) { _pool.free(point); }

  MPoint3D& operator = (const MPoint3D& point) {
    x = point.x; y = point.y; z = point.z; return *this; }

  MPoint3D& operator () (short x1, short y1, short z1) {
    x = x1; y = y1; z = z1; return *this; }

  double  distanceTo(const MPoint3D& point) const {
    return sqrt((double) SQR(x - point.x) + SQR(y - point.y) + SQR(z - point.z)); }
  double  distanceTo(short x1, short y1, short z1) const {
    return sqrt((double) SQR(x - x1) + SQR(y - y1) + SQR(z - z1)); }
  double  sqrDistanceTo(const MPoint3D& point) const {
    return (double) SQR(x - point.x) + SQR(y - point.y) + SQR(z - point.z); }
  double  sqrDistanceTo(short x1, short y1, short z1) const {
    return (double) SQR(x - x1) + SQR(y - y1) + SQR(z - z1); }

  MPoint3D  halfwayTo(const MPoint3D& point) const {
    return MPoint3D((x + point.x)/2, (y + point.y)/2, (z + point.z)/2); }

  Boolean connected4with(const MPoint3D& point) const {
    return (((x == point.x) && (abs(y - point.y) <= 1)) ||
	    ((y == point.y) && (abs(x - point.x) <= 1)) ||
	    ((z == point.z) && (abs(z - point.z) <= 1))); }

  Boolean connected8with(const MPoint3D& point) const {
    return ((abs(x - point.x) <= 1) && 
	    (abs(y - point.y) <= 1) && 
	    (abs(z - point.z) <= 1)); }

  MPoint3D& operator += (const MPoint3D& point) { 
    x += point.x; y += point.y; z += point.z; return *this; }
  MPoint3D& operator -= (const MPoint3D& point) { 
    x -= point.x; y -= point.y; z -= point.z; return *this; }
  MPoint3D  operator +  (const MPoint3D& point) const { 
    return MPoint3D(x + point.x, y + point.y, z + point.z); }
  MPoint3D  operator -  (const MPoint3D& point) const { 
    return MPoint3D(x - point.x, y - point.y, z - point.z); }
  MPoint3D  operator *  (const MPoint3D& point) const { 
    return MPoint3D(x * point.x, y * point.y, z * point.z); }
  Boolean operator < (const MPoint3D& point) const { 
    return((x < point.x) && (y < point.y) && (z < point.z)); }
  Boolean operator > (const MPoint3D& point) const { 
    return((x > point.x) && (y > point.y) && (z > point.z)); }
  Boolean operator == (const MPoint3D& point) const { 
    return((x == point.x) && (y == point.y) && (z == point.z)); }
  Boolean operator != (const MPoint3D& point) const { 
    return((x != point.x) || (y != point.y) || (z != point.z)); }
};


/*******************
 * MWorldPoint class
 *******************/

class MWorldPoint {
  static Pool<MWorldPoint> _pool;

public:
  float x;
  float y;
  float z;

  void *operator new(size_t)         { return _pool.alloc(); }
  void  operator delete(void *point) { _pool.free(point); }

  MWorldPoint(float x1 = 0.0, float y1 = 0.0, float z1 = 0.0) { 
    x = x1; y = y1; z = z1; }
  MWorldPoint(const MWorldPoint& point) { 
    x = point.x; y = point.y; z = point.z; }

  MWorldPoint& operator = (const MWorldPoint& point) {
    x = point.x; y = point.y; z = point.z; return *this; }
};

#endif
