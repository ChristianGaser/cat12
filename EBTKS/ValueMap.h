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
$RCSfile: ValueMap.h,v $
$Revision: 1.3 $
$Author: stever $
$Date: 2003/11/17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _VALUE_MAP_H
#define _VALUE_MAP_H

#include <iostream>		/* (bert) changed from iostream.h */
#include "SimpleArray.h"

/*************************
 * Abstract ValueMap class
 *************************/

class ValueMap {
public:
  // Invert the value map
  virtual ValueMap& inv() = 0;

  // Map concatenation
  virtual ValueMap& concat(const ValueMap& map) = 0;
  virtual ValueMap& operator () (const ValueMap& map) { return concat(map); }

  // Evaluate map
  virtual double operator () (double sourceValue) const = 0;
  // Reverse evaluate map
  virtual double reverse(double destValue) const = 0;

  // I/O
  virtual std::ostream& print(std::ostream&) const = 0;
};

inline std::ostream& operator<< (std::ostream& os, const ValueMap& map) 
{ 
  return map.print(os);
}

/******************
 * Linear map class
 ******************/

class LinearMap : public ValueMap {
private:
  double _factor;
  double _offset;

public:
  LinearMap(double factor = 1.0, double offset = 0.0) { 
    _factor = factor; _offset = offset;
  }
  LinearMap(const LinearMap& map) { 
    _factor = map._factor; _offset = map._offset;
  }
  LinearMap(double sourceMin, double sourceMax, double destMin, double destMax) {
    _factor = (destMax - destMin)/(sourceMax - sourceMin);
    _offset = destMin - _factor*sourceMin;
  }
  virtual ~LinearMap() { _factor = 1.0; _offset = 0.0; }

  double factor() const { return _factor; }
  double offset() const { return _offset; }

  double& factor() { return _factor; }
  double& offset() { return _offset; }

  double factor(double f) { _factor = f; return _factor; }
  double offset(double o) { _offset = o; return _offset; }

  LinearMap& operator = (const LinearMap& map) { 
    _factor = map._factor; _offset = map._offset; return *this; }

  ValueMap& inv() {
      _factor = 1.0/_factor; 
      _offset = -_offset; 
      return *this; 
  }

  ValueMap& concat(const ValueMap& map) {
      std::cerr << "LinearMap::concat() called but not implemented" << std::endl;
    return *this;
  }

  ValueMap& concat(const LinearMap& map) {
    _offset += _factor*map._offset; 
    _factor *= map._factor;
    return *this;
  }

  ValueMap& operator () (const ValueMap& map) { return concat(map); }
  double operator () (double sourceValue) const { return _offset+_factor*sourceValue;}

  LinearMap& operator () (double factor, double offset) {
    _factor = factor; 
    _offset = offset; 
    return *this; 
  }
  LinearMap& operator () (double sourceMin, double sourceMax, double destMin,
			  double destMax) {
    _factor = (destMax - destMin)/(sourceMax - sourceMin);
    _offset = destMin - _factor*sourceMin;
    return *this; 
  }

  
  double reverse(double destValue) const        { return (destValue-_offset)/_factor;}
  
  std::ostream& print(std::ostream& os) const {
    os << "(" << _factor << ", " << _offset << ")"; 
    return os;
  }

};

/********************
 * Lookup table class
 ********************/

template <class Type>
class LUT : public ValueMap {
private:
  SimpleArray<Type> _source;
  SimpleArray<Type> _dest;

public:
  LUT(unsigned length = 0); // Allocated size; actual size is zero
  LUT(const SimpleArray<Type>& source, const SimpleArray<Type>& dest);
  LUT(const LUT& map);
  virtual ~LUT();

  LUT& add(Type source, Type dest); // Add one map entry

  LUT& operator = (const LUT& map);

  ValueMap& inv();
  ValueMap& concat(const ValueMap& map);
  ValueMap& concat(const LUT& map);

  ValueMap& operator () (const ValueMap& map) { return concat(map); }
  double operator () (double sourceValue) const;
  double reverse(double destValue) const;
  
  std::ostream& print(std::ostream&) const;

private:
  void _sort();
};

/**************************
 * Mapping of other objects
 **************************/

// Map all elements through a single map
template <class Type>
SimpleArray<Type>& map(SimpleArray<Type>& array, const ValueMap&);
// Removed because of DCC resolving problems
//template <class Type>
//SimpleArray<Type>  map(const SimpleArray<Type>& array, const ValueMap& valueMap) {
//  return mapConst(array, valueMap); }
template <class Type>
SimpleArray<Type>  mapConst(const SimpleArray<Type>& array, const ValueMap& valueMap);

// Map each element through a linear map
template <class Type>
SimpleArray<Type>& map(SimpleArray<Type>& array, const Array<LinearMap>&);
// Removed because of DCC resolving problems
//template <class Type>
//SimpleArray<Type>  map(const SimpleArray<Type>& array, const Array<LinearMap>& maps) {
//  return mapConst(array, maps); }
template <class Type>
SimpleArray<Type>  mapConst(const SimpleArray<Type>& array,const Array<LinearMap>& maps);

#endif

