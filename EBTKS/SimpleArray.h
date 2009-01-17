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
$RCSfile: SimpleArray.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef SIMPLE_ARRAY_H
#define SIMPLE_ARRAY_H

#include "trivials.h"
#include "Array.h"

/********************************************************************
 * SimpleArray class
 *
 * This class is derived from Array. It extends the functionality of
 * Array with various mathematical operations; as a result, these 
 * operations must be defined for the objects stored in the array.
 * SimpleArray is therefore most suitable for simple types.
 ********************************************************************/

class ValueMap;
template <class Type> class SimpleArray;

typedef SimpleArray<char>          BoolArray;
typedef SimpleArray<int>           IntArray;
typedef SimpleArray<unsigned>      UnsignedArray;
typedef SimpleArray<double>        DblArray;

template <class Type>
struct IndexStruct {
  Type     element;
  unsigned index;

  static int compareAscending(const void *struct1, const void *struct2);
  static int compareDescending(const void *struct1, const void *struct2);
};

template <class Type>
class SimpleArray : public Array<Type> {
public:
  static int     compareAscending(const void *, const void *);
  static int     compareDescending(const void *, const void *);
  static unsigned     _rangeErrorCount;

// Manager functions
  SimpleArray(unsigned sz = DEFAULT_SIZE) : Array<Type>(sz) {}
  SimpleArray(Type value, unsigned sz) : Array<Type>(value, sz) {}
  SimpleArray(const Type *init, unsigned nElements) : Array<Type>(init, nElements) {}
  SimpleArray(const Array<Type>& array) : Array<Type>(array) {}
  SimpleArray(const SimpleArray<Type>& array) : Array<Type>(array) {}
  SimpleArray(Type minVal, double step, Type maxVal);
  virtual ~SimpleArray() {}

  // Binary I/O functions
  std::ostream&         write(std::ostream& os) const      { saveBinary(os); return os; }
  virtual std::ostream& saveBinary(std::ostream&, unsigned n = 0, unsigned start = 0) const;
  std::istream&         read(std::istream& is)             { loadBinary(is); return is; }
  virtual std::istream& loadBinary(std::istream&, unsigned n = 0, unsigned start = 0);

  // ASCII I/O functions
  std::ostream&         print(std::ostream& os) const     { saveAscii(os); return os; }
  virtual std::ostream& saveAscii(std::ostream&, unsigned n = 0, unsigned start = 0) const;
  std::istream&         scan(std::istream& is)            { loadAscii(is); return is; }
  virtual std::istream& loadAscii(std::istream&, unsigned n = 0, unsigned start = 0);

  // Matlab I/O functions
#ifdef HAVE_MATLAB
  virtual Boolean saveMatlab(const char *fileName, const char *varName = "A", 
			     const char *option = "u") const;
#endif

// Get functions
  SimpleArray operator () (unsigned nElements) const;   // Sub-array (first nElements)
  SimpleArray operator () (unsigned start, unsigned end) const; // Sub-array
  SimpleArray operator () (const BoolArray& boolArray) const;
  // Return an array containing only the elements for which boolArray is TRUE
  SimpleArray operator () (const UnsignedArray& indices) const;
  // Return an array containing only the elements with indices listed in <indices>

  Boolean  contains(Type value) const;
  Boolean  contains(Type value, unsigned start, unsigned end) const;
  Boolean  containsOnly(Type value) const;
  Boolean  containsOnly(Type value, unsigned start, unsigned end) const;
  // Total number of occurrences of object
  unsigned occurrencesOf(Type value) const {
    return occurrencesOf(value, 0, this->_size - 1); }
  // Total number of occurrences of object present in a range
  unsigned occurrencesOf(Type, unsigned start, unsigned end) const; 
  // Index of first (dir > 0) or last (dir < 0) occurrence of value
  int    indexOf(Type, int dir, unsigned start) const;
  int    indexOf(Type value, int dir = 1) const { 
    return indexOf(value, dir, (dir > 0) ? 0 : this->_size - 1); }
  int    indexOf(Type value, unsigned start) const { 
    return indexOf(value, 1, start); }

  void  removeAll(Type value); 
  // Remove all elements equal to <value>
  void  removeAllIn(Type floor, Type ceil, unsigned *N = 0);
  // Remove all elements in [floor, ceil]
  void  removeAllNotIn(Type floor, Type ceil, unsigned *nBelow=0, unsigned *nAbove=0);
  // Remove all elements not in [floor, ceil]
  UnsignedArray indicesOf(Type value) const;
  // Return an array containing the (unique) elements common with <array>
  SimpleArray common(const SimpleArray<Type>& array) const; 

  // Removes all infinite and NaN values
  SimpleArray& prune(); 
  // Fills array with uniformly distributed numbers
  SimpleArray& randuniform(double min = 0, double max = 1);
  // Fills array with normally distributed numbers
  SimpleArray& randnormal(double mean = 0, double std = 1);

  // Quicksort all elements
  virtual void qsort() { qsortAscending(); }
  virtual void qsort(int (*compare) (const void *, const void *)) { 
    ::qsort(this->_contents, this->_size, sizeof(Type), compare); }
  virtual void qsortAscending() { 
    ::qsort(this->_contents, this->_size, sizeof(Type), compareAscending); }
  virtual void qsortDescending() {
    ::qsort(this->_contents, this->_size, sizeof(Type), compareDescending); }
  virtual SimpleArray<unsigned> qsortIndexAscending() const;
  virtual SimpleArray<unsigned> qsortIndexDescending() const;

  Type   min(unsigned *index = 0) const;
  Type   max(unsigned *index = 0) const;
  void   extrema(Type *min, Type *max) const;
  Type   range(unsigned *minIndex = 0, unsigned *maxIndex = 0) const;
  double sum() const;
  double sum2() const;
  double mean() const { return sum()/double(this->_size); }
  double prod() const;
  double prod2() const;
  double var() const;
  double std() const { return ::sqrt(var()); }

  // Median functions. median() and medianVolatile() use the algorithm described
  // in Sections 8.1, 8.3, and 10.2 in Cormen, Leierson, and Rivest,
  // "Introduction to Algorithms" (aka the Big White Book).
  // The low median is returned in case the array has an even # elements.
  virtual Type median() const;
  virtual Type medianVolatile(); // Does not create a copy of the array, but changes it

  virtual Type mode(const Type binWidth) const;

  // Cumulative sum/prod2uct
  DblArray cumSum() const;
  DblArray cumProd() const;

  // Impose ceiling/floor to values
  void ceil(Type ceil);
  void floor(Type floor);

  // Boolean operations
  // Elementwise inequality/equality
  Boolean operator != (const SimpleArray&) const; 
  Boolean operator == (const SimpleArray& array) const { 
    return !((*this) != array); }

  BoolArray operator && (const SimpleArray<Type>& array) const;
  BoolArray operator || (const SimpleArray<Type>& array) const;

  BoolArray operator == (double value) const;
  BoolArray operator != (double value) const;
  BoolArray operator >= (double value) const;
  BoolArray operator >  (double value) const;
  BoolArray operator <= (double value) const;
  BoolArray operator <  (double value) const;

  BoolArray operator >= (const SimpleArray<Type>&) const;
  BoolArray operator >  (const SimpleArray<Type>&) const;
  BoolArray operator <= (const SimpleArray<Type>&) const;
  BoolArray operator <  (const SimpleArray<Type>&) const;

// Arithmetic operations
  SimpleArray  operator - () const {
    SimpleArray<Type> result(Type(0), this->_size); return result -= *this; }

  SimpleArray& operator += (Type);         
  SimpleArray  operator +  (Type value) const { 
    SimpleArray<Type> result(*this); return result += value; }
  SimpleArray& operator += (const SimpleArray&);  
  SimpleArray  operator +  (const SimpleArray& array) const { 
    SimpleArray<Type> result(*this); return result += array; }

  SimpleArray& operator -= (Type);
  SimpleArray  operator -  (Type value) const { 
    SimpleArray<Type> result(*this); return result -= value; }
  SimpleArray& operator -= (const SimpleArray&);
  SimpleArray  operator -  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result -= array; }

  SimpleArray& operator *= (double);
  SimpleArray  operator *  (double value) const { 
    SimpleArray<Type> result(*this); return result *= value; }
  SimpleArray& operator *= (const SimpleArray&);
  SimpleArray  operator *  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result *= array; }

  SimpleArray& operator /= (double value) {
    return (*this) *= 1.0/value; }
  SimpleArray  operator /  (double value) const { 
    SimpleArray<Type> result(*this); return result /= value; }
  SimpleArray& operator /= (const SimpleArray&);
  SimpleArray  operator /  (const SimpleArray& array) const {
    SimpleArray<Type> result(*this); return result /= array; }

  SimpleArray abs() const;                // Absolute value
  SimpleArray round(unsigned n = 0) const;// Round to n decimal places
  SimpleArray sqr() const;                   // Square of all elements
  SimpleArray sqrt() const;                  // Square root of all elements
  SimpleArray operator ^ (int power) const;  // Raise all elements to <power>
    // Note: all elements are cast to int first. Should be generalized!!
  SimpleArray ln() const;                    // Natural logarithm of all elements
  SimpleArray log() const;                   // Base-10 logarithm of all elements
  SimpleArray exp() const;                   // e^(all elements)
  SimpleArray exp10() const;                 // 10^(all elements)

  SimpleArray sample(unsigned maxN) const;
  SimpleArray applyElementWise(Type (*function) (Type)) const;
  SimpleArray map(const ValueMap& map) const;

protected:
  // Median support functions
  Type _randomizedSelect(int p, int r, int i);
  int  _randomizedPartition(int p, int r);
  int  _partition(int p, int r);
};

// Various non-member template functions
//
template <class Type>
unsigned size(const SimpleArray<Type>& array) { return array.size(); }

template <class Type> 
Type min(const SimpleArray<Type>& array) { return array.min(); }

template <class Type>
Type max(const SimpleArray<Type>& array) { return array.max(); }

template <class Type>
Type range(const SimpleArray<Type>& array, unsigned *minIndex = 0, 
	   unsigned *maxIndex = 0)
{ 
  return array.range(minIndex, maxIndex);
}

template <class Type>
double sum(const SimpleArray<Type>& array) { return array.sum(); }

template <class Type>
double sum2(const SimpleArray<Type>& array) { return array.sum2(); }

template <class Type>
double mean(const SimpleArray<Type>& array) { return array.mean(); }

template <class Type>
double prod(const SimpleArray<Type>& array) { return array.prod(); }

template <class Type>
double prod2(const SimpleArray<Type>& array) { return array.prod2(); }

template <class Type>
double var(const SimpleArray<Type>& array) { return array.var(); }

template <class Type>
double stdev(const SimpleArray<Type>& array) { return array.std(); }

template <class Type>
Type median(const SimpleArray<Type>& array) { return array.median(); }

template <class Type>
Type medianVolatile(SimpleArray<Type>& array) { return array.medianVolatile(); }

template <class Type>
Type mode(const SimpleArray<Type>& array, Type binWidth) { return array.mode(binWidth); }

template <class Type>
SimpleArray<Type> abs(const SimpleArray<Type>& array) { return array.abs(); }

template <class Type>
SimpleArray<Type> sqr(const SimpleArray<Type>& array) { return array.sqr(); }

template <class Type>
SimpleArray<Type> sqrt(const SimpleArray<Type>& array) { return array.sqrt(); }

template <class Type>
SimpleArray<Type> ln(const SimpleArray<Type>& array) { return array.ln(); }

template <class Type>
SimpleArray<Type> log(const SimpleArray<Type>& array) { return array.log(); }

template <class Type>
SimpleArray<Type> max(const SimpleArray<Type>& array1, 
		      const SimpleArray<Type>& array2);

template <class Type>
SimpleArray<Type> min(const SimpleArray<Type>& array1, 
		      const SimpleArray<Type>& array2);

template <class Type>
SimpleArray<double> cumSum(const SimpleArray<Type>& array) 
{ 
  return array.cumSum();
}

template <class Type>
SimpleArray<double> cumProd(const SimpleArray<Type>& array) 
{ 
  return array.cumProd();
}

template <class Type>
SimpleArray<Type> exp(const SimpleArray<Type>& array) 
{
  return ::exp(1.0) ^ array;
}

template <class Type>
SimpleArray<Type> exp10(const SimpleArray<Type>& array) { return 10^array; }

template <class Type>
SimpleArray<Type> operator ^ (double base, const SimpleArray<Type>& array);

// I/O  
template <class T> std::ostream& operator << (std::ostream& os, const SimpleArray<T>& A) 
{
  return A.saveAscii(os); 
}

template <class T> std::istream& operator >> (std::istream& is, SimpleArray<T>& A) 
{
  return A.loadAscii(is); 
}

// Type conversions
template <class Type>
SimpleArray<int> asIntArray(const SimpleArray<Type>&);

template <class Type>
SimpleArray<float> asFloatArray(const SimpleArray<Type>&);

template <class Type>
SimpleArray<double> asDblArray(const SimpleArray<Type>&);

#ifdef USE_COMPMAT
template <> 
SimpleArray<dcomplex> SimpleArray<dcomplex>::round(unsigned n) const;

template <> 
SimpleArray<dcomplex>& SimpleArray<dcomplex>::prune();

template <> 
SimpleArray<dcomplex> operator ^ (double, const SimpleArray<dcomplex>&);
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template <> 
SimpleArray<fcomplex> SimpleArray<fcomplex>::round(unsigned n) const;

template <> 
SimpleArray<fcomplex>& SimpleArray<fcomplex>::prune();

template <> 
SimpleArray<fcomplex> operator ^ (double, const SimpleArray<fcomplex>&);
#endif /* USE_FCOMPMAT */

template <class Type>
SimpleArray<Type> round(const SimpleArray<Type>& array, unsigned n = 0) 
{ 
  return array.round(n);
}

template <class Type>
void prune(SimpleArray<Type>& array)
{ 
  array.prune();
}

#endif
