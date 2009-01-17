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
$RCSfile: Array.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
/*
 * Essentially the template Array class described in Lippman's C++ Primer
 *
 * 11/12/1993
 *
 * Alex Zijdenbos
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include "trivials.h"
#include "MTypes.h"


const int DEFAULT_SIZE   = 0;
const int SIZE_INCREMENT = 32;

const int BEFORE = -1;
const int AFTER  = 1;

/******************
 * Array base class
 ******************/

template <class Type>
class Array {
// data members
  Array<Type> *_self; // Provided to circumvent the const mechanism for the iterator

protected:
  static unsigned _arrayCtr;
  static Boolean  _debug;
  static unsigned _rangeErrorCount;

  unsigned  _size;
  unsigned  _maxSize;
  Type     *_contents;

  // Iterator member
  unsigned  _itIndex;

public:
  static unsigned debug(Boolean on = TRUE);

// Manager functions
  Array (unsigned sz = DEFAULT_SIZE);
  Array (const Type&, unsigned sz);
  Array (const Type *, unsigned);
  Array (const Array&);
  virtual ~Array ();

  // Access functions
  inline Type&       operator [] (unsigned i);
  inline Type&       operator [] (int i);
  inline const Type& operator [] (unsigned i) const;
  inline const Type& operator [] (int i) const;

  inline virtual Type& getEl(unsigned i);
  inline virtual const Type& getEl(unsigned i) const;
  inline virtual const Type& getElConst(unsigned i) const;
  inline virtual void  setEl(unsigned i, Type value);

  // Iterator functions
  //virtual void  resetIterator(unsigned i = 0);       // Resets iterator to element i
  virtual void  resetIterator(unsigned i = 0) const;
  inline virtual Type& current();
  inline virtual const Type& current() const;

  // Prefix ascending iterators
  inline virtual Type& operator ++();
  inline virtual const Type&  operator ++() const;
  // Postfix ascending iterators
  inline virtual Type& operator ++(int);
  inline virtual const Type& operator ++(int) const;

  // Prefix descending iterators
  inline virtual Type& operator --();
  inline virtual const Type& operator --() const;
  // Postfix descending iterators
  inline virtual Type& operator --(int);
  inline virtual const Type& operator --(int) const;

  Array  operator () (unsigned nElements) const;           // Sub-array (first nElements)
  Array  operator () (unsigned start, unsigned end) const; // Sub-array

  Boolean operator ! () const;
  
  virtual Type&    first();
  virtual Type&    last();
  virtual unsigned size() const;
  // Direct access to contents - volatile!!
  virtual const Type *contents() const;
  virtual Type *contents();
  virtual operator const Type *() const;
  virtual operator Type *();

// Other functions
  Array&  operator = (const Array&);        // Copy
  Array&  absorb(Array&);                   // Copy (source destroyed)
  Array&  operator () (const Type *, unsigned); // Copy from C array
  // Returns C array. Uses <array> if specified; otherwise allocates a new array.
  Type   *asCarray(Type *array = 0) const;  
  Array&  append(const Type value);        // Append a value at the end
  Array&  append(const Array&);             // Concatenate two arrays
  Array&  insert(const Type&, unsigned index = 0);  // Insert an element at <index>
  Array&  insert(const Array&, unsigned index = 0); // Insert an array at <index>
  Array&  replace(const Array&, unsigned index = 0);// Replace, starting at <index>
  Type    remove(unsigned index = 0);       // Remove the element at <index>
  Type    removeLast();                     // Remove the last element
  Array&  reorder(const Array<unsigned>& indices); // Reorder the array
  Array&  shuffle();                        // Randomly re-order the array
  virtual void clear(const Type& value);
  virtual void newSize(unsigned);    // Changes size of array; keeps contents if possible
  Array&  destroy();                 // Destroys array (sets size to zero, frees memory)

  // Rotate operators
  Array& operator << (unsigned);
  Array& operator >> (unsigned);

  Array   sample(unsigned maxN) const;
  Array   applyElementWise(Type (*function) (Type)) const;

// I/O
  virtual std::ostream& print(std::ostream& os) const;

protected:
  void _grow(unsigned amount);       // Allocate more space; do not change _size
  virtual void _rangeError(unsigned& index) const;
  virtual void _notImplementedError() const;
};

template <class Type>
unsigned size(const Array<Type>& array);

template <class Type> 
std::ostream& operator << (std::ostream& os, const Array<Type>& array);

template <class Type>
Type&
Array<Type>::operator [] (unsigned i) 
{ 
  return getEl(i); 
}

template <class Type>
Type&
Array<Type>::operator [] (int i) 
{ 
  return getEl((unsigned)i); 
}

template <class Type>
const Type& 
Array<Type>::operator [] (unsigned i) const 
{ 
  return getElConst(i); 
}

template <class Type>
const Type& 
Array<Type>::operator [] (int i) const 
{ 
  return getElConst(i); 
}

//
// Access functions
//

template <class Type>
Type& 
Array<Type>::getEl(unsigned i) 
{
  if (i >= _size)_rangeError(i); 
  return _contents[i]; 
}

template <class Type>
const Type& 
Array<Type>::getEl(unsigned i) const
{
  return getElConst(i); 
}

template <class Type>
const Type& 
Array<Type>::getElConst(unsigned i) const
{ 
  if (i >= _size) _rangeError(i); 
  return _contents[i]; 
}
 
template <class Type>
void
Array<Type>::setEl(unsigned i, Type value)  
{ 
  if (i >= _size) _rangeError(i); 
  _contents[i] = value; 
}

template <class Type>
Type& 
Array<Type>::current()
{ 
  return _contents[_itIndex]; 
}

template <class Type>
const Type& 
Array<Type>::current() const  
{ 
  return _contents[_itIndex]; 
}

// Prefix ascending iterators
template <class Type>
Type& 
Array<Type>::operator ++()          
{ 
  return _contents[++_itIndex]; 
}

template <class Type>
const Type&  
Array<Type>::operator ++() const 
{
  return _contents[++_self->_itIndex]; 
}

// Postfix ascending iterators
template <class Type>
Type& 
Array<Type>::operator ++(int)       
{ 
  return _contents[_itIndex++]; 
}

template <class Type>
const Type& 
Array<Type>::operator ++(int) const 
{ 
  return _contents[_self->_itIndex++]; 
}

// Prefix descending iterators
template <class Type>
Type& 
Array<Type>::operator --()
{ 
  return _contents[--_itIndex]; 
}

template <class Type>
const Type& 
Array<Type>::operator --() const 
{ 
  return _contents[--_self->_itIndex]; 
}

// Postfix descending iterators
template <class Type>
Type& 
Array<Type>::operator --(int)
{ 
  return _contents[_itIndex--]; 
}

template <class Type>
const Type& 
Array<Type>::operator --(int) const 
{ 
  return _contents[_self->_itIndex--]; 
}

#endif
