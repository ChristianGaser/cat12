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
$RCSfile: Stack.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef STACK_H
#define STACK_H

#include "Array.h"

template <class Type>
class Stack : private Array<Type> {

  Stack(const Stack&); // Prevent copying a Stack
  void operator = (const Stack&);

public:
  Stack(unsigned size = SIZE_INCREMENT) : Array<Type>(size) { _size = 0; }
  
  void push(const Type& value) {
    if (_maxSize <= _size)
      _grow(SIZE_INCREMENT);
    _contents[_size++] = value;
  }
    
  Type pop() { return (_size)? _contents[--_size] : (Type) 0; }
  
  unsigned size() const { return _size; }
};

#endif
