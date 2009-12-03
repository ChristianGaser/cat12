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
$RCSfile: Pool.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:26 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef POOL_H
#define POOL_H

#include "OrderedCltn.h"

/**********************************************************
 *
 * Template Pool class, for optimizing the allocation
 * of small objects.
 *
 * IMPORTANT NOTE: any object allocated/freed using Pool,
 * must clean up thoroughly in its destructor since 
 * the space used for the object will be re-used (without
 * re-initialization).
 *
 * Alex Zijdenbos, 14 Sept 1995
 *
 **********************************************************/

struct _Link { 
  _Link *next; 
};

const unsigned _POOL_EXPANSION_SIZE = 512;
const unsigned _POOL_N_BLOCKS       = 50;

template <class Type>
class Pool {
  const unsigned _elementSize;
  const unsigned _expansionSize;
  _Link         *_head;
  OrderedCltn    _blocks;

  Pool(const Pool&);
  void operator = (const Pool&);

public:
  Pool(unsigned nElements = _POOL_EXPANSION_SIZE);
  ~Pool();

  void *alloc() {
    if (!_head) _grow();
    _Link *element = _head;
    _head = element->next;
    return (void *) element;
  }

  void free(void *x) {
    _Link *element = (_Link *) x;
    element->next = _head;
    _head = element;
  }

private:
  void _grow();
};

#endif
