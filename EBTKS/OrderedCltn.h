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
$RCSfile: OrderedCltn.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
/******************************************************************
 *								  *
 *		           V-ATMS V. 1.0			  *
 *								  *
 *		    Sugato Bagchi, Serdar Uckun			  *
 *     Center for Intelligent Systems, Vanderbilt University	  *
 *								  *
 *          Sponsored by: Computer Systems Laboratory		  *
 *                   Nippon Steel Corporation			  *
 *                     Kanagawa 229, JAPAN			  *
 *								  *
 *  Use and distribution for purposes other than research and 	  *
 *  education is prohibited.					  *
 *                                                                *
 *  Modified 1993/12 by Alex Zijdenbos                            *
 *    - changed private class members to protected                *
 *    - added qsort()                                             *
 *  Modified 1994/03 by Alex Zijdenbos                            *
 *    - added default constructor and attach() to ocIterator      *
 *    - changed indexOf() to return -1 if Ptr is not found        *
 *    - added check for empty collection in remove()              *
 *  Modified 1994/04 by Alex Zijdenbos                            *
 *    - Changed remove(const void *) to return 1 on success and   *
 *      0 on failure                                              *
 *    - Changed destructor to virtual                             *
 *  Modified 1994/05 by Alex Zijdenbos                            *
 *    - Added addAllFirst(const OrderedCltn&)                     *
 *  Modified 1994/08 by Alex Zijdenbos                            *
 *    - Added postfix ++ and -- operators for ocIterator          *
 *                                                                *
 ******************************************************************/

#ifndef ORDCLTN_H
#define ORDCLTN_H

#include <stdlib.h>

const unsigned EXPANSION_INCREMENT = 512;
const unsigned DEF_OC_SIZE = 512;


class OrderedCltn {
protected:
   friend class ocIterator;
   void** Contents;		// points to an array of (void *)s
   unsigned uEndIndex;		// index to the last entry in array
   unsigned uCapacity;		// capacity of the array

   void no_mem_err () const;		// displays error message and quits
   void range_err (unsigned u) const;	// displays error msg and quits
   void not_impl_err() const;           // Operation not implemented
      
public:
  OrderedCltn (unsigned uCap = DEF_OC_SIZE);	// define
  OrderedCltn (const OrderedCltn& OC);		// copy
  virtual ~OrderedCltn ();

  virtual unsigned size () const      // number of entries in the collection
      { return uEndIndex; }
  unsigned capacity () const	       // total capacity (may be increased with resize())
      { return uCapacity; }
  void reSize (unsigned uNewCap);     // increase the size of the array
  int indexOf (const void* Ptr) const; 	// returns index of Ptr if found, -1 if not found
  unsigned occurrencesOf (const void* Ptr) const; // returns the total number of occurrences of an object
  unsigned add (const void* Newptr, unsigned uIndex); 	// add entry at specified position
  unsigned add (const void* Newptr) {	// add entry at the end of the collection
    return add (Newptr, uEndIndex); }
  
  unsigned addAfter (const void* Ptr,const void *Newptr); 	// add after Ptr (if it is found)
  unsigned addBefore (const void* Ptr, const void *Newptr); 	// add before Ptr (if found)
  const OrderedCltn& addAllFirst (const OrderedCltn&);
  const OrderedCltn& addAllLast (const OrderedCltn&);
  void* remove (unsigned uIndex); 	// remove entry from specified position
  void* remove ();                      // remove the last entry
  unsigned remove (const void* Ptr);	// remove the specified entry (if it exists)
  
  void removeAll ()		// remove all entries
      { uEndIndex = 0; }
  void* at (unsigned u) const	// get the pointer at specified index
      { return (*this)[u]; }
  void* first () const		// get the first pointer in the collection
      { return (uEndIndex ? (*this)[0] : 0); }
  void* last () const		// get the last pointer in the collection
      { return (uEndIndex ? (*this)[uEndIndex-1] : 0); }
  int isEmpty () const		// is the collection empty ?
      { return (uEndIndex ? 0 : 1); }

  void qsort(int (*compare) (const void *, const void *)) { // Quicksort all elements
    ::qsort(Contents, uEndIndex, sizeof(void *), compare); }

  void*& operator [](unsigned u) const {		 // return the entry at index `u'
    if (u >= uEndIndex) range_err (u); return Contents[u]; }
  OrderedCltn& operator =  (const OrderedCltn& OC);	  // assignment
  int          operator == (const OrderedCltn& OC) const; // equality
  int          operator != (const OrderedCltn& OC) const { return !(*this == OC); }
  OrderedCltn& operator &= (const OrderedCltn& OC);	  // destructive intersection
  const OrderedCltn& operator += (const OrderedCltn& OC) { return addAllLast(OC); }
};

class ocIterator {
  const OrderedCltn* pOC;		// collection to be iterated
  unsigned uCurIndex; 		// current iteration position

public:
  ocIterator()                      { pOC = 0; uCurIndex = 0; }
  ocIterator(const OrderedCltn& OC) { pOC = &OC; uCurIndex = 0; }

  void attach(const OrderedCltn& OC) { pOC = &OC; uCurIndex = 0; }
  
  void first()       // set current iteration position to the first entry
  { uCurIndex = 0; }
  void last()        // set current iteration position to the last entry
  { uCurIndex = pOC->uEndIndex - 1; }

  void* operator ++();      // return the next entry in collection
  void* operator ++(int);   // return the current entry in collection and increment
  void* operator --();      // return the previous entry in collection
  void* operator --(int);   // return the current entry in collection and decrement
};

#endif
