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
$RCSfile: CachedArray.h,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2004/12/08 17:02:34 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _CACHED_ARRAY_H
#define _CACHED_ARRAY_H

#include "SimpleArray.h"
#include <fstream>		/* (bert) changed from fstream.h */

/********************************************************************
 * CachedArray class
 *
 * This class is derived from SimpleArray, and reduces memory use
 * by caching to disk.
 ********************************************************************/

template <class Type> class CachedArray;
template <class Type> class CacheBlock;

/*
typedef CachedArray<char>          CachedBoolArray;
typedef CachedArray<int>           CachedIntArray;
typedef CachedArray<unsigned>      CachedUnsignedArray;
typedef CachedArray<double>        CachedDblArray;
*/

template <class Type>
class CachedArray : public SimpleArray<Type> {

  static const unsigned _DEFAULT_BLOCK_SIZE;
  static const unsigned _DEFAULT_N_BLOCKS;
  static      unsigned _rangeErrorCount;

  CachedArray<Type> *_self;       // Provided to circumvent the const mechanism
  CacheBlock<Type>  *_head;       // Pointer to first cache block in linked list
  CacheBlock<Type> **_blocks;     // Pointers to cache blocks
  unsigned           _blockSize;  // # elements in each cache block
  unsigned           _nBlocks;    // Current # of cache blocks
  unsigned           _maxNblocks; // Maximum # cache blocks (= length of block table)
  std::fstream       _s;          // The stream associated with the temporary file
  // Hitrate calculation members
  unsigned           _hits;       // Cache hit counter
  unsigned           _misses;     // Cache miss counter
  // Iterator members
  unsigned           _itBlock;
  Type              *_itBlockPtr;

public:
  // Manager functions
  CachedArray(unsigned sz = 0, unsigned nBlocks = _DEFAULT_N_BLOCKS,
	      unsigned blockSize = _DEFAULT_BLOCK_SIZE);
  CachedArray(const Type *init, unsigned sz, unsigned nBlocks = _DEFAULT_N_BLOCKS,
	      unsigned blockSize = _DEFAULT_BLOCK_SIZE);
  CachedArray(const CachedArray<Type>& array);
  CachedArray(const SimpleArray<Type>& array, unsigned nBlocks = _DEFAULT_N_BLOCKS,
	      unsigned blockSize = _DEFAULT_BLOCK_SIZE);
  virtual ~CachedArray();

  // Binary I/O functions
  std::ostream& saveBinary(std::ostream&, unsigned n = 0, unsigned start = 0) const;
  std::istream& loadBinary(std::istream&, unsigned n = 0, unsigned start = 0);

  // Matlab I/O functions
#ifdef HAVE_MATLAB
  virtual Boolean saveMatlab(const char *, const char *, const char *) const {
    this->_notImplementedError(); return FALSE; }
#endif

  // Access functions
  Type& getEl(unsigned i)              { return _block(i)->getEl(i % _blockSize); }
  const Type& getEl(unsigned i) const  { return getElConst(i); }
  const Type& getElConst(unsigned i) const { 
    return _block(i)->getElConst(i % _blockSize); }
  void  setEl(unsigned i, Type value)  { _block(i)->setEl(i % _blockSize, value); }

  // Iterator functions
  // These functions do not modify the block read/write counters, as
  // the access is (supposed to be) sequential

  // Resets iterator to element i
  // g++ doesn't appear to be able to resolve the const/non-const versions.
  // see note in CachedArray.c
  //void  resetIterator(unsigned i = 0);       
  void  resetIterator(unsigned i = 0) const;
  Type& current();
  const Type& current() const;

  // Prefix ascending iterators
  Type& operator ++();
  const Type&  operator ++() const;
  // Postfix ascending iterators
  Type& operator ++(int);
  const Type&  operator ++(int) const;

  // Prefix descending iterators
  Type& operator --();
  const Type&  operator --() const;
  // Postfix descending iterators
  Type& operator --(int);
  const Type&  operator --(int) const;

  // Copy operators
  CachedArray& operator = (const CachedArray<Type>& array) {
    Array<Type>::operator = (array); return *this; }
  CachedArray& copyFromCarray(const Type *init, unsigned sz);

  // Conversions
  //operator SimpleArray<Type> () const;

  // Set functions
  void     newSize(unsigned sz);     // Resize array
  void     resetHitRate() const { 
    CachedArray<Type> *self = (CachedArray<Type>*) this; self->_hits=self->_misses=0; }

  // Not implemented
  SimpleArray<Type> operator () (unsigned) const {
    this->_notImplementedError(); return SimpleArray<Type>(0); }
  SimpleArray<Type> operator () (unsigned, unsigned) const {
    this->_notImplementedError(); return SimpleArray<Type>(0); }
  SimpleArray<Type> operator () (const BoolArray&) const {
    this->_notImplementedError(); return SimpleArray<Type>(0); }
  SimpleArray<Type> operator () (const UnsignedArray&) const {
    this->_notImplementedError(); return SimpleArray<Type>(0); }
  SimpleArray<Type> common(const SimpleArray<Type>&) const {
    this->_notImplementedError(); return SimpleArray<Type>(0); }

  Type& first()                    { return getEl(0); }
  Type& last()                     { return getEl(this->_size - 1); }
  const Type    *contents() const  { this->_notImplementedError(); return 0; }
  Type          *contents()        { this->_notImplementedError(); return 0; }
  operator const Type *() const    { this->_notImplementedError(); return 0; }
  operator       Type *()          { this->_notImplementedError(); return 0; }

  // Get functions
  double hitRate() const;

  // Quicksort functions; only straight ascending qsort is implemented
  void qsort() { qsortAscending(); }
  void qsort(int (*) (const void *, const void *)) { this->_notImplementedError(); }
  void qsortAscending();
  void qsortDescending() { this->_notImplementedError(); }
  UnsignedArray qsortIndexAscending() const { 
    this->_notImplementedError(); return UnsignedArray(0); }
  UnsignedArray qsortIndexDescending() const {
    this->_notImplementedError(); return UnsignedArray(0); }

  // Median functions. The array is first reduced, using a histogram approach,
  // to a single-block array, for which the algorithm described
  // in Sections 8.1, 8.3, and 10.2 in Cormen, Leierson, and Rivest,
  // "Introduction to Algorithms" (aka the Big White Book) is used.
  // The low median is returned in case the array has an even # elements.
  Type median() const;
  Type medianVolatile();

  CachedArray sample(unsigned maxN) const;
  CachedArray applyElementWise(Type (*function) (Type)) const;
  CachedArray map(const ValueMap& map) const;

  CachedArray ln() const;                    // Natural logarithm of all elements
  CachedArray log() const;                   // Base-10 logarithm of all elements
  CachedArray exp() const;                   // e^(all elements)
  CachedArray exp10() const;                 // 10^(all elements)

private:
  void _initialize(unsigned sz, unsigned nBlocks = _DEFAULT_N_BLOCKS,
		   unsigned blockSize = _DEFAULT_BLOCK_SIZE);
  void              _openStream();
  void              _destroy();
  CacheBlock<Type> *_block(unsigned i) const {
    if (i >= this->_size)
      this->_rangeError(i);
    CacheBlock<Type> *block = _blocks[i / _blockSize];
    return (block) ? block : _read(i / _blockSize);
  }
  void              _flush() const;              // Flush cache to disk
  CacheBlock<Type> *_read(unsigned block) const; // Reads block i
  void              _write() const;              // Write current block
  void              _revertIterator() const {    // Revert iterator (load correct block)
    _self->_itBlockPtr = _read(_self->_itBlock)->_contents; }
    
  // Median and quicksort support functions
  void _qsort(int p, int r);
  Type _histMedian(unsigned nBelow = 0, unsigned nAbove = 0);
  Type _randomizedSelect(int p, int r, int i);
  int  _randomizedPartition(int p, int r);
  int  _partition(int p, int r);
};

template <class Type>
CachedArray<Type> operator ^ (double base, const CachedArray<Type>& array);

/********************************************************************
 * CacheBlock class
 *
 * Support class for CachedArray. 
 ********************************************************************/

template <class Type>
class CacheBlock : private SimpleArray<Type> {
  friend class CachedArray<Type>;

  CacheBlock<Type> *_self;

  CacheBlock *_next;
  unsigned    _nBytes;
  Boolean     _changed;
  unsigned    _ID;
  unsigned    _nRead;
  unsigned    _nWrite;

public:  
  CacheBlock(unsigned ID, unsigned sz);
  ~CacheBlock();
  
  Type&       getEl(unsigned i) { _changed = TRUE; _nWrite++; return this->_contents[i]; }
  const Type& getEl(unsigned i) const { _self->_nRead++; return this->_contents[i]; }
  const Type& getElConst(unsigned i) const  { _self->_nRead++; return this->_contents[i]; }
  void        setEl(unsigned i, Type value) { 
    _changed = TRUE; _nWrite++; this->_contents[i]=value; }

  //operator SimpleArray<Type> () const { return SimpleArray<Type>(_contents, _size); }

  CacheBlock *addBlock(unsigned ID, unsigned sz);
  Boolean     read(std::fstream& s, unsigned ID); // Read new block; flushes existing one 
  Boolean     write(std::fstream&) const;
};

#endif
