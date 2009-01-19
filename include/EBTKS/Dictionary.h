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
$RCSfile: Dictionary.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <iostream>		/* (bert) changed from iostream.h */
#include "MTypes.h"
#include "trivials.h"

/*************************************************************
 * Dictionary class, essentially an implementation of the Map
 * class in the 2nd edition of "The C++ Programming Language"
 * by Bjarne Stroustrup, pp. 284-...
 *
 * Alex Zijdenbos, 1994/03/15
 *************************************************************/

template <class Key, class Value> class Dictionary;
template <class Key, class Value> class DictIterator;

/*****************
 * DictEntry class
 *****************/

template <class Key, class Value>
class DictEntry {
public:
  const Key key;
  Value     value;

private:
  DictEntry *_next;
  DictEntry *_previous;

  DictEntry(const Key& k, const Value& v) : key(k), value(v) {}
  ~DictEntry() { delete _next; }

friend class Dictionary<Key, Value>;
friend class DictIterator<Key, Value>;
};

/******************
 * Dictionary class
 ******************/

template <class Key, class Value>
class Dictionary {
friend class DictIterator<Key, Value>;

  DictEntry<Key, Value> *_head;
  DictEntry<Key, Value> *_current;
  DictEntry<Key, Value> *_tail;

  Key   _defaultKey;
  Value _defaultValue;

  unsigned _size;

public:
// Constructors/destructor
  Dictionary() { _initialize(); }
  Dictionary(const Key& key, const Value& value) 
  : _defaultKey(key), _defaultValue(value) { _initialize(); }
  Dictionary(const Dictionary&);
  ~Dictionary() { delete _head; }

  Dictionary& operator = (const Dictionary&);

// Get functions
  unsigned   size() const { return _size; }
  Value     *at(const Key& key);
  Value&     operator [] (const Key& key);
  Value      operator [] (const Key& key) const;
  const Key& keyOf(const Value& value);
  Boolean    containsKey(const Key& key) const;
  Boolean    containsValue(const Value& value) const;

  const DictEntry<Key, Value> *first() const { return _head; }
  const DictEntry<Key, Value> *last() const  { return _tail; }

// Entry removing functions
  void clear() { delete _head; _initialize(); }
  void removeAll() { clear(); }
  void remove(const Key& key);

// Operators
  //Dictionary<Key, Value>& operator += (const Dictionary<Key, Value>& dict);
  // Add two dictionaries together; this implies that the += operator must be 
  // defined on Value.
  Boolean operator == (const Dictionary<Key, Value>& dict) const;

private:
  void _initialize() { _size = 0; _head = 0; _current = 0; _tail = 0; }
  void _notImplementedError();
};

template <class Key, class Value>
std::ostream& operator << (std::ostream& os, const Dictionary<Key, Value>& dict);

/***************************
 * Dictionary iterator class
 ***************************/

template <class Key, class Value>
class DictIterator {
  const Dictionary<Key, Value> *_dictionary;
  DictEntry<Key, Value>        *_dictEntry;

public:
  DictIterator() { _dictionary = 0; _dictEntry = 0; }

  DictIterator(const Dictionary<Key, Value>& dictionary) {
    _dictionary = &dictionary; _dictEntry = _dictionary->_head; }

  void attach(const Dictionary<Key, Value>& dictionary) {
    _dictionary = &dictionary; _dictEntry = _dictionary->_head; }

  DictEntry<Key, Value>* first() { return _dictEntry = _dictionary->_head; }
  DictEntry<Key, Value>* last()  { return _dictEntry = _dictionary->_tail; }
  
  operator DictEntry<Key, Value>* () { return _dictEntry; }

  DictEntry<Key, Value>* operator ++();    // Prefix
  DictEntry<Key, Value>* operator ++(int); // Postfix
  DictEntry<Key, Value>* operator --();    // Prefix
  DictEntry<Key, Value>* operator --(int); // Postfix
};

#endif
