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
$RCSfile: MString.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef MSTRING_H
#define MSTRING_H

#include <string.h>
#include "SimpleArray.h"
#include "MTypes.h"

class MStringIterator;

class MString : private SimpleArray<char> {

  static unsigned _MAX_LOAD_LENGTH;

friend class MStringIterator;

public:
// Constructors/destructor
  MString(int length = 0);       // Empty, but allocated for <length> chars
  MString(const char *);         // Initialized with string
  MString(const MString&);       // Copy
  MString(char c, unsigned n);   // Initialized with n characters <c>
  ~MString();

// Get functions
  Boolean  isEmpty() const { return (_size <= 1) || !strlen(_contents); }
  unsigned length() const  { return MIN(_size - 1, strlen(_contents)); }
  char     firstChar() const   { return _contents[0]; }
  char     lastChar() const    { 
    unsigned l = length(); return (l == 0)? '\0' : _contents[length()-1]; }
  char       *string()        { return _contents; }   // Return contained char *
  const char *string() const  { return _contents; }   // Return contained char *
  Boolean  isInteger(int *value = 0) const;
  Boolean  isDouble(double *value = 0) const;
  int      indexOf(char c, int dir = 1) const { 
    return SimpleArray<char>::indexOf(c, dir); }
  int      indexOf(char c, int dir, unsigned start) const { 
    return SimpleArray<char>::indexOf(c, dir, start); }
  Boolean  contains(const MString& subString) const;
  Boolean  contains(const char *charPtr) const;
  Boolean  contains(char c) const { return SimpleArray<char>::contains(c); }
  Boolean  isPartOf(const MString& string) const { return string.contains(*this); }

  bool operator ! () const     { return (_size <= 1) || !strlen(_contents); }
  operator void *() const { return (void *) ((_size > 1) && strlen(_contents)); }
  operator char *()       { return ((_size > 1) && strlen(_contents)) ? _contents : 0; }
  operator const char *() const { return ((_size>1) && strlen(_contents)) ? _contents:0;}

// Conversion operators
  operator unsigned() const;
  operator int() const;
  operator float() const;
  operator double() const;

// Padding functions
  MString& pad(char c, int n=1);             // Pad left (n<0) or right (n>0) with <c>
  MString& pad(const MString& str, int n=1); // Pad left (n<0) or right (n>0) with <str>

// Operators
  char& operator [] (unsigned index) {
    return Array<char>::operator [] (index); }
  char  operator [] (unsigned index) const {
    return Array<char>::operator [] (index); }
  Boolean operator == (const MString& string) const {
    return (strcmp(this->string(), string.string()) == 0); }
  Boolean operator == (const char *charPtr) const {
    return (strcmp(this->string(), charPtr) == 0); }
  Boolean operator != (const MString& string) const {
    return (strcmp(this->string(), string.string()) != 0); }
  Boolean operator != (const char *charPtr) const {
    return (strcmp(this->string(), charPtr) != 0); }
  Boolean operator < (const MString& string) const;
  Boolean operator > (const MString& string) const;

  MString& operator =  (const char *);         // Copy
  MString& operator =  (const MString&);
  MString& operator =  (char);                 // Convert values
  MString& operator =  (int);
  MString& operator =  (double);
  MString& operator += (const char *);         // Concatenation
  MString& operator += (const MString&);
  MString& operator += (char);                 // Append single char
  MString& operator += (int);                  // Append values
  MString& operator += (double);               
  MString  operator +  (const char *) const;   // Nondestructive concatenation
  MString  operator +  (const MString&) const;
  MString  operator +  (char character) const {// Concatenate a single character
    MString result(*this); return result += character; }
  MString  operator +  (int value) const {     // Concatenate a value
    MString result(*this); return result += value; }
  MString  operator +  (double value) const {     // Concatenate a value
    MString result(*this); return result += value; }
  MString  operator () (unsigned) const;       // Copy part of string (to end)
  MString  operator () (unsigned, unsigned) const; 
                                               // Copy part of string (fixed length)
// Special functions
  MString& applyTemplate(const MString&, const char *separator = ".");
  MString  chop(unsigned n = 1); // Chop off and return the last n characters

friend MString operator + (const char *charPtr, const MString& mString) {
  return MString(charPtr) + mString; }
friend MString operator + (int value, const MString& mString) {
  MString result; result += value; return result += mString; }

friend std::ostream& operator << (std::ostream& os, const MString& mString);
friend std::istream& operator >> (std::istream& is, MString& mString);
};

/***********************
 * MStringIterator class
 ***********************/

class MStringIterator {
  const MString *_mString;
  unsigned       _index;
  const MString  _separator;

public:
  MStringIterator(const MString& mString, unsigned start = 0);
  MStringIterator(const MString& mString, const MString& separator, unsigned start = 0);
  ~MStringIterator();

  unsigned index() { return _index; }

  MString first();
  MString start(unsigned index);
  char    character();

  MString operator ++();    // Prefix
  MString operator ++(int) { return operator ++ (); } // Postfix == prefix

private:
  MString _nextToken();
};

#endif

