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
$RCSfile: TrainingSet.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include "SimpleArray.h"
#include "OrderedCltn.h"

/***************
 * Example class
 ***************/

struct Example {
  unsigned label;
  DblArray input;
  DblArray target;

// Constructors/destructor
  Example(unsigned nInputs, unsigned nTargets); // Empty with supplied sizes
  Example(unsigned label, const DblArray& input, const DblArray& output);
  ~Example();

// Get functions
  double *inputPtr()  { return input.contents(); }
  double *targetPtr() { return target.contents(); }

    friend std::ostream& operator << (std::ostream&, Example&);
};

/*******************
 * TrainingSet class
 *******************/

class TrainingSet : protected OrderedCltn {
  unsigned _nInputs;
  unsigned _nTargets;
  double   _minTarget;
  double   _maxTarget;

public:
// Constructors/destructor
  TrainingSet();
  TrainingSet(unsigned nExamples, unsigned nInputs, unsigned nTargets,
	      double minTarget = 0.1, double maxTarget = 0.9);
  ~TrainingSet();

// Set functions
  void set(unsigned nInputs, unsigned nTargets,
	   double minTarget = 0.1, double maxTarget = 0.9);

// Get functions
  unsigned size() const { return uEndIndex; }
  const Example& operator [] (unsigned i) const { 
    return *(Example *) OrderedCltn::operator[](i); }
  unsigned nInputs() const  { return _nInputs; }
  unsigned nTargets() const { return _nTargets; }

// Adding examples
  void add(unsigned, const double *); // Label and input values supplied
                                       // input and target values supplied
  void add(const double *, const double *); 
                                       // Add a mixture of two classes
  void add(unsigned, const double *, unsigned, const double *, double pct);
                                       // Add two-class mixtures (partial volume)
  void add(const Array<SimpleArray<unsigned> >& mixtures, unsigned nMixtures); 
                                       // Same, but using the means
  void add(const Array<SimpleArray<unsigned> >& mixtures, unsigned nMixtures, unsigned nCopies); 

// Removing examples
  void removeAll();

// Misc functions
  void     shuffle();
  std::ostream& print(std::ostream&) const;

friend class TrainingSetIterator;
};

/****************************
 * TrainingSet iterator class
 ****************************/

class TrainingSetIterator {
  const TrainingSet *_trainingSet;
  unsigned           _exampleIndex;

public:
  TrainingSetIterator(const TrainingSet& trainingSet) { 
    _trainingSet = &trainingSet; _exampleIndex = 0; }

  void first() { _exampleIndex = 0; }
  void last()  { _exampleIndex = _trainingSet->uEndIndex - 1; }

  Example *operator ++ ();    // Prefix
  Example *operator ++ (int); // Postfix
  Example *operator -- ();    // Prefix
  Example *operator -- (int); // Postfix
};

#endif
