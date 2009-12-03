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
$RCSfile: Histogram.h,v $
$Revision: 1.4 $
$Author: bert $
$Date: 2004/12/08 16:42:55 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>		/* (bert) changed from iostream.h */
#include "MTypes.h"
#include "ValueMap.h"
#include "SimpleArray.h"

class Histogram : private SimpleArray<unsigned> {
  double    _min, _max; // True extrema (i.e., not the bin centers)
  double    _binWidth;
  LinearMap _valueToBinMap;

public:
// Constructors/destructor
  // The min and max values in these constructors are bin centers,
  // not absolute extrema.
  Histogram();
  Histogram(double min, double max, unsigned nBins = 0);
  Histogram(double min, double max, double binWidth);
  Histogram(unsigned nBins, double min = 0.0, double binWidth = 1.0);
  Histogram(const Histogram&);
  virtual ~Histogram() {}

  Histogram& operator = (const Histogram&);
  Histogram& newRange(double min, double max); // Changes ranges; keeps binWidth

// Get functions
  unsigned  operator [] (unsigned i) const { 
    return SimpleArray<unsigned>::operator [] (i); }
  unsigned& operator [] (unsigned i) { 
    return SimpleArray<unsigned>::operator [] (i); }
  double    binWidth() const            { return _binWidth; }
  unsigned  nBins() const               { return _size; }
  double    binStart(unsigned i) const  { return _min + i*_binWidth; }
  DblArray  binStarts() const;
  double    binCenter(unsigned i) const { return binStart(i) + _binWidth/2; }
  DblArray  binCenters() const          { return binStarts() + _binWidth/2; }
  unsigned  n() const                   { return unsigned(sum()); }
  unsigned  count(double value) const   { return _contents[bin(value)]; }
  unsigned  bin(double value) const;
  unsigned  max(unsigned *bin = 0) const;
  double    mean() const;
  double    median(unsigned *bin, unsigned nBelow = 0, unsigned nAbove = 0) const;
  // (bert) - This median() added just to make the IRIX compiler happier.
  unsigned  median() const { return SimpleArray<unsigned>::median(); }
  double    majority(unsigned *bin = 0) const;
  double    biModalThreshold() const;
  double    varianceThreshold() const;
  double    kullbackThreshold() const;
  double    pctThreshold(double pct) const;
  double    entropy() const;
  DblArray  pdf() const;
  DblArray  cdf() const;

// Set functions
  Boolean add(double value) {
    if ((value < _min) || (value > _max))
      return FALSE;
    unsigned index = (unsigned) _valueToBinMap(value);
    if (index >= _size)
      index = _size - 1;
    _contents[index]++;
    return TRUE;
  }

// Other operators
  Histogram&  operator += (const Histogram& hist);
  LUT<double> equalize(const Histogram& hist) const;

// I/O
  std::ostream& printHeadAndTail(std::ostream& os, unsigned n = 10) const;
#ifdef HAVE_MATLAB
  Boolean  saveMatlab(const char *fileName, const char *binVarName = "X", 
		      const char *countVarName = "N", const char *option = "u") const;
#endif

// Friends
  friend std::ostream&    operator << (std::ostream& os, const Histogram& hist);
  friend SimpleArray<double> asDblArray(const Histogram& hist);
};

// (bert) - moved the following three functions from inside class definition:
//
inline DblArray pdf(const Histogram& hist) { return hist.pdf(); }

inline DblArray cdf(const Histogram& hist) { return hist.cdf(); }

inline LUT<double> equalize(const Histogram& hist1, const Histogram& hist2)
{
  return hist1.equalize(hist2);
}

template <class Type>
Histogram &
add(Histogram & hist, const SimpleArray <Type> & array)
{
  if (array.size()) {
    array.resetIterator();
    for (unsigned i = array.size(); i; i--)
      hist.add(array++);
  }

  return hist;
}

#endif
