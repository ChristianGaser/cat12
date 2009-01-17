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
$RCSfile: MatrixTest.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MATRIX_H
#define _MATRIX_H

/******************************************************************************
 * The user must use the following #defines for each matrix type that is 
 * required. Selective use of these #defines will greatly reduce mushrooming
 * of object code due to (unneccesary) template instantiations. Although most
 * compilers provide facilities to control template instantiations, this
 * approach was found more robust, for a large part due to the type conversions
 * provided by this libary.
 *   Without using the #defines, basic functionality will still be available,
 * but the typedefs and some functionality (such as type conversions) will 
 * be missing.
 *
 * Alex Zijdenbos, Aug 11, 1995.
 *
 * #define     type              typedef provided
 *-----------------------------------------------
 * USE_UCHRMAT (unsigned char)   UChrMat
 * USE_CHRMAT  (char)            ChrMat
 * USE_USHMAT  (unsigned short)  UShMat
 * USE_SHMAT   (short)           ShMat
 * USE_UINTMAT (unsigned int)    UIntMat
 * USE_INTMAT  (int)             IntMat
 * USE_FLMAT   (float)           FlMat
 * USE_DBLMAT  (double)          DblMat
 * USE_COMPMAT (complex)         CompMat
 * USE_DBLMAT  (double)          fCompMat
 *****************************************************************************/

#include <iostream.h>
#include <stream.h>
#include <math.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>

#ifdef USE_COMPMAT
#include <complex.h>
#endif

#ifndef MIN
#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))
#endif

const MATLAB = 0;
const RAW = 1;
const ASCII  = 2;

template <class Type> class Mat {
protected:
/*****************************Internal data elements***************************/
   unsigned   _rows;
   unsigned   _cols;
   unsigned   _maxrows;
   unsigned   _maxcols;
   Type	    **_el;
/******************************************************************************/
   
public:

/*******************************Constructors***********************************/
  Mat(unsigned nrows, unsigned ncols, Type value);
  //Copy constructor
  Mat(const Mat& A);
  //Calls the clear function to deallocate memory
  ~Mat();

/***********************Explicit Memory deallocation***************************/
   //May be called within the program to deallocate the dynamic memory
   //assigned to an object.  The object itself is not destroyed until
   //it goes out of scope.
   void clear();

private:
   void _allocateEl();
};

#ifdef USE_COMPMAT
class Mat<complex> {
protected:
/*****************************Internal data elements***************************/
   unsigned   _rows;
   unsigned   _cols;
   unsigned   _maxrows;
   unsigned   _maxcols;
   complex  **_el;
/******************************************************************************/
   
public:

/*******************************Constructors***********************************/
  Mat(unsigned nrows, unsigned ncols, complex value);
  //Copy constructor
  Mat(const Mat& A);
  //Calls the clear function to deallocate memory
  ~Mat();

/***********************Explicit Memory deallocation***************************/
   //May be called within the program to deallocate the dynamic memory
   //assigned to an object.  The object itself is not destroyed until
   //it goes out of scope.
   void clear();

private:
   void _allocateEl();
};
#endif

#endif
