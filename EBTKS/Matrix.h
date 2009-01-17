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
$RCSfile: Matrix.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MATRIX_H
#define _MATRIX_H

#include "SimpleArray.h"

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
 * USE_UCHRMAT  (unsigned char)   UChrMat
 * USE_CHRMAT   (char)            ChrMat
 * USE_USHMAT   (unsigned short)  UShMat
 * USE_SHMAT    (short)           ShMat
 * USE_UINTMAT  (unsigned int)    UIntMat
 * USE_INTMAT   (int)             IntMat
 * USE_FLMAT    (float)           FlMat
 * USE_DBLMAT   (double)          DblMat
 * USE_COMPMAT  (dcomplex)        CompMat
 * USE_FCOMPMAT (fcomplex)        fCompMat
 *****************************************************************************/

#include "dcomplex.h"

#ifdef USE_FCOMPMAT
  #include "fcomplex.h"
#endif

#ifdef HAVE_MATLAB
extern "C" {
  #include"mat.h" 
}
#endif

/* #include <stream.h> (bert) removed */
#include <math.h>
#include <fstream>		/* (bert) changed from fstream.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>
#include <iostream>		/* (bert) changed from iostream.h */
#include "MTypes.h"
#include "MatrixSupport.h"
#include "Histogram.h"

#ifndef MIN
#define MIN(x, y) ((x < y) ? (x) : (y))
#define MAX(x, y) ((x > y) ? (x) : (y))
#endif

const int MATLAB = 0;
const int RAW = 1;
const int ASCII  = 2;

typedef double  (*IndexFunction)(unsigned r, unsigned c);
typedef dcomplex (*ComplexIndexFunction)(unsigned r, unsigned c);

/******************** Explicit type conversions *******************/

// Forward declarations
template <class Type> class Mat;
template <class Type> class MatrixIterator;
template <class Type> class ConstMatrixIterator;

#ifdef USE_DBLMAT
//Converts Mat<Type> to Mat<double>
typedef Mat<double> DblMat;
template <class Type> Mat<double> asDblMat(const Mat<Type>&);
#endif

#ifdef USE_FLMAT
//Converts Mat<Type> to Mat<float>
typedef Mat<float> FlMat;
template <class Type> Mat<float> asFlMat(const Mat<Type>&);
#endif

#ifdef USE_INTMAT
//Converts Mat<Type> to Mat<int>
typedef Mat<int> IntMat;
template <class Type> Mat<int> asIntMat(const Mat<Type>&);
#endif

#ifdef USE_UINTMAT
//Converts Mat<Type> to Mat<unsigned int>
typedef Mat<unsigned int> UIntMat;
template <class Type> Mat<unsigned int> asUIntMat(const Mat<Type>&);
#endif

#ifdef USE_SHMAT
//Converts Mat<Type> to Mat<short>
typedef Mat<short> ShMat;
template <class Type> Mat<short> asShMat(const Mat<Type>&);
#endif

#ifdef USE_USHMAT
//Converts Mat<Type> to Mat<unsigned short>
typedef Mat<unsigned short> UShMat;
template <class Type> Mat<unsigned short> asUShMat(const Mat<Type>&);
#endif

#ifdef USE_CHRMAT
//Converts Mat<Type> to Mat<char>
typedef Mat<char> ChrMat;
template <class Type> Mat<char> asChrMat(const Mat<Type>&);
#endif

#ifdef USE_UCHRMAT
//Converts Mat<Type> to Mat<unsigned char>
typedef Mat<unsigned char> UChrMat;
template <class Type> Mat<unsigned char> asUChrMat(const Mat<Type>&);
#endif

#ifdef USE_COMPMAT
typedef Mat<dcomplex> CompMat;
template <class Type> Mat<dcomplex> asCompMat(const Mat<Type>& A);
template <class Type> Mat<dcomplex> asCompMat(const Mat<Type>& Re, const Mat<Type>& Im);
#endif

#ifdef USE_FCOMPMAT
typedef Mat<fcomplex> fCompMat;
template <class Type> Mat<fcomplex> asFcompMat(const Mat<Type>&);
template <class Type> Mat<fcomplex> asFcompMat(const Mat<Type>& Re, const Mat<Type>& Im);
#endif

// Some mathematical operations. These had to be implemented as template
// functions in order to allow different types of matrices to be operated on
// ***************************** Addition ****************************
template <class T1, class T2>
Mat<T1>& operator += (Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> operator + (const Mat<T1>& A, const Mat<T2>& B);

// ***************************** Substraction *************************
template <class T1, class T2>
Mat<T1>& operator -= (Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> operator - (const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return T -= B; }

// ***************************** Multiplication ***********************
template <class T1, class T2>
Mat<T1> operator * (const Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1>& operator *= (Mat<T1>& A, const Mat<T2>& B) {
  A = A * B; return A; }

// ***************************** Division ***** ***********************
template <class T1, class T2>
Mat<T1> operator / (const Mat<T1>& A, const Mat<T2>& B) {
  return A * inv(B); }

template <class T1, class T2>
Mat<T1>& operator /= (Mat<T1>& A, const Mat<T2>& B) {
  A = A * inv(B); return A; }

// *********************** Point multiplication ***********************
template <class T1, class T2>
Mat<T1>& pmultEquals(Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> pmult(const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return pmultEquals(T, B); }

// ************************* Point division ***************************
template <class T1, class T2>
Mat<T1>& pdivEquals(Mat<T1>& A, const Mat<T2>& B);

template <class T1, class T2>
Mat<T1> pdiv(const Mat<T1>& A, const Mat<T2>& B) {
  Mat<T1> T(A); return pdivEquals(T, B); }

/******************************************************************************
 *     MATRIX CLASS
 ******************************************************************************/

template <class Type> 
class Mat {
  friend class MatrixIterator<Type>;
  friend class ConstMatrixIterator<Type>;

protected:
/*****************************Internal data elements***************************/
   unsigned   _rows;
   unsigned   _cols;
   unsigned   _maxrows;
   unsigned   _maxcols;
   Type	    **_el;
   static unsigned _rangeErrorCount;
/******************************************************************************/
   
public:
  static Boolean flushToDisk;
  enum vector_orientation { ROW = 0, COLUMN = 1}; 

/*******************************Constructors***********************************/
   //Default constructor-try not to use this one. An object instantiated using
   //this constructor cannot be reassigned or copied to
   Mat(void) { _rows = _maxrows = 0; _cols = _maxcols = 0; _el=0; }
   
   //Copy constructor
   Mat(const Mat& A);

  // Create row (default) or column vector from SimpleArray
  // Dir spec removed because DCC can't distinghuish this from Mat(unsigned,unsigned)???
  Mat(const SimpleArray<Type>& A, char dir = Mat<Type>::ROW);
 
   //When the below constructor is used for type character, single quotes must
   // be used for the value which the matrice is to be filled with 
   Mat(unsigned nrows, unsigned ncols, Type value); 

   //Initializes all elements to zero
   Mat(unsigned nrows, unsigned ncols);

   Mat(unsigned nrows, unsigned ncols, const Type *data);
   
   //Initialize object from a matlab file
   Mat(const char *filename, int type = RAW, const char *varname = "A");

/********************************Destructor************************************/
   //Calls the clear function to deallocate memory
  virtual ~Mat();

/***********************Explicit Memory deallocation***************************/
   //May be called within the program to deallocate the dynamic memory
   //assigned to an object.  The object itself is not destroyed until
   //it goes out of scope.
   void clear();

/*******************Element wise accessing operators***************************/
   //Matrix is treated as a linear array with rows being concatenated
   //i.e.  [ 1 2 3     is viewed as [1 2 3 4 5 6]
   //        4 5 6 ]  
   //and single argument n is used to return the appropriate element
   //As with all matrix accessing operators, the first element is element 
   //zero
   //Does allow modification of the element      
   Type& operator() (unsigned n);
   
   //Same as above but does not allow modification of the actual element
   Type  operator() (unsigned n) const;

   //Returns the address of the element of the Matrix specified by row,column
   //Allows modification of the actual element
   Type& operator() (unsigned r, unsigned c);
   
   //Returns the element of the Matrix specified by row,column
   //Does not allow modification of the actual element
   Type  operator() (unsigned r, unsigned c) const;

   //Return a Matrix(Mat) object consisiting of the subsection of the original
   //Matrix(Mat) specified by r1(starting row) to r2(ending row) and c1
   //(starting column) to c2(ending column)
   Mat operator()(unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;
   Mat row(unsigned r) const                { return (*this)(r, r, 0, _cols - 1); }
   Mat rows(unsigned r1, unsigned r2) const { return (*this)(r1, r2, 0, _cols - 1); }
   Mat col(unsigned c) const                { return (*this)(0, _rows - 1, c, c); }
   Mat cols(unsigned c1, unsigned c2) const { return (*this)(0, _rows - 1, c1, c2); }   
/*****************Functions to access Internal data pointers*******************/
   //Returns a double pointer to the first element of the data
   //The primary pointer points to the first of a group of
   //secondary pointers(one for each row).  Each of the secondary pointers
   //points to the first data element for that row
   //The primary pointer that is returned cannot be modified in itself
   //However, a copy of the value of the secondary pointer that it points to
   //can be made and then changes may ne made dirrectly to the encapsulated data
   //This violates the protection provided by the class but may be 
   //necessary in certain places to increase execution speed 
   const Type **getEl() const { return (const Type **) _el; } // Be careful with this one

   //This returns a copy of the secondary pointer to the beginning of the
   //row specified by indx
   const Type* operator[] (unsigned indx) const;

/*********************Exclusive Reassignment Operators*************************/
  //This operator reassigns the calling object to a copy of the argument object
  void  operator() (const Mat&);
  
  //This operator functions identically to the one above execept that it 
  //returns a pointer to the calling object so that assignments maybe cascaded
  Mat& operator = (const Mat& A);
  
  // Functionally identical to operator = (), but directly copies the data
  // pointers from A, thus destroying it (!!). 
  Mat& absorb(Mat& A);

/**************************Matrix reshaping functions**************************/
  // Padding functions; return matrices padded with <value> all around. 
  // Generic pad function; pads to nrows x ncols and inserts at row, col.
  Mat& pad(unsigned nrows, unsigned ncols, int row, int col, Type value = 0);
  Mat  pad(unsigned nrows, unsigned ncols, int row, int col, 
	   Type value = 0) const {
    return padConst(nrows, ncols, row, col, value); }
  Mat  padConst(unsigned nrows, unsigned ncols, int row, int col, 
		Type value = 0) const {
    Mat<Type> A(*this); return A.pad(nrows, ncols, row, col, value); }
  
  // Symmetrical pad, i.e., the size of the returned matrix is (rows+2*rowpad, 
  // cols+2*colpad).
  Mat& pad(unsigned rowpad, unsigned colpad, Type value = 0) {
    return pad(_rows + 2*rowpad, _cols + 2*colpad, rowpad, colpad, value); }
  Mat  pad(unsigned rowpad, unsigned colpad, Type value = 0) const {
    return padConst(rowpad, colpad, value); }
  Mat  padConst(unsigned rowpad, unsigned colpad, Type value = 0) const {
    return padConst(_rows + 2*rowpad, _cols + 2*colpad, rowpad, colpad, value); }
  
  // Resize function. 
  // If the matrix is made smaller, the elements not included in the new matrix
  // will remain in existence (i.e., the memory occupied is not freed).
  // If the matrix is made larger, the existing elements will be preserved.
  // The old matrix elements matrix will be entered at (row, col) of the new
  // resized matrix (row and col may be negative, defaults: 0);
  Mat& resize(unsigned nrows, unsigned ncols, int row = 0, int col = 0);
  
  // This function returns the matrix dimensions to their maximum values
  // before any previous resizing had been done.  Data in the newly resized
  // portion may or may not be the same as what was previously there
  Mat& resize();
  
  //Returns a new matrix which is the side by side concatenation of the calling
  //matrix and the operand
  //No changes are made to the operand or calling matrix object
  Mat appendRight(const Mat& A) const;
  
  //Returns a new matrix which is the top and bottom concatenation of the
  //calling matrix and the operand.
  //No changes are made to the operand or calling matrix object
  Mat appendBelow(const Mat& A) const;
  
  // Swap two rows or columns
  Mat& swapRows(unsigned r1, unsigned r2);
  Mat& swapCols(unsigned c1, unsigned c2);

  //Returns the residual matrix after eliminating row and col
  Mat residual(unsigned row, unsigned col) const;
  
  // Insert a matrix into another at (row, col). If necessary, the argument
  // matrix will be clipped to the bounds of the source matrix.
  Mat& insert(const Mat<Type>& A, int row = 0, int col = 0);
  // Same, but reads the matrix from a raw file
  Mat& insert(const char *path, unsigned nrows=0, unsigned ncols=0, int row=0,int col=0);
  
/*******************Special Matrix reinitialization functions******************/
  // Fills the matrix with <value>
  Mat& fill(Type value);
  
  // Fills a sub-rectangle of the matrix with <value>
  Mat& fill(Type value, unsigned r1, unsigned r2, unsigned c1, unsigned c2);

  // Fills the matrix with <value> in a circular area, centered in the matrix.
  // When <diameter> <= 0, the circle will extend to the edge of the matrix.
  Mat& fillCircle(Type value, double diameter = 0) {
    if (diameter <= 0) diameter = MIN(_rows, _cols);
    return fillEllips(value, diameter, diameter); 
  }

  // Fills the matrix with <value> in the specified circular area, 
  // centered at (row, col). When diameter <= 0, the circle will extend to the
  // closest edge of the matrix.
  Mat& fillCircle(Type value, double row, double col, double diameter = 0) {
    if (diameter <= 0) {
      double rowDiameter = 2*MIN(row + 0.5, _rows - row - 0.5);
      double colDiameter = 2*MIN(col + 0.5, _cols - col - 0.5);
      diameter = MIN(rowDiameter, colDiameter);
    }
    return fillEllips(value, row, col, diameter, diameter); 
  }

  // Fills the matrix with <value> in a circular area, centered in the matrix.
  // When one of the diameters is <= 0, the ellips will extend to the edge of the
  // matrix.
  Mat& fillEllips(Type value, double rowDiameter = 0, double colDiameter = 0);

  // Fills the matrix with <value> in the specified elliptical area, 
  // centered at (row, col). If one of the diameters is <= 0, the ellips will
  // extend to the closest edge of the matrix
  Mat& fillEllips(Type value, double row, double col, double rowDiameter, 
		  double colDiameter);

   //Turns the calling object into a identity matrix
   Mat& eye();
   
   //Fills the calling object with random numbers
   Mat& randuniform(double min = 0, double max = 1);
   
   //Stores a random normal matrix to the calling object
   Mat& randnormal(double mean = 0, double std = 1);
      
   //Stores a hamming window to the calling object
   Mat& hamming();

   //Stores a hanning window to the calling object
   Mat& hanning();

   //Stores a blackman window to the calling object
   Mat& blackman();
   
   //This creates a diagonal matrix from the existing row or column vector objec
   //This function uses(internally) the () operator of this class(The one that 
   //takes one argument and the one that takes two arguments)      
   Mat diag() const;
   
   //This creates a toeplitz matrix using two vectors, the calling object is
   //used as the column vector and the argument as the row vector
   //With no argument, the calling object is used for both
   //This function requires use of the () operator and the appendBelow function
   Mat toeplitz(const Mat& r)const;
   Mat toeplitz() const { return  toeplitz(*this); }

/***********Functions to examine a Matrix(Mat) objects properties**************/

  // Returns non-zero if the matrix is empty
  Boolean operator ! () const { return Boolean(!_rows || !_cols || !_el); }
  operator void * () const    { return (void *) (_rows && _cols && _el); }

   //Returns the number of ACTIVE rows in a Mat object
   unsigned  getrows() const { return _rows; }

   //Returns the number of ACTIVE columns in a Mat object
   unsigned  getcols() const { return _cols; }

   //Returns the number of ACTIVE elements in a Mat object
   unsigned  nElements() const { return _rows*_cols; }

   //Returns the number of ACTUAL rows in a Mat object
   unsigned  getmaxrows() const { return _maxrows; }

   //Returns the number of ACTUAL columns in a Mat object
   unsigned  getmaxcols() const { return _maxcols; }

   //Returns 1 if the ACTIVE Matrix(Mat) is a vector else 0
   int    isvector() const       { return (_rows == 1) || (_cols == 1); }

   //Returns 1 if the ACTIVE Matrix(Mat) is a column vector else 0
   int    iscolumnvector() const { return (_cols == 1); }

   //Returns 1 if the ACTIVE Matrix(Mat) is a row vector else 0
   int    isrowvector() const    { return (_rows == 1); }
  
   //Returns the length of the larger ACTIVE dimension(row vs column) of the
   //Matrix
   unsigned    length() const         { return (_rows > _cols) ? _rows : _cols; }
  
   //Display the contents of the entire ACTIVE Matrix(Mat) object
   std::ostream& display(std::ostream& os) const;
   
   //Display a subsection of the Matrix(Mat) defined by r1 (starting row) to
   //r2(ending row) and c1(starting column) to c2(ending column)
   std::ostream& display(std::ostream& os, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;

/**********************Basic Matrix Arithmetic operators***********************/
  // Unary -
  Mat operator - () const { Mat<Type> A(_rows, _cols, Type(0)); return A -= *this; }

  // Elementwise comparison of two matrices
  int operator != (const Mat&) const;
  int operator == (const Mat& A) const { return !(operator != (A)); }
  
  //Adds a single value to each of the elements of the calling Matrix
  //(on the right hand side) and returns a pointer to this Matrix
  //Right-hand operand is modified
  Mat& operator += (dcomplex);
  
  //Same as above, except that no changes are made to either operand,
  //rather a new Matrix is constructed and returned by the function 
  Mat  operator + (dcomplex x) const      { Mat<Type> A(*this); return A += x; }
  
  //Subtracts a single value from each of the elements of the calling Matrix
  //(on the right hand side) and returns a pointer to this Matrix
  //Right-hand operand is modified
  Mat& operator -= (dcomplex x)            { return *this += (-x); }
  
  //Same as above, except that no changes are made to either operand,
  //rather a new Matrix is constructed and returned by the function
  Mat  operator -  (dcomplex x) const      { return *this +  (-x); }
  
  //Multiplys each element of the right hand side Matrix by a double 
  //on the left hand side.  Changes are propogated to the right hand
  //side Matrix
  Mat& operator *= (dcomplex);
  
  //Same as above, except that changes are not made to either operand
  Mat  operator *  (dcomplex x) const        { Mat<Type> A(*this); return A *= x; }
  
  //Divides each element of the Matrix by a complex and effects the changes in
  //the original Matrix
  Mat& operator /= (dcomplex x)              { return (*this) *=  (1.0/x); }
  
  //Same as above, except that changes are not made to either operand
  Mat  operator /  (dcomplex x) const        { return (*this) *  (1.0/x); }
  
  /********************General purpose Matrix functions**************************/
  //Returns a Matrix which is the transpose of the original matrix
  //Does not change the original matrix
  Mat t() const;
  
  //Returns a Matrix which is the reverse or 180 rotate 
  //Does not change the original matrix
  Mat rotate180() const;
  
  //Returns a Matrix which is the inverse of the original matrix
  //Does not change the original matrix. Returns an empty matrix on failure
  Mat inv() const;
  
  //Returns the hermitan transpose of the original Matrix
  //Does not change the original matrix
  Mat h() const;
  
  //Returns the value of the smallest element of the ACTIVE matrix
  //optionally returns the index (row and column) of the element
  Type min(unsigned *row = 0, unsigned *col = 0) const;
  
  //Returns the value of the largest element of the ACTIVE matrix
  //optionally returns the index (row and column) of the element
  Type max(unsigned *row = 0, unsigned *col = 0) const;
  
  // Returns the median of the matrix, only considering elements in the supplied range.
  // By default, the full range of the matrix will be used.
  Type median(Type minVal = 0, Type maxVal = 0) const;

  //Returns the sum of all the elements of the ACTIVE matrix
  double  sum() const { return real(csum()); }
  dcomplex csum() const;
  
  //Returns the sum of squares of all the elements of the ACTIVE matrix
  double  sum2() const { return real(csum2()); }
  dcomplex csum2() const;
  
  //Returns the sum of all the elements of the ACTIVE matrix divided by the
  //number of active elements
  double  mean() const  { return sum()/nElements(); }
  dcomplex cmean() const { return csum() / double(nElements()); }

  //Returns the standard deviation of the matrix elements
  double  std() const  { return ::sqrt(var()); }
  dcomplex cstd() const { return std::sqrt(cvar()); }

  //Returns the variance of the matrix elements
  double  var() const  { double  mn = mean(); return sum2()/nElements() - SQR(mn); }
  dcomplex cvar() const { dcomplex mn = cmean(); return csum2() / double(nElements()) - SQR(mn); }
  
  //Returns the norm of the matrix
  double  norm() const  { return ::sqrt(sum2()); }
  dcomplex cnorm() const { return std::sqrt(csum2()); }
  
  //Returns the trace of the matrix
  double  trace() const { return real(ctrace()); }
  dcomplex ctrace() const;
  
  //Returns the determinant of the matrix
  double  det() const { return real(cdet()); }
  dcomplex cdet() const;
  
  /**********************Other Matrix element wise functions*********************/
  //Applies a given single argument function to each of the elements
  Mat&  applyElementWise(double (*function)(double));
  Mat   applyElementWise(double (*function)(double)) const {
    return applyElementWiseConst(function); }
  Mat   applyElementWiseConst(double (*function)(double)) const {
    Mat<Type> A(*this); return A.applyElementWise(function); }

  Mat&  applyElementWiseC2D(double (*function)(dcomplex));
  Mat   applyElementWiseC2D(double (*function)(dcomplex)) const {
    return applyElementWiseConstC2D(function); }
  Mat   applyElementWiseConstC2D(double (*function)(dcomplex)) const {
    Mat<Type> A(*this); return A.applyElementWiseC2D(function); }

  Mat&  applyElementWiseC2C(dcomplex (*function)(dcomplex));
  Mat   applyElementWiseC2C(dcomplex (*function)(dcomplex)) const {
    return applyElementWiseConstC2C(function); }
  Mat   applyElementWiseConstC2C(dcomplex (*function)(dcomplex)) const {
    Mat<Type> A(*this); return A.applyElementWiseC2C(function); }
  
  Mat& applyIndexFunction(IndexFunction F);
  Mat  applyIndexFunction(IndexFunction F) const {
    return applyIndexFunctionConst(F); }
  Mat   applyIndexFunctionConst(IndexFunction F) const {
    Mat<Type> A(*this); return A.applyIndexFunction(F); }

  Mat& applyIndexFunction(ComplexIndexFunction F);
  Mat  applyIndexFunction(ComplexIndexFunction F) const {
    return applyIndexFunctionConst(F); }
  Mat   applyIndexFunctionConst(ComplexIndexFunction F) const {
    Mat<Type> A(*this); return A.applyIndexFunction(F); }

  //Takes the inverse natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& exp();
  Mat  exp() const      { return expConst(); }
  Mat  expConst() const { Mat<Type> T(*this); return T.exp(); }
  
  //Takes the natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& log();
  Mat  log() const      { return logConst(); }
  Mat  logConst() const { Mat<Type> T(*this); return T.log(); }
  
  //Takes the cosine of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat& cos();
  Mat  cos() const      { return cosConst(); }
  Mat  cosConst() const { Mat<Type> T(*this); return T.cos(); }
  
  //Takes the sine of each of the elements of the matrix
  Mat& sin();
  Mat  sin() const      { return sinConst(); }
  Mat  sinConst() const { Mat<Type> T(*this); return T.sin(); }
  
  //Returns the conjugate matrix
  Mat& conj();
  Mat  conj() const      { return conjConst(); }
  Mat  conjConst() const { Mat<Type> T(*this); return T.conj(); }
  
  //Takes the absolute value of each of the elements of the matrix
  Mat& abs();
  Mat  abs() const      { return absConst(); }
  Mat  absConst() const { Mat<Type> T(*this); return T.abs(); }
  
  //Rounds each of the elements of the matrix
  Mat& round();
  Mat  round() const      { return roundConst(); }
  Mat  roundConst() const { Mat<Type> T(*this); return T.round(); }

  //Takes the sqrt value of each of the elements of the matrix
  Mat& sqrt();
  Mat  sqrt() const      { return sqrtConst(); }
  Mat  sqrtConst() const { Mat<Type> T(*this); return T.sqrt(); }

  //Takes the power value of each of the elements of the matrix
  Mat& pow(double exponent);
  Mat  pow(double exponent) const      { return powConst(exponent); }
  Mat  powConst(double exponent) const { Mat<Type> T(*this); return T.pow(exponent); }

/************************Special purpose Matrix functions**********************/
  // Returns transpose(A)*A of the matrix. 
  Mat transposeXself() const;

   void   eig(Mat& D, Mat& V) const;

   //Returns a Matrix which is the householder transform of the calling object
   //Calling object must be a vector
   Mat house() const;

   //Returns a Matrix which is the rowhouse transform of the calling object and
   //the input column vector
   Mat rowhouse(const Mat& V) const;

   //convolve 2d Matrix //taken from red book and modified
   Mat convolv2d(const Mat& filter) const;

  //Returns a histogram for the calling object using the specified range and # bins
  Histogram histogram(double minin = 0, double maxin = 0, unsigned n = 0) const;

  Mat& histmod(const Histogram& hist1, const Histogram& hist2);
  Mat  histmod(const Histogram& hist1, const Histogram& hist2) const {
    return histmodConst(hist1, hist2); }
  Mat  histmodConst(const Histogram& hist1, const Histogram& hist2) const {
    Mat<Type> A(*this); return A.histmod(hist1, hist2); }

  // Returns an array which contains only the values in the range specified
  SimpleArray<Type> array(Type minVal = 0, Type maxVal = 0) const;

  // Returns first (topleft) element
  Type              scalar() const { return **_el; }

   //Performs a QR factorization and modifies Q and R appropriately
   //In this case, X=Q*R
   //Where X is the calling object
   void   qr(Mat& R, Mat& Q) const;

   //Performs a QR factorization using pivoting.
   //In this case X*P=Q*R
   //Where X is the calling object
   void   qr(Mat& R, Mat& Q, Mat& P) const;

   //Performs the singular value decomposition on the calling matrix
   //and modifies U, S, and V appropriately
   //In this case X = U*S*(V')
   //Where X is the calling object
   void   svd(Mat& U, Mat& S, Mat& V) const;

  // Sets all values in the matrix below minVal to minFill, and those
  // above maxVal to maxFill.
  Mat& clip(Type minVal, Type maxVal, Type minFill, Type maxFill);
  Mat  clip(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    return clipConst(minVal, maxVal, minFill, maxFill); }
  Mat  clipConst(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    Mat<Type> A(*this); return A.clip(minVal, maxVal, minFill, maxFill); }

  // Maps the values of a matrix through a value map
  Mat& map(const ValueMap& valueMap);
  Mat  map(const ValueMap& valueMap) const { return mapConst(valueMap); }
  Mat  mapConst(const ValueMap& valueMap) const { 
    Mat<Type> A(*this); return A.map(valueMap);}

  // Scales a Matrix (similar to map, but truncates the output range to be
  // [minout, maxout])
  Mat& scale(double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin = 0.0);
  Mat  scale(double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin = 0.0) const {
    return scaleConst(minout, maxout, minin, maxin); }
  Mat  scaleConst(double minout = 0.0, double maxout = 255.0, 
		  double minin = 0.0, double maxin = 0.0) const {
    Mat<Type> A(*this); return A.scale(minout, maxout, minin, maxin); }

   //Performs a linear interpolation on the calling object Matrix
   //The y (argument) Matrix recieves the results
   //Urow is the upscale factor for the rows and Drow is the downscale factor
   //for the rows.  Ucol and Dcol are analagous for the columns.
   void linearinterpolate( unsigned Urow, unsigned Drow, unsigned Ucol, unsigned Dcol,Mat& y) const;

   //Applies the filter defined by B and A to the calling object and returns the
   //resultant Matrix.
   //The filter has the form y= (B/A)*x
   //Where x is the calling object and y is the result
   Mat filter(Mat& B,Mat& A) const;

/***************************morphology functions*******************************/
   // In the following (greylevel) morphology operations, negative strel
   // values are considered NaN's and ignored.
#ifdef USE_DBLMAT
   Mat erode(const Mat<double>& strel) const;
   Mat dilate(const Mat<double>& strel) const;
   Mat open(const Mat<double>& strel) const  { return erode(strel).dilate(strel); }
   Mat close(const Mat<double>& strel) const { return dilate(strel).erode(strel); }
#endif

  //   Not converted yet
  //   void radialscale(Mat& scale, Mat& y) const;

  //This function converts a Matlab matrix object to a Mat object
  //The first argument is a pointer to a Matlab matrix object
  //This function checks the storage format of the Matlab object
  //(Full or Sparse) and does the appropriate conversions to return a Mat
  //object of the calling object's type.
#ifdef HAVE_MATLAB
  Mat matlab2Mat(Matrix *MatlabMatrix) const;
  
  //Creates a matlab type matrix and assigns it the name specified as the
  //argument.  It then returns a pointer to this object
  Matrix *mat2Matlab(const char *name) const;
#endif

/*************************Matrix file I/O functions****************************/
   //This function goes to the named mat file and searches for the given
   //variable.  Next, if possible, it reinitializes the calling object to
   //represent the same matrix as that defined by the Matlab variable.  The
   //Matlab matrix is stored in double format.  The function Matlab2Mat will
   //convert the Matlab data to the type of its calling object (this).  This may
   //mean a loss of precision.
#ifdef HAVE_MATLAB
  Boolean loadMatlab(const char *filename, const char *varname = "A");
#endif

   //Loads a binary version of a matrix from a file
   //must have previously been stroed by saveRaw
   //The calling object must have the same type as the object that was 
   //previously saved in order for this function to work correctly
  Boolean loadRaw(const char *filename, unsigned nrows = 0, unsigned ncols = 0);

   //Ascii files will always have rows and columns up top
   //Loads an Ascii version of the saved Matrix
  Boolean loadAscii(const char *filename);
   
   //Loads a matrix from either a binary, Ascii or Matlab file, determined
   //by the type flag
  Boolean load(const char *filename, int type = RAW, const char *varname = "A");
   //Saves the calling object in either binary,Ascii or Matlab style, determined
   //by a flag set in this header file
  Boolean save(const char *filename, int type = RAW, const char *varname = "A") const;

   //Saves the calling object as a Matlab readable matrix to the namedfile with
   //given variable name.  Matlab stores in double format so, the calling object
   //is first converted to a temporary Mat object of type double.
   //Matlab stores in a column based format, so the Mat<double> is transposed
   //before saving.  The third argument lets the user update a file("u") or
   //rewrite an existing file or create a new one("w") 
#ifdef HAVE_MATLAB
  Boolean saveMatlab(const char *filename, const char *varname = "A", 
		     const char *option = "u") const;
#endif

   //Saves the calling matrix in binary format
  Boolean saveRaw(const char *) const;

   //Saves the calling object in Ascii format
   //Ascii files will always save rows and columns up top 
  Boolean saveAscii(const char *) const;
  
/*****************************fft functions************************************/

  // For the following fft functions, a dimension of 0 (which is the
  // default argument) will be expanded to the next power of two, equal
  // to or above the matrix dimension. A dimension of 1 indicates that
  // the fft should not be taken along that dimension.

  // Non-const fft function
  Mat& fft(unsigned nrows = 0, unsigned ncols = 0) {
    return _fft(nrows, ncols, ::fft); }
  // Const fft function
  Mat  fft(unsigned nrows = 0, unsigned ncols = 0) const { 
    return fftConst(nrows, ncols); }
  // Explicit const fft function
  Mat  fftConst(unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat<Type> A(*this); return A.fft(nrows, ncols); }

  // Non-const ifft function
  Mat& ifft(unsigned nrows = 0, unsigned ncols = 0);
  // Const ifft function
  Mat  ifft(unsigned nrows = 0, unsigned ncols = 0) const { 
    return ifftConst(nrows, ncols); }
  // Explicit const ifft function
  Mat  ifftConst(unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat<Type> A(*this); return A.ifft(nrows, ncols); }

/****************************End of public functions***************************/

private:

/********************************Private functions*****************************/
   //Allocate the memory for a Matrix object
   //If the object already has memory allocated to it and this function is 
   //called, it first deallocates the memory of the object
   void _allocateEl();

   //This function copies the Matrix argument B to the subsection of
   //the calling object defined by r1(starting row),r2(ending row),
   //c1(starting column),c2(ending column).  
   void _modify_sub_section(unsigned r1,unsigned r2,unsigned c1,unsigned c2,const Mat& B);

  // Given the path of a raw file, attempts to infer the matrix dimensions from
  // the size of the file.
  void _checkMatrixDimensions(const char *path, unsigned& nrows, unsigned& ncols) const;

  Mat& _fft(unsigned nrows, unsigned ncols, FFTFUNC fftFunc);

/***************************End Private functions*****************************/
   
/*******************************END OF Mat CLASS******************************/
};

// Non-member template functions

template <class T>
void clear(Mat<T> A) { A.clear(); }

template <class T>  
Mat<T> pad(const Mat<T>& A, unsigned rowpad, unsigned colpad, T value = 0)
{
  return A.padConst(rowpad, colpad, value);
}

template <class T>  
Mat<T>& resize(Mat<T> A, unsigned nrows, unsigned ncols, 
	       int row = 0, int col = 0)
{
  return A.resize(nrows, ncols, row, col);
}

template <class T>
Mat<T>& resize(Mat<T> A) { return A.resize(); }

template <class T>  
Mat<T> appendRight(const Mat<T>& A, const Mat<T>& V)
{
  return A.appendRight(V);
}

template <class T>
Mat<T> appendBelow(const Mat<T>& A, const Mat<T>& V)
{
  return A.appendBelow(V);
}
  
template <class T> 
Mat<T> residual(const Mat<T>& A, unsigned row, unsigned col)
{
  return A.residual(row, col);
}
  
template <class T>
Mat<T>& fill(Mat<T> A, T value) { return A.fill(value); }

template <class T>
Mat<T>& eye(Mat<T>& A) { return A.eye(); }

template <class T>
Mat<T>& randuniform(Mat<T>& A, double min = 0, double max = 1)
{
  return A.randuniform(min, max);
}

template <class T>
Mat<T>& randnormal(Mat<T>& A, double mean = 0, double std = 1)
{ 
  return A.randnormal(mean, std);
}

template <class T>
Mat<T>& hamming(Mat<T>& A) { return A.hamming(); }

template <class T>
Mat<T>& hanning(Mat<T>& A) { return A.hanning(); }

template <class T>
Mat<T>& blackman(Mat<T>& A) { return A.blackman(); }

template <class T>
Mat<T> diag(const Mat<T>& A) { return A.diag(); }

template <class T>
Mat<T> toeplitz(const Mat<T>& c, const Mat<T>& r)
{ 
  return(c.toeplitz(r));
}

template <class T>
Mat<T> toeplitz(const Mat<T>& c) { return(c.toeplitz(c)); }

template <class T>
unsigned getrows(const Mat<T>& A) { return A.getrows(); }

template <class T>
unsigned getcols(const Mat<T>& A) { return A.getcols(); }

template <class T>
unsigned nElements(const Mat<T>& A) { return A.nElements(); }

template <class T>
unsigned getmaxrows(const Mat<T>& A) { return A.getmaxrows(); }

template <class T>
unsigned get_maxcols(const Mat<T>& A) { return A.getmaxcols(); }

template <class T>
int isvector(const Mat<T>& A) { return A.isvector(); }

template <class T>
int iscolumnvector(const Mat<T>& A) { return A.iscolumnvector(); }

template <class T>
int isrowvector(const Mat<T>& A) { return A.isrowvector(); }

template <class T>
unsigned length(const Mat<T>& A) { return A.length(); }

template <class T>
std::ostream& display(std::ostream& os, const Mat<T>& A) { return A.display(os); }

template <class T>
std::ostream& display(std::ostream& os, const Mat<T> A,unsigned r1, unsigned r2, 
		 unsigned c1, unsigned c2)
{ 
  return A.display(os, r1, r2, c1, c2); 
}
  
template <class T>
Mat<T>  t(const Mat<T>& A) { return A.t(); }

template <class T>
Mat<T>  rotate180(const Mat<T>& A) { return A.rotate180(); }

template <class T>
Mat<T>  h(const Mat<T>& A) { return A.h(); }
  
template <class T>
double sum(const Mat<T>& A) { return A.sum(); }

template <class T>
dcomplex csum(const Mat<T>& A) { return A.csum(); }

template <class T>
double sum2(const Mat<T>& A) { return A.sum2(); }

template <class T>
dcomplex csum2(const Mat<T>& A) { return A.csum2(); }

template <class T>
double mean(const Mat<T>& A) { return A.mean(); }

template <class T>
dcomplex cmean(const Mat<T>& A) { return A.cmean(); }

template <class T>
double trace(const Mat<T>& A) { return A.trace(); }

template <class T>
dcomplex ctrace(const Mat<T>& A) { return A.ctrace(); }

template <class T>
double norm(const Mat<T>& A) { return A.norm(); }

template <class T>
dcomplex cnorm(const Mat<T>& A) { return A.cnorm(); }

template <class T>
double det(const Mat<T>& A) { return A.det(); }

template <class T>
dcomplex cdet(const Mat<T>& A) { return A.cdet(); }

template <class T>
Mat<T> exp(const Mat<T>& A) { return A.exp(); }

template <class T>
Mat<T> sqrt(const Mat<T>& A) { return A.sqrt(); }

template <class T>
Mat<T> rowhouse(const Mat<T>& A, const Mat<T>& V)
{
  return A.rowhouse(V);
}

template <class T>
Mat<T> convolv2d(const Mat<T>& A, const Mat<T>& filter) 
{
  return A.convolv2d(filter);
}
  
template <class T>
void remake(Mat<T>& A, unsigned nrows, unsigned ncols)
{
  A(Mat<T>(nrows, ncols));
}
  
#ifdef USE_COMPMAT
#ifdef USE_DBLMAT
extern Mat<double> applyElementWiseC2D(const Mat<dcomplex>& A,
				       double (*function)(const dcomplex&));
extern Mat<double> arg(const Mat<dcomplex>& A);
extern Mat<double> real(const Mat<dcomplex>& A);
extern Mat<double> imag(const Mat<dcomplex>& A);
#endif
#endif
#ifdef USE_FCOMPMAT
#ifdef USE_FLMAT
extern Mat<float> applyElementWiseC2D(const Mat<fcomplex>& A,
					double (*function)(const dcomplex&));
extern Mat<float> arg(const Mat<fcomplex>& A);
extern Mat<float> real(const Mat<fcomplex>& A);
extern Mat<float> imag(const Mat<fcomplex>& A);
#endif
#endif

template <class T> 
Mat<T>& operator += (double addend, Mat<T>& A) { return A+=addend; }

template <class T> 
Mat<T>  operator +  (double addend, const Mat<T>& A) { return A+addend; }

template <class T> 
Mat<T>& operator -= (double subend, Mat<T>& A) { return A-=subend; }

template <class T> 
Mat<T>  operator -  (double subend, const Mat<T>& A) { return A-subend; }

template <class T> 
Mat<T>& operator *= (double factor,  Mat<T>& A) { return A*=factor; }

template <class T> 
Mat<T>  operator *  (double factor,  const Mat<T>& A) { return A*factor; }

template <class T> 
Mat<T>& operator /= (double factor,  Mat<T>& A) { return A/=factor; }

template <class T> 
Mat<T>  operator /  (double factor,  const Mat<T>& A) { return A/factor; }

template <class T> 
Mat<T> inv(const Mat<T>& A) { return A.inv(); }

template <class T> 
T min(const Mat<T>& A, unsigned *row = 0, unsigned *col = 0)
{ 
  return A.min(row, col);
}

template <class T> 
T max(const Mat<T>& A, unsigned *row = 0, unsigned *col = 0) 
{ 
  return A.max(row, col);
}

template <class T> T 
median(const Mat<T>& A, T minVal, T maxVal) 
{ 
  return A.median(minVal, maxVal);
}

template <class T> 
Mat<T> applyElementWise(const Mat<T>& A, double (*function)(double))
{ 
  return A.applyElementWise(function);
}

template <class T> 
Mat<T> applyElementWiseC2C(const Mat<T>& A, dcomplex (*function)(dcomplex))
{
  return A.applyElementWiseC2C(function);
}

template <class T> 
Mat<T> log(const Mat<T>& A) { return A.log(); }

template <class T> 
Mat<T> cos(const Mat<T>& A) { return A.cos(); }

template <class T>
Mat<T> sin(const Mat<T>& A) { return A.sin(); }

template <class T>
Mat<T> round(const Mat<T>& A) { return A.round(); }

template <class T>
Mat<T> abs(const Mat<T>& A) { return A.abs(); }

template <class T> 
Mat<T> pow(const Mat<T>& A, double exp) { return A.pow(exp); }

template <class T>
void eig(const Mat<T>& A, Mat<T>& D, Mat<T>& V) { A.eig(D, V); }

template <class T> 
Mat<T> house(const Mat<T>& A) { return A.house(); }

template <class T> 
Histogram histogram(const Mat<T>& A, double minin = 0, double maxin = 0, unsigned n = 0) { return A.histogram(minin, maxin, n); }

template <class T>
void qr(const Mat<T>& A, Mat<T>& R, Mat<T>& Q) { A.qr(R, Q); }

template <class T>
void svd(const Mat<T>& A, Mat<T>& U, Mat<T>& S, Mat<T>& V) { A.svd(U,S,V); }
  
template <class T>
Mat<T> clip(const Mat<T>& A, T minVal, T maxVal, T minFill, T maxFill) 
{
  return A.clipConst(minVal, maxVal, minFill, maxFill);
}

template <class T>
Mat<T> map(const Mat<T>& A, const ValueMap& valueMap) 
{ 
  return A.mapConst(valueMap); 
}

template <class T>
Mat<T> scale(const Mat<T>& A, double minout = 0.0, double maxout = 255.0, 
	     double minin = 0.0, double maxin=0.0) 
{ 
  return A.scaleConst(minout, maxout, minin, maxin);
}
  
template <class T>
void linearinterpolate(const Mat<T>& x, unsigned Urow, unsigned Drow, 
		       unsigned Ucol, unsigned Dcol, Mat<T>& y)
{
  x.linearinterpolate(Urow,Drow,Ucol,Dcol,y); 
}
  
template <class T>
Mat<T> filter(Mat<T>& B, Mat<T>& A, const Mat<T>& x) 
{ 
  return (x.filter(B, A));
}
  
#ifdef USE_DBLMAT
template <class T>
Mat<T> erode(const Mat<T>& A, const Mat<double>& strel) 
{
  return A.erode(strel);
}

template <class T>
Mat<T> dilate(const Mat<T>& A, const Mat<double>& strel) 
{
  return A.dilate(strel);
}

template <class T> 
Mat<T> open(const Mat<T>& A, const Mat<double>& strel)
{
  return A.open(strel);
}

template <class T> 
Mat<T> close(const Mat<T>& A, const Mat<double>& strel)
{
  return A.close(strel);
}
#endif // USE_DBLMAT

template <class T>
T scalar(const Mat<T>& A) { return A.scalar(); }
  
template <class T> 
void qr(const Mat<T>& A, Mat<T>& R, Mat<T>& Q, Mat<T>& P) { A.qr(R,Q,P); }
  
  
template <class T>
Boolean loadRaw(Mat<T>& A, const char *filename,
		unsigned nrows = 0, unsigned ncols = 0) 
{
  return A.loadRaw(filename, nrows, ncols);
}

template <class T>
Boolean loadAscii(Mat<T>& A, const char *filename)
{ 
  return A.loadAscii(filename);
}
  
#ifdef HAVE_MATLAB   
template <class T>
Boolean loadMatlab(Mat<T>& A, const char *filename, const char *varname = "A")
{
  return A.loadMatlab(filename, varname);
}
  
template <class T>
Boolean saveMatlab(const Mat<T>& A, const char *filename, 
		   const char *varname, const char *option = "u")
{
  return A.saveMatlab(filename, varname, option);
}
#endif // HAVE_MATLAB
  
template <class T>
Boolean saveRaw(const Mat<T>& A, const char *filename) 
{ 
  return A.saveRaw(filename);
}

template <class T>
Boolean saveAscii(const Mat<T>& A, const char *filename)
{
  return A.saveAscii(filename);
}

#ifdef USE_COMPMAT
template <class T>
Mat<dcomplex> fft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
template <class T>
Mat<dcomplex> ifft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
#endif /* USE_COMPMAT */

#ifdef USE_FCOMPMAT
template <class T>
Mat<fcomplex> ffft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
template <class T>
Mat<fcomplex> fifft(const Mat<T>& A, unsigned nrows = 0, unsigned ncols = 0);
#endif /* USE_FCOMPMAT */

template <class T> std::ostream& operator << (std::ostream&, const Mat<T>&);

/*******************************Sub classes************************************/

template <class Type>
class Ones : public Mat<Type> {
public: 
  Ones(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols, (Type) 1) {}
  Ones(const Mat<Type>& A) : Mat<Type>(A.getrows(), A.getcols(), (Type) 1) {}
  Ones(unsigned n) : Mat<Type>(n, n, (Type) 1) {}
  ~Ones() {}
};


template <class Type>
class Zeros : public Mat<Type> {
public:
  Zeros(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols) {}
  Zeros(const Mat<Type>& A) : Mat<Type>(A.getrows(), A.getcols()) {}  
  Zeros(unsigned n)  : Mat<Type>(n, n) {} 
  ~Zeros() {}
}; 


template <class Type>
class Eye : public Mat<Type> {
public:
  Eye(unsigned nrows, unsigned ncols) : Mat<Type>(nrows, ncols) { this->eye(); }
  Eye(Mat<Type>& A) :  Mat<Type>(A.getrows(), A.getcols()) { this->eye(); }
  Eye(unsigned n) : Mat<Type>(n, n)  { this->eye(); }
  ~Eye() {}
};


template <class Type>
class Randuniform : public Mat<Type> {
public:
  Randuniform(unsigned nrows, unsigned ncols, double min = 0, double max = 1) 
    : Mat<Type>(nrows, ncols) { this->randuniform(min, max); }
  Randuniform(const Mat<Type>& A, double min = 0, double max = 1) 
    : Mat<Type>(A.getrows(), A.getcols()) { this->randuniform(min, max); }  
  Randuniform(unsigned n, double min = 0, double max = 1) : Mat<Type>(n, n) { 
      this->randuniform(min, max); } 
  ~Randuniform() {}
};


template <class Type>
class Randnormal : public Mat<Type> {
public:
  Randnormal(unsigned nrows, unsigned ncols, double mean = 0, double std = 1) 
    : Mat<Type>(nrows, ncols) { this->randnormal(mean, std); }
  Randnormal(const Mat<Type>& A, double mean = 0, double std = 1) 
    : Mat<Type>(A.getrows(), A.getcols()) { this->randnormal(mean, std); }  
  Randnormal(unsigned n, double mean = 0, double std = 1)  
    : Mat<Type>(n, n) { this->randnormal(mean, std); } 
  ~Randnormal() {}
};

template <class Type>
class  Hamming: public Mat<Type> {
public:
  Hamming(unsigned n)  : Mat<Type>(n, 1) { this->hamming(); } 
  ~Hamming() {}
};

template <class Type>
class  Hanning: public Mat<Type> {
public:
   Hanning(unsigned n)  : Mat<Type>(n, 1) { this->hanning(); } 
  ~Hanning() {}
};

template <class Type>
class  Blackman: public Mat<Type> {
public:
   Blackman(unsigned n)  : Mat<Type>(n, 1) { this->blackman(); } 
  ~Blackman() {}
};


/******************************************************************************
 * MATRIX ITERATOR CLASSES
 *
 * Quick traversal of a matrix. 
 * IMPORTANT NOTE: For efficiency reasons, no range cheking is performed; I.e.,
 * the user is responsible for staying within the matrix.
 * 
 ******************************************************************************/

template<class Type>
class MatrixIterator {
  Mat<Type>& _mat;   // matrix to be iterated
  Type      *_elPtr; // Current location in matrix

public:
  MatrixIterator(Mat<Type>& mat) : _mat(mat) { first(); }

  // Reset to first element
  Type& first() { return reset(); }
  // Reset to last element
  Type& last()  { return reset(this->_rows - 1, this->_cols - 1); } 
  // Reset to element (r, c)
  Type& reset(unsigned r=0, unsigned c=0) { return *(_elPtr = (_mat._el)[r] + c); } 

  Type& operator ++()    { return *++_elPtr; } // Prefix increment
  Type& operator ++(int) { return *_elPtr++; } // Postfix increment
  Type& operator --()    { return *--_elPtr; } // Prefix decrement
  Type& operator --(int) { return *_elPtr--; } // Postfix decrement
};

template<class Type>
class ConstMatrixIterator {
  const Mat<Type>& _mat;   // matrix to be iterated
  const Type      *_elPtr; // Current location in matrix

public:
  ConstMatrixIterator(const Mat<Type>& mat) : _mat(mat) { first(); }

  // Reset to first element
  const Type& first() { return reset(); }
  // Reset to last element
  const Type& last()  { return reset(this->_rows - 1, this->_cols - 1); } 
  // Reset to element (r, c)
  const Type& reset(unsigned r=0, unsigned c=0) { return *(_elPtr = (_mat._el)[r] + c); }

  const Type& operator ++()    { return *++_elPtr; } // Prefix increment
  const Type& operator ++(int) { return *_elPtr++; } // Postfix increment
  const Type& operator --()    { return *--_elPtr; } // Prefix decrement
  const Type& operator --(int) { return *_elPtr--; } // Postfix decrement
};

template<class Type>
SimpleArray<Type> array(const Mat<Type>& A, Type minVal = 0, Type maxVal = 0);

#endif
