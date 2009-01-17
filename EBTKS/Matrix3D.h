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
$RCSfile: Matrix3D.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MATRIX3D_H
#define _MATRIX3D_H

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
 * USE_UCHRMAT  (unsigned char)   UChrMat3D
 * USE_CHRMAT   (char)            ChrMat3D
 * USE_USHMAT   (unsigned short)  UShMat3D
 * USE_SHMAT    (short)           ShMat3D
 * USE_UINTMAT  (unsigned int)    UIntMat3D
 * USE_INTMAT   (int)             IntMat3D
 * USE_FLMAT    (float)           FlMat3D
 * USE_DBLMAT   (double)          DblMat3D
 * USE_COMPMAT  (complex)         CompMat3D
 * USE_FCOMPMAT (fcomplex)        fCompMat3D
 *****************************************************************************/

#include "Matrix.h"

template <class Type> class Ones3D;
template <class Type> class Zeros3D;
template <class Type> class Randuniform3D;
template <class Type> class Randnormal3D;

typedef double  (*IndexFunction3D)(unsigned s, unsigned r, unsigned c);
typedef complex (*ComplexIndexFunction3D)(unsigned s, unsigned r, unsigned c);

/******************** Explicit type conversions *******************/

// Forward declaration needed for type conversions
template <class Type> class Mat3D;

#ifdef USE_DBLMAT
//Converts Mat3D<Type> to Mat3D<double>
typedef Mat3D<double> DblMat3D;
template <class Type> Mat3D<double> asDblMat(const Mat3D<Type>&);
#endif

#ifdef USE_FLMAT
//Converts Mat3D<Type> to Mat3D<float>
typedef Mat3D<float> FlMat3D;
template <class Type> Mat3D<float> asFlMat(const Mat3D<Type>&);
#endif

#ifdef USE_INTMAT
//Converts Mat3D<Type> to Mat3D<int>
typedef Mat3D<int> IntMat3D;
template <class Type> Mat3D<int> asIntMat(const Mat3D<Type>&);
#endif

#ifdef USE_UINTMAT
//Converts Mat3D<Type> to Mat3D<unsigned int>
typedef Mat3D<unsigned int> UIntMat3D;
template <class Type> Mat3D<unsigned int> asUIntMat(const Mat3D<Type>&);
#endif

#ifdef USE_SHMAT
//Converts Mat3D<Type> to Mat3D<short>
typedef Mat3D<short> ShMat3D;
template <class Type> Mat3D<short> asShMat(const Mat3D<Type>&);
#endif

#ifdef USE_USHMAT
//Converts Mat3D<Type> to Mat3D<unsigned short>
typedef Mat3D<unsigned short> UShMat3D;
template <class Type> Mat3D<unsigned short> asUShMat(const Mat3D<Type>&);
#endif

#ifdef USE_CHRMAT
//Converts Mat3D<Type> to Mat3D<char>
typedef Mat3D<char> ChrMat3D;
template <class Type> Mat3D<char> asChrMat(const Mat3D<Type>&);
#endif

#ifdef USE_UCHRMAT
//Converts Mat3D<Type> to Mat3D<unsigned char>
typedef Mat3D<unsigned char> UChrMat3D;
template <class Type> Mat3D<unsigned char> asUChrMat(const Mat3D<Type>&);
#endif

#ifdef USE_COMPMAT
typedef Mat3D<complex> CompMat3D;
template <class Type> Mat3D<complex> asCompMat(const Mat3D<Type>&);
template <class Type> Mat3D<complex> asCompMat(const Mat3D<Type>& Re, const Mat3D<Type>& Im);
#endif

#ifdef USE_FCOMPMAT
typedef Mat3D<fcomplex> fCompMat3D;
template <class Type> Mat3D<fcomplex> asFcompMat(const Mat3D<Type>&);
template <class Type> Mat3D<fcomplex> asFcompMat(const Mat3D<Type>& Re, const Mat3D<Type>& Im);
#endif

// Some mathematical operations. These had to be implemented as template
// functions in order to allow different types of matrices to be operated on
// ***************************** Addition ****************************
template <class T1, class T2>
Mat3D<T1>& operator += (Mat3D<T1>& A, const Mat3D<T2>& B);

template <class T1, class T2>
Mat3D<T1> operator + (const Mat3D<T1>& A, const Mat3D<T2>& B) {
  Mat3D<T1> T(A); return T += B; }

// ***************************** Substraction *************************
template <class T1, class T2>
Mat3D<T1>& operator -= (Mat3D<T1>& A, const Mat3D<T2>& B);

template <class T1, class T2>
Mat3D<T1> operator - (const Mat3D<T1>& A, const Mat3D<T2>& B) {
  Mat3D<T1> T(A); return T -= B; }

// ***************************** Multiplication ***********************
template <class T1, class T2>
Mat3D<T1> operator * (const Mat3D<T1>& A, const Mat3D<T2>& B);

template <class T1, class T2>
Mat3D<T1>& operator *= (Mat3D<T1>& A, const Mat3D<T2>& B) {
  A = A * B; return A; }

// ***************************** Division ***** ***********************
template <class T1, class T2>
Mat3D<T1> operator / (const Mat3D<T1>& A, const Mat3D<T2>& B) {
  return A * inv(B); }

template <class T1, class T2>
Mat3D<T1>& operator /= (Mat3D<T1>& A, const Mat3D<T2>& B) {
  A = A * inv(B); return A; }

// *********************** Point multiplication ***********************
template <class T1, class T2>
Mat3D<T1>& pmultEquals(Mat3D<T1>& A, const Mat3D<T2>& B);

template <class T1, class T2>
Mat3D<T1> pmult(const Mat3D<T1>& A, const Mat3D<T2>& B) {
  Mat3D<T1> T(A); return pmultEquals(T, B); }

// ************************* Point division ***************************
template <class T1, class T2>
Mat3D<T1>& pdivEquals(Mat3D<T1>& A, const Mat3D<T2>& B);

template <class T1, class T2>
Mat3D<T1> pdiv(const Mat3D<T1>& A, const Mat3D<T2>& B) {
  Mat3D<T1> T(A); return pdivEquals(T, B); }

/******************************************************************************
 *     3D MATRIX CLASS
 ******************************************************************************/

template <class Type> class Mat3D {
  static unsigned _rangeErrorCount;

protected:
  unsigned     _slis;
  unsigned     _rows;
  unsigned     _cols;
  Mat<Type>   *_slices;
  Type	    ***_el;
  
public:
  static Boolean flushToDisk;

  Mat3D() {_slis = 0; _rows = 0; _cols = 0; _el = 0; _slices = 0; }
  Mat3D(const Mat3D& A);  //copy constructor
  Mat3D(unsigned nslis, unsigned nrows, unsigned ncols, Type value);
  Mat3D(unsigned nslis, unsigned nrows, unsigned ncols);
  Mat3D(unsigned nslis, unsigned nrows, unsigned ncols, const Type *data);
  Mat3D(const char *filename, int type = RAW);
  virtual ~Mat3D();
  
  void clear();
  
  Boolean operator ! () const { 
    return Boolean(!_slis || !_rows || !_cols || !_el || !_slices); }
  operator void * () const    { 
    return (void *) (_slis && _rows && _cols && _el && _slices); }

  unsigned  getslis() const   { return _slis; }
  unsigned  getrows() const   { return _rows; }
  unsigned  getcols() const   { return _cols; }
  unsigned  nElements() const { return _slis*_rows*_cols; }

  // Be careful with these two
  const Type ***getEl() const           { return (const Type ***) _el; }
  Mat<Type>&    operator[] (unsigned s) { return _slices[s]; }

  const Mat<Type>& operator[] (unsigned s) const { return _slices[s]; }
  Mat3D&           setSlice(unsigned slice, const Mat<Type>& A);
  
  ostream& display(ostream& os) const;
  ostream& display(ostream& os, unsigned s1, unsigned s2, unsigned r1, unsigned r2, 
		   unsigned c1, unsigned c2) const;

  // Assignment operator
  Mat3D<Type>& operator = (const Mat3D& A);
  // Same, but destroys A
  Mat3D<Type>& absorb(Mat3D& A);

  Type   operator() (unsigned s, unsigned r, unsigned c) const;
  Mat3D  operator() (unsigned s1, unsigned s2, unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;
  Type   operator() (unsigned n) const;
  Type&  operator() (unsigned s, unsigned r, unsigned c);
  Type&  operator() (unsigned n);
  Mat3D& operator() (const Mat3D&);
  
  int    operator != (const Mat3D&) const;
  int    operator == (const Mat3D& A) const { return !(operator != (A)); }

  Mat3D& operator += (complex);
  Mat3D  operator +  (complex x) const       { Mat3D<Type> A(*this); return A += x; }
  Mat3D& operator -= (complex x)             { return *this += (-x); }
  Mat3D  operator -  (complex x) const       { return *this +  (-x); }
  Mat3D& operator *= (complex);
  Mat3D  operator *  (complex x) const       { Mat3D<Type> A(*this); return A *= x; }
  Mat3D& operator /= (complex x)             { return (*this) *= (1.0/x); }
  Mat3D  operator /  (complex x) const       { return (*this) *  (1.0/x); }
  
  int isvector() const       { return (_slis == 1) || (_rows == 1) || (_cols == 1); }
  int iscolumnvector() const { return (_cols == 1); }
  int isrowvector() const    { return (_rows == 1); }
  int isslicevector() const  { return (_slis == 1); }

  unsigned length() const;
  
  Type min(unsigned *sli = 0, unsigned *row = 0, unsigned *col = 0) const;
  Type max(unsigned *sli = 0, unsigned *row = 0, unsigned *col = 0) const;
  Type median(Type minVal = 0, Type maxVal = 0) const;
  
  double  sum() const { return real(csum()); }
  complex csum() const;
  double  sum2() const { return real(csum2()); }
  complex csum2() const;
  double  mean() const  { return sum()/nElements(); }
  complex cmean() const { return csum()/nElements(); }
  double  std() const  { return ::sqrt(var()); }
  complex cstd() const { return ::sqrt(cvar()); }
  double  var() const  { double  mn = mean(); return sum2()/nElements() - SQR(mn); }
  complex cvar() const { complex mn = cmean(); return csum2()/nElements() - SQR(mn); }
  double  norm() const  { return ::sqrt(sum2()); }
  complex cnorm() const { return ::sqrt(csum2()); }
  double  trace() const { return real(ctrace()); }
  complex ctrace() const;
  double  det() const { return real(cdet()); }
  complex cdet() const;
  
  /**********************Other Matrix element wise functions*********************/
  //Applies a given single argument function to each of the elements
  Mat3D&  applyElementWise(double (*function)(double));
  Mat3D   applyElementWise(double (*function)(double)) const {
    return applyElementWiseConst(function); }
  Mat3D   applyElementWiseConst(double (*function)(double)) const {
    Mat3D<Type> A(*this); return A.applyElementWise(function); }

  Mat3D&  applyElementWiseC2D(double (*function)(complex));
  Mat3D   applyElementWiseC2D(double (*function)(complex)) const {
    return applyElementWiseConstC2D(function); }
  Mat3D   applyElementWiseConstC2D(double (*function)(complex)) const {
    Mat3D<Type> A(*this); return A.applyElementWiseC2D(function); }

  Mat3D&  applyElementWiseC2C(complex (*function)(complex));
  Mat3D   applyElementWiseC2C(complex (*function)(complex)) const {
    return applyElementWiseConstC2C(function); }
  Mat3D   applyElementWiseConstC2C(complex (*function)(complex)) const {
    Mat3D<Type> A(*this); return A.applyElementWiseC2C(function); }

  Mat3D& applyIndexFunction(IndexFunction3D F);
  Mat3D& applyIndexFunction(IndexFunction3D F) const {
    return applyIndexFunctionConst(F); }
  Mat3D   applyIndexFunctionConst(IndexFunction3D F) const {
    Mat3D<Type> A(*this); return A.applyIndexFunction(F); }

  Mat3D& applyIndexFunction(ComplexIndexFunction3D F);
  Mat3D& applyIndexFunction(ComplexIndexFunction3D F) const {
    return applyComplexIndexFunctionConst(F); }
  Mat3D   applyIndexFunctionConst(ComplexIndexFunction3D F) const {
    Mat3D<Type> A(*this); return A.applyIndexFunction(F); }

  //Takes the inverse natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat3D& exp();
  Mat3D  exp() const      { return expConst(); }
  Mat3D  expConst() const { Mat3D<Type> T(*this); return T.exp(); }
  
  //Takes the natural log of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat3D& log();
  Mat3D  log() const      { return logConst(); }
  Mat3D  logConst() const { Mat3D<Type> T(*this); return T.log(); }
  
  //Takes the cosine of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat3D& cos();
  Mat3D  cos() const      { return cosConst(); }
  Mat3D  cosConst() const { Mat3D<Type> T(*this); return T.cos(); }
  
  //Takes the sine of each of the elements of the matrix
  //and returns a copy of the matrix.
  Mat3D& sin();
  Mat3D  sin() const      { return sinConst(); }
  Mat3D  sinConst() const { Mat3D<Type> T(*this); return T.sin(); }

  //Takes the absolute value of each of the elements of the matrix
  Mat3D& abs();
  Mat3D  abs() const      { return absConst(); }
  Mat3D  absConst() const { Mat3D<Type> T(*this); return T.abs(); }
  
  //Takes the absolute value of each of the elements of the matrix
  Mat3D& conj();
  Mat3D  conj() const      { return conjConst(); }
  Mat3D  conjConst() const { Mat3D<Type> T(*this); return T.conj(); }

  //Rounds each of the elements of the matrix
  Mat3D& round();
  Mat3D  round() const      { return roundConst(); }
  Mat3D  roundConst() const { Mat3D<Type> T(*this); return T.round(); }

  //Takes the sqrt value of each of the elements of the matrix
  Mat3D& sqrt();
  Mat3D  sqrt() const      { return sqrtConst(); }
  Mat3D  sqrtConst() const { Mat3D<Type> T(*this); return T.sqrt(); }

  //Takes the power value of each of the elements of the matrix
  Mat3D& pow(double exponent);
  Mat3D  pow(double exponent) const      { return powConst(exponent); }
  Mat3D  powConst(double exponent) const { Mat3D<Type> T(*this); return T.pow(exponent);}

  // Padding functions; return matrices padded with <value> all around. 
  // Generic pad function; pads to nslis x nrows x ncols and inserts at slice,row,col.
  Mat3D& pad(unsigned nslis, unsigned nrows, unsigned ncols, 
	   int slice, int row, int col, Type value = 0);
  Mat3D  pad(unsigned nslis, unsigned nrows, unsigned ncols, 
	   int slice, int row, int col, Type value = 0) const {
    return padConst(nslis, nrows, ncols, slice, row, col, value); }
  Mat3D  padConst(unsigned nslis, unsigned nrows, unsigned ncols, 
		int slice, int row, int col, Type value = 0) const {
    Mat3D<Type> A(*this); return A.pad(nslis, nrows, ncols, slice, row, col, value); }
  
  // Symmetrical pad, i.e., the size of the returned matrix is (rows+2*rowpad, 
  // cols+2*colpad).
  Mat3D& pad(unsigned slicepad, unsigned rowpad, unsigned colpad, Type value = 0) {
    return pad(_slis + 2*slicepad, _rows + 2*rowpad, _cols + 2*colpad, 
	       slicepad, rowpad, colpad, value); }
  Mat3D  pad(unsigned slicepad, unsigned rowpad, unsigned colpad, Type value = 0) const {
    return padConst(slicepad, rowpad, colpad, value); }
  Mat3D  padConst(unsigned slicepad, unsigned rowpad, unsigned colpad, Type value=0) const{
    return padConst(_slis + 2*slicepad, _rows + 2*rowpad, _cols + 2*colpad, 
		    slicepad, rowpad, colpad, value); }

  // Sets all values in the matrix below minVal to minFill, and those
  // above maxVal to maxFill.
  Mat3D& clip(Type minVal, Type maxVal, Type minFill, Type maxFill);
  Mat3D  clip(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    return clipConst(minVal, maxVal, minFill, maxFill); }
  Mat3D  clipConst(Type minVal, Type maxVal, Type minFill, Type maxFill) const { 
    Mat3D<Type> A(*this); return A.clip(minVal, maxVal, minFill, maxFill); }

  // Maps the values of a matrix through a value map
  Mat3D& map(const ValueMap& valueMap);
  Mat3D  map(const ValueMap& valueMap) const { return mapConst(valueMap); }
  Mat3D  mapConst(const ValueMap& valueMap) const { 
    Mat3D<Type> A(*this); return A.map(valueMap); }

  // Scales a Matrix (similar to map; kept for backward compatibility)
  Mat3D& scale(double minout = 0.0, double maxout = 255.0, 
	       double minin = 0.0, double maxin = 0.0);
  Mat3D  scale(double minout = 0.0, double maxout = 255.0, 
	       double minin = 0.0, double maxin = 0.0) const {
    return scaleConst(minout, maxout, minin, maxin); }
  Mat3D  scaleConst(double minout = 0.0, double maxout = 255.0, 
		    double minin = 0.0, double maxin = 0.0) const {
    Mat3D<Type> A(*this); return A.scale(minout, maxout, minin, maxin); }
  
  Boolean load(const char *filename, int type = RAW) const;
  Boolean loadRaw(const char *filename, unsigned nslis = 0, unsigned nrows = 0, 
		  unsigned ncols = 0);
  Boolean loadAscii(const char *filename);

  Boolean save(const char *filename, int type = RAW) const;
  Boolean saveRaw(const char *filename) const;
  Boolean saveAscii(const char *filename) const;
  
  Mat3D& insert(const Mat3D<Type>& A, int slice=0, int row=0, int col=0);
  Mat3D& insert(const char *path, unsigned nslis=0, unsigned nrows=0, unsigned ncols=0,
		int slice=0, int row=0, int col=0);
  
  Mat3D& fill(Type value);
  
  Mat3D& randuniform(double min = 0, double max = 1);
  Mat3D& randnormal(double mean = 0, double std = 1);
  
  Mat3D rotate180() const;

#ifdef USE_DBLMAT
  Mat3D erode(const Mat3D<double>& strel) const;
  Mat3D dilate(const Mat3D<double>& strel) const;
  Mat3D open(const Mat3D<double>& strel) const  { return erode(strel).dilate(strel); }
  Mat3D close(const Mat3D<double>& strel) const { return dilate(strel).erode(strel); }
#endif

  //Returns a histogram for the calling object using the specified range and # bins
  Histogram histogram(double minin = 0, double maxin = 0, unsigned n = 0) const;
  
  Mat3D& histmod(const Histogram& hist1, const Histogram& hist2);
  Mat3D  histmod(const Histogram& hist1, const Histogram& hist2) const {
    return histmodConst(hist1, hist2); }
  Mat3D  histmodConst(const Histogram& hist1, const Histogram& hist2) const {
    Mat3D<Type> A(*this); return A.histmod(hist1, hist2); }

  // Returns an array which contains only the values in the range specified
  SimpleArray<Type> asArray(Type minVal = 0, Type maxVal = 0) const;

/*****************************fft functions************************************/

  // For the following fft functions, a dimension of 0 (which is the
  // default argument) will be expanded to the next power of two, equal
  // to or above the matrix dimension. A dimension of 1 indicates that
  // the fft should not be taken along that dimension.

  // Non-const fft function
  Mat3D& fft(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0) {
    return _fft(nslis, nrows, ncols, ::fft); }
  // Const fft function
  Mat3D  fft(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0) const { 
    return fftConst(nslis, nrows, ncols); }
  // Explicit const fft function
  Mat3D  fftConst(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat3D<Type> A(*this); return A.fft(nslis, nrows, ncols); }
  
  // Non-const ifft function
  Mat3D& ifft(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0);
  // Const ifft function
  Mat3D  ifft(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0) const { 
    return ifftConst(nslis, nrows, ncols); }
  // Explicit const ifft function
  Mat3D  ifftConst(unsigned nslis = 0, unsigned nrows = 0, unsigned ncols = 0) const { 
    Mat3D<Type> A(*this); return A.ifft(nslis, nrows, ncols); }
  
private:
  void _allocateEl();
  void _setEl();
  void _checkMatrixDimensions(const char *path, 
			      unsigned& nslis, unsigned& nrows, unsigned& ncols) const;
  Mat3D& _fft(unsigned nslis, unsigned nrows, unsigned ncols, FFTFUNC fftFunc);

// Kept up to here for now
};

///////////////////////////////////////////////////////////////////////////
//  Various non-member template functions
///////////////////////////////////////////////////////////////////////////

template <class Type>
void clear(Mat3D<Type> A) { A.clear(); }
  
template <class Type>
Mat3D<Type>& operator += (double addend, Mat3D<Type>& A)
{ 
  return A += addend;
}

template <class Type>
Mat3D<Type>  operator +  (double addend, const Mat3D<Type>& A)
{ 
  return A + addend;
}

template <class Type>
Mat3D<Type>& operator -= (double subend, Mat3D<Type>& A)
{ 
  return A -= subend;
}

template <class Type>
Mat3D<Type>  operator -  (double subend, const Mat3D<Type>& A)
{ 
  return A - subend;
}

template <class Type>
Mat3D<Type>& operator *= (double factor, Mat3D<Type>& A)
{ 
  return A *= factor;
}
template <class Type>
Mat3D<Type>  operator *  (double factor, const Mat3D<Type>& A)
{ 
  return A * factor;
}

template <class Type>
Mat3D<Type>& operator /= (double factor, Mat3D<Type>& A)
{ 
  return A /= factor;
}

template <class Type>
Mat3D<Type>  operator /  (double factor, const Mat3D<Type>& A) 
{ 
  return A / factor;
}
  
template <class Type>
int isvector(const Mat3D<Type>& A) { return A.isvector(); }

template <class Type>
int iscolumnvector(const Mat3D<Type>& A) { return A.iscolumnvector(); }

template <class Type>
int isrowvector(const Mat3D<Type>& A) { return A.isrowvector(); }

template <class Type>
int isslicevector(const Mat3D<Type>& A) { return A.isslicevector(); }

template <class Type>
unsigned getrows(const Mat3D<Type>& A) { return A.getrows(); }

template <class Type>
unsigned getcols(const Mat3D<Type>& A) { return A.getcols(); }

template <class Type>
unsigned getslis(const Mat3D<Type>& A) { return A.getslis(); }

template <class Type>
unsigned nElements(const Mat3D<Type>& A) { return A.nElements(); }

template <class Type>
unsigned length(const Mat3D<Type>& A) { return A.length(); }

template <class Type>
Type min(const Mat3D<Type>& A, unsigned *sli = 0, unsigned *row = 0, 
	 unsigned *col = 0)
{ 
  return A.min(sli, row, col);
}

template <class Type>
Type max(const Mat3D<Type>& A, unsigned *sli = 0, unsigned *row = 0, 
	 unsigned *col = 0)
{
  return A.max(sli, row, col);
}

template <class Type>
Type median(const Mat3D<Type>& A, Type minVal = 0, Type maxVal = 0) 
{ 
  return A.median(minVal, maxVal);
}

template <class Type>
double sum(const Mat3D<Type>& A)           { return A.sum(); }

template <class Type>
complex csum(const Mat3D<Type>& A)          { return A.csum(); }

template <class Type>
double sum2(const Mat3D<Type>& A)          { return A.sum2(); }

template <class Type>
complex csum2(const Mat3D<Type>& A)  { return A.csum2(); }

template <class Type>
double mean(const Mat3D<Type>& A)    { return A.mean(); }

template <class Type>
complex cmean(const Mat3D<Type>& A)  { return A.cmean(); }

template <class Type>
double trace(const Mat3D<Type>& A)   { return A.trace(); }

template <class Type>
complex ctrace(const Mat3D<Type>& A) { return A.ctrace(); }

template <class Type>
double norm(const Mat3D<Type>& A)    { return A.norm(); }

template <class Type>
complex cnorm(const Mat3D<Type>& A)  { return A.cnorm(); }

template <class Type>
double det(const Mat3D<Type>& A)     { return A.det(); }

template <class Type>
complex cdet(const Mat3D<Type>& A)   { return A.cdet(); }

template <class Type>
Mat3D<Type> applyElementWise(const Mat3D<Type>& A, double (*function)(double))
{
  return A.applyElementWise(function);
}

template <class Type>
Mat3D<Type> applyElementWiseC2C(const Mat3D<Type>& A, 
				complex (*function)(complex))
{
  return A.applyElementWiseC2C(function);
}

#ifdef USE_COMPMAT
#ifdef USE_DBLMAT
template <class Type>
Mat3D<double> applyElementWiseC2D(const Mat3D<complex>& A,
				  double (*function)(const complex&));

template <class Type>
Mat3D<double> arg(const Mat3D<complex>& A);

template <class Type>
Mat3D<double> real(const Mat3D<complex>& A);

template <class Type>
Mat3D<double> imag(const Mat3D<complex>& A);
#endif // USE_DBLMAT
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
#ifdef USE_FLMAT
template <class Type>
Mat3D<float> applyElementWiseC2D(const Mat3D<fcomplex>& A,
				 double (*function)(const complex&));

template <class Type>
Mat3D<float> arg(const Mat3D<fcomplex>& A);

template <class Type>
Mat3D<float> real(const Mat3D<fcomplex>& A);

template <class Type>
Mat3D<float> imag(const Mat3D<fcomplex>& A);

#endif // USE_FLMAT
#endif // USE_FCOMPMAT

template <class Type>
Mat3D<Type> exp(const Mat3D<Type>& A) { return A.exp(); }

template <class Type>
Mat3D<Type> cos(const Mat3D<Type>& A) { return A.cos(); }

template <class Type>
Mat3D<Type> sin(const Mat3D<Type>& A) { return A.sin(); }

template <class Type>
Mat3D<Type> abs(const Mat3D<Type>& A) { return A.abs();}

template <class Type>
Mat3D<Type> round(const Mat3D<Type>& A) { return A.round();}
  
template <class Type>
Mat3D<Type>& pad(Mat3D<Type>& A, unsigned slicepad, unsigned rowpad, 
		 unsigned colpad, Type value = 0)
{
  return A.pad(slicepad, rowpad, colpad, value);
}

template <class Type>
Mat3D<Type> pad(const Mat3D<Type>& A, unsigned slicepad, unsigned rowpad, 
		unsigned colpad, Type value = 0)
{
  return A.padConst(slicepad, rowpad, colpad, value);
}

template <class Type>
Mat3D<Type> padConst(const Mat3D<Type>& A, unsigned slicepad, unsigned rowpad, 
		     unsigned colpad, Type value = 0)
{
  return A.padConst(slicepad, rowpad, colpad, value);
}
  
template <class Type>
Mat3D<Type> rotate180(const Mat3D<Type>& A) { return A.rotate180(); }

#ifdef HAVE_DBLMAT
template <class Type>
Mat3D<Type> erode(const Mat3D<Type>& A, const Mat3D<double>& strel)
{
  return A.erode(strel);
}

template <class Type>
Mat3D<Type> dilate(const Mat3D<Type>& A, const Mat3D<double>& strel)
{
  return A.dilate(strel);
}

template <class Type>
Mat3D<Type> open(const Mat3D<Type>& A, const Mat3D<double>& strel)
{
  return A.open(strel);
}

template <class Type>
Mat3D<Type> close(const Mat3D<Type>& A, const Mat3D<double>& strel)
{
  return A.close(strel);
}
#endif // HAVE_DBLMAT

template <class Type>
Histogram histogram(const Mat3D<Type>& A, double minin = 0, double maxin = 0,
		    unsigned n = 0)
{
  return A.histogram(minin, maxin, n);
}
  
template <class Type>
Mat3D<Type> clip(const Mat3D<Type>& A, Type minVal, Type maxVal, Type minFill,
		 Type maxFill)
{ 
  return A.clipConst(minVal, maxVal, minFill, maxFill);
}

template <class Type>
Mat3D<Type> map(const Mat3D<Type>& A, const ValueMap& valueMap)
{ 
  return A.mapConst(valueMap);
}

template <class Type>
Mat3D<Type> scale(const Mat3D<Type>& A, double minout = 0.0, 
		  double maxout = 255.0, double minin = 0.0, 
		  double maxin = 0.0)
{ 
  return A.scaleConst(minout, maxout, minin, maxin);
}
  
template <class Type>
void loadRaw(Mat3D<Type>& A, char *filename, unsigned nslis = 0, 
	       unsigned nrows = 0, unsigned ncols = 0) 
{ 
  A.loadRaw(filename, nslis, nrows, ncols);
}

template <class Type>
void loadAscii(Mat3D<Type>& A, char *filename) { A.loadAscii(filename); }
  
template <class Type>
void saveRaw(const Mat3D<Type>& A, const char *filename)
{ 
  A.saveRaw(filename);
}

template <class Type>
void saveAscii(const Mat3D<Type>& A, const char *filename)
{ 
  A.saveAscii(filename);
}

#ifdef USE_COMPMAT
template <class Type>
Mat3D<complex> fft(const Mat3D<Type>& A, unsigned nslis = 0,
		   unsigned nrows = 0, unsigned ncols = 0)
{
  return asCompMat(A).fft(nslis, nrows, ncols);
}

template <class Type>
Mat3D<complex> ifft(const Mat3D<Type>& A, unsigned nslis = 0,
		    unsigned nrows = 0, unsigned ncols = 0)
{
    return asCompMat(A).ifft(nslis, nrows, ncols); }
#endif // USE_COMPMAT

#ifdef USE_FCOMPMAT
template <class Type>
Mat3D<fcomplex> ffft(const Mat3D<Type>& A, unsigned nslis = 0,
		     unsigned nrows = 0, unsigned ncols = 0)
{
  return asFcompMat(A).fft(nslis, nrows, ncols);
}

template <class Type>
Mat3D<fcomplex> fifft(const Mat3D<Type>& A, unsigned nslis = 0,
		      unsigned nrows = 0, unsigned ncols = 0)
{
  return asFcompMat(A).ifft(nslis, nrows, ncols);
}
#endif // USE_FCOMPMAT

template <class Type> ostream& operator << (ostream&, const Mat3D<Type>&);

////////////////////////////////////////////////////////////////////////////
//
//Sub classes
//
///////////////////////////////////////////////////////////////////////////
//
//-------
//
template <class Type>
class Ones3D : public Mat3D<Type> {
public: 
  Ones3D(unsigned nslis, unsigned nrows, unsigned ncols) : Mat3D<Type>(nslis, nrows, ncols, (Type) 1) {}
  Ones3D(const Mat3D<Type>& A) : Mat3D<Type>(A.getslis(), A.getrows(), A.getcols(), (Type) 1) {}
  Ones3D(unsigned n) : Mat3D<Type>(n, n, n, (Type) 1) {}
  ~Ones3D() {}
};
//
//-------
//
template <class Type>
class Zeros3D : public Mat3D<Type> {
public:
  Zeros3D(unsigned nslis, unsigned nrows, unsigned ncols) : Mat3D<Type>(nslis, nrows, ncols) {}
  Zeros3D(const Mat3D<Type>& A) : Mat3D<Type>(A.getslis(), A.getrows(), A.getcols()) {}  
  Zeros3D(unsigned n)  : Mat3D<Type>(n, n, n) {} 
  ~Zeros3D() {} 
}; 
//
//-------
//
template <class Type>
class Randuniform3D : public Mat3D<Type> {
public:
  
  Randuniform3D(unsigned nslis, unsigned nrows, unsigned ncols) : Mat3D<Type>(nslis, nrows, ncols) {randuniform();}
  Randuniform3D(const Mat3D<Type>& A) : Mat3D<Type>(A.getslis(), A.getrows(), A.getcols()) {randuniform();}  
  Randuniform3D(unsigned n)  : Mat3D<Type>(n, n, n) {randuniform();} 
  ~Randuniform3D() {}
};
//
//-------
//
template <class Type>
class Randnormal3D : public Mat3D<Type> {
public:
  
  Randnormal3D(unsigned nslis, unsigned nrows, unsigned ncols) : Mat3D<Type>(nslis, nrows, ncols) {randnormal();}
  Randnormal3D(const Mat3D<Type>& A) : Mat3D<Type>(A.getslis(), A.getrows(), A.getcols()) {randnormal();}  
  Randnormal3D(unsigned n)  : Mat3D<Type>(n, n, n) {randnormal();} 
  ~Randnormal3D() {}
};
#endif
