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
$RCSfile: FileIO.h,v $
$Revision: 1.3 $
$Author: stever $
$Date: 2003/11/17 04:07:51 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef FILE_IO_H
#define FILE_IO_H

#include <stdio.h>
#include <iostream>		/* (bert) changed from iostream.h */
#include <fstream>		/* (bert) changed from fstream.h */
#include "Path.h"
//#include "popen.h"

#ifdef HAVE_MATLAB
extern "C"{    
  #include"mat.h" 
}
#endif

class InputFile {
  //ipopen *_ipipe;
    std::istream *_ipipe;
  
public:
// Constructors/destructor
  InputFile();
  InputFile(const Path& path) { _ipipe = 0; attach(path); }
  InputFile(std::istream *pipe)   { _ipipe = pipe; }
  ~InputFile()                { close(); }

  operator void *() const;
  operator std::istream& () { return *_ipipe; }

// Set functions
  Boolean  attach(const Path&);
  std::istream& skip(unsigned nBytes);
  Boolean  close();

// Get functions
  std::istream& stream()     { return *_ipipe; }
  Boolean operator ! () { return (!_ipipe || !*_ipipe); }

// Input operators
  std::istream& operator >> (char *val) { return *_ipipe >> val; }
  std::istream& operator >> (char& val) { return *_ipipe >> val; }
  std::istream& operator >> (short& val) { return *_ipipe >> val; }
  std::istream& operator >> (int& val) { return *_ipipe >> val; }
  std::istream& operator >> (long& val) { return *_ipipe >> val; }
  std::istream& operator >> (float& val) { return *_ipipe >> val; }
  std::istream& operator >> (double& val) { return *_ipipe >> val; }
  std::istream& operator >> (unsigned char *val) { return *_ipipe >> val; }
  std::istream& operator >> (unsigned char& val) { return *_ipipe >> val; }
  std::istream& operator >> (unsigned short& val) { return *_ipipe >> val; }
  std::istream& operator >> (unsigned int& val) { return *_ipipe >> val; }
  std::istream& operator >> (unsigned long& val) { return *_ipipe >> val; }
  std::istream& operator >> (std::streambuf *val) { return *_ipipe >> val; }
  std::istream& operator >> (std::istream& (*func)(std::istream&)) { return *_ipipe >> func; }
  std::istream& operator >> (std::ios& (*func)(std::ios&)) { return *_ipipe >> func; }
};

/******************
 * OutputFile class
 ******************/

class OutputFile : public std::ofstream {
  Path _path;
  int  _compress;

public:
  static const Boolean NO_COMPRESS, COMPRESS;

  OutputFile(const Path&, int mode = std::ios::out, int compress = COMPRESS);
  ~OutputFile();

// Set functions
  void compress(Boolean);
};

/***************************
 * C file handling functions
 ***************************/

FILE    *openFile(const Path&, const char *);
void     closeFile();

#ifdef HAVE_MATLAB
MATFile *openMatlabFile(const Path&, const char *);
void     closeMatlabFile(MATFile *);
#endif

extern int get_temp_filename(char *);

#endif
