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
$RCSfile: Path.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef PATH_H
#define PATH_H

#include "MString.h"

class Path : public MString {
  static const char *_separator;
  static const char *_imageNumberWildCard;

public:
// Constructors
  Path() : MString() {}
  Path(const Path& path) : MString(path) {}
  Path(const char *charPtr) : MString(charPtr) {}
  Path(const MString& string) : MString(string) {}
  Path(const MString&, const MString&);

// Get functions
  MString *dir() const;
  MString *file() const;

// Operators
  Boolean operator == (const Path& path) const { // Compare paths (including expansion)
    return (strcmp(expanded().string(), path.expanded().string()) == 0); }

// Special functions
  Path     expanded() const; // Expand ~ to full path
  Boolean  exists() const;   // Check if a path exists
  // Check if a path exists, possibly compressed, and return the correct extension
  Boolean  existsCompressed(MString *extension = 0) const;
  Boolean  hasCompressedExtension(MString *extension = 0) const;
  Path&    removeCompressedExtension();
  Path     removeCompressedExtension() const {
    Path strippedPath(*this); return strippedPath.removeCompressedExtension(); }
  Boolean  removeExtension(MString *extension = 0);
  Boolean isWritable() const;
  // Force template tokens onto path; retain if the template token is "*"
  void     applyTemplate(const Path&, const char *separator = ".");

// The following functions search and/or replace tokens in the Path's file part

  Boolean imageNumber(unsigned *) const; // Find the image number in the filename
  Boolean imageNumberWildCard() const {  // Search for the image number wildcard
    return _locateStringAfterSeparator(_imageNumberWildCard); }
  Boolean replaceImageNumber(unsigned);  // Replace the number after a separator
  Boolean replaceImageNumber(const char *);
  Boolean replaceImageNumberWildCard(unsigned);
  Boolean replaceImageNumberWildCard(const char *);

private:
  Boolean _locateStringAfterSeparator(const char *) const;
};

#endif
