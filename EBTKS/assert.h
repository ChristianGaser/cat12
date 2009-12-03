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
$RCSfile: assert.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
// Assert redefinition to make sure that it exits with a non-zero status code.

#include <stdio.h>
#include <stdlib.h>

#ifdef assert
#undef assert
#else
#define __ASSERT_H__
#endif

/* Warning this definition of assert is not ANSI standard */
/*  In particular it does not return (void)0 on success */
#define assert(EX)  if(!(EX)) { \
  fprintf(stderr, "Assertion failed at line %u in file %s\n", \
          __LINE__, __FILE__); ::exit(EXIT_FAILURE); }





