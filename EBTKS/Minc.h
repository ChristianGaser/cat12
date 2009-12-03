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
$RCSfile: Minc.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
// Some things to prevent clashes between the MNI volume_io library
// and X/C++

#ifndef _MIDAS_MINC_H
#define _MIDAS_MINC_H

#undef Status

#ifdef ROUND
#undef ROUND
#endif

#ifdef SIGN
#undef SIGN
#endif

#ifdef __GNUC__
// g++ barfs if these are included in an extern "C" declaration (which is the case
// here, as they are included in volume_io.h)
  #include <string.h>
#endif
extern "C" {
  #include <volume_io.h>
}

#ifdef public
#undef public
#undef private
#endif

#ifdef ROUND
#undef ROUND
#endif

#ifdef SIGN
#undef SIGN
#endif

#endif
