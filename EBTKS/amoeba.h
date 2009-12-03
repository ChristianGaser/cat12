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
$RCSfile: amoeba.h,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001/11/09 16:37:25 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef  _DEF_AMOEBA_H
#define  _DEF_AMOEBA_H

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "MTypes.h"

typedef  double    (*amoeba_function) ( void *, float [] );

typedef  struct
{
    int               n_parameters;
    float             **parameters;
    double              *values;
    amoeba_function   function;
    void              *function_data;
    double              tolerance;
    double              *sum;
    int               n_steps_no_improvement;
} amoeba_struct;

Boolean numerically_close(double  n1, double  n2, double  threshold_ratio );
double  get_function_value(amoeba_struct  *amoeba, float parameters[] );
void    initialize_amoeba(amoeba_struct     *amoeba,
			  int               n_parameters,
			  double              initial_parameters[],
			  double              parameter_deltas[],
			  amoeba_function   function,
			  void              *function_data,
			  double              tolerance );
double  get_amoeba_parameters(amoeba_struct *amoeba, double parameters[] );
void    terminate_amoeba(amoeba_struct  *amoeba );
double  try_amoeba(amoeba_struct  *amoeba,
		   double          sum[],
		   int             high,
		   double          fac );
Boolean  perform_amoeba(amoeba_struct  *amoeba );

#endif
