/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, John G. Sled, 
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
$RCSfile: splineSmooth.cc,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2005/03/08 15:55:34 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : splineSmooth.c,v
@INPUT      : 
@OUTPUT     : (none)
@RETURNS    : 
@DESCRIPTION: Tool for smoothing and extrapolating data in minc volumes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : April 21, 1996 (John G. Sled)
@MODIFIED   : Log: splineSmooth.c,v 
 * Revision 1.2  1996/04/23  13:36:58  jgsled
 * Working version.  Problems with thin plate splines have been fixed.
 *
 * Revision 1.1  1996/04/21  23:41:50  jgsled
 * Initial version of SplineSmooth tool
 * - B spline implementation appears to work
 *
@COPYRIGHT  : 1996
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[] = "$Header: /software/source/INSECT/N3/src/SplineSmooth/splineSmooth.cc,v 1.2 2005/03/08 15:55:34 bert Exp $";
#endif

#include <stdio.h>
#include <iostream>		// (bert)
using namespace std;		// (bert)
#include <math.h>
#include <EBTKS/Matrix.h>	// (bert) 
#include <EBTKS/TBSpline.h>	// (bert)
#undef ROUND
#undef SIGN

//---------------------------------------------------------------------------------
//  Implementation notes
/*
  The spline basis functions are defined in a world coordinate system aligned
  with the voxel coordinate system and sharing the same origin.

 */

//---------------------------------------------------------------------------------
// Declarations
DblMat volume_domain(double *separations, int *dims);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, double *src, 
                              const DblMat &domain,
                              int subsample, double *separations, int *dims);
void smoothVolumeLookup(TBSplineVolume *spline, double *src, int *dims);

//--------------------------------------------------------------------------------
// main program
extern "C" int splineSmooth( double *src, double lambda, double distance, int subsample, double *separations, int *dims)
     
{
  DblMat domain;   // region in world coordinates on which splines are defined

  // domain is whole volume
  domain = volume_domain(separations, dims);

  // create spline basis
  Spline *theSplines;
  
  double start[3] = { 0.0, 0.0, 0.0 };
  theSplines = new TBSplineVolume(domain, start, separations, dims,
                                      distance, lambda);
  
  // do least squares fit to data 
  fitSplinesToVolumeLookup((TBSplineVolume *)theSplines, src,
                               domain, subsample, separations, dims);
    
  // write smooth function to volume
  smoothVolumeLookup((TBSplineVolume *) theSplines, src, dims);
  
  return(0);
} 


//-----------------------------------------------------------------------------
// Supporting functions



// determine domain from size of volume in world coordinates
// Returns an 3 by 2 matrix
DblMat
volume_domain(double *separations, int *dims)
{
  DblMat domain(3,2);

  for(int i = 0; i < 3; i++)
    {
      if(separations[i] > 0) {
        domain(i,0) = -0.5*separations[i];
        domain(i,1) = (dims[i]-0.5)*separations[i];
      }
      else {
        domain(i,1) = -0.5*separations[i];
        domain(i,0) = (dims[i]-0.5)*separations[i];
      }
    }
  return domain;
}


void
fitSplinesToVolumeLookup(TBSplineVolume *spline, double* src,
                         const DblMat &domain, int subsample, double* separations, int* dims)
{
  int i,x,y,z,j,k;
  long area, vol, z_area, y_dims;
  double value;

  // only look at values within domain
  int lower[3], upper[3];
  for(i = 0; i < 3; i++)
    {
      if(separations[i] > 0) {
        lower[i] = (int) ceil(domain(i,0)/separations[i]);
        upper[i] = (int) floor(domain(i,1)/separations[i]);
      }
      else {
        upper[i] = (int) floor(domain(i,0)/separations[i]);
        lower[i] = (int) ceil(domain(i,1)/separations[i]);
      }
    }

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  for(z = lower[2]; z <= upper[2]; z += subsample) {
    z_area = z*area;
    for(y = lower[1]; y <= upper[1]; y += subsample) {
      y_dims = y*dims[0];
      for(x = lower[0]; x <= upper[0]; x += subsample)
	  {
        value = src[z_area + y_dims + x];
	    if(value > 0)
		  spline->addDataPoint(x,y,z, value);
	  }
    }
  }

  if(spline->fit() == FALSE) // fit splines to the data
    {
      cerr << "Fatal Error: Spline fit failed.\n";
      exit(3);
    }
}

void 
smoothVolumeLookup(TBSplineVolume *spline, double* src, int* dims)
{
  int x,y,z;
  double value;
  long area, vol, z_area, y_dims;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
		value = (*spline)(x,y,z); 
		src[z_area + y_dims + x] = value;
	  }
    }
  }
}
