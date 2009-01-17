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
$RCSfile: matlabSupport.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _MATLAB_SUPPORT
#define _MATLAB_SUPPORT

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

extern "C" {    
  #include"mat.h" 
}

Boolean 
saveMatlab(const char *fileName, const char *varName, const char *option,
	   unsigned nrows, unsigned ncols, const double *real, 
	   const double *imag = NULL)
{
  Boolean status = TRUE;

  if ((option[0] != 'u') && (option[0] != 'w')) {
    cerr<<"Incorrect file write option, use u for update or w to";
    cerr<<" create a new file or overwrite the existing one."<<endl;
    status = FALSE;
  }

  if (status) {
    MATFile *outFile = matOpen((char *) fileName, (char *) option);
    if (!outFile) {
      cerr << "Couldn't open file " << fileName << endl;
      status = FALSE;
    }
    else {
      // A lot of extra stuff to handle a matlab bug which leaves temporary files
      // in /var/tmp.
      
      time_t pre = time(0);
      
      if (matPutFull(outFile, (char *) varName, nrows, ncols, (double *) real, 
		     (double *) imag)) {
	cerr << "Error in writing matrix to file " << fileName << endl;
	status = FALSE;
      }
      
      time_t post = time(0);
      
      char *path = new char[50];
      assert(path);
      sprintf(path, "/var/tmp");
      
      DIR *dir = opendir(path);
      assert(dir);
      
      uid_t uid = getuid();
      
      Boolean tempFileFound = FALSE;
      struct dirent *dirEntry;
      struct stat buf;
      while (!tempFileFound && ((dirEntry = readdir(dir)) != NULL)) {
	sprintf(path, "/var/tmp/%s", dirEntry->d_name);
	stat(path, &buf);
	if ((buf.st_uid == uid) &&
	    (buf.st_atime >= pre) && (buf.st_atime <= post) &&
	    (buf.st_mtime >= pre) && (buf.st_mtime <= post) &&
	    (buf.st_ctime >= pre) && (buf.st_ctime <= post) &&
	    strstr(dirEntry->d_name, "aaa")) {
	  unlink(path);
	  tempFileFound = TRUE;
	}
      }
      
      if (matClose(outFile) == EOF) {
	cerr << "Error in closing file " << fileName << endl;
	status = FALSE;
      }
    }
  }

  return status;
}
#endif
