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
$RCSfile: popen.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef POPEN_H
#define POPEN_H

#include <unistd.h>
#include <iostream.h>

#ifndef zapeof
#define   zapeof(c) ((c)&0377)
#endif

//typedef int pid_t;		

extern int psystem(const char *);

class popenbuf : public streambuf{
public:
	popenbuf(const char *, ios::open_mode om);
	~popenbuf();

	int exit();
	void closewriteside();
private:
	enum { bufsize=1024 };

	int underflow();
	int overflow(int = EOF);

	void popen1(const char*);
	void popen2(const char*);

	char *gbuf;
	char *pbuf;
	int running;
	pid_t pid;
	ios::open_mode opm;
	int pfdin;
	int pfdout;
};

class pstreambase : public virtual ios
{
public:
	pstreambase(const char*, ios::open_mode);
	int exit();
	void closewriteside();
private:

	popenbuf buf;
};

class opopen : public pstreambase, public ostream
{
public:
	opopen(const char *);
	int exit();
};

class ipopen : public pstreambase, public istream
{
public:
	ipopen(const char *);
	int exit();
};

class iopopen : public pstreambase, public iostream
{
public:
	iopopen(const char*);
	int exit();
	void closewriteside();
};

#endif
