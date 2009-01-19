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
$RCSfile: OpTimer.h,v $
$Revision$
$Author$
$Date$
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _OP_TIMER_H
#define _OP_TIMER_H

#include <iostream>		/* (bert) changed from iostream.h */

typedef double (*TimeFunc)();

class OpTimer {
  int       _timeType;
  char      _verbose;
  char     *_operation;
  double    _start;
  unsigned int _NN;
  unsigned int _interval;
  unsigned int _i;
  std::ostream  *_os;
  TimeFunc  _time;
  
  static const char *_TIME_STRINGS[];

public:
  static const int CPU, SYS, USR;

  // Create timer for the spoecified time type (USR, CPU, or SYS).
  // <operation> is used only for reporting purposes
  // <N> is the total number of timer 'toc's (default 0 (== 1)).
  // <reportInterval> is the number of 'toc's between reporting the elapsed time.
  OpTimer(int timeType = USR, const char *operation = 0, unsigned N = 0, 
	  unsigned reportInterval = 1);
  OpTimer(const char *operation = 0, unsigned N = 0, unsigned reportInterval = 1);
  OpTimer(const OpTimer&); // Undefined; cannot copy OpTimers
  ~OpTimer();

  // Re-initialize OpTimer functions
  double operator () (int timeType = USR, const char *operation = 0, unsigned N = 0, 
		      unsigned reportInterval = 1);
  double operator () (const char *operation = 0, unsigned N = 0, 
		      unsigned reportInterval = 1);

  // Attribute setting functions
  void verbose(char on)          { _verbose = on; }
  void outputStream(std::ostream& os) { _os = &os; }
  void timeType(int timeType);
  void timeFunction(TimeFunc F)  { _time = F; }

  // Reset timer with new <N> and <reportInterval>
  double tic(unsigned N = 0, unsigned reportInterval = 1);
  // <i> timer 'toc's have expired; report time if needed.
  double toc(unsigned i = 1);

private:
  static double _CPUtime();
  static double _SYStime();
  static double _USRtime();
  void          _newOperation(const char *operation);
  std::ostream&      _printTime(double sec) const;
};

#endif
