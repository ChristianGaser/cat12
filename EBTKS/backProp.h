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
$RCSfile: backProp.h,v $
$Revision: 1.4 $
$Author: stever $
$Date: 2003/11/17 04:07:52 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef BACKPROP_H
#define BACKPROP_H

#include <iostream>		/* (bert) changed from iostream.h */
#include <stdio.h>
#include "Array.h"
#include "TrainingSet.h"

typedef void (ErrorMonitor)(unsigned cycle, double error);

typedef struct {
  double out;
  double delta;
  double bias;
  double dBias;
} BPNODE;

typedef struct {
  double value;
  double delta;
} WEIGHT;

class BP_ANN {
  BPNODE  **_node;
  WEIGHT  **_weight;
  unsigned  _nLayers;
  unsigned *_nNodesInLayer;
  unsigned *_nWeightsForLayer;
  unsigned  _nInputNodes, _nOutputNodes; // Provided for speed only
  double    _learningRate, _momentum, _temperature;
  DblArray  _outputError;
  DblArray  _lut;
  double    _lutStep;
  unsigned  _nCycles, _cycleCtr;
  unsigned  _nSamples, _sampleCtr;
  double    _maxError, _maxDerror;
  Boolean   _stopRequest;
  Boolean   _verbose;
  unsigned  _shuffleInterval;
  Boolean   _sigmoidAtOutputLayer;
  Boolean   _softmaxAtOutputLayer;

// Consts and default values
  static const int   _LUT_LENGTH;
  static const long  _SEED;
  static const int   _N_CYCLES;
  static const float _LEARNING_RATE, _MOMENTUM, _TEMPERATURE;
  static const float _MAX_ERROR, _MAX_D_ERROR;
  static const unsigned _SHUFFLE_INTERVAL;

public:
  static void default_error_monitor(unsigned cycle, double error) {
      std::cout << "Cycle: " << cycle << " error: " << error << std::endl;
  }

// Constructors/destructor
  BP_ANN(const UnsignedArray&, Boolean verbose = FALSE);
  BP_ANN(std::istream&, Boolean verbose = FALSE);
  ~BP_ANN();

// Get functions
  unsigned nInputNodes() const   { return _nInputNodes; }
  unsigned nOutputNodes() const  { return _nOutputNodes; }
  UnsignedArray topology() const { return UnsignedArray(_nNodesInLayer, _nLayers); }
  unsigned nCycles() const       { return _nCycles; }

// Set functions
  Boolean nInputNodes(unsigned);  // Set #input nodes. Re-initializes network!
  Boolean nOutputNodes(unsigned); // Set #output nodes. Re-initializes network!
  Boolean topology(const UnsignedArray&); // Set topology. Re-initializes network!

  void randomize(long seed = 0);  // Randomize network values; a seed value of 0 will
                                  // initialize with srand48(time())
  void setDefaults();
  void set(double learningRate, double momentum, double temperature,
	   unsigned nCycles, double maxError, double maxDerror,
	   unsigned shuffleInterval, Boolean sigmoidAtOutputLayer = TRUE,
	   Boolean softmaxAtOutputLayer = FALSE);
  void shuffleInterval(unsigned I)         { _shuffleInterval = I; }
  //  void errorMonitor(ErrorMonitor F = NULL) { _errorMonitor = F; }

// Operations
  // Full training, based on a TrainingSet
  int      train(TrainingSet&, ErrorMonitor F = NULL);
  // The following functions should be used for incremental training
  // (w/o TraiingSet)
  int      initTraining(unsigned nSamples);
  int      train(const double *input, const double *target, ErrorMonitor F = NULL);
  // Stop training
  void     stop(Boolean on = TRUE) { _stopRequest = on; }
  // Evaluate output for a single input
  void     evaluate(const double *input, double *output);
  // Evaluate output and classify (select max output)
  unsigned classify(const double *input, double *output = 0);

// I/O
  int      load(std::istream&);
  int      save(std::ostream&, Boolean includeContents = TRUE);
  std::ostream& printWeights(std::ostream&);
  std::ostream& printNodes(std::ostream&);

private:
  void _create(const UnsignedArray&);
  void _destroy();
  void _createLut(double T);
  void _adjustWeights();
  void _forward(const double *);
  void _calculateDeltas(const double *);
};

#endif
