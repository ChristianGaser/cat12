#!/usr/bin/env make -f
#
# $Id$

include Makefile.var

CC=gcc

all: AmapMex.$(SUF) ornlmMex.$(SUF) AmapMexNu.$(SUF) KmeansMex.$(SUF)

AmapMex.$(SUF): AmapMex.c Amap.o KmeansProper.o MrfPrior.o Pve.o
	$(MEX) AmapMex.c Amap.o KmeansProper.o MrfPrior.o Pve.o $(MEXEND)

AmapMexNu.$(SUF): AmapMexNu.c Amap.o Kmeans.o MrfPrior.o Pve.o SplineSmooth.o
	$(MEX) AmapMexNu.c Amap.o Kmeans.o MrfPrior.o Pve.o SplineSmooth.o $(MEXEND) ./$(EXT)/libEBTKS.a $(MEXEND)

ornlmMex.$(SUF): ornlmMex.c ornlm.o
	$(MEX) ornlmMex.c ornlm.o $(MEXEND)

KmeansMex.$(SUF): KmeansMex.c KmeansProper.o
	$(MEX) KmeansMex.c KmeansProper.o $(MEXEND)

PveAmapMex.$(SUF): PveAmap.o Amap.o MrfPrior.o Pve.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o SplineSmooth.o $(OBS2)
	$(MEX) PveAmapMex.c PveAmap.o Amap.o MrfPrior.o Pve.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o SplineSmooth.o -I. ./$(EXT)/libEBTKS.a $(MEXEND)

%.o : %.c
	$(CC) -fPIC -c -O2 $< $(MEXEND)

SplineSmooth.o : SplineSmooth.cc
	$(CXX) -fPIC -O2 -I. -c $< $(MEXEND)

%.$(SUF) : %.c %.cc
	$(MEX)  $< $(MEXEND)

clean: 
	$(DEL) *.o AmapMex.$(SUF).a

-include Makefile.vbm

