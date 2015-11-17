#!/usr/bin/env make -f
#
# $Id$

include Makefile.var

CC=gcc

all: AmapMex.$(SUF) ornlmMex.$(SUF) KmeansMex.$(SUF)

AmapMex.$(SUF): AmapMex.c Amap.o Kmeans.o MrfPrior.o Pve.o vollib.o
	$(MEX) AmapMex.c Amap.o Kmeans.o MrfPrior.o Pve.o vollib.o $(MEXEND)

ornlmMex.$(SUF): ornlmMex.c ornlm.o
	$(MEX) ornlmMex.c ornlm.o $(MEXEND)

KmeansMex.$(SUF): KmeansMex.c Kmeans.o
	$(MEX) KmeansMex.c Kmeans.o vollib.o $(MEXEND)

%.o : %.c
	$(CC) -fPIC -c -O2 $< $(MEXEND)

%.$(SUF) : %.c %.cc
	$(MEX)  $< $(MEXEND)

clean: 
	$(DEL) *.o *.$(SUF).a

-include Makefile.cat

