#!/usr/bin/env make -f
#
# $Id$

include Makefile.var

OBS = PveAmap.o Amap.o MrfPrior.o Pve5.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o splineSmooth.o

PveAmapMex.$(SUF): PveAmapMex.c $(OBS)
#	$(MEX) PveAmapMex.c PveAmap.$(SUF).a -lEBTKS -L./$(EXT) -I./ $(MEXEND)
	$(MEX) PveAmapMex.c $(OBS) ./$(EXT)/libEBTKS.a $(MEXEND)

%.o : %.c
	$(CC) -c $< $(MEXEND)

splineSmooth.o : splineSmooth.cc
	$(CXX) -I./ -c $< $(MEXEND)

%.$(SUF) : %.c %.cc
	$(MEX)  $< $(MEXEND)

clean: 
	$(DEL) $(OBS) PveAmap.$(SUF).a

-include Makefile.vbm

