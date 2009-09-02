#!/usr/bin/env make -f
#
# $Id$

include Makefile.var

OBS = PveAmap.o Amap.o MrfPrior.o Pve5.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o splineSmooth.o
OBS2 = Amap.o MrfPrior.o Pve6.o

AmapMex.$(SUF): AmapMex.c $(OBS2)
	$(MEX) AmapMex.c $(OBS2) $(MEXEND)

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
	$(DEL) $(OBS) $(OBS2) PveAmap.$(SUF).a

-include Makefile.vbm

