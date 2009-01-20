#!/usr/bin/env make -f
#
# $Id$

include Makefile.var

OBS = PveAmap.o Amap.o MrfPrior.o Pve5.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o splineSmooth.o

PveAmapMex.$(SUF): PveAmapMex.c PveAmap.$(SUF).a
#	$(MEX) PveAmapMex.c PveAmap.$(SUF).a -lEBTKS -L./$(EXT) -I./ $(MEXEND)
	$(MEX) PveAmapMex.c PveAmap.$(SUF).a ./$(EXT)/libEBTKS.a $(MEXEND)

#archive: PveAmap.$(SUF).a

PveAmap.$(SUF).a: $(OBS)
	$(DEL) $@
	$(AR) $@ $(OBS)

%.o : %.c
	$(CC) -c $< $(MEXEND)

%.o : %.cc
	$(CC) -I./ -c $< $(MEXEND)

%.$(SUF) : %.c %.cc
	$(MEX)  $< $(MEXEND)

clean: 
	$(DEL) $(OBS) PveAmapMex.$(SUF) PveAmap.$(SUF).a

-include Makefile.vbm

