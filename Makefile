#!/usr/bin/env make -f

include Makefile.var

OBS = PveAmap.o Amap.o MrfPrior.o Pve5.o Kmeans.o WarpPriors.o Bayes.o optimizer3d.o diffeo3d.o splineSmooth.o

PveAmapMex.$(SUF): PveAmapMex.c PveAmap.$(SUF).a
	$(MEX) PveAmapMex.c PveAmap.$(SUF).a -lEBTKS -L$(EXT) -I./ $(MEXEND)

archive: PveAmap.$(SUF).a

PveAmap.$(SUF).a: $(OBS)
	$(DEL) $@
	$(AR) $@ $(OBS)

%.o : %.c %.cc
	$(MEX) -c $< $(MEXEND)

%.$(SUF) : %.c %.cc
	$(MEX) $< $(MEXEND)

-include Makefile.vbm

