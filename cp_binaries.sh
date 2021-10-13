#!/bin/sh

rm ~/work/c/CAT/build-*/Progs/*.o

# check for new arm64
if [ `uname -m` == "arm64" ]; then
	cd ../CAT.maca64
else
	cd CAT.glnx86
	for i in CAT*; do cp ~/work/c/CAT/build-x86_64-pc-linux/Progs/${i} .; done
	
	cd ../CAT.w32
	for i in CAT*; do cp ~/work/c/CAT/build-i586-mingw32/Progs/${i} .; done
	chmod a+x *.exe

	cd ../CAT.maci64
fi

for i in CAT*; do cp ~/work/c/CAT/build-native/Progs/${i} .; done
cd ..
