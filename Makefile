# Personal Makefile variables
#
# $Id$

OLDVERSION="CAT12.1"
NEWVERSION="CAT12.2"
REVISION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm12/toolbox/cat12
TARGET2=/Volumes/UltraMax/spm12/toolbox/cat12

STARGET_HOST=dbm.neuro.uni-jena.de
STARGET_HTDOCS=${STARGET_HOST}:/Applications/xampp/htdocs/
STARGET_FOLDER=/Applications/xampp/htdocs/cat12
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=Contents.* cat_*.m spm_cat12.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m kmeans3D.m
C_FILES=Amap.[ch] ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c* cat_*.mex* vollib.c genus0.[ch] tricases.h
MISC_FILES=CAT12-Manual.pdf CHANGES.txt INSTALL.txt templates_1.50mm html templates_surfaces templates_surfaces_32k atlases_surfaces atlases_surfaces_32k cat12.* CAT.* distribute_to_server.sh cat_*.sh

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=cat12_r$(REVISION).zip

install: 
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@mkdir ${TARGET}
	-@cp -R ${FILES} ${TARGET}

install2:
	-@echo install2
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}
	-@mkdir ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}

help:
	-@echo Available commands:
	-@echo install zip scp scp_manual update cp_binaries

update:
	-@svn update
	-@echo '% Computational Anatomy Toolbox' > Contents.m
	-@echo '% Version' ${REVISION}' ('${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@cp Contents.m Contents.txt
	-@echo '% Computational Anatomy Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_cat12.m
	-@cat html/cat.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > html/cat.html

zip: update
	-@echo zip
	-@test ! -d cat12 || rm -r cat12
	-@mkdir cat12
	-@cp -rp ${FILES} cat12
	-@zip ${ZIPFILE} -rm cat12

scp: zip
	-@echo scp to http://dbm.neuro.uni-jena.de/cat12/${ZIPFILE}
	-@scp -P 2222 CHANGES.txt CAT12-Manual.pdf ${ZIPFILE} ${STARGET}
	-@test ! -d cat12-html || rm -r cat12-html
	-@cp -R html cat12-html
	-@perl -p -i -e "s/\','-browser'\);//g" cat12-html/*.html
	-@perl -p -i -e "s/\','-browser'\)//g" cat12-html/*.html
	-@perl -p -i -e "s/matlab:web\(\'//g" cat12-html/*.html
	-@cp cat12-html/cat.html cat12-html/index.html
	-@scp -r -P 2222 cat12-html ${STARGET_HTDOCS}/
	-@bash -c "ssh ${STARGET_HOST} ln -Fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/cat12_latest.zip"
	
scp_manual:
	-@echo scp CAT12-Manual.pdf to http://dbm.neuro.uni-jena.de/cat12
	-@scp -P 2222 CAT12-Manual.pdf ${STARGET}

cp_binaries: 
	-@echo copy binaries
	-@test ! -f ~/work/c/CAT/build-*/Progs/*.o || rm ~/work/c/CAT/build-*/Progs/*.o
	-@for i in CAT.glnx86/CAT*; do cp ~/work/c/CAT/build-x86_64-pc-linux/Progs/`basename $${i}` CAT.glnx86/ ; done
	-@for i in CAT.w32/CAT*; do cp ~/work/c/CAT/build-i586-mingw32/Progs/`basename $${i}` CAT.w32/ ; done
	-@for i in CAT.maci64/CAT*; do cp ~/work/c/CAT/build-native/Progs/`basename $${i}` CAT.maci64/ ; done
