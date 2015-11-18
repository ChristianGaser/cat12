# Personal Makefile variables
#
# $Id$

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm12/toolbox/cat12
TARGET2=/Volumes/UltraMax/spm12/toolbox/cat12

STARGET=dbm.neuro.uni-jena.de:/Applications/xampp/htdocs/cat12

MATLAB_FILES=Contents.m cat_*.m spm_cat12.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m
C_FILES=Amap.[ch] cat_*mex.m ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c cat_*.mex* vollib.c genus0.[ch] tricases.h
MISC_FILES=CAT12-Manual.pdf CHANGES.txt INSTALL.txt templates_1.50mm html templates_surfaces cat12.* CAT.* distribute_to_server.sh cat_*.sh

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=cat12_r$(VERSION).zip

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
	-@echo install zip scp update

update:
	-@svn update
	-@echo '% Computational Anatomy Toolbox' > Contents.m
	-@echo '% Version ' ${VERSION} ' (CAT12) ' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@echo '% Computational Anatomy Toolbox' > INSTALL.txt
	-@echo '% Version ' ${VERSION} ' (CAT12) ' ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt
	-@cat html/cat.txt | sed -e 's/RELNUMBER/r'${VERSION}'/g' -e 's/DATE/'${DATE}'/g' > html/cat.html

zip: update
	-@echo zip
	-@test ! -d cat12 || rm -r cat12
	-@mkdir cat12
	-@cp -rp ${FILES} cat12
	-@zip ${ZIPFILE} -rm cat12

scp: zip
	-@echo scp to http://dbm.neuro.uni-jena.de/cat12/${ZIPFILE}
	-@scp -P 2222 CHANGES.txt CAT12-Manual.pdf ${ZIPFILE} ${STARGET}
