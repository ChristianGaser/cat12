# Personal Makefile variables
#
# $Id$

OLDVERSION="CAT12.6"
NEWVERSION="CAT12.7"
REVISION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm12/toolbox/cat12
TARGET2=/Volumes/UltraMax/spm12/toolbox/cat12

STARGET_HOST=141.35.69.218
STARGET_HTDOCS=${STARGET_HOST}:/volume1/web/
STARGET_FOLDER=/volume1/web/cat12
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=Contents.* cat_*.m spm_cat12.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m kmeans3D.m cat_run*
C_FILES=Amap.[ch] ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c* cat_*.mex* vollib.c genus0.[ch] tricases.h spm_diffeo_old.mex*
MISC_FILES=CAT12-Manual.pdf CHANGES.txt INSTALL.txt standalone templates_1.50mm html templates_surfaces templates_surfaces_32k atlases_surfaces atlases_surfaces_32k cat12.* CAT.* distribute_to_server.sh cat_*.sh

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=cat12_r${REVISION}.zip

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
	-@echo install zip scp scp_manual update cp_binaries archive check_pipeline

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
	-@echo scp to http://${STARGET_HOST}/cat12/${ZIPFILE}
	-@scp -P 2222 CHANGES.txt CAT12-Manual.pdf ${ZIPFILE} ${STARGET}
	-@test ! -d cat12-html || rm -r cat12-html
	-@cp -R html cat12-html
	-@perl -p -i -e "s/\','-browser'\);//g" cat12-html/*.html
	-@perl -p -i -e "s/\','-browser'\)//g" cat12-html/*.html
	-@perl -p -i -e "s/matlab:web\(\'//g" cat12-html/*.html
	-@cp cat12-html/cat.html cat12-html/index.html
	-@scp -r -P 2222 cat12-html ${STARGET_HTDOCS}/
	-@bash -c "ssh -p 2222 ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/cat12_latest.zip"
	
scp_manual:
	-@echo scp CAT12-Manual.pdf to http://${STARGET_HOST}/cat12
	-@scp -P 2222 CAT12-Manual.pdf ${STARGET}

cp_binaries: 
	-@echo copy binaries
	-@test ! -f ~/work/c/CAT/build-*/Progs/*.o || rm ~/work/c/CAT/build-*/Progs/*.o
	-@for i in CAT.glnx86/CAT*; do cp ~/work/c/CAT/build-x86_64-pc-linux/Progs/`basename $${i}` CAT.glnx86/ ; done
	-@for i in CAT.w32/CAT*; do cp ~/work/c/CAT/build-i586-mingw32/Progs/`basename $${i}` CAT.w32/ ; done
	-@for i in CAT.maci64/CAT*; do cp ~/work/c/CAT/build-native/Progs/`basename $${i}` CAT.maci64/ ; done

archive:
	-@echo available archives to install
	-@ls cat12_r*zip
	-@test ! -d cat12 || rm -rf cat12
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@read -p "Type release number (3 or 4 digits), followed by [ENTER]:" ver; unzip cat12_r$${ver}.zip; cp -R cat12 ${TARGET}
	
check_pipeline: update install
	-@echo Check pipeline
	-@ls ./check_pipeline.sh -s ~/spm/spm12 -bg 8 -f /Volumes/UltraMax/check_pipeline/check_pipeline_files.txt
	-@echo Please finally call post-processing with the resp. pid: check_pipeline.sh -p pid
	
