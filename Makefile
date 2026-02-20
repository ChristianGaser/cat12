# Personal Makefile variables
#
# $Id$

.PHONY: clean install zip docs update cp_binaries archive check_pipeline release standalone

OLDVERSION="CAT26.0.rc1"
NEWVERSION="CAT26.0.rc1"
REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`
VERSION=`echo ${NEWVERSION} | sed -e 's/CAT//g'`
MAINVERSION=`echo ${NEWVERSION} | sed -e 's/CAT//g' | cut -f1 -d'.'`

ZIPFOLDER=/Users/gaser/matlab/cat12

TARGET=/Users/gaser/spm/spm12/toolbox/CAT
TARGET2=/Volumes/UltraMax/spm12/toolbox/CAT
TARGET4=/Users/gaser/spm/spm-octave/toolbox/CAT

PRECOMPILED=/Users/gaser/matlab/Matlab_standalone

STARGET3_HOST=paris.biomag.uni-jena.de
STARGET3_FOLDER=/home/gaser/spm12/toolbox/CAT
STARGET3=${STARGET3_HOST}:${STARGET3_FOLDER}

MATLAB_FILES=Contents.* cat_*.m spm_CAT.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m cat_run* compile.m
C_FILES=Amap.[ch] ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c* cat_*.mex* vollib.c genus0.[ch] tricases.h spm_diffeo.* tfceMex_pthread.*
MISC_FILES=README.md CHANGES.txt INSTALL.txt doc standalone templates_MNI152NLin2009cAsym templates_surfaces templates_surfaces_32k atlases_surfaces atlases_surfaces_32k cat12.* CAT.* distribute_to_server.sh cat_*.sh  cat_long_main*txt glass_brain.mat

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=CAT_r${REVISION}.zip

# remove .DS_Store files and correct file permissions
clean:
	-@find . -type f -name .DS_Store -exec rm {} \;
	-@chmod -R a+r,g+w,o-w .
	-@find . -type f \( -name "*.sh" -o -name "*.mex*" \) -exec chmod a+x {} \;
	-@find . -type f \( -name "*.c" -o -name "*.c??" -o -name "*.m" -o -name "*.gii" -o -name "*.nii" -o -name "*.txt" -o -name "*.html" -o -name "*.annot" \) -exec chmod a-x {} \;

# prepare txt file for deployed versions
copy_longmode:
	-@cp -R cat_long_main.m cat_long_main.txt

# install
install: copy_longmode
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}/*
	-@mkdir -p ${TARGET}
	-@cp -R ${FILES} ${TARGET}/
	-@gzip -df ${TARGET}/*/*.nii.gz

#install on UltraMax
install2: copy_longmode
	-@echo install on ${TARGET2}
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}/*
	-@mkdir -p ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}/
	-@gzip -d ${TARGET2}/*/*.nii.gz

#install on paris
install3: copy_longmode
	-@echo install on ${STARGET3}
	-@bash -c "ssh ${STARGET3_HOST} 'test ! -d ${STARGET3_FOLDER} || rm -rf ${STARGET3_FOLDER}/*'"
	-@bash -c "ssh ${STARGET3_HOST} 'test -d ${STARGET3_FOLDER} || mkdir ${STARGET3_FOLDER}'"
	-@scp -r ${FILES} ${STARGET3}/
	-@gzip -d ${TARGET3}/*/*.nii.gz

# install for octave
install4: copy_longmode
	-@echo install
	-@test ! -d ${TARGET4} || rm -rf ${TARGET4}/*
	-@mkdir -p ${TARGET4}
	-@cp -R ${FILES} ${TARGET4}/
	-@gzip -d ${TARGET4}/*/*.nii.gz

# print available commands
help:
	-@echo Available commands:
	-@echo clean install zip docs update cp_binaries archive check_pipeline release standalone

#make html documentation
docs:
	-@cat doc/cat.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > doc/cat.html
	-@test ! -d ../cat12-help || cp -R doc/* ../cat12-help/

# update version numbers
update: docs copy_longmode
	-@git fetch
	-@echo '% Computational Anatomy Toolbox' > Contents.m
	-@echo '% Version' ${REVISION}' ('${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@cp Contents.m Contents.txt
	-@echo '% Computational Anatomy Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt
	-@test ! -d ../tfce/ || cp cat_spm_results_ui.m ../tfce/
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_CAT.m
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" cat_batch_bids.sh
	-@cp cat12.m cat${MAINVERSION}.m
	-@chmod a+x CAT.*/*

# zip release
zip: update clean
	-@echo zip
	-@test ! -d CAT || rm -r CAT
	-@mkdir CAT
	-@cp -rp ${FILES} CAT
	-@gzip -d CAT/templates_MNI152NLin2009cAsym/*.nii.gz
	-@bash update_revision.sh
	-@zip ${ZIPFOLDER}/${ZIPFILE} -rm CAT

# copy binaries after cross-compiling
cp_binaries: 
	-@echo copy binaries
	-@for file in ~/Dropbox/GitHub/CAT-Surface/build-*/Progs/*.o; do \
			[ -f "$$file" ] && rm "$$file"; \
	done
	-@for i in CAT.glnx86/CAT*; do cp ~/Dropbox/GitHub/CAT-Surface/build-x86_64-pc-linux/Progs/`basename $${i}` CAT.glnx86/ ; done
	-@for i in CAT.w32/CAT*; do cp ~/Dropbox/GitHub/CAT-Surface/build-x86_64-w64-mingw32/Progs/`basename $${i}` CAT.w32/ ; done
	-@for i in CAT.maci64/CAT*; do cp ~/Dropbox/GitHub/CAT-Surface/build-native/Progs/`basename $${i}` CAT.maci64/ ; done
	-@for i in CAT.maca64/CAT*; do cp ~/Dropbox/GitHub/CAT-Surface/build-native-arm64/Progs/`basename $${i}` CAT.maca64/ ; done

# print check list for releasing
release: checklist
checklist:
	@printf '%s\n' \
	'' \
	'Checklist for testing CAT12 in order to release' \
	'-----------------------------------------------' \
	'1. Check Pipeline' \
	'   make check_pipeline' \
	'   check_pipeline.sh -p pid' \
	'   check_all_matrix.sh' \
	'   check_pipeline_ROIs.m -> check render views check_r*matrix.png and histograms' \
	'   check_pipeline_homogeneity.m -> check sample homogeneity' \
	'' \
	'2. Check Batches and Dependencies' \
	'   cd check_pipeline' \
	'   batch_volume_pipeline.m' \
	'   batch_surface_pipeline.m' \
	'' \
	'3. Check Expert Mode' \
	'   cat12('\''expert'\'')' \
	'   CAT12 GUI Segment' \
	'' \
	'4. Check Simple Preprocessing' \
	'   SPM->Tools->CAT12->CAT12 Simple Preprocessing' \
	'' \
	'5. Create and Check Standalone Versions' \
	'   https://github.com/ChristianGaser/cat12/actions/workflows/install_test_standalone.yml' \
	'' \
	'6. Check Skull-Stripping if Pipeline Changed' \
	'   cat12_all.m in /Volumes/UltraMax/validate_skullstripping_withT12' \
	'   calc_kappa_c0_SPM12_T12.m' \
	'' \
	'7. Make fork from new version!' \
	'' \
	'8. Check thickness phantom'

# print help for standalone
standalone:
	@printf '%s\n' \
	'' \
	'Use CAT action to compile and test standalone versions' \
	'https://github.com/ChristianGaser/cat12/actions/workflows/install_test_standalone.yml'

# rescue from archives
archive:
	-@echo available archives to install
	-@ls cat12_r*zip
	-@test ! -d cat12 || rm -rf cat12
	-@test ! -d ${TARGET} || rm -rf ${TARGET}
	-@read -p "Type release number (3 or 4 digits), followed by [ENTER]:" ver; unzip cat12_r$${ver}.zip; cp -R cat12 ${TARGET}

# run check pipeline
check_pipeline: update install
	-@echo Check pipeline
	-@cd /Volumes/UltraMax/check_pipeline/
	-@/Users/gaser/Dropbox/GitHub/cat12/check_pipeline.sh -d /Volumes/UltraMax/check_pipeline/ -s ~/spm/spm12 -bg 8 -f /Volumes/UltraMax/check_pipeline/check_pipeline_files.txt -l /Volumes/UltraMax/check_pipeline/check_pipeline_files_long.txt
	-@echo Please finally call post-processing with the resp. pid: check_pipeline.sh -p pid