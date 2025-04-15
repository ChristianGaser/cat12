# Personal Makefile variables
#
# $Id$

OLDVERSION="CAT12.8.2"
NEWVERSION="CAT12.9"
REVISION=`git rev-list --count HEAD`
DATE=`git log --date short |grep "Date:"|head -1|cut -f2 -d':'|sed -e s'/ //g'`
VERSION=`echo ${NEWVERSION} | sed -e 's/CAT//g'`

ZIPFOLDER=/Users/gaser/matlab/cat12

TARGET=/Users/gaser/spm/spm12/toolbox/cat12
TARGET2=/Volumes/UltraMax/spm12/toolbox/cat12
TARGET4=/Users/gaser/spm/spm-octave/toolbox/cat12

PRECOMPILED=/Users/gaser/matlab/Matlab_standalone

STARGET_HOST=141.35.69.218
STARGET_HTDOCS=${STARGET_HOST}:/volume1/web/
STARGET_FOLDER=/volume1/web/cat12
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

STARGET3_HOST=paris.biomag.uni-jena.de
STARGET3_FOLDER=/home/gaser/spm12/toolbox/cat12
STARGET3=${STARGET3_HOST}:${STARGET3_FOLDER}

MATLAB_FILES=Contents.* cat_*.m spm_cat12.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m cat_run* compile.m
C_FILES=Amap.[ch] ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c* cat_*.mex* vollib.c genus0.[ch] tricases.h spm_diffeo.* tfceMex_pthread.*
MISC_FILES=README.md CHANGES.txt INSTALL.txt doc standalone templates_MNI152NLin2009cAsym templates_surfaces templates_surfaces_32k atlases_surfaces atlases_surfaces_32k cat12.* CAT.* distribute_to_server.sh cat_*.sh  cat_long_main*txt glass_brain.mat

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=cat12_r${REVISION}.zip

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
	-@gzip -d ${TARGET}/*/*.nii.gz

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
	-@echo clean install zip scp scp_precompile scp_standalone docs update cp_binaries archive check_pipeline checklist precompile standalone

#make html documentation
docs:
	-@cat doc/cat.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > doc/cat.html
	-@cp -R doc/* ../cat12-help/

# update version numbers
update: docs copy_longmode
	-@git fetch
	-@git tag -f v${VERSION} -m "Release version ${VERSION}"
	-@echo '% Computational Anatomy Toolbox' > Contents.m
	-@echo '% Version' ${REVISION}' ('${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@cp Contents.m Contents.txt
	-@echo '% Computational Anatomy Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt
	-@cp cat_spm_results_ui.m ../tfce/
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_cat12.m
	-@chmod a+x CAT.*/*

# zip release
zip: update clean
	-@echo zip
	-@test ! -d cat12 || rm -r cat12
	-@mkdir cat12
	-@cp -rp ${FILES} cat12
	-@gzip -d cat12/templates_MNI152NLin2009cAsym/*.nii.gz
	-@bash update_revision.sh
	-@zip ${ZIPFOLDER}/${ZIPFILE} -rm cat12

# scp release
scp: docs zip
	-@echo scp to http://${STARGET_HOST}/cat12/${ZIPFILE}
	-@scp -O -P ${PORT} CHANGES.txt ${ZIPFOLDER}/${ZIPFILE} ${STARGET}
	-@bash -c "ssh -p ${PORT} ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/cat12_latest.zip"
	-@open http://${STARGET_HOST}/cat12/

# scp deployed versions
scp_precompile:
	-@echo scp_precompile
	-@find ${PRECOMPILED} -type f -name .DS_Store -exec rm {} \;
	-@chmod -R a+r,go-w ${PRECOMPILED}
	-@find ${PRECOMPILED} -type f \( -name "*.sh" -o -name "spm25" \) -exec chmod a+x {} \;
	# removed Mac_arm64 that is not working yet
	-@for i in Mac Linux Win MAC_arm64; do \
	   mkdir -p ${NEWVERSION}_R2023b_MCR_$${i} ;\
	   ln -s ${PRECOMPILED}/MCR_$${i}/*spm25* ${PRECOMPILED}/MCR_$${i}/readme.txt ${PRECOMPILED}/MCR_$${i}/MCR_v232.webloc ${NEWVERSION}_R2023b_MCR_$${i}/ ;\
	   cp -r standalone ${NEWVERSION}_R2023b_MCR_$${i}/ ;\
	   cp -r standalone ${PRECOMPILED}/MCR_$${i}/ ;\
	   zip ${ZIPFOLDER}/${NEWVERSION}_R2023b_MCR_$${i}.zip -r ${NEWVERSION}_R2023b_MCR_$${i} ; \
	   scp -O -P ${PORT} ${ZIPFOLDER}/${NEWVERSION}_R2023b_MCR_$${i}.zip ${STARGET}; \
	   bash -c "ssh -p ${PORT} ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${NEWVERSION}_R2023b_MCR_$${i}.zip ${STARGET_FOLDER}/cat12_latest_R2023b_MCR_$${i}.zip"; \
	done
	-@rm -r ${NEWVERSION}_R2023b_MCR*	
	-@echo Please keep in mind to change ../enigma-cat12/index.html
	-@see ../enigma-cat12/index.html
	

# scp deployed versions
scp_standalone: scp_precompile

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
checklist:
	-@echo    
	-@echo Checklist for testing CAT12 in order to release
	-@echo -----------------------------------------------
	-@echo 1. Check Test Data
	-@echo    "cd  ~/matlab/vbm8/test/maci64 && calc_all.sh"
	-@echo    mv surf surf_rXXXX
	-@echo    render_surf.sh surf_rXXXX -range 0 6
	-@echo    
	-@echo 2. Check Pipeline
	-@echo    make check_pipeline
	-@echo    check_pipeline.sh -p pid
	-@echo    check_all_matrix.sh
	-@echo    check_pipeline_ROIs.m	-> check render views check_r*matrix.png and histograms
	-@echo    check_pipeline_homogeneity.m	-> check sample homogeneity
	-@echo    
	-@echo 3. Check Batches and Dependencies
	-@echo    cd check_pipeline
	-@echo    batch_volume_pipeline.m
	-@echo    batch_surface_pipeline.m
	-@echo    
	-@echo 4. Check Expert Mode
	-@echo    "cat12('expert')"
	-@echo    CAT12 GUI Segment
	-@echo    
	-@echo 5. Check Simple Preprocessing
	-@echo    SPM->Tools->CAT12->CAT12 Simple Preprocessing
	-@echo    
	-@echo 6. Check Precompiled Versions
	-@echo    ${PRECOMPILED}/MCR_Mac/standalone/cat_standalone.sh 
	-@echo    
	-@echo 7. Check Skull-Stripping
	-@echo    cat12_all.m in /Volumes/UltraMax/validate_skullstripping_withT12
	-@echo    calc_kappa_c0_SPM12_T12.m
	-@echo    
	-@echo 8. Check Windows 11 + Ubuntu 18.04
	-@echo    UTM.app
	-@echo    CAT12 GUI Segment
	-@echo    
	-@echo 9. Check old SPM12 version on UltraMax
	-@echo    
	-@echo 10.Make fork from new version!
	-@echo    
	-@echo 11.Check thickness phantom 

# print help for precompiling
precompile:
	-@echo    
	-@echo Checklist for precompiling CAT12
	-@echo -----------------------------------------------
	-@echo    spm12_R2023b
	-@echo    spm fmri
	-@echo    cd spm12/config
	-@echo    spm_make_standalone
	-@echo    "Ubuntu 18.04 (run under UTM) : mv /Users/gaser/spm/standalone/spm25.ctf ${PRECOMPILED}/MCR_Linux/"
	-@echo    "Windows 11 ARM64(run under UTM): mv /Users/gaser/spm/standalone/spm25.[ce][tx][fe] ${PRECOMPILED}/MCR_Win/"
	-@echo    "MacOS INTEL (use /Applications/MATLAB_R2023b_intel.app/bin/matlab): rm -rf ${PRECOMPILED}/MCR_Mac/spm25.app; mv /Users/gaser/spm/standalone/spm25.app ${PRECOMPILED}/MCR_Mac/"
	-@echo    "MacOS ARM64: rm -rf ${PRECOMPILED}/MCR_Mac_arm64/spm25.app; mv /Users/gaser/spm/standalone/spm25.app ${PRECOMPILED}/MCR_Mac_arm64/"
	-@echo    

# print help for standalone
standalone: precompile

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
