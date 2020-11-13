# Personal Makefile variables
#
# $Id$

OLDVERSION="CAT12.7-RC2"
NEWVERSION="CAT12.7"
REVISION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm12/toolbox/cat12
TARGET2=/Volumes/UltraMax/spm12/toolbox/cat12
TARGET3=paris.biomag.uni-jena.de:/home/gaser/spm12/toolbox/cat12

PRECOMPILED=/Users/gaser/install/Matlab/Matlab_R2017b
CAT12=/Users/gaser/matlab/cat12

STARGET_HOST=141.35.69.218
STARGET_HTDOCS=${STARGET_HOST}:/volume1/web/
STARGET_FOLDER=/volume1/web/cat12
STARGET=${STARGET_HOST}:${STARGET_FOLDER}

MATLAB_FILES=Contents.* cat_*.m spm_cat12.m tbx_cfg_cat.m sliderPanel.m slice_overlay.m kmeans3D.m cat_run* compile.m
C_FILES=Amap.[ch] ornlm_float.c sanlm_float.c MrfPrior.c Pve.c Kmeans.c cat_*.c* cat_*.mex* vollib.c genus0.[ch] tricases.h spm_diffeo_old.mex*
MISC_FILES=CAT12-Manual.pdf CHANGES.txt INSTALL.txt standalone templates_volumes html templates_surfaces templates_surfaces_32k atlases_surfaces atlases_surfaces_32k cat12.* CAT.* distribute_to_server.sh cat_*.sh  cat_long_main*txt

FILES=${MATLAB_FILES} ${C_FILES} ${MISC_FILES}

ZIPFILE=cat12_r${REVISION}.zip

# remove .DS_Store files and correct file permissions
clean:
	-@find . -type f -name .DS_Store -exec rm {} \;
	-@chmod -R a+r .
	-@find . -type f \( -name "*.c" -o -name "*.c??" -o -name "*.m" -o -name "*.gii" -o -name "*.nii" -o -name "*.txt" -o -name "*.html" -o -name "*.annot" \) -exec chmod a-x {} \;
# 	-@svn propset svn:executable OFF *.c *.c?? *.m */*.gii */*.nii *.txt *.html */*.c */*.c?? */*.m */*.annot */*.txt */*.html

# prepare txt file for deployed versions
copy_longmode:
	-@cp -R cat_long_main.m cat_long_main.txt

# install
install: copy_longmode
	-@echo install
	-@test ! -d ${TARGET} || rm -rf ${TARGET}/*
	-@mkdir -p ${TARGET}
	-@cp -R ${FILES} ${TARGET}/

#install on UltraMax
install2: copy_longmode
	-@echo install2
	-@test ! -d ${TARGET2} || rm -rf ${TARGET2}/*
	-@mkdir -p ${TARGET2}
	-@cp -R ${FILES} ${TARGET2}/

#install on paris
install3: copy_longmode
	-@echo install3
	-@scp -r ${FILES} ${TARGET3}/

# print available commands
help:
	-@echo Available commands:
	-@echo clean install zip scp scp_manual scp_precompile doc update cp_binaries archive check_pipeline checklist precompile

#make html documentation
doc:
	-@cat html/cat.txt | sed -e 's/VERSION/'${NEWVERSION}'/g' -e 's/RELNUMBER/r'${REVISION}'/g' -e 's/DATE/'${DATE}'/g' > html/cat.html
	-@test ! -d cat12-html || rm -r cat12-html
	-@cp -R html cat12-html
	-@perl -p -i -e "s/\','-browser'\);//g" cat12-html/*.html
	-@perl -p -i -e "s/\','-browser'\)//g" cat12-html/*.html
	-@perl -p -i -e "s/matlab:web\(\'//g" cat12-html/*.html
	-@cp cat12-html/cat.html cat12-html/index.html

# update version numners
update: doc copy_longmode
	-@svn update
	-@echo '% Computational Anatomy Toolbox' > Contents.m
	-@echo '% Version' ${REVISION}' ('${NEWVERSION}')' ${DATE} >> Contents.m
	-@cat Contents_info.txt >> Contents.m
	-@cp Contents.m Contents.txt
	-@echo '% Computational Anatomy Toolbox' > INSTALL.txt
	-@echo '% Version ' ${REVISION} ${NEWVERSION} ${DATE} >> INSTALL.txt
	-@cat INSTALL_info.txt >> INSTALL.txt
	-@perl -p -i -e "s/${OLDVERSION}/${NEWVERSION}/g" spm_cat12.m

# zip release
zip: update clean
	-@echo zip
	-@test ! -d cat12 || rm -r cat12
	-@mkdir cat12
	-@cp -rp ${FILES} cat12
	-@zip ${ZIPFILE} -rm cat12

# scp release
scp: html zip
	-@echo scp to http://${STARGET_HOST}/cat12/${ZIPFILE}
	-@scp -P 2222 CHANGES.txt CAT12-Manual.pdf ${ZIPFILE} ${STARGET}
	-@scp -r -P 2222 cat12-html ${STARGET_HTDOCS}/
	-@bash -c "ssh -p 2222 ${STARGET_HOST} ln -fs ${STARGET_FOLDER}/${ZIPFILE} ${STARGET_FOLDER}/cat12_latest.zip"

# scp manual
scp_manual:
	-@echo scp CAT12-Manual.pdf to http://${STARGET}
	-@scp -P 2222 CAT12-Manual.pdf ${STARGET}

# scp deployed versions
scp_precompile:
	-@echo scp_precompile
	-@for i in Linux Mac Win; do \
	    mkdir -p MCR_$${i} ;\
	    ln -s ${PRECOMPILED}/MCR_$${i}/*spm12* ${PRECOMPILED}/MCR_$${i}/readme.txt ${PRECOMPILED}/MCR_$${i}/MCR_v93.webloc MCR_$${i}/ ;\
	    cp -r standalone MCR_$${i}/ ;\
	    zip cat12_latest_R2017b_MCR_$${i}.zip -r MCR_$${i} ;\
	  done
	-@scp -P 2222 cat12_latest_R2017b_MCR* ${STARGET}
	-@rm -r cat12_latest_R2017b_MCR* MCR_*

# copy binaries after cross-compiling
cp_binaries: 
	-@echo copy binaries
	-@test ! -f ~/work/c/CAT/build-*/Progs/*.o || rm ~/work/c/CAT/build-*/Progs/*.o
	-@for i in CAT.glnx86/CAT*; do cp ~/work/c/CAT/build-x86_64-pc-linux/Progs/`basename $${i}` CAT.glnx86/ ; done
	-@for i in CAT.w32/CAT*; do cp ~/work/c/CAT/build-i586-mingw32/Progs/`basename $${i}` CAT.w32/ ; done
	-@for i in CAT.maci64/CAT*; do cp ~/work/c/CAT/build-native/Progs/`basename $${i}` CAT.maci64/ ; done

# print check list for releasing
checklist:
	-@echo    
	-@echo Checklist for testing CAT12 in order to release
	-@echo -----------------------------------------------
	-@echo 1. Check Test data
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
	-@echo 5. Check Previous Release
	-@echo    "cat12('expert')"
	-@echo    CAT12 GUI Segment CAT12.X
	-@echo    
	-@echo 6. Check Skull-Stripping
	-@echo    cat12_all.m in /Volumes/UltraMax/validate_skullstripping_withT12
	-@echo    calc_kappa_c0_SPM12_T12.m
	-@echo    
	-@echo 7. Check Windows + Linux
	-@echo    VirtualBox.app
	-@echo    CAT12 GUI Segment
	-@echo    
	-@echo 8. Check old SPM12 version on UltraMax
	-@echo    
	-@echo 9. Check thickness phantom 

# print help for precompiling
precompile:
	-@echo    
	-@echo Checklist for precompiling CAT12
	-@echo -----------------------------------------------
	-@echo    spm fmri
	-@echo    cd spm12/config
	-@echo    spm_make_standalone
	-@echo    "Ubuntu 14.10: mv  /Users/gaser/spm/standalone/spm12.ctf /Users/gaser/install/Matlab/Matlab_R2017b/MCR_Linux/"
	-@echo    "Windows 10: mv  /Users/gaser/spm/standalone/spm12.exe /Users/gaser/install/Matlab/Matlab_R2017b/MCR_Win/"
	-@echo    "Mac OS: rm -rf /Users/gaser/install/Matlab/Matlab_R2017b/MCR_Mac/spm12.app; mv /Users/gaser/spm/standalone/spm12.app /Users/gaser/install/Matlab/Matlab_R2017b/MCR_Mac/"
	-@echo    

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
	-@./check_pipeline.sh -d /Volumes/UltraMax/check_pipeline/ -s ~/spm/spm12 -bg 8 -f /Volumes/UltraMax/check_pipeline/check_pipeline_files.txt
	-@echo Please finally call post-processing with the resp. pid: check_pipeline.sh -p pid
