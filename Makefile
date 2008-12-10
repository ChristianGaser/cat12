VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`

TARGET=/Users/gaser/spm/spm8b/toolbox/vbm8

STARGET=141.35.200.101:/Applications/xampp/htdocs/

FILES=cg_config_vbm8.m Amap.c AmapMex.* MrfPrior.c cg_vbm8_run.m cg_vbm8_write.m cg_check_sample_sd.m cg_showslice_all.m cg_spmT2x.m cg_vbm8_tools.m cg_vbm8_debug.m cg_morph_vol.m cg_cleanup_gwc.m spm_vbm8.m vbm8.man brainmask_LPBA40.nii

ZIPFILE=vbm8_r$(VERSION).zip

install: 
	-@echo install
	-@test ! -d ${TARGET} || rm -r ${TARGET}
	-@mkdir ${TARGET}
	-@cp -r ${FILES} ${TARGET}

help:
	-@echo Available commands:
	-@echo install zip scp upgrade

zip:
	-@echo zip
	-@test ! -d vbm8 || rm -r vbm8
	-@cp -rp ${TARGET} .
	-@zip ${ZIPFILE} -rm vbm8

scp: zip
	-@echo scp
	-@cp ${ZIPFILE} vbm8_latest.zip
	-@scp -pr vbm8_latest.zip ${ZIPFILE} ${STARGET}
