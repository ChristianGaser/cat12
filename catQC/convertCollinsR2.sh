TC=$(find mnc -name "*.mnc")



# convert GT files
if [ 0 == 1 ]
then
  echo convert ground truth data
  GT=$(find mnc_gt -name "*.mnc.gz"); gunzip $GT; GT=$(find mnc_gt -name "*.mnc"); 
  for F in $GT
  do
    PP=$(dirname $F); 
    FF=$(basename $F .mnc);
    mri_convert $PP/$FF.mnc $PP/$FF.nii > /dev/null 2>&1
  done
fi



# resultdirnames
HDIR=$(pwd)/results;
FULL_RDIR="$HDIR/BWPC_full";          # main result with all test cases                - removed
NOISE_RDIR="$HDIR/BWPC_noise";         # smaller resultdir with only the important test 
BIAS_RDIR="$HDIR/BWPC_bias";           # smaller resultdir with only the important test 
#RES_RDIR="$HDIR/BWPC_res";            # smaller resultdir with only the important test - download required
FIELD_RDIR="$HDIR/BWPC_field";         # smaller resultdir with only the important test 
RESR_RDIR="$HDIR/BWPC_resr";           # smaller resultdir with only the important test 
RESI_RDIR="$HDIR/BWPC_resi";           # smaller resultdir with only the important test 
NIR_RDIR="$HDIR/BWPC_NIR";             # smaller resultdir with only the important test 
CON_RDIR="$HDIR/BWPC_con";             # smaller resultdir with only the important test - further downloads required
WC_RDIR="$HDIR/BWPC_WC";               # smaller resultdir with only the important test 
MS_RDIR="$HDIR/BWPM";                  # smaller resultdir with only the important test - further downloads required T2 PD ...
GT_RDIR="$HDIR/BWPC_gt";               # groud thruth data                              
REWRITE=1;

# further downloads
# - other non T1-T2-PD modalities with pn3 and rf040pA
# - MS PD and T2 with pn3 and rf040pA (and pn7 and rf080pA)
# - other resolution of T1, T2, and PD with pn3 and rf040pA 
 
#init_matlab="addpath '/Users/dahnke/Neuroimaging/SPM12Rbeta' '/Users/dahnke/Neuroimaging/SPM12Rbeta/toolbox/vbm12'";

# create resultdirs
for TT in T1 T2 PD
do
  if [ ! -d "$FULL_RDIR/$TT" ];  then mkdir -p "$FULL_RDIR/$TT";   else if [ "$REWRITE" == "1" ]; then find $FULL_RDIR/$TT  -name "*.nii" -delete; fi
  fi
  if [ ! -d "$NOISE_RDIR/$TT" ]; then mkdir -p "$NOISE_RDIR/$TT";  else if [ $REWRITE ]; then find $NOISE_RDIR/$TT -name "*.nii" -delete; fi 
  fi
  if [ ! -d "$BIAS_RDIR/$TT" ];  then mkdir -p "$BIAS_RDIR/$TT";   else if [ $REWRITE ]; then find $BIAS_RDIR/$TT  -name "*.nii" -delete; fi
  fi
  #if [ ! -d "$RES_RDIR/$TT" ];   then mkdir -p "$RES_RDIR/$TT";    else if [ $REWRITE ]; then find $RES_RDIR/$TT   -name "*.nii" -delete; fi 
  #fi
  if [ ! -d "$FIELD_RDIR/$TT" ]; then mkdir -p "$FIELD_RDIR/$TT";  else if [ $REWRITE ]; then find $FIELD_RDIR/$TT  -name "*.nii" -delete; fi
  fi
  if [ ! -d "$RESR_RDIR/$TT" ];  then mkdir -p "$RESR_RDIR/$TT";   else if [ $REWRITE ]; then find $RESR_RDIR/$TT  -name "*.nii" -delete; fi
  fi
  if [ ! -d "$RESI_RDIR/$TT" ];  then mkdir -p "$RESI_RDIR/$TT";   else if [ $REWRITE ]; then find $RESI_RDIR/$TT  -name "*.nii" -delete; fi
  fi
  if [ ! -d "$NIR_RDIR/$TT" ];   then mkdir -p "$NIR_RDIR/$TT";    else if [ $REWRITE ]; then find $NIR_RDIR/$TT   -name "*.nii" -delete; fi
  fi
  if [ ! -d "$WC_RDIR/$TT" ];    then mkdir -p "$WC_RDIR/$TT";     else if [ $REWRITE ]; then find $WC_RDIR/$TT    -name "*.nii" -delete; fi
  fi
  if [ ! -d "$MS_RDIR/$TT" ];    then mkdir -p "$MS_RDIR/$TT";     else if [ $REWRITE ]; then find $MS_RDIR/$TT    -name "*.nii" -delete; fi
  fi
  if [ ! -d "$GT_RDIR/$TT" ];    then mkdir -p "$GT_RDIR/$TT";     else if [ $REWRITE ]; then find $GT_RDIR/$TT    -name "*.nii" -delete; fi
  fi
done
if [ ! -d "$CON_RDIR" ];  then mkdir -p "$CON_RDIR";   else if [ $REWRITE ]; then find $CON_RDIR  -name "*.nii" -delete; fi
fi

# GT-copy


#RES="1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75 3.00"
#RES="1.00 2.00 3.00"
#RES2="1.00 2.00"
i=0;
for T in $TC
do 
 (( i++ ))

  F=$(dirname  $T)/$(basename $T .mnc);
  FF=$(basename $T .mnc);
  
  # find txt file & check if correct number -> diplay success or error
  if [ -a $F.txt ]
  then
    S=$(grep $FF $F.txt); 
    if [ "$S" == "" ]
    then printf "%04d - $F: ERROR wrong content in txt-file\n" $i; continue; 
    else printf "%04d - $F " $i
    fi
  else
    printf "%04d -$F: ERROR no txt-file\n" $i; continue;
  fi  

  # read parameter
  ECHO_TIME1=$(grep echo_times $F.txt | cut -d = -f 2- | cut -d , -f 1 | tr -d " " );
  ECHO_TIME2=$(grep echo_times $F.txt | cut -d = -f 2- | cut -d , -f 2 | tr -d " " );
  FLIP_ANGLE=$(grep flip_angle $F.txt| cut -d = -f 2- | tr -d " " );
  INU_FIELD=$(grep inu_field $F.txt | cut -d = -f 2- | tr -d " " );
  NO_OF_ECHOES=$(grep no_of_echoes $F.txt | cut -d = -f 2- | tr -d " " );
  PERCENT_INU=$(grep percent_inu $F.txt | cut -d = -f 2- | tr -d " " );
  PERCENT_NOISE=$(grep percent_noise $F.txt | cut -d = -f 2- | tr -d " " );
  PHANTOM=$(grep phantom $F.txt | cut -d = -f 2- | tr -d " " );
  RANDOM_SEED=$(grep "random_seed" $F.txt | cut -d = -f 2- | tr -d " " );
  REFERENCE_TISSUE=$(grep "reference_tissue" $F.txt | cut -d = -f 2- | tr -d " " );
  SCAN_TECHNIQUE=$(grep "scan_technique" $F.txt | cut -d = -f 2- | tr -d " " );
  SLICE_THICKNESS=$(grep "slice_thickness" $F.txt | cut -d = -f 2- | tr -d " " );
  SLICE_RES=$(echo "scale=0;$SLICE_THICKNESS*100" | bc | cut -d . -f 1);
  TI=$(grep "ti = " $F.txt | cut -d = -f 2- | tr -d " " );
  TR=$(grep "tr =" $F.txt | cut -d = -f 2- | tr -d " " );
  TR=$(grep "tr =" $F.txt | cut -d = -f 2- | tr -d " " );

  case $SCAN_TECHNIQUE in
    DSE_LATE)  TT=T2;;
    DSE_EARLY) TT=PD;;
    SFLASH)    TT=T1;;
    FLASH)     TT=FL;;
    FISP)      TT=FI;;
    IR)        TT=IR;;
    SE)        TT=SE;;
    CEFAST)    TT=CE;;
    *)         printf "> non-default scan_technique $SCAN_TECHNIQUE \n"; continue;;
  esac
  
  case $PHANTOM in
    normal)    TG=HC;;
    msles1)    TG=MS1;;
    msles2)    TG=MS2;;
    msles3)    TG=MS3;;
    *)         printf "> non-default scan_technique $TG \n"; continue;;
  esac
      
  # create filename
  APERCENT_INU=$(echo $PERCENT_INU | tr -d -);
  case $PERCENT_INU in -*) SPERCENT_INU="n";; *) SPERCENT_INU="p";; esac
  if [ "$APERCENT_INU" == "0" ]; then SPERCENT_INU="0"; INU_FIELD="0"; fi
  if [ \( "$TT" == "T1" -a "$FLIP_ANGLE" == "30" -a "$TR" == "18" -a "$TI" == "" \) -o \( "$TT" == "T2" -a "$FLIP_ANGLE" == "90" -a "$TR" == "3300" -a "$TI" == "" \) -o \( "$TT" == "PD" -a "$FLIP_ANGLE" == "90" -a "$TR" == "3300" -a "$TI" == "" \) ]
  then 
    # default case
    FN=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_vx100x100x%03d" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $SLICE_RES);
  else
  # spacial case
    FN=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_vx100x100x%03d_fa%03d_tr%03d_ti%03d" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $SLICE_RES $FLIP_ANGLE $TR $TI);  
    if [ ! -d "$CON_RDIR" ]; then mkdir -p "$CON_RDIR"; fi
    
 #   if [ "$TG" == "HC" -a  "$SLICE_THICKNESS" == "1" -a \( "$SPERCENT_INU" == "p" -o "$SPERCENT_INU" == "0" \)  -a "$INU_FIELD" == "A" -a "$APERCENT_INU" == "40" -a "$PERCENT_NOISE" == "3" ] 
 #   then
      # convert to nii
      mri_convert $F.mnc $CON_RDIR/$FN.nii > /dev/null 2>&1
 #   fi
    
    printf "> $FN done\n";   
    continue     
  fi
  
   # convert to nii
   mri_convert $F.mnc $FN.nii > /dev/null 2>&1

  FD=$(pwd);
  
    # /Applications/MATLAB_R2013b.app/bin/matlab -nodesktop -nodisplay -r "addpath(fullfile(spm('dir'),'toolbox','vbm12')); cd(fullfile(spm('dir'),'toolbox','vbm12','private')); vbm_tst_reduceRes(fullfile('$FD','$FN.nii'),[2 2 2;1 1 2;2 2 1]); exit" #> /dev/null 2>&1
    NFN=$(find *.nii  -depth 0); mv $NFN $FULL_RDIR/$TT/
 continue
  
  # cp for NIR (noise inhomogeneity resolution) directory ... 64 cases
  #     N   C   R
  # N   5   2   2   =  20 
  # I   2   6   2   =  24 
  # R   2   2   5   =  20 
  #                    64 * 3 Felder = 196
  if [ "$TG" == "HC" -a  "$SLICE_THICKNESS" == "1" -a \( "$SPERCENT_INU" == "p" -o "$SPERCENT_INU" == "0" \) -a \( "$APERCENT_INU" -le "20" -o "$PERCENT_NOISE" -ge "1" \) ] #-a \( "$INU_FIELD" == "A" -o "$INU_FIELD" == "0" \) 
  then
    /Applications/MATLAB_R2013b.app/bin/matlab -nodesktop -nodisplay -r "addpath(fullfile(spm('dir'),'toolbox','vbm12')); cd(fullfile(spm('dir'),'toolbox','vbm12','private')); vbm_tst_reduceRes(fullfile('$FD','$FN.nii'),[2 2 2;1 1 2;2 2 1]); exit" #> /dev/null 2>&1
    # if [ "$PERCENT_NOISE" == "3" -a  "$APERCENT_INU" == "40" ]  
    # then 
    #   printf ">>";
    #   /Applications/MATLAB_R2013b.app/bin/matlab -nodesktop -nodisplay  -r "vbm_tst_reduceRes('$FN.nii'); exit" > /dev/null 2>&1
    # fi
    # find *i.nii -depth 0 -delete
    #IFN=$(find *i.nii -depth 0); rm $IFN;  #mv $IFN $NIR_RDIR/$TT/
    RFN=$(find *r.nii -depth 0); mv $RFN $NIR_RDIR/$TT/
    NFN=$(find *.nii  -depth 0); cp $NFN $NIR_RDIR/$TT/
  fi

 
   # cp for CON directory - 6 cases
  if [ "$TG" == "HC" -a "$APERCENT_INU" == "40" -a "$SPERCENT_INU" == "p" -a "$INU_FIELD" == "A" -a  "$SLICE_THICKNESS" == "1" -a \( "TT" == "T2" -o "TT" == "PD" \) ]
  then
    FN2=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_vx100x100x%03d_fa%03d_tr%03d_ti%03d" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $SLICE_RES $FLIP_ANGLE $TR $TI);  
    cp $FN.nii $CON_RDIR/$FN2.nii; 
  fi
  
  # cp for noise directory - 6 cases
  if [ "$TG" == "HC" -a "$APERCENT_INU" == "40" -a "$SPERCENT_INU" == "p" -a "$INU_FIELD" == "A" -a  "$SLICE_THICKNESS" == "1" ]
  then
    cp $FN.nii $NOISE_RDIR/$TT/; 
  fi

  # cp for bias directory - 6 cases
  if [ "$TG" == "HC" -a  "$PERCENT_NOISE" == "3" -a \( "$SPERCENT_INU" == "p" -o "$SPERCENT_INU" == "0" \) -a \( "$INU_FIELD" == "A" -o "$INU_FIELD" == "0" \) -a  "$SLICE_THICKNESS" == "1" ]
  then
    cp $FN.nii $BIAS_RDIR/$TT/; 
  fi
  
  # cp for bias field directory - 6 cases (pA,nA,pB,nB,pC,nC)
  if [ "$TG" == "HC" -a  "$PERCENT_NOISE" == "3" -a "$APERCENT_INU" == "40" -a  "$SLICE_THICKNESS" == "1" ]
  then
    cp $FN.nii $FIELD_RDIR/$TT/; 
  fi

  # cp for res directory - 5 cases (1,3,5,7 mm slice thickness)
  if [ "$TG" == "HC" -a  "$PERCENT_NOISE" == "3" -a  "$APERCENT_INU" == "40" -a "$SPERCENT_INU" == "p" -a "$INU_FIELD" == "A" ]
  then
    cp $FN.nii $BIAS_RDIR/$TT/; 
  fi
    
  # copy for ground truth directory - 1 case
  if [ "$TG" == "HC" -a  "$PERCENT_NOISE" == "0" -a "$APERCENT_INU" == "0" -a "$SLICE_THICKNESS" == "1" ]; 
  then 
    cp $FN.nii $GT_RDIR/$TT/; 
  fi
  
  # copy for ground truth directory - 1 case
  if [ "$TT" == "T1w" -a "$TG" == "HC" -a  "$PERCENT_NOISE" == "3" -a "$APERCENT_INU" == "40" -a "$SPERCENT_INU" == "p" -a "$INU_FIELD" == "A" -a "$SLICE_THICKNESS" == "1" ]; 
  then 
    cp $FN.nii $CON_RDIR/$TT/; 
  fi
  
  # copy for worst case directory - 6 cases
  if [ "$TG" == "HC" -a  "$PERCENT_NOISE" == "9" -a "$APERCENT_INU" == "100" -a "$SLICE_THICKNESS" == "1" ]; 
  then 
    cp $FN.nii $WC_RDIR/$TT/; 
  fi
  
  # create own reduced resolutions
  if [ "$TG" == "HC" -a "$PERCENT_NOISE" == "3" -a  "$APERCENT_INU" == "40" -a "$SPERCENT_INU" == "p" -a "$INU_FIELD" == "A"  -a "$SLICE_THICKNESS" == "1" ] 
  then 
    /Applications/MATLAB_R2013b.app/bin/matlab -nodesktop -nodisplay  -r "addpath(fullfile(spm('dir'),'toolbox','vbm12')); cd(fullfile(spm('dir'),'toolbox','vbm12','private')); vbm_tst_reduceRes(fullfile('$FD','$FN.nii'));; exit" > /dev/null 2>&1
    IFN=$(find *i.nii -depth 0); mv $IFN $RESI_RDIR/$TT/
    RFN=$(find *r.nii -depth 0); mv $RFN $RESR_RDIR/$TT/
    NFN=$(find *.nii  -depth 0);         cp $NFN $RESR_RDIR/$TT/
  
  # another way to reduce image - more standard, but wit increased noise and bias because of the resampling rather than PVE reduction 
    if [ "0" == "1" ]
    then 
      echo simple resampling  
      for R in $RES
      do 
      # isotropic
        RX=$(echo $R*100 | bc | tr . ,); RY=$(echo $R*100 | bc | tr . ,); RZ=$(echo $R*100 | bc | tr . ,); 
        FN2=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_res%03.0fx%03.0fx%03.0f" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $RX $RY $RZ);
        caret_command -volume-resample "$FN.nii" "$RESR_RDIR/$TT/$FN2.nii"  $R $R $R INTERP_CUBIC
        caret_command -volume-resample "$RESR_RDIR/$TT/$FN2.nii"  "$RESI_RDIR/$TT/${FN2}I.nii" 1 1 1 INTERP_CUBIC > /dev/null 2>&1
      # sliceresolution
        RX=$(echo $R*100 | bc | tr . ,); RY=$(echo $R*100 | bc | tr . ,); RZ=$(echo 1*100 | bc | tr . ,);
        FN2=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_res%03.0fx%03.0fx%03.0f" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $RX $RY $RZ);
        caret_command -volume-resample "$FN.nii" "$RESR_RDIR/$TT/$FN2.nii"  $R $R 1 INTERP_CUBIC
        caret_command -volume-resample "$RESR_RDIR/$TT/$FN2.nii"  "$RESI_RDIR/$TT/${FN2}I.nii"  1 1 1 INTERP_CUBIC > /dev/null 2>&1
      # slicethickness
        RX=$(echo 1*100 | bc | tr . ,); RY=$(echo 1*100 | bc | tr . ,); RZ=$(echo $R*100 | bc | tr . ,);
        FN2=$(printf "BWPC_%s_%s_pn%01d_rf%03d%s%s_res%03.0fx%03.0fx%03.0f" $TG $TT $PERCENT_NOISE $APERCENT_INU $SPERCENT_INU $INU_FIELD $RX $RY $RZ);
        caret_command -volume-resample "$FN.nii" "$RESR_RDIR/$TT/$FN2.nii"  1 1 $R INTERP_CUBIC
        caret_command -volume-resample "$RESR_RDIR/$TT/$FN2.nii" "$RESI_RDIR/$TT/${FN2}I.nii"  1 1 1 INTERP_CUBIC > /dev/null 2>&1
      done
    fi
  fi


  # move to final directory - HC or MS
  if [ "$TG" == "HC" ] 
  then
    rm $FN.nii
    #mv $FN.nii $FULL_RDIR/$TT
  else
    mv $FN.nii $MS_RDIR/$TT/; 
  fi
  printf "> $FN done\n";    

done
