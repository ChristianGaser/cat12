#! /bin/bash
# ______________________________________________________________________
#
# Christian Gaser, Robert Dahnke
# Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
# Departments of Neurology and Psychiatry
# Jena University Hospital
# ______________________________________________________________________
# $Id$

########################################################
# global parameters
########################################################
version='cat_standalone.sh $Id$'

cwd=$(dirname "$0")

if [ ! -n "$SPMROOT" ]; then
  SPMROOT=$(dirname "${cwd}")
fi

# get cat12 dir
ARCH=`uname`
if [ "$ARCH" == "Darwin" ]; then
  cat12_dir="${SPMROOT}/spm12.app/Contents/MacOS/spm12/toolbox/cat12" 
else
  cat12_dir="your_folder/spm12/toolbox/cat12" 
fi

# default values
standalone=1
expert=0
fg=0
matlab=matlab # you can use other matlab versions by changing the matlab parameter

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_files
  run_cat

  exit 0
}

########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg
  count=0
  paras=

  if [ $# -lt 1 ]; then
    help
    exit 1
  fi
    
  while [ $# -gt 0 ]; do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2`"
    paras="$paras $optname $optarg"
    case "$1" in
        --batch* | -b*)
            exit_if_empty "$optname" "$optarg"
            BATCHFILE=$optarg
            shift
            ;;
        --mcr* | -m*)
            exit_if_empty "$optname" "$optarg"
            MCRROOT=$optarg
            shift
            ;;
        --spm* | -s*)
            exit_if_empty "$optname" "$optarg"
            SPMROOT=$optarg
            shift
            ;;
        --arg1* | -a1*)
            exit_if_empty "$optname" "$optarg"
            ARG1=$optarg
            shift
            ;;
        --arg2* | -a2*)
            exit_if_empty "$optname" "$optarg"
            ARG2=$optarg
            shift
            ;;
        --arg3* | -a3*)
            exit_if_empty "$optname" "$optarg"
            ARG3=$optarg
            shift
            ;;
        --add* | -a*)
            exit_if_empty "$optname" "$optarg"
            add_to_batch="$optarg"
            shift
            ;;
        --fg* | -fg*)
            exit_if_empty "$optname" "$optarg"
            fg=1
            ;;
        --no-s* | -ns*)
            exit_if_empty "$optname" "$optarg"
            standalone=0
            ;;
        --e* | -e*)
            exit_if_empty "$optname" "$optarg"
            expert=1
            ;;
        -h* | --h* | -v | --version | -V)
            help
            exit 1
            ;;
        -*)
            echo "`basename $0`: ERROR: Unrecognized option \"$1\"" >&2
            ;;
        *)
            ARRAY[$count]=$1
            ((count++))
            ;;
    esac
    shift
  done
}

########################################################
# check arguments
########################################################

exit_if_empty ()
{
  local desc val

  desc="$1"
  shift
  val="$*"

}

########################################################
# check files
########################################################

check_files ()
{
  
  # check for SPM parameter
  if [ ! -n "$SPMROOT" ]; then
    echo "No SPM folder defined."
    help
    exit 1  
  fi
  
  # check this only for standalone version
  if [ $standalone == 1 ]; then
    # check for MCR parameter
    if [ ! -n "$MCRROOT" ]; then
      echo "No MCR folder defined."
      help
      exit 1  
    fi
    
    # check for MCR folder
    if [ ! -d "$MCRROOT" ]; then
      echo "No MCR folder found."
      help
      exit 1  
    fi

    # check for SPM folder
    if [ ! -f "$SPMROOT/run_spm12.sh" ]; then
      echo "File $SPMROOT/run_spm12.sh not found found."
      help
      exit 1  
    fi
  else
    # we use the same flag as for MCRROOT, but here for Matlab command
    if [ -n "$MCRROOT" ]; then
      eval "matlab=\"$MCRROOT\";"
    fi
    
    found=`which "${matlab}" 2>/dev/null`
    if [ ! -n "$found" ]; then
      echo $matlab not found.
      exit 1
    fi
  fi
  
  # check for batch file
  if [ ! -n "$BATCHFILE" ]; then
    echo "No batch file defined."
    help
    exit 1  
  fi

  # check for files
  i=0
  while [ "$i" -lt "$count" ]
  do
    if [ ! -f "${ARRAY[$i]}" ]; then
      if [ ! -L "${ARRAY[$i]}" ]; then
        echo ERROR: File ${ARRAY[$i]} not found
        help
        exit 1
      fi
    fi
    ((i++))
  done

}

########################################################
# run CAT
########################################################

run_cat ()
{
  
  # if no files are given expect that file name is defined
  # in batch file and execute that file
  if [ "$count" -eq "0" ] && [ $standalone == 1 ] ; then
    eval "\"${SPMROOT}/run_spm12.sh\"" $MCRROOT "batch" $BATCHFILE
    exit 0
  fi
  
  if [ "$count" -eq "0" ]  ; then
    c2=1
    c3=2
    c4=3
  else
    c2=2
    c3=3
    c4=4
  fi

  # create temporary batch file
  TMP=/tmp/cat_$$.m

  # copy everything except rows with UNDEFINED to temp file 
  grep -v "<UNDEFINED>" $BATCHFILE > $TMP
  
  if [ ! "$count" -eq "0" ]  ; then
    # extract parameter name of data structure (1st occurance of "<UNDEFINED>")
    data=`grep -m 1 "<UNDEFINED>" $BATCHFILE | cut -f1 -d'='| sed -e 's,%,,'`

    # surface data need an additional curly bracket
    if grep -q -e "\.datalong" $BATCHFILE ; then
      echo "$data = {{" >> $TMP
    else
      echo "$data = {" >> $TMP
    fi
  fi
  
  # extract parameter name of optional argument(s) (additional occurances of "<UNDEFINED>")
  if [ -n "$ARG1" ]; then # ARG1 defined?
    param1=`grep -m $c2 "<UNDEFINED>" $BATCHFILE | tail -n 1 | cut -f1 -d'=' | sed -e 's,%,,'`
    # extract parameter name of optional argument (3rd occurance of "<UNDEFINED>")
    if [ -n "$ARG2" ]; then # ARG2 defined?
      param2=`grep -m $c3 "<UNDEFINED>" $BATCHFILE | tail -n 1 | cut -f1 -d'=' | sed -e 's,%,,'`
      # extract parameter name of optional argument (4th occurance of "<UNDEFINED>")
      if [ -n "$ARG3" ]; then # ARG3 defined?
        param3=`grep -m $c4 "<UNDEFINED>" $BATCHFILE | tail -n 1 | cut -f1 -d'=' | sed -e 's,%,,'`
      fi
    fi
  fi
    
  i=0
  ARG_LIST=""
  while [ "$i" -lt "$count" ]; do

    # check whether absolute or relative names were given
    if [ ! -f ${ARRAY[$i]} ];  then
      if [ -f ${PWD}/${ARRAY[$i]} ]; then
        FILE=${PWD}/${ARRAY[$i]}
      fi
    else
      FILE=${ARRAY[$i]}
    fi

    # add file list
    echo "'${FILE}'" >> $TMP
              
    ((i++))
  done

  if [ ! "$count" -eq "0" ]  ; then
    # surface data need an additional curly bracket
    if grep -q -e "\.datalong" $BATCHFILE ; then
      echo "     }};" >> $TMP
    else
      echo "     };" >> $TMP
    fi
  fi
    
  if [ -n "$ARG1" ]; then # ARG1 defined?
    echo "$param1 = $ARG1 ;" >> $TMP
    if [ -n "$ARG2" ]; then # ARG2 defined?
      echo "$param2 = $ARG2 ;" >> $TMP
      if [ -n "$ARG3" ]; then # ARG3 defined?
        echo "$param3 = $ARG3 ;" >> $TMP
      fi
    fi
  fi
  
  # add optional lines to batch file
  if [ -n "$add_to_batch" ]; then
    echo "${add_to_batch}" >> ${TMP}
  fi

  if [ $standalone == 1 ]; then
    eval "\"${SPMROOT}/run_spm12.sh\"" $MCRROOT "batch" $TMP
    rm $TMP
    exit 0
  else
    DIR=$(dirname "${TMP}")
    
    # we have to check where spm.m is found
    if [ ! -f ${SPMROOT}/spm.m ]; then
      SPMROOT=$(dirname "${SPMROOT}")
      SPMROOT=$(dirname "${SPMROOT}")
    fi
    
    BATCHFILE=$(basename "${TMP}")
    BATCHFILE=$(echo "${BATCHFILE}" | cut -f1 -d'.')
    cat12_dir="${SPMROOT}/toolbox/cat12" 
    export MATLABPATH=${SPMROOT}:${cat12_dir}:${DIR}
    eval "COMMAND=\"$BATCHFILE\";"
    
    if [ $expert == 1 ]; then
      COMMAND="spm; spm_get_defaults; cat_get_defaults; global defaults cat matlabbatch; cat12('expert'); $COMMAND;spm_jobman('run',matlabbatch); exit;";
    else
      COMMAND="spm; spm_get_defaults; cat_get_defaults; global defaults cat matlabbatch;$COMMAND;spm_jobman('run',matlabbatch); exit;";
    fi
    if [ "$fg" -eq 0 ]; then
      nohup nice ${matlab} -nodisplay -nosplash -r "$COMMAND"  2>&1 &
    else
      nohup nice  ${matlab} -nodisplay -nosplash -r "$COMMAND" 2>&1
    fi
  fi
}


########################################################
# help
########################################################

help ()
{

  # do not use a single dot
  if [ "$SPMROOT" == "." ]; then
    SPMROOT="SPMROOT"
  fi

cat <<__EOM__

USAGE:
   cat_standalone.sh filename(s) [-s spm_standalone_folder] [-m mcr_folder] [-b batch_file] 
                                 [-a1 additional_argument1] [-a2 additional_argument2]
                                 [-a add_to_batch] [-ns -e -fg]
   
   -s  <DIR>   | --spm <DIR>     SPM12 folder of standalone version (can be also defined by SPMROOT)
                                 Often you don't need to define that option because the SPM12 folder
                                 is automatically found.
   -m  <DIR>   | --mcr <DIR>     Matlab Compiler Runtime (MCR) folder (can be also defined by MCRROOT)
   -b  <FILE>  | --batch <FILE>  batch file to execute
   -a1 <STRING>| --arg1 <STRING> 1st additional argument (otherwise use defaults or batch)
   -a2 <STRING>| --arg2 <STRING> 2nd additional argument (otherwise use defaults or batch)
   -a3 <STRING>| --arg3 <STRING> 3rd additional argument (otherwise use defaults or batch)
   -a  <STRING>| --add  <STRING> add option to batch file
   
   options for calling standard Matlab mode
   -ns         | --no-standalone call the standard Matlab version instead of the standalone version
   -e          | --expert        call CAT12 in expert mode (needed for using standalone batches!)
   -fg         | --fg            do not run matlab process in background
   -m  <FILE>  | --m <FILE>      Matlab command (matlab version) (default $matlab)

   The first occurance of the parameter "<UNDEFINED>" in the batch file will be replaced by the
   list of input files. You can use the existing batch files in this folder or create your own batch 
   file with the SPM12 batch editor and leave the data field undefined. Please note that for creating
   your own batch file CAT12 has to be called in expert mode because the CAT12 standalone installation 
   will only run in expert mode to allow more options.
   See cat_standalone_segment.m for an example. 
   
   You can also define one or two optional arguments to change other parameters that are indicated by 
   "<UNDEFINED>" in the batch file. Please take care of the order of the "<UNDEFINED>" fields in the 
   batch file! If no additional arguments are defined the default values are used.
   Also, you must use multiple quotes if the argument is a string (e.g. " 'your_string' ").
   
   If you use a computer cluster it is recommended to use the batch files to only process one data set 
   and use a job or queue tool to call the (single) jobs on the cluster.
   
   Although this tool is mainly intended for calling scripts for the standalone version of Matlab without 
   a Matlab license, you can also use it to call the standard version of Matlab. If you have a Matlab
   license, this has the advantage that you can use the same scripts as for the standalone version, but 
   you can run CAT12 without a GUI. Please note that standalone batches must be called in CAT12 expert
   mode. Of course, you can also create and use your own batches that you use regulary in CAT12 or
   SPM12. With this script you canrun these batches in headless mode without any display.
   
   PURPOSE:
   Command line call of (CAT12) batch files for SPM12 standalone installation

EXAMPLES
   -----------------------------------------------------------------------------------------------
   Dicom Import
     -a1 directory structure
     -a2 output directory
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_dicom2nii.m *.dcm \ 
       -a1 " 'patid_date' " -a2 "{'converted'}"
   Import DICOM files *.dcm and save converted nifti files in directory "converted" with structure 
   ./<PatientID>/<StudyDate-StudyTime>/<ProtocollName>
   Other options for directory structure are:
     'flat'       No directory hierarchy
     'series'     ./<ProtocollName>
     'patid_date' ./<PatientID>/<StudyDate-StudyTime>/<ProtocollName>
     'patid'      ./<PatientID>/<ProtocollName>
     'date_time'  ./<StudyDate-StudyTime>/<ProtocollName>
   Please note the multiple quotes for parameter a1.

   -----------------------------------------------------------------------------------------------
   De-facing
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_deface.m sTRIO*.nii
   Apply de-facing to sTRIO*.nii and save the files prefixed by "anon_".

   -----------------------------------------------------------------------------------------------
   Segmentation
     -a1 TPM
     -a2 Shooting template
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_segment.m sTRIO0001.nii
   Preprocess (segment) the single file sTRIO0001.nii using the default CAT12 preprocessing batch. 
   SPM12 standalone version is located in $SPMROOT and Matlab Compiler Runtime in
   /Applications/MATLAB/MATLAB_Runtime/v93.

   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_segment.m sTRIO000*.nii.gz \ 
       -a1 " '${cat12_dir}/templates_MNI152NLin2009cAsym/TPM_Age11.5.nii' " \ 
       -a2 " '${cat12_dir}/templates_MNI152NLin2009cAsym/Template_0_GS1mm.nii' "
   Unzip and preprocess (segment) the files sTRIO0001.nii.gz using the default CAT12 preprocessing 
   batch, but use the children TPM provided with CAT12 and a 1mm Shooting template (not provided 
   with CAT12). Please note that zipped file can only be handled with this standalone batch and also
   note the multiple quotes for parameter a1 and a2.

   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_segment.m sTRIO0001.nii \ 
       -a "matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;"
   Preprocess (segment) the single file sTRIO0001.nii using the default CAT12 preprocessing batch, 
   but skip surface estimation.

   -----------------------------------------------------------------------------------------------
   Longitudinal segmentation
     -a1 longitudinal model (0 - developmental; 1 - plasticity/learning; 2 - aging; 3 - save models 1 and 2)
     -a2 TPM
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_segment_long.m sTRIO000*.nii \ 
       -a1 "2"
   Preprocess (segment) the files sTRIO000*.nii with the longitudinal pipeline optimized for 
   detecting aging/developmental effects. In order to choose the longitudinal model optimized for 
   detecting small changes due to plasticity/learning change the a1 parameter to "1".

   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_segment_long.m sTRIO000*.nii \ 
       -a1 "1" -a2 " '${cat12_dir}/templates_MNI152NLin2009cAsym/TPM_Age11.5.nii' "
   Preprocess (segment) the files sTRIO000*.nii with the longitudinal pipeline optimized for 
   detecting plasticity/learning effects and use the children TPM provided with CAT12.
   Please note the multiple quotes for parameter a2.

   -----------------------------------------------------------------------------------------------
   Segmentation (simple mode)
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_simple.m sTRIO0001.nii
   Process the single file sTRIO0001.nii using the simple processing batch. 

   -----------------------------------------------------------------------------------------------
   Resample & smooth surfaces
     -a1 smoothing filter size surface values
     -a2 use 32k mesh from HCP (or 164k mesh from Freesurfer)
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_resample.m surf/lh.thickness.sTRIO0001 \ 
       -a1 "12" -a2 "1" 
   Resample and smooth the single thickness file lh.thickness.sTRIO0001 with 12mm and save the 
   resampled mesh as 32k mesh (HCP conform mesh). Only the left surface file has to be defined.
   The right hemisphere is processed automatically.

   -----------------------------------------------------------------------------------------------
   Smoothing
     -a1 smoothing filter size
     -a2 prepending string for smoothed file (e.g. 's6')
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_smooth.m mri/sTRIO*nii \ 
       -a1 "[6 6 6]" -a2 " 's6' "
   Smooth the volume files sTRIO*nii with 6mm and prepend the string "s6" to the smoothed files.
   Please note the multiple quotes for parameter a2.

   -----------------------------------------------------------------------------------------------
   Estimate and save quality measures for volumes or surfaces
     -a1 csv output filename
     -a2 enable global scaling with TIV (only for volumes meaningful)
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_get_quality.m mri/mwp1sTRIO*nii \ 
       -a1 " 'Quality_measures.csv' " -a2 "1"
   Estimate sample homogeneity (after preprocessing) using mean z-scores with global scaling with TIV 
   for the files mwp1sTRIO*nii and save quality measures in Quality_measures.csv for external analysis. 
   Processing of surface meshes is also supported.
   Please note the multiple quotes for parameter a1.

   -----------------------------------------------------------------------------------------------
   Estimate and save weighted overall image quality
     -a1 csv output filename
     -a2 enable global scaling with TIV (only for volumes meaningful)
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_get_IQR.m report/cat_*.xml \ 
       -a1 " 'IQR.txt' "
   Estimate weighted overall image quality (before preprocessing) using xml-files in report folder 
   and save IQR measures in IQR.txt for external analysis.
   Please note the multiple quotes for parameter a1.

   -----------------------------------------------------------------------------------------------
   Estimate mean/volume inside ROI
     -a1 output-file string
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_get_ROI_values.m label/catROI_*.xml \ 
       -a1 " 'ROI' " 
   Save mean volume values in mL (e.g. GM volume) or the mean surface values (e.g. thickness) for 
   all data catROI_*.xml in a csv-file. The csv-file is named "ROI_" followed by the atlas name
   and the name of the measure (e.g. Vgm).
   Please note the multiple quotes for parameter a1.
 
   -----------------------------------------------------------------------------------------------
   Estimate total intra-cranial volume (TIV) or all global tissue volumes (in ml)
     -a1 output filename
     -a2 save TIV only
     -a3 add filenames
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_getTIV.m report/cat_*.xml \ 
       -a1 " 'TIV.txt' " -a2 "1" -a3 "1"
   Estimate TIV only and save file names and values for each data set in TIV.txt.
   The parameter a3 allows to add file names to 1st column:
     0 - save values only; 1 - add file name; 2 - add folder and file name
   Please note the multiple quotes for parameter a1.

   -----------------------------------------------------------------------------------------------
   TFCE statistical estimation
     -a1 contrast number
     -a2 number of permutations
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -m /Applications/MATLAB/MATLAB_Runtime/v93 \ 
       -b ${cwd}/cat_standalone_tfce.m SPM.mat \ 
       -a1 "2" -a2 "20000"
   Call estimation of TFCE statistics for the given SPM.mat file for contrast number 2 with 20000 
   permutations.

   -----------------------------------------------------------------------------------------------
   Calling standard Matlab mode
   -----------------------------------------------------------------------------------------------
   cat_standalone.sh -ns -e -m matlab_R2017b \ 
       -b ${cwd}/cat_standalone_segment.m sTRIO0001.nii \ 
       -a "matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;"
   Preprocess (segment) the single file sTRIO0001.nii in standard Matlab mode (using CAT12 expert
   mode) and the default CAT12 preprocessing batch, but skip surface estimation. As Matlab command
   matlab_R2017b is used.

   -----------------------------------------------------------------------------------------------
   Parallelization
   -----------------------------------------------------------------------------------------------
   cat_parallelize.sh -p 8 -l /tmp \ 
       -c "cat_standalone.sh  -m /Applications/MATLAB/MATLAB_Runtime/v93 -b ${cwd}/cat_standalone_segment.m" sTRIO*.nii
   Parallelize CAT12 preprocessing by splitting all sTRIO*.nii files into 8 jobs 
   (processes) and save log file in /tmp folder. 

   The parameters SPMROOT and MCRROOT have to be defined (exported) to skip the use of the flags -s -m.

INPUT:
   nifti files or surface data

OUTPUT:
   processed images and optionally surfaces according to settings in cat_standalone_*.m

USED FUNCTIONS:
   cat_parallelize.sh
   SPM12 standalone version (compiled)
   CAT12 toolbox (compiled within SPM12 if installed)
   MATLAB Compiler Runtime R2017b (Version 9.3)

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
