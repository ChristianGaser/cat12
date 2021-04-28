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

matlab=matlab # you can use other matlab versions by changing the matlab parameter
cwd=$(dirname "$0")
cat12_dir=$cwd
spm12=$(dirname "$cwd")
spm12=$(dirname "$spm12")
LOGDIR=$PWD
CPUINFO=/proc/cpuinfo
ARCH=`uname`
time=`date "+%Y%b%d_%H%M"`
nicelevel=0
defaults_tmp=/tmp/defaults$$.m
no_mwp=0
no_surf=0
rp=0
TEST=0
fg=0

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
  get_no_of_cpus
  modifiy_defaults
  run_cat12

  exit 0
}


########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg
  count=0

  if [ $# -lt 1 ]; then
    help
    exit 1
  fi
  
  while [ $# -gt 0 ]; do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    paras="$paras $optname $optarg"
    case "$1" in
      --matlab* | -m*)
        exit_if_empty "$optname" "$optarg"
        matlab=$optarg
        shift
        ;;
      --defaults* | -d*)
        exit_if_empty "$optname" "$optarg"
        defaults=$optarg
        shift
        ;;
      --nprocesses* | -np*)
        exit_if_empty "$optname" "$optarg"
        NUMBER_OF_JOBS="-$optarg"
        shift
        ;; 
      --processes* | -p*)
        exit_if_empty "$optname" "$optarg"
        NUMBER_OF_JOBS=$optarg
        shift
        ;;
      --logdir* | -l*)
        exit_if_empty "$optname" "$optarg"
        LOGDIR=$optarg
        if [ ! -d $LOGDIR ]; then
          mkdir -p $LOGDIR
        fi
        shift
        ;;
      --no-mwp* | -nm*)
        exit_if_empty "$optname" "$optarg"
        no_mwp=1
        ;;
      --no-surf* | -ns*)
        exit_if_empty "$optname" "$optarg"
        no_surf=1
        ;;
      --rp* | -r*)
        exit_if_empty "$optname" "$optarg"
        rp=1
        ;;
      --nojvm | -nj)
        exit_if_empty "$optname" "$optarg"
        nojvm=" -nojvm "
        ;;
      --no-overwrite* | -no*)
        exit_if_empty "$optname" "$optarg"
        no_overwrite=$optarg
        shift
        ;;
      --n* | -n* | --nice* | -nice*)
        exit_if_empty "$optname" "$optarg"
        nicelevel=$optarg
        shift
        ;;
      --add* | -a*)
        exit_if_empty "$optname" "$optarg"
        add_to_defaults="$optarg"
        shift
        ;;
      --fg* | -fg*)
        exit_if_empty "$optname" "$optarg"
        fg=1
        ;;
      --files* | -f*)
        exit_if_empty "$optname" "$optarg"
        listfile=$optarg
        shift
        list=$(< $listfile);
        for F in $list; do
          ARRAY[$count]=$F
          ((count++))
        done
        ;;
      --s* | -s* | --shell* | -shell*)
        exit_if_empty "$optname" "$optarg"
        shellcommand=$optarg
        shift
        ;; 
      --tpm* | -tpm*)
        exit_if_empty "$optname" "$optarg"
        tpm=$optarg
        shift
        ;;
      --test* | -t*)
        TEST=1
        ;;
      --c* | -c* | --command* | -command*)
        exit_if_empty "$optname" "$optarg"
        matlabcommand=$optarg
        shift
        ;; 
      -h | --help | -v | --version | -V)
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

  if [ ! -n "$val" ]; then
    echo 'ERROR: No argument given with \"$desc\" command line argument!' >&2
  exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  
  if [ "$no_surf" -eq 1 ] && [ "$no_mwp" -eq 1 ] && [ "$rp" -eq 0 ]; then
    echo 'WARNING: You have deselected all outputs! Only the p0-image is saved.' >&2
  fi

  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]; then
    echo 'ERROR: No files given!' >&2
    help
    exit 1
  fi

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]; do
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
# get # of cpus
########################################################
# modified code from
# PPSS, the Parallel Processing Shell Script
# 
# Copyright (c) 2009, Louwrentius
# All rights reserved.

get_no_of_cpus () {

  if [ ! -n "$NUMBER_OF_JOBS" ]; then
    if [ "$ARCH" == "Linux" ]; then
      NUMBER_OF_PROC=`grep ^processor $CPUINFO | wc -l`
    elif [ "$ARCH" == "Darwin" ]; then
      NUMBER_OF_PROC=`sysctl -a hw | grep -w logicalcpu | awk '{ print $2 }'`
    elif [ "$ARCH" == "FreeBSD" ]; then
      NUMBER_OF_PROC=`sysctl hw.ncpu | awk '{ print $2 }'`
    else
      NUMBER_OF_PROC=`grep ^processor $CPUINFO | wc -l`
    fi
  
    if [ ! -n "$NUMBER_OF_PROC" ]; then
      echo "$FUNCNAME ERROR - number of CPUs not obtained. Use -p to define number of processes."
      exit 1
    fi
  
    # use all processors if not defined otherwise
    if [ ! -n "$NUMBER_OF_JOBS" ]; then
      NUMBER_OF_JOBS=$NUMBER_OF_PROC
    fi
    
    if [ $NUMBER_OF_JOBS -le -1 ]; then
      NUMBER_OF_JOBS=$(echo "$NUMBER_OF_PROC + $NUMBER_OF_JOBS" | bc)
    if [ "$NUMBER_OF_JOBS" -lt 1 ]; then
      NUMBER_OF_JOBS=1
    fi
    fi
    if [ "$NUMBER_OF_JOBS" -gt "$NUMBER_OF_PROC" ]; then
      NUMBER_OF_JOBS=$NUMBER_OF_PROC
    fi
    echo "Found $NUMBER_OF_PROC processors. Use $NUMBER_OF_JOBS."
    echo
  fi
}

########################################################
# modify defaults
########################################################

modifiy_defaults ()
{

  pwd=$PWD

  # argument empty?
  if [ -n "${defaults}" ]; then
    # check whether absolute or relative names were given
    if [ -f "${pwd}/${defaults}" ]; then
      defaults="${pwd}/${defaults}"
    fi

    # check whether defaults file exist
    if [ ! -f "${defaults}" ]; then
      echo Default file "$defaults" not found.
      exit
    fi
  else
    defaults=${cat12_dir}/cat_defaults.m
  fi

  cp ${defaults} ${defaults_tmp}
  
  if [ "$no_surf" -eq 1 ]; then
    echo "cat.output.surface = 0;" >> ${defaults_tmp}
  else
    echo "cat.output.surface = 1;" >> ${defaults_tmp}
  fi

  if [ "$no_mwp" -eq 1 ]; then
    echo "cat.output.GM.mod      = 0;" >> ${defaults_tmp}
    echo "cat.output.WM.mod      = 0;" >> ${defaults_tmp}
    echo "cat.output.ROI         = 0;" >> ${defaults_tmp}
    echo "cat.output.bias.warped = 0;" >> ${defaults_tmp}
    echo "cat.output.warps       = [0 0];" >> ${defaults_tmp}
  else
    echo "cat.output.GM.mod = 1;" >> ${defaults_tmp}
    echo "cat.output.WM.mod = 1;" >> ${defaults_tmp}
    echo "cat.output.ROI    = 1;" >> ${defaults_tmp}
  fi

  if [ "$rp" -eq 1 ]; then
    echo "cat.output.GM.dartel = 2;" >> ${defaults_tmp}
    echo "cat.output.WM.dartel = 2;" >> ${defaults_tmp}
  else
    echo "cat.output.GM.dartel = 0;" >> ${defaults_tmp}
    echo "cat.output.WM.dartel = 0;" >> ${defaults_tmp}
  fi
  
  if [ -n "$tpm" ]; then
    # check whether absolute or relative tpm was given
    if [ -f "${pwd}/${tpm}" ]; then
      tpm="${pwd}/${tpm}"
    fi
    echo "cat.opts.tpm = {'${tpm}'};" >> ${defaults_tmp}
  fi

  if [ -n "$add_to_defaults" ]; then
    echo "${add_to_defaults}" >> ${defaults_tmp}
  fi
}

########################################################
# run cat12
########################################################

run_cat12 ()
{
  pwd=$PWD
  
  # we have to go into the toolbox folder to find matlab files
  cd $cwd
  
  if [ ! -n "${LOGDIR}" ]; then
    LOGDIR=$(dirname "${ARRAY[0]}")
  fi
  
  export MATLABPATH=${spm12}

  SIZE_OF_ARRAY="${#ARRAY[@]}"

  i=0
  j=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]; do

    # check whether absolute or relative names were given
    if [ ! -f ${ARRAY[$i]} ]; then
      if [ -f "${pwd}/${ARRAY[$i]}" ]; then
        FILE="${pwd}/${ARRAY[$i]}"
      fi
    else
      FILE=${ARRAY[$i]}
    fi

    # replace white spaces
    FILE=$(echo "$FILE" | sed -e 's/ /\\ /g')
    
    # check whether processed files exist if no-overwrite flag is used
    if [ -n "${no_overwrite}" ]; then
      dn=$(dirname "$FILE")
      bn=$(basename "$FILE" |cut -f1 -d'.')
      processed=`eval ls "${dn}/${no_overwrite}${bn}*" 2>/dev/null`
    fi

    if [ ! -n "${processed}" ]; then
      if [ ! -n "${ARG_LIST[$i]}" ]; then
        ARRAY2[$j]="$FILE"
      else
        ARRAY2[$j]="${ARRAY2[$i]} $FILE"
      fi      
      ((j++))
    else
      echo Skip processing of ${FILE}
    fi
    ((i++))
  done

  SIZE_OF_ARRAY="${#ARRAY2[@]}"
  BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_JOBS ))
  ARG_LIST="" 
  
  # split files and prepare tmp-file with filenames
  TMP=/tmp/cat_$$
  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]; do
    count=$((10000* $i / $BLOCK ))

    FILE=${ARRAY2[$i]}
    
    ARG_LIST[$count]="${ARG_LIST[$count]} $FILE"
    
    if [ "$TEST" -eq 0  ]; then
      echo ${FILE} >> ${TMP}${count}
    else
      echo ${FILE}
    fi
    ((i++))
  done
  vbmlog="${LOGDIR}/cat_${HOSTNAME}_${time}"

  i=0
  while [ "$i" -lt "$NUMBER_OF_JOBS" ]; do
    if [ -n "${ARG_LIST[$i]}" ] && [ "$TEST" -eq 0 ]; then
      j=$(($i+1))
      if [ ! -n "$matlabcommand" ]; then
        COMMAND="cat_batch_cat('${TMP}${i}','${defaults_tmp}')"
      else
        for F in ${ARG_LIST[$i]} ; do 
          CFILES=$CFILES","\'$F\';
        done
        CFILES=$(echo $CFILES | cut -c 2-);
        matlabcommand2=$matlabcommand
        matlabcommand2=$(echo $matlabcommand2 |sed 's/CFILES/$CFILES/g');
        eval "COMMAND=\"$matlabcommand2\";"
        COMMAND="try, spm; spm_get_defaults; cat_get_defaults; global defaults cat matlabbatch; $COMMAND; catch caterr, sprintf('\n%s\nVBM Preprocessing error: %s:\n%s\n', repmat('-',1,72),caterr.identifier,caterr.message,repmat('-',1,72)); for si=1:numel(caterr.stack), cat_io_cprintf('err',sprintf('%5d - %s\n',caterr.stack(si).line,caterr.stack(si).name)); end; cat_io_cprintf('err',sprintf('%s\\n',repmat('-',1,72))); exit; end; fprintf('VBM batch processing done.'); exit;";
      fi
      SHCOMMAND="$shellcommand ${ARG_LIST[$i]}"       
            
      echo Calculate
      for F in ${ARG_LIST[$i]}; do echo $F; done
      # File Output
      echo ---------------------------------- >> "${vbmlog}_${j}.log"
      date                  >> "${vbmlog}_${j}.log"
      echo ---------------------------------- >> "${vbmlog}_${j}.log"
      echo                  >> "${vbmlog}_${j}.log"
      echo Calling string of this batch:   >> "${vbmlog}_${j}.log"
      echo " $0 $paras"           >> "${vbmlog}_${j}.log"
      echo                  >> "${vbmlog}_${j}.log"
      echo MATLAB command of this batch:   >> "${vbmlog}_${j}.log"
      echo " $COMMAND"            >> "${vbmlog}_${j}.log"
      echo                  >> "${vbmlog}_${j}.log"
      echo Shell command of this batch:    >> "${vbmlog}_${j}.log"
      echo " $SHCOMMAND"           >> "${vbmlog}_${j}.log"
      echo                  >> "${vbmlog}_${j}.log"
      
      if [ ! -n "$shellcommand" ]; then
        # do nohup in background or not
        if [ "$fg" -eq 0 ]; then
          nohup nice -n $nicelevel ${matlab} -nodisplay "$nojvm" -nosplash -r "$COMMAND" >> "${vbmlog}_${j}.log" 2>&1 &
        else
          nohup nice -n $nicelevel ${matlab} -nodisplay "$nojvm" -nosplash -r "$COMMAND" >> "${vbmlog}_${j}.log" 2>&1
        fi
      else
        # do nohup in background or not
        if [ "$fg" -eq 0 ]; then
          nohup nice -n $nicelevel $SHCOMMAND >> "${vbmlog}_${j}.log" 2>&1 &
        else
          nohup nice -n $nicelevel $SHCOMMAND >> "${vbmlog}_${j}.log" 2>&1
        fi
      fi
      echo Check "${vbmlog}_${j}.log" for logging information
      echo
    fi
    ((i++))
  done

  exit 0
}

########################################################
# check if matlab exist
########################################################

check_matlab ()
{
  found=`which "${matlab}" 2>/dev/null`
  if [ ! -n "$found" ]; then
    echo $matlab not found.
    exit 1
  fi
}

########################################################
# help
########################################################

help ()
{
 
get_no_of_cpus
cat <<__EOM__

USAGE:
 cat_batch_cat.sh filenames|filepattern [-m matlab_command] [-d default_file] [-l log_folder] 
                    [-p number_of_processes] [-p number_of_processes] [-tpm TPM-file] 
                    [-ns] [-nm] [-rp] [-no output_pattern] [-n nicelevel] 
                    [-s shell_command -f files_for_shell] [-c matlab_command] 
                    [-a add_to_defaults] [-t] [-fg] [-noj]
 
  -m  <FILE>  | --matlab <FILE>         matlab command (matlab version) (default $matlab)
  -d  <FILE>  | --defaults <FILE>       optional default file (default ${cat12_dir}/cat_defaults.m)
  -l  <FILE>  | --logdir                directory for log-file (default $LOGDIR)
  -p  <NUMBER>| --processes <NUMBER>    number of parallel jobs (=number of processors)
                                       (default $NUMBER_OF_JOBS)
  -np <NUMBER>| --nprocesses <NUMBER>   set number of jobs by number_of_processors - number_of_processes
                                       (=number of free processors)
  -tpm <FILE> | --tpm <FILE>            define own TPM
  -a  <STRING>| --add                   add option to default file
  -ns         | --no-surf               skip surface and thickness estimation
  -nm         | --no-mwp                skip estimating modulated and warped segmentations and ROI measures
  -rp         | --rp                    additionally estimate affine registered segmentations
  -no <STRING>| --no-overwrite <STRING> do not overwrite existing results
  -n  <NUMBER>| --nice <NUMBER>         nice level (default 0)
  -s  <STRING>| --shell <STRING>        shell command to call other shell scripts
  -f  <FILES> | --files <FILES>         files to process with shell command
  -c  <STRING>| --command <STRING>      alternative matlab function that can be called such as the SANLM-filter
  -t          | --test                  do not call command, but print files to be processed
  -fg         | --fg                    do not run matlab process in background
  -nj         | --nojvm                 supress call of jvm using the -nojvm flag
 
 Only one filename or pattern is allowed. This can be either a single file or a pattern
 with wildcards to process multiple files. 

PURPOSE:
 Command line call of CAT12 segmentation

EXAMPLE
 cat_batch_cat.sh ${spm12}/canonical/single_subj_T1.nii
   This command will process only the single file single_subj_T1.nii. 
 
 cat_batch_cat.sh ${spm12}/canonical/single_subj_T1.nii -d your_cat_defaults_file.m
   This command will process only the single file single_subj_T1.nii. The defaults defined
   in your_cat_defaults_file.m is used instead of cat_defaults.m.

 cat_batch_cat.sh ${spm12}/canonical/*152*.nii
   Using wildcards all files containing the term "152" are processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii and $NUMBER_OF_JOBS parallel
   jobs are used.

 cat_batch_cat.sh ${spm12}/canonical/*152*.nii -no "mri/mwp1"
   Using wildcards all files containing the term "152" are processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii and $NUMBER_OF_JOBS parallel
   jobs are used. If processed files "mwp1*" in the subfolder "mri" are
   found the processing will be skipped.

 cat_batch_cat.sh ${spm12}/canonical/*152*.nii -p 2 -m /usr/local/bin/matlab7
   Using wildcards all files containing the term "152" are processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii and 2 parallel jobs
   jobs are used. As matlab-command /usr/local/bin/matlab7 is used.
 
 cat_batch_cat.sh ${spm12}/canonical/single_subj_T1.nii -ns -nm -rp -a "cat.extopts.WMHC = 3;"
   This command will process only the single file single_subj_T1.nii with the defaults in cat_defaults.m and
   the additional option for handling WMHs as separate class. No surfaces and modulated and warped segmentations
   are estimated. Only the affine registered segmentations are saved.
 
 cat_batch_cat.sh ${spm12}/canonical/single_subj_T1.nii -tpm ${cat12_dir}/templates_volumes/TPM_Age11.5.nii
   This command will process only the single file single_subj_T1.nii with the defaults in cat_defaults.m
   and the children template that is provided with cat12.

 cat_batch_cat.sh -p 2 -c "cat_vol_sanlm(struct('data',char(CFILES),'prefix','sanlm_'))" /Volumes/4TBWD/raw-cg/r[12][0-9][0-9][0-9]*.nii
   This command will call the SANLM-filter using the given files, that have to be indicated with CFILES
   as first argument. As prefix 'sanlm_' is used.
 

INPUT:
 analyze or nifti files

OUTPUT:
 segmented images according to settings in cat_defaults.m
 ${LOGDIR}/cat_${HOSTNAME}_${time}.log for log information

USED FUNCTIONS:
 cat_batch_cat.m
 CAT12 toolbox
 SPM12

SETTINGS
 matlab command: $matlab
  
This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
