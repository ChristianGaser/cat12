#! /bin/bash
# Call CAT12 longitudinal pipeline from shell
# ______________________________________________________________________
#
# Christian Gaser, Robert Dahnke
# Structural Brain Mapping Group (https://neuro-jena.github.io)
# Departments of Neurology and Psychiatry
# Jena University Hospital
# ______________________________________________________________________
# $Id$

########################################################
# global parameters
########################################################

matlab=matlab   # you can use other matlab versions by changing the matlab parameter
cwd=$(dirname "$0")
cat12_dir=$cwd
spm12=$(dirname "$cwd")
spm12=$(dirname "$spm12")
LOGDIR=$PWD
export_dartel=0
output_surface=1
long_model=1
time=`date "+%Y%b%d_%H%M"`
defaults_tmp=/tmp/defaults$$.m
fg=0
bids=0
bids_folder=

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
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
      --defaults* | -d*)
        exit_if_empty "$optname" "$optarg"
        defaults=$optarg
        shift
        ;;
      --model* | -model*)
        exit_if_empty "$optname" "$optarg"
        long_model=$optarg
        shift
        ;;
      --matlab* | -m*)
        exit_if_empty "$optname" "$optarg"
        matlab=$optarg
        shift
        ;;
      --export-dartel* | -e*)
        export_dartel=1
        ;;
      --no-surf* | -ns*)
        output_surface=0
        ;;
      --nojvm | -nj*)
        nojvm=" -nojvm "
        ;;
      --fg* | -fg*)
        fg=1
        ;;
      --bids_folder* | --bids-folder* | -bf*)
        exit_if_empty "$optname" "$optarg"
        bids_folder=$optarg
        shift
        ;;
      --b* | -b*)
        exit_if_empty "$optname" "$optarg"
        bids=1
        ;;
      --logdir* | -log*)
        exit_if_empty "$optname" "$optarg"
        LOGDIR=$optarg
        if [ ! -d $LOGDIR ]; then
           mkdir -p $LOGDIR
        fi
        shift
        ;;
      --large* | -l*)
        long_model=2
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
 
  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]; then
    echo 'ERROR: No files given!' >&2
    help
    exit 1
  fi

  if [ "$SIZE_OF_ARRAY" -lt 2 ]; then
    echo 'ERROR: You have to define at least two files for longitudinal processing!' >&2
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

  # modifiy defaults if needed
  cp ${defaults} ${defaults_tmp}

  if [ -n "$bids_folder" ]; then
    echo "cat.extopts.bids_folder = '${bids_folder}';" >> ${defaults_tmp}
    echo "cat.extopts.bids_yes = 1;" >> ${defaults_tmp}
  fi
  
  if [ "$bids" -eq 1 ]; then
    echo "cat.extopts.bids_yes = 1;" >> ${defaults_tmp}
  fi

}

########################################################
# run cat12 long pipeline
########################################################

run_cat12 ()
{
  pwd=$PWD
  
  # we have to go into toolbox folder to find matlab files
  cd $cwd
  
  if [ ! -n "${LOGDIR}" ]; then
    LOGDIR=$(dirname "${ARRAY[0]}")
  fi
  
  # we have to add current path if cat_batch_cat.sh was called from relative path
  if [ -d ${pwd}/${spm12} ]; then
    spm12=${pwd}/${spm12}
  fi

  export MATLABPATH=$spm12

  SIZE_OF_ARRAY="${#ARRAY[@]}"
  
  TMP=/tmp/cat_$$
  i=0
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
    
    if [ ! -n "${ARG_LIST}" ]; then
      ARG_LIST="$FILE"
    else
      ARG_LIST="${ARG_LIST} $FILE"
    fi
    ((i++))
  done
  
  echo ${ARG_LIST} >> ${TMP}

  time=`date "+%Y%b%d_%H%M"`
  vbmlog=${LOGDIR}/cat_${HOSTNAME}_${time}.log
  echo Check $vbmlog for logging information
  echo

  COMMAND="cat_batch_long('${TMP}','${output_surface}','${long_model}','${defaults_tmp}','${export_dartel}')"
  echo Running ${ARG_LIST}
  echo > $vbmlog
  echo ---------------------------------- >> $vbmlog
  date >> $vbmlog
  echo ---------------------------------- >> $vbmlog
  echo >> $vbmlog
  echo $0 $ARG_LIST >> $vbmlog
  echo >> $vbmlog

  if [ "$fg" -eq 0 ]; then
    nohup ${matlab} "$nojvm" -nodisplay -nosplash -r "$COMMAND" >> $vbmlog 2>&1 &
  else
    nohup ${matlab} "$nojvm" -nodisplay -nosplash -r "$COMMAND" >> $vbmlog 2>&1
  fi
  
  exit 0
}

########################################################
# check if matlab exist
########################################################

check_matlab ()
{
  found=`which ${matlab} 2>/dev/null`
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
cat <<__EOM__

USAGE:
  cat_batch_long.sh filenames|filepattern [-d default_file] [-m matlabcommand] 
                      [-log logdir] [-ns] [-large] [-model longmodel] [-e] [-nj] 
  
  -m <FILE>   | --matlab  <FILE>        matlab command (default $matlab)
  -d <FILE>   | --defaults <FILE>       optional default file (default ${cat12_dir}/cat_defaults.m)
  -log <FILE> | --logdir                directory for log-file (default $LOGDIR)
  -fg         | --fg                    do not run matlab process in background
  -ns         | --no-surf               disable surface and thickness estimation
  -e          | --export-dartel         export affine registered segmentations for Dartel
  -large      | --large                 use longitudinal model for detecting large changes (i.e. ageing or development)
                                        This option is only thought for compatibility with older scripts. Do not use that option together with the model flag. 
  -nj         | --nojvm                 supress call of jvm using the -nojvm flag
  -model      | --model                 longitudinal model:
                                          0 - detect large changes with brain/head growth (i.e. developmental effects)
                                          1 - detect small changes (i.e. due to plasticity)
                                          2 - detect large changes (i.e. ageing or development)
                                          3 - save results for both models 1 and 2
  -b          | --bids                  use default BIDS path (i.e. '../derivatives/CAT12.x_rxxxx')
  -bf <STRING>| --bids_folder <STRING>  define BIDS path

  Processing is only supported for one subject.
  Optionally you can set the matlab command with the "-m" option. As default no display
  is used (via the -nodisplay option in matlab). However sometimes the batch file needs
  a graphical output and the display should be enabled with the option "-d".

PURPOSE:
  Command line call of longitudinal segmentation pipeline

EXAMPLE
  cat_batch_long.sh all_files*.nii -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab
  This command will process all given files in the longitudinal pipeline. As matlab command 
  /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab will be used.
  
INPUT:
  filenames

OUTPUT:
  ${LOGDIR}/spm_${HOSTNAME}_${time}.log for log information

USED FUNCTIONS:
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
