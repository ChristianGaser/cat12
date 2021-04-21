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

matlab=matlab   # you can use other matlab versions by changing the matlab parameter
cwd=$(dirname "$0")
cat12_dir=$cwd
spm12=$(dirname "$cwd")
spm12=$(dirname "$spm12")
LOGDIR=$PWD
output_surface=1
large_changes=0
fg=0

########################################################
# run main
########################################################

main ()
{
 parse_args ${1+"$@"}
 check_matlab
 check_files
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
      --large* | -large*)
        exit_if_empty "$optname" "$optarg"
        large_changes=1
        ;;
      --no-surf* | -ns*)
        exit_if_empty "$optname" "$optarg"
        output_surface=0
        ;;
      --nojvm | -nj*)
        exit_if_empty "$optname" "$optarg"
        nojvm=" -nojvm "
        ;;
      --fg* | -fg*)
        exit_if_empty "$optname" "$optarg"
        fg=1
        ;;
      --logdir* | -log*)
        exit_if_empty "$optname" "$optarg"
        LOGDIR=$optarg
        if [ ! -d $LOGDIR ]; then
           mkdir -p $LOGDIR
        fi
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
# run cat12 long pipeline
########################################################

run_cat12 ()
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

  # we have to go into toolbox folder to find matlab files
  cd $cwd
  
  if [ ! -n "${LOGDIR}" ]; then
    LOGDIR=$(dirname "${ARRAY[0]}")
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

  COMMAND="cat_batch_long('${TMP}','${output_surface}','${large_changes}','${defaults}')"
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
                      [-log logdir] [-ns] [-l] [-nj] 
  
  -m <FILE>   | --matlab  <FILE> matlab command (default $matlab)
  -d <FILE>   | --default <FILE> optional default file (default ${cat12_dir}/cat_defaults.m)
  -log <FILE> | --logdir         directory for log-file (default $LOGDIR)
  -fg         | --fg             do not run matlab process in background
  -ns         | --no-surf        disable surface and thickness estimation
  -large      | --large          use longitudinal model for detecting large changes (e.g. ageing or development)
  -nj         | --nojvm          supress call of jvm using the -nojvm flag

  Processing is omly supported for one subject.
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
