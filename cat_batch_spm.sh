#! /bin/bash
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

cwd=$(dirname "$0")
matlab=matlab     # you can use other matlab versions by changing the matlab parameter
display=0         # use nodisplay option for matlab or not
LOGDIR=$PWD
spm12=$(dirname "$cwd")
spm12=$(dirname "$spm12")

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  run_batch

  exit 0
}


########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg

  if [ $# -lt 1 ]; then
    help
    exit 1
  fi

  while [ $# -gt 0 ]; do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    case "$1" in
      --matlab* | -m*)
        exit_if_empty "$optname" "$optarg"
        matlab=$optarg
        shift
        ;;
      --display* | -d*)
        display=1
        ;;
      --nojvm | -nj*)
        exit_if_empty "$optname" "$optarg"
        nojvm=" -nojvm "
        ;;
      --logdir* | -l*)
        exit_if_empty "$optname" "$optarg"
        LOGDIR=$optarg
        if [ ! -d $LOGDIR ] 
        then
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
        file="$1"
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
# run batch
########################################################

run_batch ()
{
  pwd=$PWD
  
  # we have to go into toolbox folder to find matlab files
  cd $cwd

  if [ ! -n "${LOGDIR}" ]; then
    LOGDIR=$(dirname "${ARRAY[0]}")
  fi

  # add current folder to matlabfile if file was not found
  if [ ! -f $file ]; then
    file=${pwd}/$file
  fi

  if [ ! -f $file ]; then
    echo File $file does not exist.
    exit 0
  fi

  dname=$(dirname "$file")
  file=$(basename "$file")
  
  if [ ! `echo "$file" | cut -f2 -d'.'` == "m" ]; then
    echo File "$file" is not a matlab script.
    exit 0
  fi

  # we have to add current path if cat_batch_cat.sh was called from relative path
  if [ -d ${pwd}/${spm12} ]; then
    spm12=${pwd}/${spm12}
  fi

  export MATLABPATH=$spm12:$dname
  
  time=`date "+%Y%b%d_%H%M"`
    spmlog=${LOGDIR}/spm_${HOSTNAME}_${time}.log
  echo Check $spmlog for logging information
  echo
    
  file=`echo $file| sed -e 's/\.m//g'`

  X="cat_batch_spm('${file}')"
  echo Running $file
  echo > $spmlog
  echo ---------------------------------- >> $spmlog
  date >> $spmlog
  echo ---------------------------------- >> $spmlog
  echo >> $spmlog
  echo $0 $file >> $spmlog
  echo >> $spmlog
  if [ $display == 0 ]; then
    nohup ${matlab} -nodisplay "$nojvm" -nosplash -r $X >> $spmlog 2>&1 &
  else
    nohup ${matlab} -nosplash -r $X >> $spmlog 2>&1 &
  fi
  exit 0
}

########################################################
# check if matlab exist
########################################################

check_matlab ()
{
  found=`which ${matlab} 2>/dev/null`
  if [ ! -n "$found" ];then
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
   cat_batch_spm.sh batchfile.m [-d] [-m matlabcommand]
   
  -m  <FILE>  | --matlab <FILE>         matlab command (matlab version) (default $matlab)
  -d  <FILE>  | --display <FILE>        use display option in matlab in case that batch file needs graphical output
  -l  <FILE>  | --logdir                directory for log-file (default $LOGDIR)
  -nj         | --nojvm                 supress call of jvm using the -nojvm flag

   Only one batch filename is allowed. Optionally you can set the matlab command 
   with the "-m" option. As default no display is used (via the -nodisplay option 
   in matlab). However sometimes the batch file needs a graphical output and the 
   display should be enabled with the option "-d".

PURPOSE:
   Command line call of SPM12 batch files

EXAMPLE
   cat_batch_spm.sh test_batch.m -m /usr/local/bin/matlab7
   This command will process the batch file test_batch.m. As matlab command 
   /usr/local/bin/matlab7 will be used.
   
INPUT:
   batch file saved as matlab-script or mat-file

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
