#! /bin/bash
# Call SPM12 batch jobs from shell
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
mpath="'"$(dirname "$1 ")"'"
ivar=

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
    optarg2="`echo $3 | sed 's,^[^=]*=,,'`"
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
      -p | --path)
        exit_if_empty "$optname" "$optarg"
        mpath="$mpath,'$optarg'"
        shift
        ;;
      -f | --files)
      # read all files until nothing remains or the next command is coming and
      # pack everything into a matlab cellstr
      	  exit_if_empty2 "$optname" "$optarg" "$optarg2"
      	  ivar="$ivar,'$optarg',{'$optarg2'"
      	  runon=1
      	  shift
      	  shift
		  while [ $# -gt 0 ] && [ $runon -gt 0 ]; do
			case "$2" in
				-*)
					runon=0
					;;
				*) 
					ivar="$ivar;'$2'"
					shift
					;;
			esac
		  done
		  ivar=$ivar"}"
		  ;;
      -i | --var | --ivar)
        exit_if_empty2 "$optname" "$optarg" "$optarg2"
        ivar="$ivar,'$optarg',$3"
        shift
        shift
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
exit_if_empty2 ()
{
  local desc var val

  desc="$1"
  shift
  var="$2"
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

  # prepare matlab code with additional path 
  X="addpath($mpath); cat_batch_spm('${file}'${ivar})"
  
  echo SPM command:
  echo "  "$X; 
  echo

  echo Running $file
  echo > $spmlog
  echo ---------------------------------- >> $spmlog
  date >> $spmlog
  echo ---------------------------------- >> $spmlog
  echo >> $spmlog
  echo $0 $file >> $spmlog
  echo >> $spmlog
  if [ $display == 0 ]; then
    #nohup 
    ${matlab} -nodisplay "$nojvm" -nosplash -r "$X" #>> $spmlog 2>&1 &
  else
    nohup ${matlab} -nosplash -r $X #>> $spmlog 2>&1 &
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
  -d  <FILE>  | --display <FILE>        use display option in matlab in case that 
                                        batch file needs graphical output
  -l  <FILE>  | --logdir                directory for log-file (default $LOGDIR)
  -nj         | --nojvm                 supress call of jvm using the -nojvm flag
  -p  <PATH>  | --mpath <PATH>          add directory to matlab path (the path of
                                        the called file is added automaticly)
  -f  <varname> <FILES> | -files <varname> 
                                        Create a cellstr by a set of files,e.g.,
                                         -f data path1/file1.ext path2/file2.ext 
                                        will create a matlab variable 
                                          data = {'path1/file1.ext';'path1/file1.ext'};                                       
  -i  <varname> "<variable>"  |  --ivar <varname> "<variable>" 
                                        Create a variable "varname" to be used in  
                                        the matlabbatch to set up variables, e.g.,
                                        to specify directories, files, or options.  
                                        Examples: 
                                         (1) To create a simple matrix:
                                             -i mat   "[0 8 1 5]"                                       
                                         (2) To create a structure with subfields:
                                             -i opt   "struct('flag1',1,'verb',0,'files',{'file.ext})"
                                         (3) To create a cell array with filenames:
                                             -i files "{'path1/fname1.ext','path2/fname2.ext'}"         

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
   
   cat_batch_spm.sh /Users/cat/matlabbatches/batch_prepana_smoothandmore.m 
     -m /Applications/MATLAB_R2020a.app/bin/matlab 
     -f pdirs /Volumes/drive/MRData/ADNI/derivatives/CAT12.8.1/sub-ADNI002S0955
              /Volumes/drive/MRData/ADNI/derivatives/CAT12.8.1/sub-ADNI002S0954
     -i vdata "{'mwp1'}" -i smoothing "[4 8]"  
   
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
