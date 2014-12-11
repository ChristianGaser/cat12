#! /bin/sh

########################################################
# global parameters
########################################################
# $Id: cg_vbm_longitudinal_batch.sh 653 2014-11-21 13:13:47Z gaser $

matlab=matlab     # you can use other matlab versions by changing the matlab parameter
display=0         # use nodisplay option for matlab or not
LOGDIR=`PWD`

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
  run_batch

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

  while [ $# -gt 0 ]
  do
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

  if [ -z "$val" ]
  then
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
  if [ "$SIZE_OF_ARRAY" -eq 0 ]
  then
      echo 'ERROR: No files given!' >&2
      help
      exit 1
  fi

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]
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
# run batch
########################################################

run_batch ()
{
    cwd=`dirname $0`
    pwd=$PWD
    
    # we have to go into toolbox folder to find matlab files
    cd $cwd
    
    spm12=`dirname $cwd`
    spm12=`dirname $spm12`

    if [ "${LOGDIR}" == "" ]; then
        LOGDIR=`dirname ${ARRAY[0]}`
    fi
    
    export MATLABPATH=$spm12

    SIZE_OF_ARRAY="${#ARRAY[@]}"
	
    TMP=/tmp/vbm_$$
    i=0
    ARG_LIST=""
    while [ "$i" -lt "$SIZE_OF_ARRAY" ]
    do
        # check wether absolute or relative names were given
        if [ ! -f ${ARRAY[$i]} ];  then
            FILE=${pwd}/${ARRAY[$i]}
        else
            FILE=${ARRAY[$i]}
        fi
        if [ -z "${ARG_LIST}" ]; then
            ARG_LIST="$FILE"
        else
            ARG_LIST="${ARG_LIST} $FILE"
        fi
        ((i++))
    done
    
    echo ${ARG_LIST} >> ${TMP}

	time=`date "+%Y%b%d_%H%M"`
    vbmlog=${LOGDIR}/vbm_${HOSTNAME}_${time}.log
	echo Check $vbmlog for logging information
	echo
		
	X="cg_vbm_longitudinal_batch('${TMP}')"
	echo Running $file
	echo > $vbmlog
	echo ---------------------------------- >> $vbmlog
	date >> $vbmlog
	echo ---------------------------------- >> $vbmlog
	echo >> $vbmlog
	echo $0 $file >> $vbmlog
	echo >> $vbmlog

	if [ $display == 0 ]; then
		nohup ${matlab} -nodisplay -nosplash -r "$X" >> $vbmlog 2>&1 &
	else
		nohup ${matlab} -nosplash -r $X >> $vbmlog 2>&1 &
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
   cg_vbm_longitudinal_batch.sh file1.nii file2.nii ... filex.nii [-d] [-m matlabcommand]
   
   -d   use display option in matlab in case that batch file needs graphical output
   -m   matlab command

   Only one batch filename is allowed. Optionally you can set the matlab command 
   with the "-m" option. As default no display is used (via the -nodisplay option 
   in matlab). However sometimes the batch file needs a graphical output and the 
   display should be enabled with the option "-d".

PURPOSE:
   Command line call of SPM12 batch files

EXAMPLE
   cg_vbm_longitudinal_batch.sh all_files*.nii -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab
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
