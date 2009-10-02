#! /bin/sh

########################################################
# global parameters
########################################################
# $Id$

spm8=~/spm/spm8	# this parameter has to be set to your spm8 directory
matlab=matlab	# you can use other matlab versions by changing the matlab parameter
writeonly=0
CPUINFO=/proc/cpuinfo
ARCH=`uname`

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_matlab
  check_files
  get_no_of_cpus
  run_vbm

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
		--m* | -m*)
			exit_if_empty "$optname" "$optarg"
			matlab=$optarg
			shift
			;;
		--s* | -s*)
			exit_if_empty "$optname" "$optarg"
			spm8=$optarg
			shift
			;;
		--p* | -p*)
			exit_if_empty "$optname" "$optarg"
			NUMBER_OF_JOBS=$optarg
			shift
			;;
		--w* | -w*)
			writeonly=1
			;;
		-h | --help | -v | --version | -V)
			help
			exit 1
			;;
		-*)
			echo "`basename $0`: ERROR: Unrecognized option \"$1\"" >&2
			;;
		*)
            ARRAY[$count]="$1"
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
      echo ERROR: File ${ARRAY[$i]} not found
  	  help
      exit 1
    fi
    ((i++))
  done

}

########################################################
# get # of cpus
########################################################

get_no_of_cpus () {

  if [ -z "$NUMBER_OF_JOBS" ];
  then
    if [ "$ARCH" == "Linux" ]
    then
      NUMBER_OF_JOBS=`grep ^processor $CPUINFO | wc -l`

    elif [ "$ARCH" == "Darwin" ]
    then
      NUMBER_OF_JOBS=`sysctl -a hw | grep -w logicalcpu | awk '{ print $2 }'`

    elif [ "$ARCH" == "FreeBSD" ]
    then
      NUMBER_OF_JOBS=`sysctl hw.ncpu | awk '{ print $2 }'`

    else
      NUMBER_OF_JOBS=`grep ^processor $CPUINFO | wc -l`

    fi
    echo "Found $NUMBER_OF_JOBS processors."

    if [ -z "$NUMBER_OF_JOBS" ]
    then
        echo "$FUNCNAME ERROR - number of CPUs not obtained. Use -p to define number of processes."
        exit 1
    fi
  fi
}

########################################################
# run vbm tool
########################################################

run_vbm ()
{
	cwd=`dirname $0`
	pwd=$PWD
	
    SIZE_OF_ARRAY="${#ARRAY[@]}"
    BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_JOBS ))

    # split files
    i=0
    while [ "$i" -lt "$SIZE_OF_ARRAY" ]
    do
        count=$((10000* $i / $BLOCK ))
        
        # check wheter absolute or relative names were given
        if [ ! -f ${ARRAY[$i]} ];  then
            FILE=${pwd}/${ARRAY[$i]}
        else
            FILE=${ARRAY[$i]}
        fi
        if [ -z "${ARG_LIST[$count]}" ]; then
            ARG_LIST[$count]="$FILE"
        else
            ARG_LIST[$count]="${ARG_LIST[$count]} $FILE"
        fi
        ((i++))
    done

	# we have to go into toolbox folder to find matlab files
	cd $cwd

	if [ $# -eq 2 ]; then
		if [ ! -d $spm8 ]; then
			spm8=${pwd}/$spm8
	    fi
		if [ ! -d $spm8 ]; then
			echo Directory $spm8 does not exist.
			exit 0
		fi
	fi

	export MATLABPATH=$MATLABPATH:${spm8}/toolbox/vbm8:$spm8
	
	time=`date "+%Y%b%d_%H%M"`
	vbmlog=${pwd}/vbm8_${time}.log
	echo Check $vbmlog for logging information
	echo
	
    echo >> $vbmlog
    i=0
    while [ "$i" -lt "$NUMBER_OF_JOBS" ]
    do
        if [ ! "${ARG_LIST[$i]}" == "" ]; then
            j=$(($i+1))
            COMMAND="cg_vbm8_batch('${ARG_LIST[$i]}',${writeonly})"
            echo Calculate ${ARG_LIST[$i]}
            echo ---------------------------------- >> $vbmlog
            date >> $vbmlog
            echo ---------------------------------- >> $vbmlog
            echo >> $vbmlog
            echo $0 $file >> $vbmlog
            echo >> $vbmlog
            nohup ${matlab} -nodisplay -nojvm -nosplash -r $COMMAND >> $vbmlog 2>&1 &
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
   cg_vbm8_batch.sh filename|filepattern [-s spm8-path] [-m matlabcommand] [-w]
   
   -m   matlab command
   -s   spm8 directory
   -p   number of parallel jobs (=number of processors)
   -w
   Only one filename or pattern is allowed. This can be either a single file or a pattern
   with wildcards to process multiple files. For the latter case you have to (single) quote
   the pattern. Optionally you can set the spm8 directory with the "-s" option and the
   matlab command with the "-m" option and force to write already estimated segmentations with
   the "-w" option.

PURPOSE:
   Command line call of VBM8 segmentation

EXAMPLE
   cg_vbm8_batch.sh spm/spm8/canonical/single_subj_T1.nii -s ~/spm/spm8
   This command will process only the single file single_subj_T1.nii. The spm8 directory
   is set to ~/spm/spm8.
   
   cg_vbm8_batch.sh 'spm/spm8/canonical/*152*.nii'
   Using the quotes and wildcards all files containing the term "152" will
   be processed. In this case these are the files avg152PD.nii, avg152T1.nii,
   and avg152T2.nii.

   cg_vbm8_batch.sh 'spm/spm8/canonical/*152*.nii' -m /usr/local/bin/matlab7
   Using the quotes and wildcards all files containing the term "152" will
   be processed. In this case these are the files avg152PD.nii, avg152T1.nii,
   and avg152T2.nii. As matlab command /usr/local/bin/matlab7 will be used.

INPUT:
   analyze/nifti files

OUTPUT:
   segmented images according to settings in cg_vbm8_defaults.m
   vbm8_log_$time.txt for log information

USED FUNCTIONS:
   cg_vbm8_batch.m
   VBM8 toolbox
   SPM8

SETTINGS
   spm8 path: $spm8
   matlab command: $matlab
   
This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
