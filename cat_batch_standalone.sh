#! /bin/bash

# $Id: cat_batch_standalone.sh 1361 2018-08-31 08:44:57Z gaser $

########################################################
# run main
########################################################

main ()
{
  cwd=`dirname $0`
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
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
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

  if [ -z "$val" ]; then
    echo 'ERROR: No argument given with \"$desc\" command line argument!' >&2
    exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  
  # check for MCR parameter
  if [ ! -n "$MCRROOT" ]; then
    echo "No MCR folder given."
    help
    exit 1  
  fi
  
  # check for MCR folder
  if [ ! -d "$MCRROOT" ]; then
    echo "No MCR folder found."
    help
    exit 1  
  fi

  # check for SPM parameter
  if [ ! -n "$SPMROOT" ]; then
    echo "No SPM folder given."
    help
    exit 1  
  fi
  
  # check for SPM folder
  if [ ! -f "$SPMROOT/run_spm12.sh" ]; then
    echo "File $SPMROOT/run_spm12.sh not found found."
    help
    exit 1  
  fi

  # check for batch file
  if [ ! -n "$BATCHFILE" ]; then
    echo "No batch file given."
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
  if [ "$count" -eq "0" ]; then
    
    eval "\"${SPMROOT}/run_spm12.sh\"" $MCRROOT "batch" $BATCHFILE
    exit 0
  fi
  
  # create temporary batch file
	TMP=/tmp/cat_$$.m

  # copy everything except rows with UNDEFINED to temp file	
	grep -v "UNDEFINED" $BATCHFILE > $TMP

	# extract parameter name of data structure
	data=`grep "UNDEFINED" $BATCHFILE | grep -e "\.data" -e "\.cdata" -e "\.mov" | cut -f1 -d'='`

  # surface data need an additional curly bracket
  if grep -q "\.data_surf" $BATCHFILE ; then
    echo "$data = {{" >> $TMP
  else
    echo "$data = {" >> $TMP
	fi
	
  i=0
  ARG_LIST=""
  while [ "$i" -lt "$count" ]; do

    # check wether absolute or relative names were given
    if [ ! -f ${ARRAY[$i]} ];  then
      if [ -f ${pwd}/${ARRAY[$i]} ]; then
        FILE=${pwd}/${ARRAY[$i]}
      fi
    else
      FILE=${ARRAY[$i]}
    fi

    # add file list
		echo "'${FILE}'" >> $TMP
		          
    ((i++))
  done

  # surface data need an additional curly bracket
  if grep -q "\.data_surf" $BATCHFILE ; then
	  echo "     }};" >> $TMP
	else
	  echo "     };" >> $TMP
	fi
	eval "\"${SPMROOT}/run_spm12.sh\"" $MCRROOT "batch" $TMP
#	rm $TMP

  exit 0
}


########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
   cat_batch_standalone.sh filename(s) [-s spm_standalone_folder] [-m mcr_folder] [-b batch_file]
   
   -s   SPM12 folder of standalone version (can be also defined by SPMROOT)
   -m   Matlab Compiler Runtime (MCR) folder (can be also defined by MCRROOT)
   -b   batch file
   
PURPOSE:
   Command line call of CAT12 segmentation for SPM12 standalone installation

EXAMPLE
   cat_batch_standalone.sh -b ${cwd}/cat_batch_standalone.m sTRIO0001.nii
   Preprocess the single file sTRIO0001.nii using the default CAT12 
   preprocessing batch. 

   cat_batch_standalone.sh -b ${cwd}/cat_batch_simple_standalone.m sTRIO0001.nii
   Process the single file sTRIO0001.nii using the simple processing batch. 

   cat_batch_standalone.sh  -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_batch_standalone.m sTRIO0001
   Preprocess the single file sTRIO0001.nii using the simple processing batch
   and the SPM12 standalone version in ~/spm/standalone and Matlab Compiler Runtime in
   /Applications/MATLAB/MATLAB_Runtime/v93/
   
   cat_batch_standalone.sh -b ${cwd}/cat_batch_resample_standalone.m lh.thickness.sTRIO0001
   Resample and smooth the single thickness file sTRIO0001.nii. Only the left surface file has to be defined.

   cat_parallelize.sh -p 8 -l /tmp -c "cat_batch_standalone.sh  -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_batch_standalone.m" sTRIO*.nii
   Parallelize CAT12 preprocessing by splitting all sTRIO*.nii files into 8 jobs 
   (processes) and save log file in /tmp folder. 

   The parameter SPMROOT and MCRROOT have to be defined if no additional flags -s -m are used.      

INPUT:
   analyze or nifti files

OUTPUT:
   segmented images and optionally surfaces according to settings in cat_batch_standalone.m or
   cat_batch_simple_standalone.m.

USED FUNCTIONS:
   cat_batch_standalone.m
   cat_batch_simple_standalone.m
   cat_batch_resample_standalone.m
   cat_parallelize.sh
   SPM12 standalone version (compiled)
   CAT12 toolbox (compiled within SPM12)
   MATLAB Compiler Runtime R2017b (Version 9.3)

This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
