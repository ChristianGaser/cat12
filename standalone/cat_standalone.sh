#! /bin/bash

# $Id$

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
	grep -v "<UNDEFINED>" $BATCHFILE > $TMP
  
	# extract parameter name of data structure (1st occurance of "<UNDEFINED>")
	data=`grep -m 1 "<UNDEFINED>" $BATCHFILE | cut -f1 -d'='`

	# extract parameter name of optional argument(s) (additional occurances of "<UNDEFINED>")
  if [ -n "$ARG1" ]; then # ARG1 defined?
	  param1=`grep -m 2 "<UNDEFINED>" $BATCHFILE | tail -n 1 | cut -f1 -d'='`
	  # extract parameter name of optional argument (3rd occurance of "<UNDEFINED>")
		if [ -n "$ARG2" ]; then # ARG2 defined?
			param2=`grep -m 3 "<UNDEFINED>" $BATCHFILE | tail -n 1 | cut -f1 -d'='`
		fi
  fi
  
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
	
  if [ -n "$ARG1" ]; then # ARG1 defined?
    echo "$param1 = $ARG1 ;" >> $TMP
    if [ -n "$ARG2" ]; then # ARG2 defined?
      echo "$param2 = $ARG2 ;" >> $TMP
    fi
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
   cat_standalone.sh filename(s) [-s spm_standalone_folder] [-m mcr_folder] [-b batch_file] 
                                 [-a1 additional_argument1] [-a2 additional_argument2]
   
   -s   SPM12 folder of standalone version (can be also defined by SPMROOT)
   -m   Matlab Compiler Runtime (MCR) folder (can be also defined by MCRROOT)
   -b   batch file to execute
   -a1  1st additional argument 
   -a2  2nd additional argument 
   
   The first occurance of the parameter "<UNDEFINED>" in the batch file will be replaced by the
   list of input files. You can use the existing batch files in this folder or create your own batch 
   file with the SPM12 batch editor and leave the data field undefined. Please note that for creating
   your own batch file CAT12 has to be called in expert mode because the CAT12 standalone installation 
   will only run in expert mode to allow more options.
   See cat_standalone_segment.txt for an example. 
   
   You can also define an one or two optional arguments to change other parameters that are indicated by "<UNDEFINED>"
   in the batch file. Please take care of the order of the "<UNDEFINED>" fields in the batch file.
   If you use a computer cluster it is recommended to use the batch files to only process one data set 
   and use a job or queue tool to provide the single jobs to the cluster.
   
PURPOSE:
   Command line call of (CAT12) batch files for SPM12 standalone installation

EXAMPLE
   cat_standalone.sh -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_standalone_segment.txt sTRIO0001.nii
   Preprocess (segment) the single file sTRIO0001.nii using the default CAT12 preprocessing batch. 
   SPM12 standalone version is located in ~/spm/standalone and Matlab Compiler Runtime in
   /Applications/MATLAB/MATLAB_Runtime/v93/.

   cat_standalone.sh -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_standalone_simple.txt sTRIO0001.nii
   Process the single file sTRIO0001.nii using the simple processing batch. 

   cat_standalone.sh -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_standalone_resample.txt -a1 "12" -a2 "1" lh.thickness.sTRIO0001
   Resample and smooth the single thickness file lh.thickness.sTRIO0001 with 12mm and save the resampled mesh as 32k mesh from HCP. 
   Only the left surface file has to be defined.

   cat_standalone.sh -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_standalone_smooth.txt -a1 "[6 6 6]" -a2 "'s6'" sTRIO*nii
   Smooth the volume files sTRIO*nii with 6mm and prepend the string "s6" to smoothed files.

   cat_parallelize.sh -p 8 -l /tmp -c "cat_standalone.sh  -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b ${cwd}/cat_standalone_segment.txt" sTRIO*.nii
   Parallelize CAT12 preprocessing by splitting all sTRIO*.nii files into 8 jobs 
   (processes) and save log file in /tmp folder. 

   The parameters SPMROOT and MCRROOT have to be defined (exported) to skip the use of the flags -s -m.

INPUT:
   nifti files or surface data

OUTPUT:
   segmented images and optionally surfaces according to settings in cat_standalone_segment.txt or
   cat_standalone_simple.txt.

USED FUNCTIONS:
   cat_standalone_segment.txt
   cat_standalone_simple.txt
   cat_standalone_resample.txt
   cat_parallelize.sh
   SPM12 standalone version (compiled)
   CAT12 toolbox (compiled within SPM12 if installed)
   MATLAB Compiler Runtime R2017b (Version 9.3)

This script was written by Christian Gaser (christian.gaser@uni-jena.de).

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}
