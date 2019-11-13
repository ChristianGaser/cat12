#! /bin/sh

version='cat_parallelize.sh $Id$'

########################################################
# global parameters
########################################################

COMMAND=""
TEST=""
CPUINFO=/proc/cpuinfo
ARCH=`uname`
LOGDIR=$PWD
time=`date "+%Y%b%d_%H%M"`

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  check_files
  get_no_of_cpus
  parallelize

  exit 0
}

########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg
  count=0
  while [ $# -gt 0 ]
  do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="$2"
    case "$1" in
        --command* | -c*)
            exit_if_empty "$optname" "$optarg"
            COMMAND=$optarg
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
            shift
            ;;
        --test* | -t*)
            TEST=1
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
    echo ERROR: "No argument given with \"$desc\" command line argument!" >&2
    exit 1
  fi
}

########################################################
# check files
########################################################

check_files ()
{
  if [ -z "$COMMAND" ];
  then
    echo "$FUNCNAME ERROR - no command defined."
      help
    exit 1
  fi
  
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
    if [ ! -f "${ARRAY[$i]}" ] && [ ! -d "${ARRAY[$i]}" ]; then
      echo ERROR: File or directory ${ARRAY[$i]} not found
      help
      exit 1
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
# run parallelize
########################################################

parallelize ()
{
  SIZE_OF_ARRAY="${#ARRAY[@]}"
  BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_JOBS ))

  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]
  do
    count=$((10000* $i / $BLOCK ))
    if [ -z "${ARG_LIST[$count]}" ]; then
      ARG_LIST[$count]="${ARRAY[$i]}"
    else
      ARG_LIST[$count]="${ARG_LIST[$count]} ${ARRAY[$i]}"
    fi
    ((i++))
  done

  time=`date "+%Y%b%d_%H%M"`
  log=${LOGDIR}/parallelize_${HOSTNAME}_${time}.log
  if [ "${TEST}" == "" ]; then
    echo Check $log for logging information
    echo > $log
    echo
  fi
    
  i=0
  while [ "$i" -lt "$NUMBER_OF_JOBS" ]
  do
    if [ ! "${ARG_LIST[$i]}" == "" ]; then
      j=$(($i+1))
      echo job ${j}/"$NUMBER_OF_JOBS":
      echo $COMMAND ${ARG_LIST[$i]}
      if [ "${TEST}" == "" ]; then
          echo job ${j}/"$NUMBER_OF_JOBS": $COMMAND ${ARG_LIST[$i]} >> $log
          nohup bash -c "for k in ${ARG_LIST[$i]}; do $COMMAND \$k; done" >> $log 2>&1 &
      fi
    fi
    ((i++))
  done

}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
  cat_parallelize.sh [-p number_of_processes] [-l log_folder] [-t] -c command_to_parallelize filename|filepattern
  
   -p   number of parallel jobs (=number of processors)
   -c   command that should be parallelized
   -t   do not call command, but print files to be processed
   -l   directory where log-file will be saved

   Only one filename or pattern is allowed. This can be either a single file or a pattern
   with wildcards to process multiple files. Optionally you can set number of processes,
   that are automatically set to the number of processors as default.

PURPOSE:
   Parallelize a job or command

OUTPUT:
   parallelize_${HOSTNAME}_${time}.log with current data and time in name as log-file

EXAMPLE
   cat_parallelize.sh -c "niismooth -v -fwhm 8" sTRIO*.nii
     Parallelize smoothing with fwhm of 8mm for all files sTRIO*.nii. Use
     verbose mode to see diagnostic output.
   
   cat_parallelize.sh -c gunzip *.zip
     Parallelize unzipping of all zip-files in current folder. 
     
   cat_parallelize.sh -p 8 -l /tmp -c "cat_standalone.sh  -s ~/spm/standalone/ -m /Applications/MATLAB/MATLAB_Runtime/v93/ -b cat_standalone_segment.m" sTRIO*.nii
     Parallelize CAT12 preprocessing by splitting all sTRIO*.nii files into 8 jobs 
     (processes) and save log-file in /tmp folder. 

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}

