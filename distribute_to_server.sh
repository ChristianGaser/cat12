#! /bin/sh

########################################################
# global parameters
########################################################
version='distribute_to_server.sh $Id: distribute_to_server.sh 212 2013-08-12 08:57:32Z gaser $'

COMMAND=""
SERVER=localhost
PATTERN=""
DIR=""

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}
  
  check_files
  distribute

  exit 0
}

########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg files_to_calculate all_files
  count=0
  files_to_calculate=0
  all_files=0
  while [ $# -gt 0 ]
  do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    case "$1" in
        --command* | -c*)
            exit_if_empty "$optname" "$optarg"
            COMMAND=$optarg
            shift
            ;;
        --server* | -s*)
            exit_if_empty "$optname" "$optarg"
            SERVER=$optarg
            shift
            ;;
        --pattern* | -p*)
            exit_if_empty "$optname" "$optarg"
            PATTERN=$optarg
            shift
            ;;
        --dir* | -d*)
            exit_if_empty "$optname" "$optarg"
            DIR=$optarg
            shift

  if [ -z "$PATTERN" ]; then
    echo Pattern have to be defined.
    exit 0
  fi

list=`find $DIR -name "*.[in][mi][gi]" \! -name "*wrp[0-3]*.nii"  \! -name "*wp[0-3]*.nii" \! -name "wm*.nii"   \! -name "wrm*.nii"  \! -name "bf*.nii"  \! -name "p[0-3]*.nii"  \! -name "iy_*.nii"  \! -name "y_*.nii"  \! -name "rp[0-3]*.nii"  \! -name "._*.nii"`

for i in ${list} ; do
  # change extension to .nii and remove leading "./"
  name=`echo $i|sed -e 's/.img/.nii/' -e 's/\.\///g'`
  bname=`basename $name|cut -f1 -d'.'`
  dname=`dirname $name` 
  for j in ${dname}/${PATTERN}${bname}.nii ; do
    if [ ! -f "$j" ]; then
      ARRAY[$count]=$i
      ((count++))
    fi
  done
done

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

  if [ "$count" == "0" ]; then
    echo All files are already processed.
    exit 0
  fi
  
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
# run distribute
########################################################

distribute ()
{

  NUMBER_OF_SERVERS=0
  for k in ${SERVER}; do
    ((NUMBER_OF_SERVERS++))
  done

  SIZE_OF_ARRAY="${#ARRAY[@]}"
  BLOCK=$((10000* $SIZE_OF_ARRAY / $NUMBER_OF_SERVERS ))

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
    
  i=0
  for x in ${SERVER};
  do
    if [ ! "${ARG_LIST[$i]}" == "" ]; then
      j=$(($i+1))
      echo job ${j}/"$NUMBER_OF_SERVERS":
      echo $COMMAND ${ARG_LIST[$i]}
      if [ "$x" == "localhost" ]; then
        $COMMAND ${ARG_LIST[$i]}
      else
        bash -c "ssh ${x} $COMMAND ${ARG_LIST[$i]}"
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
  distribute_to_server.sh [-s server] -c command_to_distribute_to_server.sh filename|filepattern
  
   -c   command that should be parallelized
   -s   server list

   Only one filename or pattern is allowed. This can be either a single file or a pattern
   with wildcards to process multiple files. Optionally you can set number of processes,
   that are automatically set to the number of processors as default.

PURPOSE:
   distribute_to_server.sh a job or command

OUTPUT:

EXAMPLE
   distribute_to_server.sh -c "niismooth -v -fwhm 8" sTRIO*.nii
   smoothing with fwhm of 8mm for all files sTRIO*.nii. Use verbose mode to see diagnostic output.
   
   distribute_to_server.sh -s "141.35.68.68 141.35.68.72 141.35.68.73 141.35.68.74 141.35.68.75" -c "/Volumes/UltraMax/spm12b/toolbox/vbm12/cg_vbm_batch.sh -p 8 -d /Volumes/UltraMax/cg_vbm_defaults_p0123.m -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab" /Volumes/UltraMax/SVE.LPBA40.testdata/S*.nii
   VBM12 batch for all files in /Volumes/UltraMax/SVE.LPBA40.testdata/S*.nii with 8 parallel jobs and optional default file 

   distribute_to_server.sh -s "141.35.68.68 141.35.68.73 141.35.68.74 141.35.68.75" -c "/Volumes/UltraMax/spm12b/toolbox/vbm12/cg_vbm_batch.sh -p 8 -w -d /Volumes/UltraMax/cg_vbm_defaults_m0wrp12.m -m /Volumes/UltraMax/MATLAB_R2010b.app/bin/matlab" -d /Volumes/UltraMax/SVE.LPBA40.testdata -p m0wrp1
   VBM12 batch with 8 parallel jobs and optional default file using "Write alreayd estimated segmentations" as option. Only those files in /Volumes/UltraMax/SVE.LPBA40.testdata/ are processed where no prepended m0wrp1 pattern can be found. All other files are skipped.

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}

