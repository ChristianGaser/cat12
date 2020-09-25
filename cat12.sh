#! /bin/bash

# $Id$

########################################################
# global parameters
########################################################
version='cat12.sh $Id$'

cat12_dir=`dirname $0`
defaults=""
matlab=matlab
no_mwp=0
rp=0
no_surf=0
bg_flag=" -fg -p 1"
tpm=
add_to_defaults=
defaults_tmp=/tmp/defaults$$.m

# add full path if necessary
if [ -d ${PWD}/${cat12_dir} ]; then
  cat12_dir=${PWD}/`echo $cat12_dir | sed -e 's/\.\///g'`
fi

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
        --bg* | -b*)
            exit_if_empty "$optname" "$optarg"
            bg_flag=" -p "$optarg
            bg=1
            shift
            ;;
        --tpm* | -t*)
            exit_if_empty "$optname" "$optarg"
            tpm=$optarg
            shift
            ;;
        --add* | -a*)
            exit_if_empty "$optname" "$optarg"
            add_to_defaults="$optarg"
            shift
            ;;
        --no-mwp* | -nm*)
            exit_if_empty "$optname" "$optarg"
            no_mwp=1
            ;;
        --no-surf* | -ns*)
            exit_if_empty "$optname" "$optarg"
            no_surf=1
            ;;
        --rp* | -r*)
            exit_if_empty "$optname" "$optarg"
            rp=1
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
  if [ $no_surf -eq 1 ] && [ $no_mwp -eq 1 ] && [ $rp -eq 0 ]; then
    echo 'WARNING: You have deselected all outputs! Only the p0-image will be saved.' >&2
  fi
  
  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]; then
    echo 'ERROR: No files given!' >&2
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

  # argument empty?
  if [ ! "${defaults}" == "" ]; then
      # check wether absolute or relative names were given
      if [ ! -f ${defaults} -a -f ${pwd}/${defaults} ]; then
          defaults=${pwd}/${defaults}
      fi

      # check whether defaults file exist
      if [ ! -f ${defaults} ];  then
          echo Default file $defaults not found.
          exit
      fi
  else
    defaults=${cat12_dir}/cat_defaults.m
  fi

  cp ${defaults} ${defaults_tmp}
  
  if [ $no_surf -eq 1 ]; then
    echo "cat.output.surface = 0;" >> ${defaults_tmp}
  else
    echo "cat.output.surface = 1;" >> ${defaults_tmp}
  fi

  if [ $no_mwp -eq 1 ]; then
    echo "cat.output.GM.mod      = 0;" >> ${defaults_tmp}
    echo "cat.output.WM.mod      = 0;" >> ${defaults_tmp}
    echo "cat.output.ROI         = 0;" >> ${defaults_tmp}
    echo "cat.output.bias.warped = 0;" >> ${defaults_tmp}

  else
    echo "cat.output.GM.mod = 1;" >> ${defaults_tmp}
    echo "cat.output.WM.mod = 1;" >> ${defaults_tmp}
    echo "cat.output.ROI    = 1;" >> ${defaults_tmp}
  fi

  if [ $rp -eq 1 ]; then
    echo "cat.output.GM.dartel = 2;" >> ${defaults_tmp}
    echo "cat.output.WM.dartel = 2;" >> ${defaults_tmp}
  else
    echo "cat.output.GM.dartel = 0;" >> ${defaults_tmp}
    echo "cat.output.WM.dartel = 0;" >> ${defaults_tmp}
  fi
  
  if [ ! -z "$tpm" ]; then
    echo "cat.opts.tpm = ${tpm};" >> ${defaults_tmp}
  fi

  if [ ! -z "$add_to_defaults" ]; then
    echo "${add_to_defaults}" >> ${defaults_tmp}
  fi
}

########################################################
# run cat12
########################################################

run_cat12 ()
{
  ${cat12_dir}/cat_batch_cat.sh -m ${matlab} -d ${defaults_tmp} ${bg_flag} ${ARRAY[@]}
}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
  cat12.sh filename|filepattern [-m matlab_command] [-d default_file] [-b number_of_processes] [-a add_to_defaults] [-ns] [-nm] [-rp] [-tpm TPM-file]
  
   --matlab   | -m      matlab command (matlab version)
   --bg       | -b      run cat12.sh in the background with given number of processes
   --defaults | -d      optional default file
   --add      | -a      add option to default file
   --no-surf  | -ns     skip surface and thickness estimatation
   --no-mwp   | -nm     skip estimating modulated and warped segmentations and ROI measures
   --rp       | -r      additionally estimate affine registered segmentations

  Only one filename or pattern is allowed. This can be either a single file or a pattern
  with wildcards to process multiple files. 

PURPOSE:
   Wrapper to call cat_batch_cat.sh with predefined set of defaults that can be changed

INPUT:
   nifti files

OUTPUT:
   segmented images and surfaces according to defined settings

EXAMPLE:
   cat12.sh -ns -nm -rp spm/spm12/canonical/*152*.nii -m /usr/local/bin/matlab
   Using wildcards all files containing the term "152" will be processed. In this case these 
   are the files avg152PD.nii, avg152T1.nii, and avg152T2.nii. Only the affine registered segmentations are saved.
   As matlab-command /usr/local/bin/matlab will be used.
   
   cat12.sh spm/spm12/canonical/single_subj_T1.nii
   This command will process only the single file single_subj_T1.nii with the defaults in cat_defaults.m

   cat12.sh spm/spm12/canonical/single_subj_T1.nii -a "cat.extopts.WMHC = 3;"
   This command will process only the single file single_subj_T1.nii with the defaults bin cat_defaults.m and
   the additional option for handling WMHs as separate class. 

   cat12.sh spm/spm12/canonical/single_subj_T1.nii -tpm ${cat12_dir}/templates_volumes/TPM_Age11.5.nii
   This command will process only the single file single_subj_T1.nii with the defaults bin cat_defaults.m
   and the children template that is provided with cat12.

USED FUNCTIONS:
   CAT12 toolbox
   SPM12

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}


