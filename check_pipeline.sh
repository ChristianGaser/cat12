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
version='check_pipeline.sh $Id$'

spm12_tmp=/tmp/spm12_$$
calc_tmp=/tmp/calc$$
proc_dir=$PWD
bg_flag=" -fg -p 1"
bg_flag_long=" -fg"
bg=0
postprocess_only=0
volumes_only=0
scp_target="dbm.neuro.uni-jena.de:/volume1/web/check_pipeline/"

########################################################
# run main
########################################################

main ()
{
  parse_args ${1+"$@"}

  if [ $postprocess_only -eq 0 ]; then
    copy_files
    get_release
    run_pipeline
  fi

  # don't run postprocess if check_pipeline is running in the background
  if [  $bg -eq 0 ]; then
    postprocess
  fi
    
  exit 0
}

########################################################
# check arguments and files
########################################################

parse_args ()
{
  local optname optarg files_to_calculate all_files
  count=0
  count_long=0
  files_to_calculate=0
  all_files=0
  while [ $# -gt 0 ]
  do
    optname="`echo $1 | sed 's,=.*,,'`"
    optarg="`echo $2 | sed 's,^[^=]*=,,'`"
    case "$1" in
      --release* | -r*)
          exit_if_empty "$optname" "$optarg"
          release=$optarg
          shift
          ;;
      --spm* | -s*)
          exit_if_empty "$optname" "$optarg"
          spm12_dir=$optarg
          shift
          ;;
      --dir* | -d*)
          exit_if_empty "$optname" "$optarg"
          proc_dir=$optarg
          shift
          ;;
      --file* | -f*)
          exit_if_empty "$optname" "$optarg"
          listfile=$optarg
          shift
          list=$(< $listfile);
          for F in $list; do
            ARRAY[$count]=$F
            ((count++))
          done
          ;;
      --long* | -l*)
          exit_if_empty "$optname" "$optarg"
          listfile=$optarg
          shift
          list_long=$(< $listfile);
          for F in $list_long; do
            ARRAY_LONG[$count_long]=$F
            ((count_long++))
          done
          ;;
      --bg* | -b*)
          exit_if_empty "$optname" "$optarg"
          bg_flag=" -p "$optarg
          bg_flag_long=""
          bg=1
          shift
          ;;
      --post* | -p*)
          exit_if_empty "$optname" "$optarg"
          postprocess_only=$optarg
          shift
          ;;
      --no-surf* | -ns*)
          exit_if_empty "$optname" "$optarg"
          volumes_only=1
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
# check and copy files
########################################################

copy_files ()
{
  
  if [ ! -n "$spm12_dir" ]; then
    echo "SPM12 directory is undefined!"
  fi

  SIZE_OF_ARRAY=${#ARRAY[@]}
  SIZE_OF_ARRAY_LONG=${#ARRAY_LONG[@]}

  old_dir=$PWD

  if [ "$SIZE_OF_ARRAY" -eq 0 ] && [ "$SIZE_OF_ARRAY_LONG" -eq 0 ]; then
    echo 'ERROR: No files given!' >&2
    help
    exit 1
  fi

  if [ "$SIZE_OF_ARRAY_LONG" -gt 2 ]; then
    echo 'ERROR: Only two longitudinal scans are currently supported!' >&2
    help
    exit 1
  fi

  if [ "$SIZE_OF_ARRAY" -gt 0 ]; then

    mkdir -p $calc_tmp
    cd $calc_tmp
    
    i=0
    while [ "$i" -lt "$SIZE_OF_ARRAY" ]; do
      if [ ! -f "${ARRAY[$i]}" ]; then
        if [ ! -L "${ARRAY[$i]}" ]; then
          if curl --output /dev/null --silent --head --fail "${ARRAY[$i]}"; then
            curl -O ${ARRAY[$i]}
          else
            echo File or url ${ARRAY[$i]} does not exist
            exit
          fi
        fi
      fi
      cp ${ARRAY[$i]} ${calc_tmp}/
      ((i++))
    done
    
    cd $old_dir
  fi

  if [ "$SIZE_OF_ARRAY_LONG" -gt 0 ]; then
    mkdir -p ${calc_tmp}/long
  
    cd ${calc_tmp}/long

    i=0
    while [ "$i" -lt "$SIZE_OF_ARRAY_LONG" ]; do
      if [ ! -f "${ARRAY_LONG[$i]}" ]; then
        if [ ! -L "${ARRAY_LONG[$i]}" ]; then
          if curl --output /dev/null --silent --head --fail "${ARRAY_LONG[$i]}"; then
            curl -O ${ARRAY_LONG[$i]}
          else
            echo File or url ${ARRAY_LONG[$i]} does not exist
            exit
          fi
        fi
      fi
      cp ${ARRAY_LONG[$i]} ${calc_tmp}/long/
      ((i++))
    done
    cd $old_dir
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

  if [ ! -n "$val" ]
  then
    echo ERROR: "No argument given with \"$desc\" command line argument!" >&2
    exit 1
  fi
}

########################################################
# run get_release
########################################################

get_release ()
{
  if [ ! -d "${spm12_dir}" ]; then
    echo Directory $spm12_dir does not exist!
    exit 1
  fi
  
  # copy current spm12 installation to tmp folder
  cp -r $spm12_dir $spm12_tmp
  
  if [ ! -n "$release" ]; then
    echo "Use current release."
  else
    # remove old cat12 folder
    rm -r ${spm12_tmp}/toolbox/cat12
        
    # check whether it's a file
    if [ -f "$release" ]; then
      unzip -q $release -d ${spm12_tmp}/toolbox/    
    else # or a web-address
      cat12_tmp=/tmp/cat12$$.zip
      curl -o $cat12_tmp $release 
      unzip -q $cat12_tmp -d ${spm12_tmp}/toolbox/
      rm $cat12_tmp
    fi
  fi

  # allow execution of mexmaci64 file on MAC
  if [ "$ARCH" == "Darwin" ]; then
    echo "Please login as admin to allow execution of mex files on MAC OS"
    sudo xattr -r -d com.apple.quarantine $spm12_tmp
    sudo find $spm12_tmp -name \*.mexmaci64 -exec spctl --add {} \;  
  fi
}

########################################################
# run run_pipeline
########################################################

run_pipeline ()
{
  
  echo PID is $$
  
  # set ROI output and surface output
  if [ $volumes_only -eq 0 ]; then
    echo "cat.output.surface = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  fi
  echo "cat.output.ROI = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  echo "cat.extopts.ignoreErrors = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  
  # run cat12 in foreground with all files in tmp folder
  if [ "$SIZE_OF_ARRAY" -gt 0 ]; then
    if [ -f "${spm12_tmp}/toolbox/cat12/cat_batch_cat.sh" ]; then
      ${spm12_tmp}/toolbox/cat12/cat_batch_cat.sh -l ${proc_dir} ${bg_flag} ${calc_tmp}/*.[in][mi][gi] 
    else
      ${spm12_tmp}/toolbox/cat12/cat_batch_vbm.sh -l ${proc_dir} ${bg_flag} ${calc_tmp}/*.[in][mi][gi] 
    fi
  fi
  
  if [ "$SIZE_OF_ARRAY_LONG" -gt 0 ]; then
    large=`grep "\-large" ${spm12_tmp}/toolbox/cat12/cat_batch_long.sh`
    # call "-large" option only if available for that release
    if [ -n "$large" ]; then
      ${spm12_tmp}/toolbox/cat12/cat_batch_long.sh -large ${bg_flag_long} ${calc_tmp}/long/*.[in][mi][gi]
    else
      ${spm12_tmp}/toolbox/cat12/cat_batch_long.sh ${bg_flag_long} ${calc_tmp}/long/*.[in][mi][gi] 
    fi
  fi
  
}

########################################################
# run postprocess
########################################################

postprocess ()
{

  # if postprocess_only > 0 we assume that this is a pid
  if [ $postprocess_only -gt 0 ]; then
    pid="$postprocess_only"
    if [ -d /tmp/calc${pid} ]; then
      spm12_tmp=/tmp/spm12${pid}
      calc_tmp=/tmp/calc${pid}
    else
      echo Please check process ID. Directory /tmp/calc${pid} was not found.
      exit 1
    fi
  fi

  # if postprocess_only < 0 we assume that this is the release number
  if [ $postprocess_only -lt 0 ]; then
    release=`echo $postprocess_only| cut -f2 -d'-'`
    if [ ! -d ${proc_dir}/check_r${release} ]; then
      echo Please check release number. Directory ${proc_dir}/check_r${release} was not found.
      exit 1
    fi
    calc_tmp=${proc_dir}/check_r${release}
  fi
  

  # check whether xml files are found
  for folder in "" "/long"; do
    calc_tmp2=${calc_tmp}/${folder}
    
    # remove all average files from long mode
    rm ${calc_tmp}/long/*/*avg* 2>/dev/null
    
    tmp=`ls ${calc_tmp2}/report/cat_*xml 2>/dev/null`
    if [ -n "$tmp" ]; then
      for i in ${calc_tmp2}/report/cat_*xml; do
        revision_cat=`grep revision_cat ${i}| cut -f2 -d">"|cut -f1 -d"<"`
        if [ ! -n "$revision_cat" ]; then
          revision_cat=`grep version_cat ${i}| cut -f2 -d">"|cut -f1 -d"<"`
        fi
        label=${calc_tmp2}/label/catROI_`basename $i| sed -e 's/cat_//g'`
        labels=${calc_tmp2}/label/catROIs_`basename $i| sed -e 's/cat_//g'`
        report=${calc_tmp2}/report/cat_`basename $i| sed -e 's/cat_//g'`
        subj=`basename $i | sed -e 's/\.xml//g' -e 's/cat_//g'`
  
        echo Finalize $subj with revision $revision_cat
        
        # get current csv files from dbm server
        scp -q -P $PORT ${scp_target}/${subj}*csv .
  
        # grep for vol_TIV and vol_abs_CGW and update csv file
        # check first for keywords and print next 5 lines
        vol_TIV=`grep -A5 "<vol_abs_CGW" $report |grep vol_TIV |cut -f2 -d">"|cut -f1 -d"<"`
        vol_CGW=`grep "<vol_abs_CGW" $report | sed -e 's/\ /,/g'|cut -f2 -d"["|cut -f1 -d"]"|cut -f1-4 -d','`
        if [ -n "$vol_TIV" ] && [ -n "$vol_CGW" ]; then
          # add entry to csv file and sort and only keep unique lines
          echo "${revision_cat},${vol_TIV},${vol_CGW}" >> ${subj}_vol.csv
          cat ${subj}_vol.csv |sort -r|uniq > tmp$$
          mv tmp$$ ${proc_dir}/${subj}_vol.csv
        fi
  
        # grep for Vgm and update csv file
        # check first for keyword neuromorphometrics and print next 200 lines
        Vgm=`grep -A200 "<neuromorphometrics" $label |grep Vgm | sed -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
        if [ -n "$Vgm" ]; then
          # add entry to csv file and sort and only keep unique lines
          echo "${revision_cat},${Vgm}" >> ${subj}_Vgm.csv
          cat ${subj}_Vgm.csv |sort -r|uniq > tmp$$
          mv tmp$$ ${proc_dir}/${subj}_Vgm.csv
        fi
        
        # grep for Vcsf and update csv file
        # check first for keyword neuromorphometrics and print next 200 lines
        Vcsf=`grep -A200 "<neuromorphometrics" $label |grep Vcsf | sed -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
        if [ -n "$Vcsf" ]; then
          # add entry to csv file and sort and only keep unique lines
          echo "${revision_cat},${Vcsf}" >> ${subj}_Vcsf.csv
          cat ${subj}_Vcsf.csv |sort -r|uniq > tmp$$
          mv tmp$$ ${proc_dir}/${subj}_Vcsf.csv
        fi
  
        # grep for thickness and update csv file
        # check first for keyword neuromorphometrics and print next 200 lines
        thickness=`grep -A200 "<aparc_DK40" $labels |grep thickness | sed -e 's/\ /,/g' -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
        thickness=`grep -A200 "<aparc_DK40" $labels |grep thickness | sed -e 's/;/,/g' -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
        if [ -n "$thickness" ]; then
          # add entry to csv file and sort and only keep unique lines
          echo "${revision_cat},${thickness}" >> ${subj}_thickness.csv
          cat ${subj}_thickness.csv |sort -r|uniq > tmp$$
          mv tmp$$ ${proc_dir}/${subj}_thickness.csv
        fi
  
        # scp updated csv files to dbm server
        scp -q -P $PORT *.csv ${scp_target}
  
      done
    fi
  done
     
  # rename tmp-folder and zip and scp them
  if [ -n "$revision_cat" ]; then
    if [ ! $postprocess_only -lt 0 ]; then
      if [ -d ${proc_dir}/check_r${revision_cat} ]; then
        mv ${calc_tmp}/*/ ${proc_dir}/check_r${revision_cat}/
      else
        mv ${calc_tmp} ${proc_dir}/check_r${revision_cat}
      fi
    fi
    
    # prepare renderview if tool is found and surface processing is enabled
    if ([ ! -z `which render_surf.sh` ] && [ ! -z `which CAT_View` ]) && [ $volumes_only -eq 0 ]; then
      mkdir -p ${proc_dir}/check_r${revision_cat}/surf
      ln -s ${proc_dir}/check_r${revision_cat}/long/surf/* ${proc_dir}/check_r${revision_cat}/surf/ >/dev/null 2>&1
      render_surf.sh -range 0 6 ${proc_dir}/check_r${revision_cat}/surf
      mv check_r${revision_cat}*.png ${proc_dir}/ >/dev/null 2>&1
      scp -q -P $PORT ${proc_dir}/check_r${revision_cat}*.png $scp_target
    else
      echo "You need render_surf.sh and image_matrix.sh for preparing render view."
    fi

    # delete original files, WM and normalized T1 images
    if [ ! $postprocess_only -lt 0 ]; then
      rm ${proc_dir}/check_r${revision_cat}/*.[in][mi][gi] 2>/dev/null
      rm ${proc_dir}/check_r${revision_cat}/mri/*mwp2* 2>/dev/null
      rm ${proc_dir}/check_r${revision_cat}/mri/wm* 2>/dev/null

      if [ -d ${proc_dir}/check_r${revision_cat}/long ]; then
        rm ${proc_dir}/check_r${revision_cat}/long/*.[in][mi][gi] 2>/dev/null
        rm ${proc_dir}/check_r${revision_cat}/long/mri/mwp2* 2>/dev/null
        rm ${proc_dir}/check_r${revision_cat}/long/mri/wm* 2>/dev/null
      fi
      
      zip -q ${proc_dir}/check_r${revision_cat}.zip -r ${proc_dir}/check_r${revision_cat}
      scp -q -P $PORT ${proc_dir}/check_r${revision_cat}.zip $scp_target
    fi
  fi
  
}

########################################################
# help
########################################################

help ()
{
cat <<__EOM__

USAGE:
  check_pipeline.sh -s spm12_folder [-p process_id] [-r cat12_zip_file] [-d proc_folder] [-b number_of_processes] [-ns] [-f file_list | filenames]
  
   --release <FILE|URL> | -r <FILE|URL> zip-file of CAT12 release that will be used for checking. If no zip-file is given, then the current CAT12 release will be used.
   --spm <DIR>          | -s <DIR>      folder of spm12 installation
   --dir <DIR>          | -d <DIR>      folder for writing processed files and results (default $proc_folder)
   --file <FILE>        | -f <FILE>     list with file names for checking
   --post <STRING>      | -p <STRING>   post-process given pid
   --bg                 | -b            run check_pipeline.sh in the background
   --no-surf            | -ns           run check_pipeline.sh without surface processing

   All given files will be processed using either the current CAT12 version or the defined CAT12 release with the "-r" flag.
   For the latter case the zip-file can be defined as local file or as url-address. During processing temporary folder are 
   saved in /tmp and modulated gray matter images and surfaces are calculated for each file. Finally, the post-processing will 
   be applied. During post-processing the temporary folders are renamed according to the found release number (check_rxxxx).
   Finally, the csv-files, the zipped check-folders, and the render views are transfered to ${scp_target}.
   If you run check_pipeline.sh in the background the post-processing has to be called manually using the "-p" flag if the
   processing for all data is finished.
    

PURPOSE:
   check_pipeline.sh a cat12 release

INPUT:
   nifti files

OUTPUT:
   The following files are saved in the current folder:
     csv-files with ROI volumes
     render views
     check_r* folders with release name

EXAMPLE:
   check_pipeline.sh -s ~/spm/spm12 -r http://www.neuro.uni-jena.de/cat12/cat12_r1318.zip file1 file2 file3 file4 file5
   Run check_pipeline.sh with given 5 files and use CAT12 release 1318 from dbm-server and the defined SPM12 folder.
   
   check_pipeline.sh -s ~/spm/spm12 -bg 8 -f /Volumes/UltraMax/check_pipeline/check_pipeline_files.txt -l /Volumes/UltraMax/check_pipeline/check_pipeline_files_long.txt
   check_pipeline.sh -p pid
   Run check_pipeline.sh in the background with 8 processes with given file list and use current CAT12 release and the defined SPM12 folder.
   Finally, post-processing can be called with the given pid if processing is finished.

   check_pipeline.sh -p -xxxx
   Run post-processing again for release xxxx. This command is helpful if you have to re-calculated single data set because of a previous crash.
   Please note the "-" before the release number to indicate that this is not a pid but a release number.

USED FUNCTIONS:
   CAT12 toolbox
   SPM12
   render_surf.sh
   image_matrix.sh

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}

