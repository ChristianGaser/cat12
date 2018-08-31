 #! /bin/bash

# $Id$

########################################################
# global parameters
########################################################
version='check_pipeline.sh $Id$'

release=""
spm12_dir=""
spm12_tmp=/tmp/spm12$$
calc_tmp=/tmp/calc$$
file_list=""
bg_flag=" -fg -p 1"
bg=0
postprocess_only=0
volumes_only=0
scp_target="dbm.neuro.uni-jena.de:/Applications/xampp/htdocs/check_pipeline/"

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
        --bg* | -b*)
            exit_if_empty "$optname" "$optarg"
            bg_flag=" -p "$optarg
            bg=1
            shift
            ;;
        --post* | -p*)
            exit_if_empty "$optname" "$optarg"
            postprocess_only=$optarg
            shift
            ;;
        --vol* | -v*)
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
  
  if [ -z "$spm12_dir" ]; then
    echo "SPM12 directory is undefined!"
  fi

  SIZE_OF_ARRAY="${#ARRAY[@]}"
  if [ "$SIZE_OF_ARRAY" -eq 0 ]; then
      echo 'ERROR: No files given!' >&2
      help
      exit 1
  fi

  mkdir $calc_tmp

  old_dir=$PWD
  cd $calc_tmp
  
  i=0
  while [ "$i" -lt "$SIZE_OF_ARRAY" ]; do
    if [ ! -f "${ARRAY[$i]}" ]; then
      if [ ! -L "${ARRAY[$i]}" ]; then
        if curl --output /dev/null --silent --head --fail "${ARRAY[$i]}"
        then
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
  
  if [ -z "$release" ]; then
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
  
}

########################################################
# run run_pipeline
########################################################

run_pipeline ()
{
  
  # set ROI output and surface output
  if [ $volumes_only -eq 0 ]; then
    echo "cat.output.surface = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  fi
  echo "cat.output.ROI = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  echo "cat.extopts.ignoreErrors = 1;" >> ${spm12_tmp}/toolbox/cat12/cat_defaults.m
  
  # run cat12 in foreground with all files in tmp folder
  if [ -f "${spm12_tmp}/toolbox/cat12/cat_batch_cat.sh" ]; then
    ${spm12_tmp}/toolbox/cat12/cat_batch_cat.sh ${bg_flag} ${calc_tmp}/*.[in][mi][gi] 
  else
    ${spm12_tmp}/toolbox/cat12/cat_batch_vbm.sh ${bg_flag} ${calc_tmp}/*.[in][mi][gi] 
  fi
  
}

########################################################
# run postprocess
########################################################

postprocess ()
{

  # if postprocess_only > 0 we assume that this is the pid
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
    if [ ! -d check_r${release} ]; then
      echo Please check release number. Directory check_r${release} was not found.
      exit 1
    fi
    calc_tmp=check_r${release}
  fi
  

  # check whether xml files are found
  tmp=`ls ${calc_tmp}/report/cat_*xml 2>/dev/null`
  if [ ! -z "$tmp" ]; then
    for i in ${calc_tmp}/report/cat_*xml; do
      revision_cat=`grep revision_cat ${i}| cut -f2 -d">"|cut -f1 -d"<"`
      if [ -z "$revision_cat" ]; then
        revision_cat=`grep version_cat ${i}| cut -f2 -d">"|cut -f1 -d"<"`
      fi
      label=${calc_tmp}/label/catROI_`basename $i| sed -e 's/cat_//g'`
      subj=`basename $i | sed -e 's/\.xml//g' -e 's/cat_//g'`

      echo Finalize $subj with revision $revision_cat
      
      # get current csv files from dbm server
      scp -q -P 2222 ${scp_target}/${subj}*csv .

      # grep for Vgm and update csv file
      # check first for keyword neuromorphometrics and print next 200 lines
      Vgm=`grep -A200 "<neuromorphometrics" $label |grep Vgm | sed -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
      if [ ! -z "$Vgm" ]; then
        # add entry to csv file and sort and only keep unique lines
        echo "${revision_cat},${Vgm}" >> ${subj}_Vgm.csv
        cat ${subj}_Vgm.csv |sort -r|uniq > tmp$$
        mv tmp$$ ${subj}_Vgm.csv
      fi
      
      # grep for Vcsf and update csv file
      # check first for keyword neuromorphometrics and print next 200 lines
      Vcsf=`grep -A200 "<neuromorphometrics" $label |grep Vcsf | sed -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
      if [ ! -z "$Vcsf" ]; then
        # add entry to csv file and sort and only keep unique lines
        echo "${revision_cat},${Vcsf}" >> ${subj}_Vcsf.csv
        cat ${subj}_Vcsf.csv |sort -r|uniq > tmp$$
        mv tmp$$ ${subj}_Vcsf.csv
      fi

      # grep for ct and update csv file
      # check first for keyword neuromorphometrics and print next 200 lines
      ct=`grep -A200 "<neuromorphometrics" $label |grep ct | sed -e 's/;/,/g'|cut -f2 -d"["|cut -f1 -d"]"`
      if [ ! -z "$ct" ]; then
        # add entry to csv file and sort and only keep unique lines
        echo "${revision_cat},${ct}" >> ${subj}_ct.csv
        cat ${subj}_ct.csv |sort -r|uniq > tmp$$
        mv tmp$$ ${subj}_ct.csv
      fi

      # scp updated csv files to dbm server
      scp -q -P 2222 *.csv ${scp_target}

    done
        
    # rename tmp-folder and zip and scp them
    if [ ! -z "$revision_cat" ]; then
      if [ ! $postprocess_only -lt 0 ]; then
        if [ -d check_r${revision_cat} ]; then
          mv ${calc_tmp}/*/ check_r${revision_cat}/
        else
          mv ${calc_tmp} check_r${revision_cat}
        fi
        if [ ! -d spm12_r${revision_cat} ]; then
          mv ${spm12_tmp} spm12_r${revision_cat}
        fi
      fi
      
      # prepare renderview if tool is found and surface processing is enabled
      if [ ! -z `which CAT_View_Render_Matrix_ui` ] & [ $volumes_only -eq 0 ]; then
        CAT_View_Render_Matrix_ui check_r${revision_cat}/surf
      else
        echo "You need CAT_View_Render_Matrix_ui and image_matrix.sh for preparing render view."
      fi

      # delete original files, WM and normalized T1 images
      if [ ! $postprocess_only -lt 0 ]; then
        rm check_r${revision_cat}/*.[in][mi][gi]
        rm check_r${revision_cat}/mri/mwp2*
        rm check_r${revision_cat}/mri/wm*

        zip -q check_r${revision_cat}.zip -r check_r${revision_cat}
        scp -q -P 2222 check_r${revision_cat}.zip
      fi
      scp -q -P 2222 check_r*png $scp_target
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
  check_pipeline.sh -s spm12_folder [-p process_id] [-r cat12_zip_file] [-b number_of_processes] [-v] [-f file_list | filenames]
  
   -r   zip-file of CAT12 release that will be used for checking. If no zip-file is given, then the current CAT12 release will be used.
   -s   folder of spm12 installation
   -f   list with file names for checking
   -p   post-process given pid
   -b   run check_pipeline.sh in the background
   -v   run check_pipeline.sh without surface processing

   All given files will be processed using either the current CAT12 version or the defined CAT12 release with the "-r" flag.
   For the latter case the zip-file can be defined as local file or as url-address. During processing temporary folder are 
   saved in /tmp and modulated gray matter images and surfaces are calculated for each file. Finally, the post-processing will 
   be applied. During post-processing the temporary folders are renamed according to the found release number (check_rxxxx).
   Finally, the csv-files, the zipped check-folders, and the render views are transfered to ${scp_target}.
   If you run check_pipeline.sh in the background the post-processing has to be called manually using the "-b" flag if the
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
   check_pipeline.sh -s ~/spm/spm12 -r http://www.neuro.uni-jena.de/cat12/cat12_r1318.zip -f /Volumes/UltraMax/check_pipeline/check_pipeline_files.txt
   Run check_pipeline.sh with given file list and use CAT12 release 1318 from dbm-server and the defined SPM12 folder.
   
   check_pipeline.sh -s ~/spm/spm12 -bg 8 file1 file2 file3 file4 file5
   check_pipeline.sh -p pid
   Run check_pipeline.sh in the background with 8 processes with given 5 files and use current CAT12 release and the defined SPM12 folder.
   Finally, post-processing can be called with the given pid if processing is finished.

   check_pipeline.sh -p -xxxx
   Run post-processing again for release xxxx. This command is helpful if you have to re-calculated single data set because of a previous crash.
   Please note the "-" before the release number to indicate that this is not a pid but a release number.

USED FUNCTIONS:
   CAT12 toolbox
   SPM12
   CAT_View_Render_Matrix_ui
   image_matrix.sh

This script was written by Christian Gaser (christian.gaser@uni-jena.de).
This is ${version}.

__EOM__
}

########################################################
# call main program
########################################################

main ${1+"$@"}

