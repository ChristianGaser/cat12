#! /bin/sh

########################################################
# global parameters
########################################################
# $Id$

spm8=~/spm/spm8b  # this parameter has to be set to your spm8 directory
matlab=matlab     # you can use other matlab versions by changing the matlab parameter
display=0         # use nodisplay option for matlab or not

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
		--d* | -d*)
			display=1
			;;
		-h | --help | -v | --version | -V)
			help
			exit 1
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

  if [ -z "$val" ]
  then
	echo 'ERROR: No argument given with \"$desc\" command line argument!' >&2
	exit 1
  fi
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

	if [ $# -eq 2 ]; then
		if [ ! -d $spm8 ]; then
			spm8=${pwd}/$spm8_
	    fi
		if [ ! -d $spm8 ]; then
			echo Directory $spm8 does not exist.
			exit 0
		fi
	fi

	if [ ! -f $file ]; then
		echo File $file does not exist.
		exit 0
	fi

	dname=`dirname $file`
	file=`basename $file`
	
	if [ ! `echo $file | cut -f2 -d'.'` == "m" ]; then
		echo File $file is not a matlab script.
		exit 0
	fi

	export MATLABPATH=$MATLABPATH:${spm8}/toolbox/vbm8:$spm8:$dname
	
	time=`date "+%Y%b%d_%H%M"`
	vbmlog=${pwd}/spm8_${time}.log
	echo Check $vbmlog for logging information
	echo
		
	file=`echo $file| sed -e 's/\.m//g'`

	X="cg_spm8_batch('${file}')"
	echo Running $file
	echo > $vbmlog
	echo ---------------------------------- >> $vbmlog
	date >> $vbmlog
	echo ---------------------------------- >> $vbmlog
	echo >> $vbmlog
	echo $0 $file >> $vbmlog
	echo >> $vbmlog
	if [ $display == 0 ]; then
		${matlab} -nodisplay -nojvm -nosplash -r $X -logfile $vbmlog &
	else
		${matlab} -nojvm -nosplash -r $X -logfile $vbmlog &
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
   cg_spm8_batch.sh batchfile.m [-d] [-s spm8-path] [-m matlabcommand]
   
   -d   use display option in matlab in case that batch file needs graphical output
   -m   matlab command
   -s   spm8 directory

   Only one batch filename is allowed. Optionally you can set the spm8 
   directory with the "-s" option and the matlab command with the "-m" option.
   As default no display is used (via the -nodisplay option in matlab). However
   sometimes the batch file needs a graphical output and the display should
   be enabled with the option "-d".

PURPOSE:
   Command line call of SPM8 batch files

EXAMPLE
   cg_spm8_batch.sh test_batch.m -s ~/spm/spm8 -m /usr/local/bin/matlab7
   This command will process the batch file test_batch.m. The spm8 directory
   is set to ~/spm/spm8. As matlab command /usr/local/bin/matlab7 will be used.
   
INPUT:
   batch file saved as matlab-script or mat-file

OUTPUT:
   spm8_log_$time.txt for log information

USED FUNCTIONS:
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
