#!/bin/bash
# from a run list, prints the list of runs that are on disk.
# written by Nathan Kelley-Hoskins Aug 2013

NOTFLAG=false # flag for if the -n flag was used
HELPFLAG=false # if true, print help text and exit
#echo "INP:'`basename $0`' '$1' '$2' '$3'"

ISPIPEFILE=`readlink /dev/fd/0` # check to see if input is from terminal, or from a pipe
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} ]] ; then # its a terminal (not a pipe)
	if [ "$#" -eq "1" ] ; then # format is "exe <fname>"
		RUNFILE=$1
	elif [ "$#" -eq "2" ] ; then # format is "exe -flag <fname>"
		if [ "$1" = "-n" ] ; then
			NOTFLAG=true
			RUNFILE=$2
		else
			echo " Error: `basename $0` doesn't understand flag $1.  Only acceptable flag is -n"
			HELPFLAG=true
		fi
	else
		echo "needs at least one argument"
		HELPFLAG=true
	fi
else # it is a pipe
	if [ "$#" -eq "0" ] ; then # format is " cat runlist.dat | exe "
		RUNFILE=$1
	elif [ "$#" -eq "1" ] ; then # format is " cat runlist.dat | exe -flags"
		if [ "$1" = "-n" ] ; then
			NOTFLAG=true
			RUNFILE=$2
		else
			echo " Error: `basename $0` doesn't understand flag $1.  Only acceptable flag is -n"
			HELPFLAG=true
		fi
	else
		echo "needs at least one argument"
		HELPFLAG=true
	fi
fi

if $HELPFLAG ; then
	echo
	echo "Prints the run numbers that ARE stored on disk." ; echo
	echo " $ `basename $0` <file of runs>" ; echo
	echo "Or, prints the run numbers that are NOT stored on disk" ; echo
	echo " $ `basename $0` -n <file of runs>" ; echo
	echo " $ cat <file of runs> | `basename $0`" ; echo
	echo " $ cat <file of runs> | `basename $0` -n" ; echo
	exit
fi
	
#echo "NOTFLAG:$NOTFLAG"

# list of run_id's to read in
#RUNFILE=$1
if [ ! -e $RUNFILE ] ; then
	echo "File '$RUNFILE' could not be found, sorry."
	exit	
fi
RUNLISTTMP=`cat $RUNFILE`
RUNLIST=$(echo "$RUNLISTTMP" | grep -oP "^\d+$" )
if [ -z "$RUNLIST" ] ; then
  >&2 echo "Error, RUNLIST.whichRunsAreOnDisk.sh : input file/pipe contains no runs, exiting..." 
  exit 1
fi
#echo "RUNLIST:$RUNLIST"
#echo "Files not on disk:"
    
# find the veritas db url
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
if [ ! -n "$MYSQLDB" ] ; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit
fi 


# mysql login info
MYSQL="mysql -u readonly -h $MYSQLDB -A"

# generate list of runs to ask for ( run_id = RUNID[1] OR run_id = RUNID[2] etc)
COUNT=0
SUB=""
for ARUN in $RUNLIST ; do
	if (( $ARUN > 0 )); then
		if [[ "$COUNT" -eq 0 ]] ; then
			SUB="run_id = $ARUN"
		else 
			SUB="$SUB OR run_id = $ARUN"
		fi
		COUNT=$((COUNT+1))
	fi
done
#echo "SUB:"
#echo "$SUB"

# search through mysql result rows, where each row's elements
# are assigned to RUNID and RUNDATE
while read -r RUNID RUNDATE ; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		
		# decode the date tag
		read YY MM DD HH MI SE <<< ${RUNDATE//[-:]/ }
		#echo "  YEARMONTHDAY:$YY$MM$DD"
		
		# generate the filename
		TARGFILE="$VERITAS_DATA_DIR/data/d$YY$MM$DD/$RUNID.cvbf"
		
		# test to see if the file exists
		#echo "  Does file exist: $TARGFILE"
		if [ -e $TARGFILE ] ; then # file exists
			if ! $NOTFLAG ; then # $NOTFLAG is false, and we should print the runnumber
				echo "$RUNID"
			fi
		else # file does not exist
			if $NOTFLAG ; then # $NOTFLAG is true, and we should print the runnumber
				echo "$RUNID"
			fi
		fi
	fi
# This is where the MYSQL command is executed, with the list of requested runs
# You have to do it this way, because using a pipe | calls the command in a
# subshell, and that prevents variables from being saved within the 'while' loop
# http://stackoverflow.com/questions/14585045/is-it-possible-to-avoid-pipes-when-reading-from-mysql-in-bash
done < <($MYSQL -e "USE VERITAS ; SELECT run_id, data_start_time FROM tblRun_Info WHERE $SUB")

exit
