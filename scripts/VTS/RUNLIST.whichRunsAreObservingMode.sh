#!/bin/bash
# from a run list, prints the list of runs that were taken in a specific observing mode (observing, obsLowHV, obsFilter,...)
#

ISPIPEFILE=`readlink /dev/fd/0` # check to see if input is from terminal, or from a pipe
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} ]] ; then # its a terminal (not a pipe)
	if ! [ $# -eq 2 ] ; then # the human didn't add any arguments, and we must tell them so
		echo
		echo "From a runlist or pipe, prints the run numbers that have been taken in a particular observing mode."
		echo "Usage: "
		echo " $ `basename $0` <mode> <file of runs>" 
		echo "or"
		echo " $ cat myrunlist.dat | `basename $0` <mode>" ; echo
		echo "   <mode> = observing for regular runs"
		echo "            obsLowHV for runs taken with reduced HV"
		echo "            obsFilter for runs taken with UV filters"
		echo "            \"obs[LF]*\" for runs taken with reduced HV or UV filters"
		echo "Full list of run types: observing, chargeInjection, laser, pedestal, biasCurve, roadlaser, flasher, highlow, test, other, obsFilter, obsLowHV, rasterScan"
		echo
		exit
	fi
fi

# list of run_ids to read in
RUNFILE=$2
if [ ! -e $RUNFILE ] ; then
	echoerr "File $RUNFILE could not be found in $PWD , sorry."
	exit 1
fi
RUNLIST=`cat $RUNFILE`
#echo "RUNLIST:$RUNLIST"

MODE="$1"

# get database url from parameter file
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' "$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter" | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
    
if [ ! -n "$MYSQLDB" ] ; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit
#else
#    echo "MYSQLDB: $MYSQLDB"
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
while read -r RUNID RUNMODE ; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		
	        #echo "$RUNID :: $RUNWOBANGLE"
	        if [ "$RUNMODE" ] && [[ "$RUNMODE" == $MODE ]] ; then
	            echo "$RUNID" 
	        fi
	fi
# This is where the MYSQL command is executed, with the list of requested runs
# You have to do it this way, because using a pipe | calls the command in a
# subshell, and that prevents variables from being saved within the 'while' loop
# http://stackoverflow.com/questions/14585045/is-it-possible-to-avoid-pipes-when-reading-from-mysql-in-bash
done < <($MYSQL -e "USE VERITAS ; SELECT run_id, run_type FROM tblRun_Info WHERE $SUB")

exit
