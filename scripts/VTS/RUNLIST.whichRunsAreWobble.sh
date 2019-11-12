#!/bin/bash
# from a run list, prints the list of runs that were taken in a specific atmosphere, summer(22) or winter(21)
# written by Nathan Kelley-Hoskins Sept 2013

CONORM="\e[0m"
CORED='\e[1;31m'

ISPIPEFILE=`readlink /dev/fd/0` # check to see if input is from terminal, or from a pipe
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} ]] ; then # its a terminal (not a pipe)
	if ! [ $# -eq 2 ] ; then # the human didn't add any arguments, and we must tell them so
		echo
		echo "From a runlist or pipe, prints the run numbers that are of a particular wobble or wobbles."
		echo " $ `basename $0` [nsew] <file of runs>" ; echo
		echo "Print list of only north wobble runs:"
		echo " $ `basename $0` n myrunlist.dat" ; echo
		echo "Print list of only south and east wobble runs:"
		echo " $ `basename $0` se myrunlist.dat" ; echo
		echo "Works with pipes : " 
		echo " $ cat myrunlist.dat | `basename $0` w" ; echo
		exit
	fi
fi

function echoerr(){ echo -e "${CORED}$@${CONORM}" 1>&2; } #for spitting out error text

# list of run_id's to read in
RUNFILE=$2
if [ ! -e $RUNFILE ] ; then
	echoerr "File $RUNFILE could not be found in $PWD , sorry."
	exit 1
fi
RUNLIST=`cat $RUNFILE`
#echo "RUNLIST:$RUNLIST"

NORTFLAG=false
SOUTFLAG=false
EASTFLAG=false
WESTFLAG=false

LOWARG=`echo "$1" | tr '[A-Z]' '[a-z]'` # make all uppercase letters in arg 1 lowercase, for easier handling
#echo "\$LOWARG: '$LOWARG'"
if [[ "$LOWARG" == *n* ]] ; then NORTFLAG=true ; fi
if [[ "$LOWARG" == *s* ]] ; then SOUTFLAG=true ; fi
if [[ "$LOWARG" == *e* ]] ; then EASTFLAG=true ; fi
if [[ "$LOWARG" == *w* ]] ; then WESTFLAG=true ; fi

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
while read -r RUNID RUNWOBANGLE ; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		
        #echo "$RUNID :: $RUNWOBANGLE"
        if $NORTFLAG ; then
            if [ "$RUNWOBANGLE" -eq "0"   ] ; then echo "$RUNID" ; fi
        fi
        if $EASTFLAG ; then
            if [ "$RUNWOBANGLE" -eq "90"  ] ; then echo "$RUNID" ; fi
        fi
        if $SOUTFLAG ; then
            if [ "$RUNWOBANGLE" -eq "180" ] ; then echo "$RUNID" ; fi
        fi
        if $WESTFLAG ; then
            if [ "$RUNWOBANGLE" -eq "270" ] ; then echo "$RUNID" ; fi
        fi
		
	fi
# This is where the MYSQL command is executed, with the list of requested runs
# You have to do it this way, because using a pipe | calls the command in a
# subshell, and that prevents variables from being saved within the 'while' loop
# http://stackoverflow.com/questions/14585045/is-it-possible-to-avoid-pipes-when-reading-from-mysql-in-bash
done < <($MYSQL -e "USE VERITAS ; SELECT run_id, offset_angle FROM tblRun_Info WHERE $SUB")

exit
