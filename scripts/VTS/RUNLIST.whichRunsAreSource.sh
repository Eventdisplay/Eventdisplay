#!/bin/bash
# from a run list, prints the list of runs that target a specific source
# written by Nathan Kelley-Hoskins Sept 2013

# variables for coloring terminal output
CONORM="\e[0m"
CORED='\e[1;31m'

NOTFLAG=false # flag for if the -n flag was used
HELPFLAG=false # if true, print help text and exit

# check to see if input is from terminal, or from a pipe
ISPIPEFILE=`readlink /dev/fd/0`
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} ]] ; then # its a terminal (not a pipe)
	if [ "$#" -eq "1" ] ; then
		NOTFLAG=true
		SOURCE=""
		RUNFILE="$1"
	elif [ "$#" -eq "2" ] ; then # format is "exe src <fname>"
		SOURCE=$1
		RUNFILE=$2
	elif [ "$#" -eq "3" ] ; then # format is "exe -flag src <fname>"
		if [ "$1" = "-n" ] ; then
			NOTFLAG=true
			SOURCE="$2"
			RUNFILE="$3"
		else
			echo " Error: `basename $0` doesn't understand flag '$1'.  Only acceptable flag is -n"
			HELPFLAG=true
		fi
	else
		#echo "needs at least one argument"
		HELPFLAG=true
	fi
else # it is a pipe
	if [ "$#" -eq "0" ] ; then
		NOTFLAG=true
		SOURCE=""
		RUNFILE="$1"
	elif [ "$#" -eq "1" ] ; then # format is " cat runlist.dat | exe src"
		SOURCE=$1
		RUNFILE=$2
	elif [ "$#" -eq "2" ] ; then # format is " cat runlist.dat | exe -flag src"
		if [ "$1" = "-n" ] ; then
			NOTFLAG=true
			SOURCE=$2
			RUNFILE=$3
		else
			echo " Error: `basename $0` doesn't understand flag '$1'.  Only acceptable flag is -n"
			HELPFLAG=true
		fi
	else
		#echo "needs at least one argument"
		HELPFLAG=true
	fi
fi

if $HELPFLAG ; then
    echo
    echo "Prints the source names for a list of run numbers,"
	echo "  OR, only print runs that do or do not target a particular source." ; echo
    echo "$ `basename $0` [-n] [source name] <file of runs>" ; echo
	echo "Required arguments:"
	echo "  <file of runs> : The runlist file, containing the list of runs, 1 runnumber per line"
	echo "                   Can also take input from a piped runlist" ; echo
	echo "Optional arguments:"
	echo "  [source name]  : The name of the source, as stored in VERITAS.tblObserving_Sources"
	echo "  -n             : Print all run numbers that do *not* target the given source"
	echo "                   (only has effect if [source name] is specified)" ; echo
	echo "Examples:" ; echo
	echo "  Print source names for each run:"
	echo "  $ `basename $0` myrunlist.dat"
	echo "  43744 PSRJ2229+6114"
	echo "  43745 PSRJ2229+6114"
	echo "  47146 Tycho"
	echo "  47147 Tycho"
	echo "  47478 Tycho" ; echo
    echo "  Only print the Tycho runs from a runlist:"
    echo "  $ `basename $0` \"Tycho\" myrunlist.dat"
	echo "  47146"
	echo "  47147"
	echo "  47478" ; echo
    echo "  Only print the PSRJ2229+6114 runs from a runlist:"
    echo "  $ `basename $0` \"PSRJ2229+6114\" myrunlist.dat"
	echo "  43744"
	echo "  43745" ; echo
	echo "  Print list of non-Tycho runs, along with their source name:"
	echo "  $ `basename $0` -n \"Tycho\" myrunlist.dat"
	echo "  43744 PSRJ2229+6114"
	echo "  43745 PSRJ2229+6114" ; echo
    echo "Works with pipes : " 
    echo "  $ cat myrunlist.dat | `basename $0` \"Tycho\""
	echo "  47146"
	echo "  47147"
	echo "  47478" ; echo
    exit
fi

function echoerr(){ echo "$@" 1>&2; } #for spitting out error text

# list of run_id's to read in
if [ ! -e $RUNFILE ] ; then
	echo "File '$RUNFILE' could not be found in $PWD , sorry."
	exit	
fi
RUNLIST=`cat $RUNFILE`
#echo "RUNLIST:$RUNLIST"

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
while read -r RUNID RUNTYPE RUNSOURCE; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		
		if [[ "$RUNSOURCE" == "$SOURCE" ]] || [[ "$RUNTYPE" == "$SOURCE" ]] ; then
			if ! $NOTFLAG ; then
				echo "$RUNID"
			fi
		else
			if $NOTFLAG ; then
				if [[ "$RUNSOURCE" == "NOSOURCE" ]]; then
					echo "$RUNID $RUNTYPE" 
				else
					echo "$RUNID $RUNSOURCE"
				fi
			fi
		fi
		
	fi
# This is where the MYSQL command is executed, with the list of requested runs
# You have to do it this way, because using a pipe | calls the command in a
# subshell, and that prevents variables from being saved within the 'while' loop
# http://stackoverflow.com/questions/14585045/is-it-possible-to-avoid-pipes-when-reading-from-mysql-in-bash
done < <($MYSQL -e "USE VERITAS ; SELECT run_id, run_type, source_id FROM tblRun_Info WHERE $SUB ")

exit
