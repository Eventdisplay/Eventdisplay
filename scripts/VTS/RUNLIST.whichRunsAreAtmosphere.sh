#!/bin/bash
# from a run list, prints the list of runs that were taken in a specific atmosphere, summer(22) or winter(21)
# written by Nathan Kelley-Hoskins Sept 2013

# check to see if input is from terminal, or from a pipe
ISPIPEFILE=`readlink /dev/fd/0`
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} && $# < 2 ]]; then # its a terminal (not a pipe)
    echo
    echo "From a runlist or pipe, prints the run numbers that are of a particular atmosphere."
    echo " $ `basename $0` [w|21|s|22] <file of runs>" ; echo
    echo "w = 21 = winter, s = 22 = summer" ; echo
    echo "Print list of summer runs:"
    echo " $ `basename $0` s myrunlist.dat" ; echo
    echo "Print list of winter runs:"
    echo " $ `basename $0` 21 myrunlist.dat" ; echo
    echo "Works with pipes : " 
    echo " $ cat myrunlist.dat | `basename $0` w" ; echo
    echo "Summer/winter transition dates taken from $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/VERITAS.Epochs.runparameter"
    exit
fi

function echoerr(){ echo "$@" 1>&2; } #for spitting out error text

# list of run_id's to read in
RUNFILE=$2
if [ ! -e $RUNFILE ] ; then
	echo "File $RUNFILE could not be found in $PWD , sorry."
	exit	
fi
RUNLIST=`cat $RUNFILE`
#echo "RUNLIST:$RUNLIST"

SUMMFLAG=false
WINTFLAG=false
LOWARG=`echo "$1" | tr '[A-Z]' '[a-z]'` # make all uppercase letters in arg 1 lowercase, for easier handling
#echo "\$LOWARG: '$LOWARG'"
if [[ "$LOWARG" == *w* ]] || [[ "$LOWARG" == "21" ]] ; then
	WINTFLAG=true
fi
if [[ "$LOWARG" == *s* ]] || [[ "$LOWARG" == "22" ]] ; then
	SUMMFLAG=true
fi
#echo "\$WINTFLAG:$WINTFLAG    \$SUMMFLAG:$SUMMFLAG"
if $WINTFLAG && $SUMMFLAG ; then
	echo "recognized both summer[s] and winter[w] in arg '$1', can only specify one or the other"
	exit
elif ! $WINTFLAG && ! $SUMMFLAG ; then
	echo "Need to specifiy an atmosphere: Argument 1 '$1' needs to be either 'w' or 's'."
	exit
fi

# how should we calculate which atmosphere to use?
METHOD="useparamfile"

function IsWinter {
    local date="$1"
    local month=${date:4:2}
	
	# use boundaries from param file
	if [[ "$METHOD" == "useparamfile" ]] ; then
		
		# epoch file to load
		local PARAMFILE="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/VERITAS.Epochs.runparameter"
		
		# get only lines that start with '*'
		local ATMOTHRESH=$( cat $PARAMFILE | grep -P "^\s??\*" | grep "ATMOSPHERE" )
		#echo "$ATMOTHRESH"
		
		# flag for if we found the atmo
		local FOUNDATMO=false
		
		# other vars
		local ATMOCODE=""
		local MINDATE=""
		local MAXDATE=""
		(IFS='
'
		for line in $ATMOTHRESH ; do
			#echoerr "line:$line"
			ATMOCODE=$( echo "$line" | awk '{ print $3 }' | grep -oP "\d+" )
			MINDATE=$(  echo "$line" | awk '{ print $4 }' | tr -d '-' | grep -oP "\d+" )
			MAXDATE=$(  echo "$line" | awk '{ print $5 }' | tr -d '-' | grep -oP "\d+" )
			#echoerr "  ATMOCODE:$ATMOCODE"
			#echoerr "  MINDATE: $MINDATE"
			#echoerr "  MAXDATE: $MAXDATE"
			if (( "$date" >= "$MINDATE" )) && (( "$date" <= "$MAXDATE" )) ; then
				
				# winter
				if [[ "$ATMOCODE" == "21" ]] ; then
					echo 1 	
					#echoerr "$date - winter!"
					FOUNDATMO=true
					break
				# summer
				elif [[ "$ATMOCODE" == "22" ]] ; then
					echo 2
					#echoerr "$date - summer!"
					FOUNDATMO=true
					break
				fi
			fi
		done 
		)
		# 3 = did not find valid atmo range
		if [ ! $FOUNDATMO ] ; then
			echo 3
		fi
		
	fi
	
	# use hardcoded boundaries
	if [[ "$METHOD" == "hardcoded" ]] ; then
		if  [ "$date" -gt "20071026" ] && [ "$date" -lt "20080420" ] ||
			[ "$date" -gt "20081113" ] && [ "$date" -lt "20090509" ] ||
			[ "$date" -gt "20091102" ] && [ "$date" -lt "20100428" ] ||
			[ "$date" -gt "20101023" ] && [ "$date" -lt "20110418" ] ||
			[ "$date" -gt "20111110" ] && [ "$date" -lt "20120506" ] ||
			[ "$date" -gt "20121029" ] && [ "$date" -lt "20130425" ] ; then
			#echo true
			echo 1  # winter
		elif [ "$date" -ge "20130425" -o "$date" -le "20071026" ] ; then
			# don't have specific dates for summer/winter boundary, so we will generalize to the months
			# may through october inclusive is 'summer'
			if   [ "$month" -ge 5 -a "$month" -le 10 ] ; then
				echo 2 # summer
			# november through april inclusive is 'winter'
			elif [ "$month" -le 4 -o "$month" -ge 11 ] ; then
				echo 1 # winter
			else
				echo 3 # unassignable
			fi
		else
			#echo false
			echo 2 # summer
		fi
	fi
}

function badAtmosphere {
	echoerr "Error, in 'RUNLIST.whichRunsAreAtmosphere.sh', run $2 is too new and could not have its atmosphere assigned on date $1, exiting..."
	exit 1
}

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
while read -r RUNID RUNDATE ; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		
		# decode the date tag
		read YY MM DD HH MI SE <<< ${RUNDATE//[-:]/ }
		#echo "  YEARMONTHDAY:$YY$MM$DD"

		# get the atmosphere code
        STATUSFLAG=`IsWinter "$YY$MM$DD"`
		#echo "$RUNID '$STATUSFLAG'"
        
		# did the user ask for summer runs?
		if $SUMMFLAG ; then
			if   [[ "$STATUSFLAG" == "2" ]] ; then echo "$RUNID"
			elif [[ "$STATUSFLAG" == "3" ]] ; then badAtmosphere "$YY$MM$DD" "$RUNID"
			fi
		# did the user ask for winter runs?
		elif $WINTFLAG ; then
			if   [[ "$STATUSFLAG" == "1" ]] ; then echo "$RUNID"
			elif [[ "$STATUSFLAG" == "3" ]] ; then badAtmosphere "$YY$MM$DD" "$RUNID"
			fi
		fi
		
	fi
# This is where the MYSQL command is executed, with the list of requested runs
# You have to do it this way, because using a pipe | calls the command in a
# subshell, and that prevents variables from being saved within the 'while' loop
# http://stackoverflow.com/questions/14585045/is-it-possible-to-avoid-pipes-when-reading-from-mysql-in-bash
done < <($MYSQL -e "USE VERITAS ; SELECT run_id, data_start_time FROM tblRun_Info WHERE $SUB")

exit
