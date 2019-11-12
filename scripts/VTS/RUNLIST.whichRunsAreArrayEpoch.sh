#!/bin/bash
# from a run list, prints the list of runs that are considered V4 runs, before T1 was moved
# written by Nathan Kelley-Hoskins Aug 2013

#echo "\$#:$#   \$1:$1   \$2:$2   \$3:$3   \$4:$4"

ISPIPEFILE=`readlink /dev/fd/0` # check to see if input is from terminal, or from a pipe
#echo "\$ISPIPEFILE: '$ISPIPEFILE'"
if [[ "$ISPIPEFILE" =~ ^/dev/pts/[0-9]{1,2} ]] ; then # its a terminal (not a pipe)
	if ! [ $# -eq 2 ] ; then # the human didn't add any arguments, and we must tell them so
		echo "Prints the run numbers that are of the specific run versions runs."
		echo "  for just V4 runs, do"
		echo "    $ `basename $0` 4 <file of runs>"
		echo "  for just V5 runs, do"
		echo "    $ `basename $0` 5 <file of runs>"
		echo "  to print all V5 and V6 runs, do"
		echo "    $ `basename $0` 56 <file of runs>"
		exit
	fi
fi

INPUTEPOCH="$1"

# list of run_id's to read in
RUNFILE=$2
if [ ! -e $RUNFILE ] ; then
	echo "File $RUNFILE could not be found in $PWD , sorry."
	exit	
fi
RUNLIST=`cat $RUNFILE`
#echo "RUNLIST:$RUNLIST"

# how should we get the array epochs?
# only change this if you know what you're doing!
METHOD="useparamfile"

#Paramfile method
if [[ "$METHOD" == "useparamfile" ]] ; then
	
	# epoch file to load
	PARAMFILE="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/VERITAS.Epochs.runparameter"
	
	# get only lines that start with '*'
	EPOCHTHRESH=$( cat $PARAMFILE | grep -P "^\s??\*" | grep "EPOCH" | grep -P "V\d" )
	#echo "$EPOCHTHRESH"
	
	# find out what are the smallest and largest epochs to work with
	# so we don't loop over V1, V2, V3.... V10, V11, etc
	MINEPOCH=$( echo "$EPOCHTHRESH" | awk '{ print $3 }' | grep -oP "\d" | awk '{ if(min==""){min=$1}; if($1<min){min=$1};} END {print min }' )
	MAXEPOCH=$( echo "$EPOCHTHRESH" | awk '{ print $3 }' | grep -oP "\d" | awk '{ if(max==""){max=$1}; if($1>max){max=$1};} END {print max }' )
	#echo "EPOCHTHRESH:"
	#echo "$EPOCHTHRESH"
	#echo "MINEPOCH:$MINEPOCH"
	#echo "MAXEPOCH:$MAXEPOCH"
	
	# loop over runs in runlist
	for run in ${RUNLIST[@]} ; do
		
		# loop through all epochs between min and max
		for epoch in $(seq $MINEPOCH $MAXEPOCH) ; do
			
			# check to see if the user wants each epoch
			if [[ $INPUTEPOCH == *$epoch* ]] ; then
				#echo "  testing for $epoch"
				
				# find out run boundaries for this
				MINRUN=$( echo "$EPOCHTHRESH" | grep -P "V$epoch" | awk '{ print $4 }' | grep -oP "\d+" )
				MAXRUN=$( echo "$EPOCHTHRESH" | grep -P "V$epoch" | awk '{ print $5 }' | grep -oP "\d+" )
				#echo "RUN:   '$run'"
				#echo "MINRUN:'$MINRUN'"
				#echo "MAXRUN:'$MAXRUN'"
				if (( "$run" <= "$MAXRUN" )) && (( "$run" >= "$MINRUN" )) ; then
					echo "$run"
					break # break out of epoch loop, but not the runlist loop
				fi
			
			fi 

		done

	done
	
fi

# Hardcoded method
# old, as of 20140616
if [[ "$METHOD" == "hardcoded" ]] ; then

	# if we find 4, 5, or 6 in $1, then set appropriate flags
	V4FLAG=false
	V5FLAG=false
	V6FLAG=false
	if [[ "$INPUTEPOCH" == *4* ]] ; then V4FLAG=true ; fi
	if [[ "$INPUTEPOCH" == *5* ]] ; then V5FLAG=true ; fi
	if [[ "$INPUTEPOCH" == *6* ]] ; then V6FLAG=true ; fi

	# first run of V5 : 46642 
	# first run of V6 : 63373  
	for i in ${RUNLIST[@]} ; do
		if $V4FLAG ; then
			if [ "$i" -le "46641" ] ; then 
				echo "$i"
			fi
		fi
		if $V5FLAG ; then
			if [ "$i" -ge "46642" -a "$i" -le "63372" ] ; then 
				echo "$i"
			fi
		fi
		if $V6FLAG ; then
			if [ "$i" -ge "63373" ] ; then 
				echo "$i"
			fi
		fi
	done

fi

exit
