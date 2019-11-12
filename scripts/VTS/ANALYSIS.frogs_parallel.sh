#!/bin/bash
# script to run eventdisplay analysis with FROGS

# qsub parameters
h_cpu=00:29:30; h_vmem=2000M; tmpdir_size=10G

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP data analysis: evndisp FROGS analysis for a simple run list

ANALYSIS.evndisp_frogs.sh <runlist> [output directory] [mscw directory] [runparameter file] [pedestals] [VPM]

required parameters:

    <runlist>               simple run list with one run number per line

optional parameters:

    [output directory]    directory where output ROOT files will be stored
    
    [mscw directory]      directory which contains mscw_energy files
    
    [runparameter file]   file with integration window size and reconstruction cuts/methods, expected in $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/

                          Default: EVNDISP.reconstruction.runparameter (long sumwindow -> for use with CARE IRFs; DISP disabled )

                          other options:

                          EVNDISP.reconstruction.runparameter.DISP              (long sumwindow -> for use with CARE IRFs;
                                                                                 DISP enabled, use RecID 1 in later stages to access it)

    [pedestals]           calculate pedestals (default); set to 0 to skip

    [VPM]                 set to 0 to switch off (default is on)
	
	[eventsplit]          break the input file up into chunks at most [eventsplit] events long,
						  will need to be remerged with ANALYSIS.frogs_merge.sh

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RLIST=$1
[[ "$2" ]] && ODIR=$2       || ODIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$3" ]] && MSCWDIR=$3    || MSCWDIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$4" ]] && ACUTS=$4      || ACUTS=EVNDISP.reconstruction.runparameter
[[ "$4" ]] && CALIB=$5      || CALIB=5
[[ "$5" ]] && VPM=$6        || VPM=1
[[ "$6" ]] && EVENTSPLIT=$7 || EVENTSPLIT=0

# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
RUNNUMS=`cat $RLIST`

# Log file directory
DATE=`date +"%y%m%d"`

LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/FROGS"
#LOGDIR="frogslogs"
LOGDIR="/afs/ifh.de/group/cta/scratch/nkelhos/fastfrogs/logs" # HARDCODE

echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

# output directory
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.frogs_parallel_sub"

TIMETAG=`date +"%s"`

NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

joblist=()

# loop over all files in files loop to figure out how many jobs to submit
for RUN in $RUNNUMS; do

	DATAEVENTS=$( $EVNDISPSYS/bin/frogsGetNEvents $MSCWDIR/$RUN.mscw.root "data" | tail -n 1 | awk '{ print $2 }' ) 
	echo "DATAEVENTS:$DATAEVENTS"
	
	if [ $EVENTSPLIT == 0 ] ; then
		maxjobs=1
	else
		maxjobs=$( python -c "from math import ceil; print int(ceil( $DATAEVENTS / ( 1.0 * $EVENTSPLIT ) ) )" )
	fi

	MERGELIST="$ODIR/$RUN.mergelist"
	rm -f "$MERGELIST"
	
	echo "$RUN: $maxjobs"
	ijob=0
	debugbreak=0
	for ijob in `seq 1 $maxjobs` ; do
		
		# stop early for debugging, 
		#debugbreak=$((debugbreak+1))
		#if [ $debugbreak -gt 3 ] ; then continue ; fi
		#echo "  job:$ijob"
		firstevent=$( python -c "print ( $ijob - 1 ) * $EVENTSPLIT " ) 
		#echo "     first=$firstevent"
		#echo "     neven=$EVENTSPLIT"
		OUTPUTFILE="$ODIR/$RUN-$ijob.extracted.root"
		joblist+=("jobiter=$ijob:run=$RUN:inputfile=$MSCWDIR/$RUN.mscw.root:startevent=$firstevent:nevents=$EVENTSPLIT:outputfile=$OUTPUTFILE:")
		echo "$OUTPUTFILE" >> "$MERGELIST"
	done
	
done

echo
echo "${#joblist[@]} jobs:"
for job in ${joblist[@]} ; do
	echo "$job"
done

#for RUN in $RUNNUMS; do
for job in ${joblist[@]} ; do
		
    echo "jobline '$job'"

	echo
	RUN=$( echo "$job"     | grep -oP "run=\d+\:"     | grep -oP "\d+" ) 
	JOBITER=$( echo "$job" | grep -oP "jobiter=\d+\:" | grep -oP "\d+" ) 

	#echo " RERUN:$RUN"
    FSCRIPT="$LOGDIR/EVN.data-$RUN-$JOBITER"
	
	
    sed -e "s|RUNFILE|$RUN|"             \
        -e "s|CALIBRATIONOPTION|$CALIB|" \
        -e "s|OUTPUTDIRECTORY|$ODIR|"    \
        -e "s|MSCWDIRECTORY|$MSCWDIR|"   \
        -e "s|USEVPMPOINTING|$VPM|"      \
		-e "s|SETSETSTRING|$job|" $SUBSCRIPT.sh > $FSCRIPT.sh

    chmod u+x $FSCRIPT.sh
    echo "  JOBITER:$JOBITER = SCRIPT:'$FSCRIPT.sh'"

    # run locally or on cluster
    SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
    SUBC=`eval "echo \"$SUBC\""`
    if [[ $SUBC == *"ERROR"* ]]; then
        echo $SUBC
        exit
    fi
    echo "  command:'$SUBC $FSCRIPT.sh'"
    if [[ $SUBC == *qsub* ]]; then
		
        JOBID=`$SUBC $FSCRIPT.sh`
		
        # account for -terse changing the job number format
        if [[ $SUBC != *-terse* ]] ; then
            echo "without -terse!"      # need to match VVVVVVVV  8539483  and 3843483.1-4:2
            JOBID=$( echo "$JOBID" | grep -oP "Your job [0-9.-:]+" | awk '{ print $3 }' )
        fi
    
        echo "RUN $RUN JOBID $JOBID"
        echo "RUN $RUN SCRIPT $FSCRIPT.sh"
        if [[ $SUBC != */dev/null* ]] ; then
            echo "RUN $RUN OLOG $FSCRIPT.sh.o$JOBID"
            echo "RUN $RUN ELOG $FSCRIPT.sh.e$JOBID"
        fi
    fi
        
	if [[ $SUBC == *parallel* ]]; then
	echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.$TIMETAG.dat
	fi

	if [[ "$SUBC" == *simple* ]] ; then
		"$FSCRIPT.sh" |& tee "$FSCRIPT.log"
	fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.$TIMETAG.dat | $SUBC
fi

exit

