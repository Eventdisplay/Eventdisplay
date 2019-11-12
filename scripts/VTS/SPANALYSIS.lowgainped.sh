#!/bin/bash
# script to run eventdisplay analysis for VTS data

# qsub parameters
h_cpu=11:59:00; h_vmem=2000M; tmpdir_size=25G

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
Low gain pedestal analysis: submit jobs from a simple run list

SPANALYSIS.lowgainped.sh <runlist> <start of sumwindow> [calib directory] [sum window] [teltoana] [nevents]

required parameters:

    <runlist>               simple run list with one run number per line
    
    <start of sumwindow>    default 100 for calib runs with 128 samples, 40 for calib runs with 64 samples.

optional parameters:
    
    [calib directory]       directory where output files will be stored. Default $VERITAS_EVNDISP_AUX_DIR/ .

    [sum window]            calibration sum window. Default 20.

    [teltoana]              restrict telescope combination to be analyzed:
                            e.g.: teltoana=123 (for tel. 1,2,3), 234, ...
                            default is 1234 (all telescopes)

    [nevents] 		    number of events to be used, default all

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`

# create extra stdout for duplication of command output
# look for ">&5" below
exec 5>&1

# Parse command line arguments
RLIST=$1
SUMSTART=$2
[[ "$3" ]] && ODIR=$3 || ODIR="$VERITAS_EVNDISP_AUX_DIR"
mkdir -p $ODIR
[[ "$4" ]] && SUMWINDOW=$4   || SUMWINNDOW=20
[[ "$5" ]] && TELTOANA=$5 || TELTOANA=1234
[[ "$6" ]] && NEVENTS=$6 || NEVENTS=-1

# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
FILES=`cat $RLIST`

# Output directory for error/output
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.LGAINPED"
mkdir -p $LOGDIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/SPANALYSIS.lowgainped_sub"


NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

#########################################
# loop over all files in files loop
for AFILE in $FILES
do
    echo "Now starting run $AFILE"
    FSCRIPT="$LOGDIR/EVN.lped-$AFILE"

    sed -e "s|RUNFILE|$AFILE|"              \
        -e "s|OUTPUTDIRECTORY|$ODIR|"       \
        -e "s|TELTOANACOMB|$TELTOANA|"      \
	-e "s|NNNN|$NEVENTS|" 		    \
	-e  "s|CALIBFIRST|$SUMSTART|" 	    \
	-e  "s|CALIBSUMWINDOW|$SUMWINDOW|"  \
	$SUBSCRIPT.sh > $FSCRIPT.sh

    chmod u+x $FSCRIPT.sh
    echo $FSCRIPT.sh
	# output selected input during submission:
	echo "calibrationsumwindow $SUMWINDOW, start at $SUMSTART, $NEVENTS events "
	echo "Output: $ODIR"

	if [[ $TELTOANA == "1234" ]]; then
	echo "Analyzed telescopes: $TELTOANA (default, all telescopes)"
	else
	echo "Analyzed telescopes: $TELTOANA"
	fi 

    # run locally or on cluster
    SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
    SUBC=`eval "echo \"$SUBC\""`
    if [[ $SUBC == *"ERROR"* ]]; then
        echo $SUBC
        exit
    fi
    echo $SUBC
    if [[ $SUBC == *qsub* ]]; then
        JOBID=`$SUBC $FSCRIPT.sh`
        # account for -terse changing the job number format
        if [[ $SUBC != *-terse* ]] ; then
            echo "without -terse!"      # need to match VVVVVVVV  8539483  and 3843483.1-4:2
            JOBID=$( echo "$JOBID" | grep -oP "Your job [0-9.-:]+" | awk '{ print $3 }' )
        fi
        
        echo "RUN $AFILE JOBID $JOBID"
        echo "RUN $AFILE SCRIPT $FSCRIPT.sh"
        if [[ $SUBC != */dev/null* ]] ; then
            echo "RUN $AFILE OLOG $FSCRIPT.sh.o$JOBID"
            echo "RUN $AFILE ELOG $FSCRIPT.sh.e$JOBID"
        fi
    elif [[ $SUBC == *parallel* ]]; then
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
        echo "RUN $AFILE OLOG $FSCRIPT.log"
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.dat | $SUBC
fi

exit
