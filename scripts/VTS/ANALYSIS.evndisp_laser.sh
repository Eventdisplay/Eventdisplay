#!/bin/bash
# script run eventdisplay laser analysis with a queue system

# qsub parameters
h_cpu=11:29:00; h_vmem=2000M; tmpdir_size=5G

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP laser analysis: submit jobs from a simple run list

ANALYSIS.evndisp_laser.sh <runlist> [teltoana]

required parameters:

    <runlist>           simple run list with one run number per line
                        runlist should contain laser run numbers

                        example for run list:
                        48626
                        58453
                        61429

    [teltoana]          restrict telescope combination to be analyzed:
                        e.g.: teltoana=123 (for tel. 1,2,3), 234, ...
                        default is 1234 (all telescopes)

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

[[ "$2" ]] && TELTOANA=$2 || TELTOANA=1234

# Read runlist
if [[ ! -f "$RLIST" ]] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
RUNNUMS=`cat "$RLIST"`
echo "Laser files to analyze:"
echo "$RUNNUMS"

# Output directory for error/output from batch system
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANADATA"
echo "Log files will be written to: $LOGDIR"
mkdir -p "$LOGDIR"

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.evndisp_laser_sub"

#########################################
# loop over all files in files loop
for RUN in $RUNNUMS; do
    # check if laser file exists
    DFILE=`find -L "$VERITAS_DATA_DIR/data/" -name "$RUN.cvbf"`
    if [ -z "$DFILE" ]; then
        echo "Error: laser vbf file not found for run $RUN"
    else
        echo "Now starting laser run $RUN"

        # output selected input during submission:
        if [[ $TELTOANA == "1234" ]]; then
        echo "Analyzed telescopes: $TELTOANA (default, all telescopes)"
        else
        echo "Analyzed telescopes: $TELTOANA"
        fi
        
        FSCRIPT="$LOGDIR/EVN.laser-$RUN"

        sed -e "s|RUNFILE|$RUN|" \
            -e "s|TELTOANACOMB|$TELTOANA|" \
            -e "s|LOGDIRECTORY|$LOGDIR|" "$SUBSCRIPT.sh" > "$FSCRIPT.sh"

        chmod u+x "$FSCRIPT.sh"
        echo "$FSCRIPT.sh"

        # run locally or on cluster
        SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
        SUBC=`eval "echo \"$SUBC\""`
        if [[ $SUBC == *"ERROR"* ]]; then
            echo "$SUBC"
            exit
        fi
        if [[ $SUBC == *qsub* ]]; then
            JOBID=`$SUBC $FSCRIPT.sh`
            
            # account for -terse changing the job number format
            if [[ $SUBC != *-terse* ]] ; then
                echo "without -terse!"      # need to match VVVVVVVV  8539483  and 3843483.1-4:2
                JOBID=$( echo "$JOBID" | grep -oP "Your job [0-9.-:]+" | awk '{ print $3 }' )
            fi
            
            echo "RUN $RUN: JOBID $JOBID"
        elif [[ $SUBC == *parallel* ]]; then
            echo "$FSCRIPT.sh &> $FSCRIPT.log" >> "$LOGDIR/runscripts.dat"
        elif [[ "$SUBC" == *simple* ]] ; then
	    "$FSCRIPT.sh" |& tee "$FSCRIPT.log"	
        fi
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat "$LOGDIR/runscripts.dat" | "$SUBC"
fi

exit
