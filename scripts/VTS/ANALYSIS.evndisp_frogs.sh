#!/bin/bash
# script to run eventdisplay analysis with FROGS

# qsub parameters
h_cpu=48:00:00; h_vmem=2000M; tmpdir_size=10G

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
[[ "$2" ]] && ODIR=$2    || ODIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$3" ]] && MSCWDIR=$3 || MSCWDIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$4" ]] && ACUTS=$4   || ACUTS=EVNDISP.reconstruction.runparameter
[[ "$4" ]] && CALIB=$5   || CALIB=5
[[ "$5" ]] && VPM=$6     || VPM=1

# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
RUNNUMS=`cat $RLIST`

# Log file directory
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/FROGS"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

# output directory
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.evndisp_frogs_sub"

TIMETAG=`date +"%s"`

NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

# loop over all files in files loop
for RUN in $RUNNUMS; do
    echo "Now starting run $RUN"
    FSCRIPT="$LOGDIR/EVN.data-$RUN"

    # get run array epoch using a run info function
    #EPOCH=`getRunArrayVersion $RUN`

    sed -e "s|RUNFILE|$RUN|"             \
        -e "s|CALIBRATIONOPTION|$CALIB|" \
        -e "s|OUTPUTDIRECTORY|$ODIR|"    \
        -e "s|MSCWDIRECTORY|$MSCWDIR|"   \
        -e "s|USEVPMPOINTING|$VPM|" $SUBSCRIPT.sh > $FSCRIPT.sh

    chmod u+x $FSCRIPT.sh
    echo $FSCRIPT.sh

    # run locally or on cluster
    SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
    SUBC=`eval "echo \"$SUBC\""`
    if [[ $SUBC == *"ERROR"* ]]; then
        echo $SUBC
        exit
    fi
    if [[ $SUBC == *qsub* ]]; then
        JOBID=`$SUBC $FSCRIPT.sh`
    fi
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

