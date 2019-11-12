#!/bin/bash
# script to make DST files from a raw data file

# qsub parameters
h_cpu=11:29:00; h_vmem=4000M; tmpdir_size=40G

if [[ $# < 2 ]]; then
# begin help message
echo "
EVNDISP DST maker: submit jobs from a simple run list

UTILITY.make_DST.sh <runlist> <summation window> [pedestal calculation] [LMULT] [output dir] [run parameter file]

required parameters:

    <runlist>               simple run list with one run number per line
    
    <summation window>      FADC trace summation window (in samples)

optional parameters:

    [pedestal calculation]  flag to specify if pedestals/tzeros are calculated
                            (default = 1 = on)

    [LMULT]		    Make DSTs for LMULT calculation (use all events, special low gain calib file) 
			    Default: 0 (off), only the first 5000 events are analysed.

    [output dir] 	    Output directory
			    Default: \$VERITAS_USER_DATA_DIR/analysis/EVD400_DST/\$SUMW/

    [runparameter file]     file with integration window size and reconstruction cuts/methods, expected in $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/
			    Default: EVNDISP.reconstruction.LMULT.SWXX.runparameter, sed -e \"s/XX/\$SUMW/g\"
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
SUMW=$2
[[ "$3" ]] && PED=$3 || PED="1"
[[ "$4" ]] && LMULT=$4 || LMULT="0"
[[ "$5" ]] && ODIR=$5 || ODIR="$VERITAS_USER_DATA_DIR/analysis/EVD400_DST/$SUMW/"
[[ "$6" ]] && RUNPFILE=$6 

# Read runlist
if [[ ! -f "$RLIST" ]]; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
RUNNUMS=`cat $RLIST`

# run scripts are written into this directory
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANADATA"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR


# check for the existence of a run parameter file corresponding
# to the given summation window; if it exists, remove it and replace it
# with the LOWGAIN runparameter file, modified for the given sum window
RUNPDIR="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/"

if [[ -z $RUNPFILE ]]; then

	RUNPFILE="EVNDISP.reconstruction.LMULT.SW$SUMW.runparameter"

	if [ ! -e "$RUNPDIR/$RUNPFILE" ]; then
    		sed -e "s|XX|$SUMW|g" $RUNPDIR/EVNDISP.reconstruction.LMULT.SWXX.runparameter > $RUNPDIR/$RUNPFILE
	fi
fi


# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.make_DST_sub"

#########################################
# loop over all files in files loop
for RUN in $RUNNUMS; do
    echo "Now starting run $RUN"
    FSCRIPT="$LOGDIR/EVN.DST-$RUN-$SUMW"

    sed -e "s|RUNFILE|$RUN|" \
        -e "s|PEDESTALS|$PED|" \
        -e "s|LLLOWGAIN|$LMULT|" \
        -e "s|OUTPUTDIR|$ODIR|" \
        -e "s|RRRRPFILE|$RUNPFILE|" \
        -e "s|SUMWINDOW|$SUMW|" $SUBSCRIPT.sh > $FSCRIPT.sh

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
		echo "RUN $RUN: JOBID $JOBID"
    elif [[ $SUBC == *parallel* ]]; then
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
    elif [[ $SUBC == *simple* ]]; then
        $FSCRIPT.sh |& tee $FSCRIPT.log
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.dat | $SUBC
fi

exit
