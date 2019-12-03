#!/bin/bash
# script to analyse VTS data files with lookup tables

# qsub parameters
h_cpu=00:29:00; h_vmem=2000M; tmpdir_size=4G

# EventDisplay version
"$EVNDISPSYS"/bin/mscw_energy--version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`"$EVNDISPSYS"/bin/mscw_energy --version | tr -d .`
else
    EDVERSION="g500"
fi

if [ $# -lt 2 ]; then
# begin help message
echo "
MSCW_ENERGY data analysis: submit jobs from a simple run list

ANALYSIS.mscw_energy.sh <table file> <runlist> [evndisp directory] [Rec ID] [output directory]

required parameters:

    <table file>            mscw_energy lookup table file. Expected in $VERITAS_EVNDISP_AUX_DIR/Tables/ .
			    For example:
				table-v480-auxv01-CARE_June1425-ATM61-V6-DISP.root
				table-v480-auxv01-CARE_June1425-ATMX-VX-DISP.root
				table-v480-auxv01-GRISU-SW6-ATMX-VX-GEO.root
			    ATMX and VX will be replaced by the proper atmosphere/detector configuration.			

    <runlist>               simple run list with one run number per line.    
    
optional parameters:
    
    [evndisp directory]     directory containing evndisp output ROOT files.
			    Default: $VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/

    [Rec ID]                reconstruction ID. Default 0.
                            (see EVNDISP.reconstruction.runparameter)
    
    [output directory]      directory where mscw.root files are written
                            default: <evndisp directory>
                            If RecID given: <evndisp directory>/RecID#

--------------------------------------------------------------------------------
"

#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# create extra stdout for duplication of command output
# look for ">&5" below
exec 5>&1

# Parse command line arguments
TABFILE=$1
TABFILE=${TABFILE%%.root}.root
RLIST=$2
[[ "$3" ]] && INPUTDIR=$3 || INPUTDIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$4" ]] && ID=$4 || ID=0
if [[ "$5" ]]; then
    ODIR=$5
elif [[ "$4" ]]; then
    # user passed a Rec ID value
    ODIR=$INPUTDIR/RecID$ID
else
    ODIR=$INPUTDIR
fi

# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
FILES=`cat "$RLIST"`

NRUNS=`cat "$RLIST" | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

# make output directory if it doesn't exist
mkdir -p "$ODIR"
echo -e "Output files will be written to:\n $ODIR"

# run scripts are written into this directory
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/MSCW.ANADATA"
mkdir -p "$LOGDIR"
echo -e "Log files will be written to:\n $LOGDIR"

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.mscw_energy_sub"

TIMETAG=`date +"%s"`

#########################################
# loop over all files in files loop
for AFILE in $FILES
do
    BFILE="$INPUTDIR/$AFILE.root"
    echo "Now analysing $BFILE (ID=$ID)"

     # get array epoch, atmosphere and telescope combination for this run
    RUNINFO=`"$EVNDISPSYS"/bin/printRunParameter "$INPUTDIR/$AFILE.root" -runinfo`
    EPOCH=`echo "$RUNINFO" | awk '{print $(1)}'`
    ATMO=`echo "$RUNINFO" | awk '{print $(2)}'`
    if [[ $ATMO == *error* ]]; then
       echo "error finding atmosphere; skipping run $AFILE"
       continue
    fi
    echo "RUN $AFILE at epoch $EPOCH and atmosphere $ATMO"
    TABFILERUN=${TABFILE/VX/${EPOCH}}
    TABFILERUN=${TABFILERUN/ATMX/ATM${ATMO}}

    # Check that table file exists
    if [[ "$TABFILERUN" == `basename "$TABFILERUN"` ]]; then
        TABFILERUN="$VERITAS_EVNDISP_AUX_DIR/Tables/$TABFILERUN"
    fi
    if [ ! -f "$TABFILERUN" ]; then
        echo "Error, table file '$TABFILERUN' not found, skipping run..."
        continue
    else
        echo "Table file: $TABFILERUN" 
        echo "Table file: $TABFILERUN" 
    fi

    FSCRIPT="$LOGDIR/MSCW.data-ID$ID-$AFILE-$(date +%s)"
    rm -f "$FSCRIPT.sh"

    sed -e "s|TABLEFILE|$TABFILE|" \
        -e "s|RECONSTRUCTIONID|$ID|" \
        -e "s|OUTPUTDIRECTORY|$ODIR|" \
        -e "s|EVNDISPFILE|$BFILE|" "$SUBSCRIPT.sh" > "$FSCRIPT.sh"

    chmod u+x "$FSCRIPT.sh"
    echo "$FSCRIPT.sh"

    # run locally or on cluster
    SUBC=`"$EVNDISPSYS"/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
    SUBC=`eval "echo \"$SUBC\""`
    if [[ $SUBC == *"ERROR"* ]]; then
        echo "$SUBC"
        exit
    fi
    echo "Submission command: $SUBC"
    if [[ $SUBC == *qsub* ]]; then
        JOBID=`$SUBC "$FSCRIPT.sh"`
        
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
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> "$LOGDIR/runscripts.$TIMETAG.dat"
        echo "RUN $AFILE OLOG $FSCRIPT.log"
    elif [[ "$SUBC" == *simple* ]] ; then
	"$FSCRIPT.sh" |& tee "$FSCRIPT.log"
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat "$LOGDIR/runscripts.$TIMETAG.dat" | "$SUBC"
fi

exit
