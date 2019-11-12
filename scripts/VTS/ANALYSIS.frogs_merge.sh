#!/bin/bash
# script to run eventdisplay analysis for VTS data

# qsub parameters
h_cpu=00:10:00; h_vmem=2000M; tmpdir_size=25G

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
merge Frogs Parallel datafiles by submiting jobs from a simple run list

ANALYSIS.frogs_merge.sh <runlist> <input directory>

required parameters:

    <runlist>              simple run list with one run number per line
    
	<input directory>      directory containing the #####.mergelist files
			   
will save the resulting '#####.frogs.root' to the <input directory>
	 

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
RLIST=$1
ODIR=$2

# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
FILES=`cat $RLIST`

# Output directory for error/output
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANADATA"
LOGDIR="/afs/ifh.de/group/cta/scratch/nkelhos/fastfrogs/logs/"
mkdir -p $LOGDIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.frogs_merge_sub"

SECONDS=`date +"%s"`

NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

#########################################
# loop over all files in files loop
for AFILE in $FILES
do
	echo
    echo "Now starting run $AFILE"
    FSCRIPT="$LOGDIR/EVN.data-$AFILE"
	
	MERGELIST=$ODIR/$AFILE.mergelist
	OUTNAME=$ODIR/$AFILE.frogs.root

    sed -e "s|SUBRUN|$AFILE|"           \
		-e "s|SUBOUTDIR|$ODIR|"         \
        -e "s|SUBMERGELIST|$MERGELIST|" \
		-e "s|SUBOUTPUTFILE|$OUTNAME|"  \
        $SUBSCRIPT.sh > $FSCRIPT.sh

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
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.$SECONDS.dat
        echo "RUN $AFILE OLOG $FSCRIPT.log"
    elif [[ "$SUBC" == *simple* ]] ; then
	"$FSCRIPT.sh" |& tee "$FSCRIPT.log"	
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.$SECONDS.dat | sort -u | $SUBC
fi

exit
