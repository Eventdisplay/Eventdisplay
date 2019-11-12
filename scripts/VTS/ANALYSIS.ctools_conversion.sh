#!/bin/bash
# script to run eventdisplay analysis for VTS data

# qsub parameters
h_cpu=29:55:00; h_vmem=2000M; tmpdir_size=5G

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
CTOOLS data conversion: submit jobs from a simple run list

ANALYSIS.ctools_conversion.sh <runlist> <input directory> <output directory> [max chunk length]

required parameters:

    <runlist>              simple run list with one run number per line

    <input directory>      directory where input ROOT files are stored.

    <output directory>     directory where output FITS files will be saved.
    
optional parameters:

    [max chunk length]     maximum allowable length of a chunk
                           each run will be divided up into pieces by time cuts, 
                           then each piece will be broken up into chunks.
                           30s = 30 second-long-chunks
                           5m  = 5 minute-long-chunks
                           2h  = 2 hour-long-chunks
                           A chunk can be shorter than this number, but will never
                           be longer.
                           Default: 5m

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
IDIR=$2
ODIR=$3
[[ "$4" ]] && CHUNKLEN=$4 || CHUNKLEN="5m"


# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
# clean up runlist, VV remove leading spaces   VV remove empty lines   VV only get 5-6 digit numbers
FILES=`cat $RLIST | sed 's/^ *//'            | sed '/^$/d'           | grep -oP "^\d{5,6}"`
#echo "FILES:'$FILES'"

# check input dir exists
if [ ! -d "$IDIR" ] ; then
  echo "Error, input directory '$IDIR' doesn't appear to be a directory, exiting..."
  exit 1
fi

# check that the user has ctools and gammalib
$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.init_gammalib_and_ctools.sh

# check that we can access veripy.py, error out if we can't.
VERSTAT=$( python -c "import veripy ; print veripy" 2>&1 )
if ! echo $VERSTAT | grep -q "<module 'veripy'" ; then
  echo "Error, unable to locate python module 'veripy'."
  echo "Please make sure veripy.py is in \$PYTHONPATH before attempting"
  echo "  to use this script."
  exit 1
fi

# check that we can access the backgrounds fits file
BACKARCHIVE="$VERITAS_EVNDISP_AUX_DIR/RadialAcceptances/ctools_background_archive.fits"
if [ ! -e "$BACKARCHIVE" ] ; then
  echo "Error, unable to locate backgrounds archive '$BACKARCHIVE', exiting..."
  exit 1
fi

# create output dir if needed
mkdir -p $ODIR

# Output directory for error/output
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.CTOOLSDATA"
mkdir -p $LOGDIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.ctools_conversion_sub"

TIMETAG=`date +"%s"`

NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"

for AFILE in $FILES ; do
    echo
    echo "Now starting run $AFILE"
    FSCRIPT="$LOGDIR/EVN.ctools-$AFILE"
  
    sed -e "s|INPUTANASUMROOTFILE|$IDIR/$AFILE.anasum.root|" \
        -e "s|OUTPUTFITSDIRECTORY|$ODIR|" \
        -e "s|MAXIMUMCHUNKLENGTH|$CHUNKLEN|"  \
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
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.$TIMETAG.dat
        echo "RUN $AFILE OLOG $FSCRIPT.log"
    elif [[ "$SUBC" == *simple* ]] ; then
        "$FSCRIPT.sh" |& tee "$FSCRIPT.log"	
    fi

done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.$TIMETAG.dat | sort -u | $SUBC
fi

exit 0

