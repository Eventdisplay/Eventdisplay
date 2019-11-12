#!/bin/bash
# script to analyse data files with anasum

# qsub parameters
h_cpu=8:00:00; h_vmem=4000M; tmpdir_size=10G

if [[ $# < 3 ]]; then
# begin help message
echo "
ANASUM data analysis: submit jobs from an anasum run list

ANALYSIS.anasum.sh <anasum run list> <output directory> <output file name>
 [run parameter file] [mscw directory]

required parameters:

    <anasum run list>       full anasum run list
    
    <output directory>      anasum output files are written to this directory
    
    <output file name>      name of combined anasum file
                            (written to same location as anasum files)
    
optional parameters:

    [run parameter file]    anasum run parameter file (located in 
                            \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/;
                            default is ANASUM.runparameter)

    [mscw directory]        directory containing the mscw.root files

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
FLIST=$1
ODIR=$2
ONAME=$3
ONAME=${ONAME%%.root}
[[ "$4" ]] && RUNP=$4  || RUNP="ANASUM.runparameter"
[[ "$5" ]] && INDIR=$5 || INDIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/RecID0"

# Check that run list exists
if [[ ! -f "$FLIST" ]]; then
    echo "Error, anasum runlist $FLIST not found, exiting..."
    exit 1
fi

# Check that run parameter file exists
if [[ "$RUNP" == `basename $RUNP` ]]; then
    RUNP="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/$RUNP"
fi
if [[ ! -f "$RUNP" ]]; then
    echo "Error, anasum run parameter file not found, exiting..."
    exit 1
fi

# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/ANASUM.ANADATA"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

# output directory
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.anasum_sub"

TIMETAG=`date +"%s"`

FSCRIPT="$LOGDIR/ANA.$ONAME-$(date +%s)"
sed -e "s|FILELIST|$FLIST|" \
    -e "s|DATADIR|$INDIR|"  \
    -e "s|OUTDIR|$ODIR|"    \
    -e "s|OUTNAME|$ONAME|"  \
    -e "s|RUNPARAM|$RUNP|" $SUBSCRIPT.sh > $FSCRIPT.sh

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
    
elif [[ $SUBC == *parallel* ]]; then
    echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.$TIMETAG.dat
    cat $LOGDIR/runscripts.$TIMETAG.dat | $SUBC

elif [[ $SUBC == *simple* ]]; then
	"$FSCRIPT.sh" |& tee $FSCRIPT.log
fi

elif [[ $SUBC == *simple* ]]; then
	"$FSCRIPT.sh" |& tee $FSCRIPT.log
fi

exit
