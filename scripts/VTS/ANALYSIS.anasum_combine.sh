#!/bin/bash
# script to combine anasum files processed in parallel mode

if [[ $# < 3 ]]; then
# begin help message
echo "
ANASUM parallel data analysis: combine parallel-processed anasum runs

ANALYSIS.anasum_combine.sh <anasum run list> <anasum directory> <output file name> [run parameter file]

required parameters:

    <anasum run list>       full anasum run list
                            (with effective areas, file names, etc.)
        
    <anasum directory>      input directory containing anasum root files
        
    <output file name>      name of combined anasum file
                            (written to same location as anasum files)
        
optional parameters:

    [run parameter file]    anasum run parameter file (located in 
                            \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/;
                            default is ANASUM.runparameter)

IMPORTANT! Run only after all ANALYSIS.anasum_parallel.sh jobs have finished!

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RUNLIST=$1
DDIR=$2
OUTFILE=$3
OUTFILE=${OUTFILE%%.root}
[[ "$4" ]] && RUNP=$4 || RUNP="ANASUM.runparameter"

# Check that run list exists
if [ ! -f "$RUNLIST" ]; then
    echo "Error, anasum runlist $RUNLIST not found, exiting..."
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

$EVNDISPSYS/bin/anasum -i 1 -l $RUNLIST -d $DDIR -f $RUNP -o $DDIR/$OUTFILE.root 2>&1 | tee $DDIR/$OUTFILE.log

exit
