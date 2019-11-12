#!/bin/bash
# analyze data files with eventdisplay

## HARD-CODED VALUES
# run options (will be overwritten in DST mode)
# OPT="-nevents=2000"
## END OF HARD-CODED VALUES

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: evndisp single data file analysis

SPANALYSIS.evndisp.sh <run number> [teltoana] [readcalibdb] [DST]

required parameters:

    <run number>            perform evndisp analysis for this run

optional parameters:

    [teltoana]              use 1 for T1 only, 13 for T1 and T3 only, etc.
                            (default telescope combination is 1234)

    [readcalibdb]           set to 0 to switch off (default is on)

    [DST]                   output format is DST (default OFF=0)

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RUNNUM=$1
[[ "$2" ]] && TELTOANA=$2 || TELTOANA="1234"
if [[ $TELTOANA == "-1" ]]; then
    TELTOANA="1234"
fi
[[ "$3" ]] && CALDB=$3   || CALDB="1"
[[ "$4" ]] && DSTMODE=$4 || DSTMODE="0"

# Check if source vbf file exists
SF=`find -L $VERITAS_DATA_DIR/data -name "$RUNNUM.cvbf"`
if [[ ${#SF} = 0 ]]; then
    echo "ERROR: VERITAS source file $RUNNUM.cvbf not found in $VERITAS_DATA_DIR/data/"
    exit 1
fi

# Set run options
OPT="-runnumber=$RUNNUM -teltoana=$TELTOANA"
if [[ $CALDB -eq 1 ]]; then
    OPT="$OPT -readCalibDB"
else
    OPT="$OPT -calibrationfile calibrationlist.dat"
    OPT="$OPT -nocalibnoproblem"
fi

# DST production
if [[ DSTMODE -ne 0 ]]; then
    OPT="$OPT -runmode=4 -dstfile $VERITAS_DATA_DIR/eventdisplay_output/$RUNNUM.dst.root"
fi

# Run evndisp
echo "$EVNDISPSYS/bin/evndisp $OPT"
$EVNDISPSYS/bin/evndisp $OPT

exit
