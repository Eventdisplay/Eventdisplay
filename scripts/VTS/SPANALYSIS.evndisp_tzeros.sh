#!/bin/bash
# calculate mean tzeros

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: calculate mean tzeros from a data file

VTS.EVNDISP.analyse_tzeros <run number> [teltoana] [readcalibdb]

required parameters:

    <run number>            calculate tzeros for this run

optional parameters:

    [teltoana]              use 1 for T1 only, 13 for T1 and T3 only, etc.
                            (default telescope combination is 1234)

    [readcalibdb]           set to 0 to switch off (default is on)

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

# Check if source vbf file exists
#SF=`find -L $VERITAS_DATA_DIR/data -name "$RUNNUM.cvbf"`
#if [[ ${#SF} = 0 ]]; then
#    echo "ERROR: VERITAS source file $RUNNUM.cvbf not found in $VERITAS_DATA_DIR/data/"
#    exit 1
#fi

if [[ $CALDB == "1" ]]; then
    OPT="$OPT -readCalibDB"
else
    OPT="$OPT -nocalibnoproblem"
fi
echo $OPT

# run options
OPT="-runmode=7 -runnumber=$RUNNUM -teltoana=$TELTOANA -sumwindowaveragetime=6 $OPT -nevents=1000 -calibrationdirectory $VERITAS_USER_DATA_DIR"

# Run evndisp
echo "$EVNDISPSYS/bin/evndisp $OPT"
$EVNDISPSYS/bin/evndisp $OPT

exit
