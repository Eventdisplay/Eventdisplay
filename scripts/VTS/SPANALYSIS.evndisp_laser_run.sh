#!/bin/bash
# analyze a laser/flasher run (calculate pedestals, calculate gains)

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: analyse a laser/flasher file and calculate
gains and time offsets

SPANALYSIS.evndisp_laser_run.sh <run number> [teltoana] [lasermin] [gain]

required parameters:

    <run number>            calculate gains and time offsets for this run

optional parameters:

    [teltoana]              use 1 for T1 only, 13 for T1 and T3 only, etc.
                            (default telescope combination is 1234)
    
    [lasermin]              minimum image size to identify as laser event
                            (default = 50000; check results in display)
                            
    [gain]                  calculate gains for high (0) or low (1) gain
                            (default = 0 = high gain; EXPERT USER ONLY)

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
[[ "$3" ]] && LASERMIN=$3 || LASERMIN=50000
[[ "$4" ]] && RUNMODE=$4  || RUNMODE=2
if [[ $RUNMODE != 1 ]]; then
    RUNMODE=2   # high-gain mode
else
    RUNMODE=5   # low gain mode
fi
CALIBDIR="$VERITAS_USER_DATA_DIR/"

# Check if source vbf file exists
SF=`find -L $VERITAS_DATA_DIR/data -name "$RUNNUM.cvbf"`
if [ ${#SF} = 0 ]; then
    echo "ERROR: VERITAS source (VBF laser/flasher) file $RUNNUM.cvbf not found in $VERITAS_DATA_DIR/data/"
    exit 1
fi

# Run options
OPT="-runmode=$RUNMODE -runnumber=$RUNNUM -lasermin=$LASERMIN -calibrationsumwindow=18 -calibrationsumfirst=2 -reconstructionparameter EVNDISP.reconstruction.SW18_noDoublePass.runparameter -calibrationdirectory $CALIBDIR -writeextracalibtree"

# calculate pedestals (for high gain only)
if [[ $RUNMODE == 2 ]]; then
    echo "Calculating pedestals for run $RUNNUM"
    $EVNDISPSYS/scripts/VTS/SPANALYSIS.evndisp_pedestal_events.sh $RUNNUM
fi

# calculate gains, looping over all telescopes
TELTOANA=`echo $TELTOANA | fold -w1`
for i in $TELTOANA; do
    echo "Calculating gains for run $RUNNUM, telescope $i"
    $EVNDISPSYS/bin/evndisp -teltoana=$i $OPT
done

exit
