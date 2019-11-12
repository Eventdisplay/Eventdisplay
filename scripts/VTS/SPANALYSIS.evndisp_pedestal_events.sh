#!/bin/bash
# calculate pedestals

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: calculate pedestals from a data file
                                  (for high and low gain)

SPANALYSIS.evndisp_pedestal_events.sh <run number> [samples] [calibdir]

required parameters:

    <run number>        calculate pedestals for this run

optional parameters:

    [lowgain]		for low gain pedestal calculation. Set to either 64 or 128 samples for low gain mode,
                        0 for high gain mode (expert user only!)

    [calibdir]		eg. \$VERITAS_USER_DATA_DIR/LowGainCalibration for low gain calibration.
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
[[ "$2" ]] && SAMPLES=$2 || SAMPLES="0"
[[ "$3" ]] && CALIBDIR=$3 || CALIBDIR="$VERITAS_USER_DATA_DIR/"

OPT=" -calibrationsumwindow=20 -calibrationdirectory $CALIBDIR " 

#high gain mode
if [[ $SAMPLES == "0" ]]; then
	OPT="$OPT -runmode=1 -calibrationsumfirst=0 "

#low gain mode
else 
	OPT="$OPT -runmode=6 -reconstructionparameter EVNDISP.reconstruction.LGCalibration.runparameter "

	if [[ $SAMPLES == "64" ]]; then
		OPT="$OPT -calibrationsumfirst=40 "
	elif  [[ $SAMPLES == "128" ]]; then
		OPT="$OPT -calibrationsumfirst=100 "
	elif [[ $SAMPLES == "100" ]]; then
		OPT="$OPT -calibrationsumfirst=75 "
	else 
		echo "Invalid number of samples given, please fix the script."
		exit 1
	fi
fi

# Run evndisp
echo "$EVNDISPSYS/bin/evndisp -runnumber=$RUNNUM $OPT "
$EVNDISPSYS/bin/evndisp -runnumber=$RUNNUM $OPT 

exit
