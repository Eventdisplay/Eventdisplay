#!/bin/bash
# read a calibration file produced by write_analysis_scripts.pl

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: calculate gains from a calibration file

SPANALYSIS.evndisp_laser_runs_from_calibfile.sh <calib file>

required parameters:

    <calib file>         calibration file produced by write_analysis_scripts.pl

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
CALIBFILE=$1

# locations of vbf files
DDIR="$VERITAS_DATA_DIR/data/"

FILES=`grep LASER $CALIBFILE | awk '{print $2"_"$3}'`

for AFILE in $FILES
do
    RUN=${AFILE:0:5}
    DTEL=${AFILE:6}
    echo $AFILE $RUN $DTEL
    DFILE=`find -L $DDIR -name "$RUN.cvbf"`

    if [[ -f $DFILE ]]; then
        $EVNDISPSYS/scripts/VTS/SPANALYSIS.evndisp_laser_run.sh $DTEL $DFILE
    else
	    echo "Missing laser/flasher file $DFILE, please download it"
    fi
done

exit
