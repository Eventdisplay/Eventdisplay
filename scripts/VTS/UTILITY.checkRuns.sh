#!/bin/bash
#
# UTILITY.checkRuns.sh
# Small script that quickly check the ED and mscw log files for errors and warnings.
# Also prints out the ped values - quick check for missing flasher info.
#
# Gareth
#
# NOTE: THIS SCRIPT DOESN'T REALLY WORK YET (SG)
# 140528: 	I made it working with some adaptions to the ANALYSIS.mscw_energy.sh
#			input style. Moritz


if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
    echo "UTILITY.checkRuns.sh <runlist> [ED version] [RecID] [MSCW output directory]"
    exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RUNLIST=$1
[[ "$2" ]] && EDVERSION=$2 || EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`



OUTPUTDIR=${VERITAS_USER_DATA_DIR}/analysis/Results/EVD-${EDVERSION}

if [[ "$4" ]]; then
    MSCWDIR=$4
elif [[ "$3" ]]; then
	ID=$3
    # user passed a Rec ID value
    MSCWDIR=$OUTPUTDIR/RecID$ID
else
    MSCWDIR=$OUTPUTDIR/MSCW-output
fi

FILE=`cat $RUNLIST`
for AFIL in $FILE
do
    LOGFILE="$OUTPUTDIR/${AFIL}.log"
    MSCWFILE="$MSCWDIR/${AFIL}.mscw.log"

    if [ -n "$3" ]; then
    	echo "run" ${AFIL} ":"
        echo $MSCWFILE
        grep -iH error $MSCWFILE 
        grep -iH warn  $MSCWFILE
        grep -H ped $MSCWFILE
        echo
    else
    	echo "run" ${AFIL} ":"
        grep -iH exit $LOGFILE
        grep -iH error $LOGFILE
        #grep -iH warn  $LOGFILE
        grep -iH 'Final checks on result file (seems to be OK)' $LOGFILE
        echo
    fi
done

exit
