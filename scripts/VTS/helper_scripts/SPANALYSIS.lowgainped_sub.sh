#!/bin/bash
# script to analyse VTS raw files (VBF) with eventdisplay

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RUN=RUNFILE
CALDIR=OUTPUTDIRECTORY
TELTOANA=TELTOANACOMB
NEVENTS=NNNN
CALIBRATIONSUMFIRST=CALIBFIRST
CALIBRATIONSUMWINDOW=CALIBSUMWINDOW

# temporary (scratch) directory
if [[ -n $TMPDIR ]]; then
    TEMPDIR=$TMPDIR/$RUN
else
    TEMPDIR="$VERITAS_USER_DATA_DIR/TMPDIR"
fi
echo "Scratch dir: $TEMPDIR"
mkdir -p $TEMPDIR

LOGDIR="$TEMPDIR"


# eventdisplay reconstruction parameter
ACUTS="EVNDISP.reconstruction.LGCalibration.runparameter"

OPT=" -calibrationsumwindow=$CALIBRATIONSUMWINDOW -calibrationsumfirst=$CALIBRATIONSUMFIRST -nevents=$NEVENTS -calibrationdirectory $TEMPDIR -teltoana=$TELTOANA "

#########################################
# pedestal calculation
rm -f $LOGDIR/$RUN.ped.log
$EVNDISPSYS/bin/evndisp -runmode=6 -runnumber=$RUN -reconstructionparameter $ACUTS $OPT &> $LOGDIR/$RUN.ped.log
echo "$EVNDISPSYS/bin/evndisp -runmode=6 -runnumber=$RUN -reconstructionparameter $ACUTS $OPT "

for ((i=1; i<5; i++)) 
do
	mkdir -p $CALDIR/Calibration/Tel_$i/
	mv $TEMPDIR/Calibration/Tel_$i/$RUN.lped $CALDIR/Calibration/Tel_$i/
	mv $TEMPDIR/Calibration/Tel_$i/$RUN.lped.root $CALDIR/Calibration/Tel_$i/
	mv $LOGDIR/$RUN.ped.log $CALDIR/
done	

exit
