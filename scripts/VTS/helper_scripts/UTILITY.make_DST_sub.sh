#!/bin/bash
# script to make DSTs

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RUN=RUNFILE
PED=PEDESTALS
SUMW=SUMWINDOW
LMULT=LLLOWGAIN

# temporary (scratch) directory
if [[ -n $TMPDIR ]]; then
    TEMPDIR=$TMPDIR/$RUN
else
    TEMPDIR="$VERITAS_USER_DATA_DIR/TMPDIR"
fi
echo "Temporary directory: $TEMPDIR"
mkdir -p $TEMPDIR

# output data files are written to this directory
ODIR=OUTPUTDIR
mkdir -p $ODIR

# output log files are written to this directory
LOGDIR="$ODIR"
mkdir -p $LOGDIR

# eventdisplay reconstruction parameter
ACUTS=RRRRPFILE

#########################################
# pedestal and tzero calculation. No tzeros needed for lmult. 
if [[ $PED == "1" ]]; then
    rm -f $LOGDIR/$RUN.ped.log
    $EVNDISPSYS/bin/evndisp -runnumber=$RUN -runmode=1 &> $LOGDIR/$RUN.ped.log

    if [[ $LMULT == "0" ]] ; then    
	rm -f $LOGDIR/$RUN.tzero.log
    	$EVNDISPSYS/bin/evndisp -runnumber=$RUN -runmode=7 -nocalibnoproblem &> $LOGDIR/$RUN.tzero.log
    fi
fi


#########################################
# Other options. hilo runs should have all events analysed.
if [[ $LMULT == "0" ]] ; then
	OPT=" -nevents=5000 "
else
	OPT=" -lowgaincalibrationfile calibrationlist.LowGainForCalibration.dat "
fi

#########################################
# run eventdisplay
rm -f $LOGDIR/$RUN.log
$EVNDISPSYS/bin/evndisp -runnumber=$RUN -runmode=4 -nocalibnoproblem $OPT -reconstructionparameter $ACUTS -dstfile $TEMPDIR/$RUN.DST.root &> $LOGDIR/$RUN.DST.log

# move data file from temp dir to data dir
cp -f -v $TEMPDIR/$RUN.DST.root $ODIR/$RUN.DST.root
rm -f $TEMPDIR/$RUN.DST.root

exit
