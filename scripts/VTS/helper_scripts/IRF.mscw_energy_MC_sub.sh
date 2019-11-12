#!/bin/bash
# script to analyse MC files with lookup tables

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
INDIR=INPUTDIR
ODIR=OUTPUTDIR
TABFILE=TABLEFILE
ZA=ZENITHANGLE
NOISE=NOISELEVEL
WOBBLE=WOBBLEOFFSET
NROOTFILES=NFILES
ANAMETHOD=ANALYSISMETHOD
RECID="RECONSTRUCTIONID"

# output directory
OSUBDIR="$ODIR/MSCW_RECID$RECID"
mkdir -p $OSUBDIR
chmod g+w $OSUBDIR
echo "Output directory for data products: " $OSUBDIR

# file names
OFILE="${ZA}deg_${WOBBLE}wob_NOISE${NOISE}"

# temporary directory
if [[ -n "$TMPDIR" ]]; then 
    DDIR="$TMPDIR/MSCW_${ZA}deg_${WOBBLE}deg_NOISE${NOISE}_ID${RECID}"
else
    DDIR="/tmp/MSCW_${ZA}deg_${WOBBLE}deg_NOISE${NOISE}_ID${RECID}"
fi
mkdir -p $DDIR
echo "Temporary directory: $DDIR"

# mscw_energy command line options
MOPT="-noNoTrigger -nomctree -writeReconstructedEventsOnly=1 -arrayrecid=$RECID -tablefile $TABFILE"
MOPT="$MOPT -useMedian=0 -distance_energyCuts=1.3"
echo "MSCW options: $MOPT"

if [[ $ANAMETHOD == "MODEL3D" ]]; then
    MODEL3D="-model3d"
    echo "using Model3D for direction and core values"
fi


# run mscw_energy
if [[ $NROOTFILES == 1 ]]; then
      rm -f $OSUBDIR/$OFILE.log
      inputfilename="$INDIR/*[0-9]*.root"
      outputfilename="$DDIR/$OFILE.mscw.root"
      logfile="$OSUBDIR/$OFILE.mscw.log"
elif	[[ $NROOTFILES > 1 ]]; then
      ITER=$((SGE_TASK_ID - 1))
      rm -f $OSUBDIR/${OFILE}_$ITER.log
      inputfilename="$INDIR/*[0-9]_$ITER.root"
      outputfilename="$DDIR/${OFILE}_$ITER.mscw.root"
      logfile="$OSUBDIR/${OFILE}_$ITER.mscw.log"
fi
$EVNDISPSYS/bin/mscw_energy $MOPT $MODEL3D -inputfile $inputfilename -outputfile $outputfilename -noise=$NOISE &> $logfile

# cp results file back to data directory and clean up
outputbasename=$( basename $outputfilename )
cp -f -v $outputfilename $OSUBDIR/$outputbasename
rm -f $outputfilename
rmdir $DDIR
chmod g+w $OSUBDIR/$outputbasename
chmod g+w $logfile

exit
