#!/bin/bash
# script to run energy3d on batch

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
SIMTYPE=SIMULATION
TELSETUP=TELESTEUP
SUMWINT=SUMERWINTER
NOICE=OLJUD
INDIR=INPUTDIR
OUTDIR=OUTPUTDIR

#################################
# temporary directory
if [[ -n "$TMPDIR" ]]; then 
    DDIR="$TMPDIR/energy3d_reference"
else
    DDIR="/tmp/energy3d_reference"
fi
mkdir -p $DDIR
echo "Temporary directory: $DDIR"

if [ $SIMTYPE = "CARE" ]; then
  POSTERUM="1200"
elif [ $SIMTYPE = "GRISU" ]; then
  POSTERUM="6500"
fi
OUTFILENAME="energy3d_referencetables_V$TELSETUP""_ATM$SUMWINT""_NOISE$NOICE""_$SIMTYPE"
for AFILE in 00 20 30 35 40 45 50 55 60 65; do
  SUBDIR="V$TELSETUP""_ATM$SUMWINT""_gamma/ze$AFILE""deg_offset0.5deg_NSB$NOICE""MHz"
  CURGETRNAME="9$TELSETUP""$POSTERUM.root"
  CURBATCHNAME="ze$AFILE""deg_offset0.5deg_NSB$NOICE""MHz_9$TELSETUP""$POSTERUM.root"
  if [ ! -e "$INDIR/$SUBDIR/$CURGETRNAME" ] ; then
    echo "Error, unable to find 3D-modeled file: '$INDIR/$SUBDIR/$CURGETRNAME', exiting..."
    exit 1
  fi
  cp "$INDIR/$SUBDIR/$CURGETRNAME" "$DDIR/$CURBATCHNAME"
  echo "moved $INDIR/$SUBDIR/$CURGETRNAME to $DDIR/$CURBATCHNAME"
done

###############################################
# run eventdisplay
###############################################
# run options
if [ ! -e "$EVNDISPSYS/bin/create_energy3d_referencetables" ] ; then
  echo "Error, unable to find the analysistool!"
  exit 1
fi
echo "$EVNDISPSYS/bin/create_energy3d_referencetables $SIMTYPE $TELSETUP $SUMWINT $NOICE $DDIR 2>&1 > $OUTDIR/$OUTFILENAME.log"
      $EVNDISPSYS/bin/create_energy3d_referencetables $SIMTYPE $TELSETUP $SUMWINT $NOICE $DDIR 2>&1 > $OUTDIR/$OUTFILENAME.log

# copying analysed temp files to outputfolder and removing the temp files
cp -f -v $DDIR/$OUTFILENAME.root $OUTDIR/$OUTFILENAME.root
chmod g+w $OUTDIR/$OUTFILENAME.root
chmod g+w $OUTDIR/$OUTFILENAME.log
for AFILE in 00 20 30 35 40 45 50 55 60 65; do
  rm -f -v $DDIR/V$TELSETUP""_ATM$SUMWINT""_gamma/ze$AFILE""deg_offset0.5deg_NSB$NOICE""MHz/9$TELSETUP""$POSTERUM.root
done

echo "energy3d output root file written to $OUTDIR/$OUTFILENAME.root"
echo "energy3d log file written to $OUTDIR/$OUTFILENAME.log"
exit
