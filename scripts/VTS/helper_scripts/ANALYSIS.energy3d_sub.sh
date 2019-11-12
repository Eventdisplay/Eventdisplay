#!/bin/bash
# script to run energy3d on batch

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RSFIL=REALSIM
RFIL=RNUM
NFIL=NOICELEVELS
DFIL=INDATADIR
REFDIR=RDIR
OUTDIR=OUTPUTDIR
TVFIL=TELV
SWFIL=SUMW
CGFIL=CAREGRISU

#################################

#tempout=/afs/ifh.de/group/that/work-cta/VERITAS/EVNDISP/trunk/jklö.txt
#echo "jfdkslfjsklfjkl" > $tempout

# temporary directory
if [[ -n "$TMPDIR" ]]; then 
    DDIR="$TMPDIR/Energy3D"
else
    DDIR="/tmp/Energy3D"
fi
mkdir -p $DDIR
echo "Temporary directory: $DDIR" 

# get biascorrection file  
BIASCORECTIONFILE=energy3d_biascorrections.txt
if [ ! -e "$REFDIR/$BIASCORECTIONFILE" ] ; then
   echo "Error, unable to find biascorrection file: '$REFDIR/$BIASCORECTIONFILE', exiting..."
   exit 1
fi
cp "$REFDIR/$BIASCORECTIONFILE" $DDIR/$BIASCORECTIONFILE
echo "moved $REFDIR/$BIASCORECTIONFILE to $DDIR/"
BIASCORR=$DDIR/$BIASCORECTIONFILE

# get referencetable
REFFILE="energy3d_referencetables_V${TVFIL}_ATM${SWFIL}_NOISE${NFIL}_${CGFIL}.root"
if [ ! -e "$REFDIR/$REFFILE" ] ; then
   echo "Error, unable to find Reference file: '$REFDIR/$REFFILE', exiting..."
   exit 1
fi
cp $REFDIR/$REFFILE $DDIR
echo "moved $REFDIR/$REFFILE to $DDIR"

# create RUNLIST in temporary directory
TARGETRUN=$DDIR/targetfile.txt

THEREST="$TVFIL $SWFIL $CGFIL"

if [ $RSFIL == "REAL" ] ; then

  # file the tempdir runlist
  echo $RFIL > $TARGETRUN
  echo 
  echo "$RFIL > $TARGETRUN"
  echo 

  # get the root file for analysis
  if [ ! -e "$VERITAS_USER_DATA_DIR/analysis/$DFIL/$RFIL.root" ] ; then
     echo "Error, unable to find 3D-modeled file: '$VERITAS_USER_DATA_DIR/analysis/$DFIL/$RFIL.root', exiting..."
     exit 1
  fi
  cp "$VERITAS_USER_DATA_DIR/analysis/$DFIL/$RFIL.root" $DDIR/
  echo "moved $VERITAS_USER_DATA_DIR/analysis/$DFIL/$RFIL.root to $DDIR/"

  THEFILE=$RFIL

else #if simulation

  # file the tempdir runlist
  WORKFILE=${RFIL}_${NFIL}_${DFIL}_${SWFIL}
  echo $WORKFILE > $TARGETRUN
  echo 
  echo "$WORKFILE > $TARGETRUN"
  echo 
  
  # get the root file for analysis
  if [ $CGFIL == "CARE" ]; then
     DATADIR="$VERITAS_USER_DATA_DIR/analysis/CARE/v500-dev08/CARE_June1425/V${TVFIL}_ATM${SWFIL}_gamma/ze${DFIL}deg_offset0.5deg_NSB${NFIL}MHz"
  else
     DATADIR="$VERITAS_USER_DATA_DIR/analysis/v444/GRISU/V${TVFIL}_ATM${SWFIL}_gamma/ze${DFIL}deg_offset0.5deg_NSB${NFIL}MHz"
  fi
  if [ ! -e "$DATADIR/$RFIL.root" ] ; then
    echo "Error, unable to find 3D-modeled file: '$DATADIR/$RFIL.root', exiting..."
    exit 1
  fi
  cp $DATADIR/$RFIL.root $DDIR/$WORKFILE.root
  echo "moved $DATADIR/$RFIL.root to $DDIR/$WORKFILE.root"

  THEFILE=$WORKFILE
fi

###############################################
# run eventdisplay
###############################################
# run options
if [ ! -e "$EVNDISPSYS/bin/energy3d" ] ; then
  echo "Error, unable to find the analysistool!"
  exit 1
fi
echo "$EVNDISPSYS/bin/energy3d TARGETRUN $TARGETRUN DDIR $DDIR THEREST $THEREST BIASCORR $BIASCORR 2>&1 > $OUTDIR/$THEFILE.energy3d.log"
      $EVNDISPSYS/bin/energy3d $TARGETRUN $DDIR $THEREST $BIASCORR 2>&1 > $OUTDIR/$THEFILE.energy3d.log

# copying analysed temp files to outputfolder and removing the temp files
cp -f -v $DDIR/$THEFILE.root $OUTDIR/$THEFILE.energy3d.root
chmod g+w $OUTDIR/$THEFILE.energy3d.root
chmod g+w $OUTDIR/$THEFILE.energy3d.log
rm -f -v $DDIR/$THEFILE.root
rm -f -v $DDIR/$REFFILE
rm -f -v $DDIR/$BIASCORECTIONFILE

echo "energy3d output root file written to $OUTDIR/$THEFILE.energy3d.root"
echo "energy3d log file written to $OUTDIR/$THEFILE.energy3d.log"
exit







