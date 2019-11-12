#!/bin/bash
# script to analyse files with anasum

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# set up gammalib and ctools
source $EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.init_gammalib_and_ctools.sh

INFILE=INPUTANASUMROOTFILE
ODIR=OUTPUTFITSDIRECTORY
MAXLEN=MAXIMUMCHUNKLENGTH

echo "INFILE:'$INFILE'"
echo "ODIR:'$ODIR'"
echo "MAXLEN:'$MAXLEN'"

RUNNUM=$( basename $INFILE | grep -oP "^\d+" )

# temporary (scratch) directory
if [[ -n $TMPDIR ]]; then
    export TEMPDIR=$TMPDIR/$RUN
else
    export TEMPDIR="$VERITAS_USER_DATA_DIR/TMPDIR"
fi
echo "Scratch dir: $TEMPDIR"
mkdir -p $TEMPDIR
cd $TEMPDIR/

# copy the input file to our local temporary directory
if [ ! -e "$INFILE" ] ; then
  echo "Error, ANALYSIS.ctools_conversion_sub.sh could not locate input file \$INFILE='$INFILE', exiting..."
  exit 1
fi
INFILEBASE=$( basename $INFILE )
LOCALINPUTFILE="$TEMPDIR/$INFILEBASE"
cp $INFILE $LOCALINPUTFILE

# local output directory
LOCALOUTPUTDIR=$TEMPDIR/output
mkdir -p $LOCALOUTPUTDIR

# local intermediate directory
LOCALINTERMEDIATEDIR=$TEMPDIR/intermediate
mkdir -p $LOCALINTERMEDIATEDIR

# figure out which effective area file we need (with the -e option), then copy it to the local temporary directory
# THIS DOES NOT ACTUALLY DO THE CONVERSION!
# writeCTAEventListFromAnasum automatically looks locally first for the effective area file, 
# then looks in $VERITAS_EVNDISP_AUX_DIR
REQEFFFILENAME=$( $EVNDISPSYS/bin/writeCTAEventListFromAnasum -i $LOCALINPUTFILE -o $LOCALOUTPUTDIR -c "$MAXLEN" -e | grep -P "^requires" | awk '{ print $2 }' )
FULLEFFFILENAME="$VERITAS_EVNDISP_AUX_DIR/EffectiveAreas/$REQEFFFILENAME"
echo "To properly perform the conversion, we require the effective area file '$FULLEFFFILENAME'..."
if [ ! -e "$FULLEFFFILENAME" ] ; then
  echo "Error, ANALYSIS.ctools_conversion_sub.sh could not locate the needed effective area file \$FULLEFFFILENAME='$FULLEFFFILENAME', exiting..."
  exit 1
fi
LOCALEFF=$TEMPDIR/$REQEFFFILENAME
cp $FULLEFFFILENAME $LOCALEFF

# quick ls before we run the command
echo "ls -lh"
ls -lh
echo

# DO CONVERSION HERE
# store files in an intermediate directory
LOGBASE="VR$RUNNUM.ctools.log"
echo "COMMAND: '$EVNDISPSYS/bin/writeCTAEventListFromAnasum -i $LOCALINPUTFILE -o $LOCALINTERMEDIATEDIR -c $MAXLEN'"
$EVNDISPSYS/bin/writeCTAEventListFromAnasum -i $LOCALINPUTFILE -o $LOCALINTERMEDIATEDIR -c "$MAXLEN" > $LOCALOUTPUTDIR/$LOGBASE

# loop over intermediate fits files
# (they have bad backgrounds, we must replace them with better ones)
echo >> $LOCALOUTPUTDIR/$LOGBASE
echo "now replacing backgrounds with better ones" >> $LOCALOUTPUTDIR/$LOGBASE
INTERFITS=$( ls -1 $LOCALINTERMEDIATEDIR/*.fits )
for interfile in $INTERFITS ; do
  interbase=$( basename $interfile )
  
  echo "" >> $LOCALOUTPUTDIR/$LOGBASE
  echo "$VERIPY/scripts/background.improve.py : $interbase" >> $LOCALOUTPUTDIR/$LOGBASE

  # REPLACE BACKGROUND HERE
  echo "COMMAND: '$VERIPY/scripts/background.improve.py $LOCALINTERMEDIATEDIR/$interbase $LOCALOUTPUTDIR/$interbase 2>&1 >>  $LOCALOUTPUTDIR/$LOGBASE'"
  $VERIPY/scripts/background.improve.py $LOCALINTERMEDIATEDIR/$interbase $LOCALOUTPUTDIR/$interbase 2>&1 >>  $LOCALOUTPUTDIR/$LOGBASE

done

# copy chunklist from intermediate dir to the local output dir
cp $LOCALINTERMEDIATEDIR/*.chunklist $LOCALOUTPUTDIR/

# another ls afterwards
echo ; echo "ls -lh"
ls -lh
echo ; echo "ls -lh $LOCALOUTPUTDIR/"
ls -lh $LOCALOUTPUTDIR
echo

# copy output files to our output dir
cp $LOCALOUTPUTDIR/* $ODIR/
echo "RUN $RUNNUM CTOOLSLOG $ODIR/$LOGBASE"
echo "RUN $RUNNUM CTOOLSMETA $ODIR/VR$RUNNUM.chunklist"

exit
