#!/bin/bash
# script to analyse VTS raw files (VBF) with eventdisplay

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RUN=SUBRUN
MERGELIST=SUBMERGELIST
OUTPUTFILE=SUBOUTPUTFILE
ODIR=SUBOUTDIR
LOGDIR="$ODIR"

# temporary (scratch) directory
if [[ -n $TMPDIR ]]; then
    TEMPDIR=$TMPDIR/$RUN
else
    TEMPDIR="$VERITAS_USER_DATA_DIR/TMPDIR"
fi
echo "Scratch dir: $TEMPDIR"
mkdir -p $TEMPDIR

# copy needed files to scratch dir
SCRATCHMERGELIST=$TEMPDIR/scratchmerge
MERGEFILES=$( cat $MERGELIST ) 
for file in ${MERGEFILES[@]} ; do
	cp $file $TMPDIR
	echo "$TMPDIR/`basename $file`" >> "$SCRATCHMERGELIST"
done

SCRATCHOUT=$TEMPDIR/merged.root

echo
echo "just before we merge, here's whats in $TEMPDIR:"
ls -l $TEMPDIR/
echo

$EVNDISPSYS/bin/frogsMergeDatafile $SCRATCHMERGELIST $SCRATCHOUT 2>&1

echo
echo "just after we merge, here's whats in $TEMPDIR:"
ls -l $TEMPDIR/
echo

rm -rf $OUTPUTFILE
cp -f $SCRATCHOUT $OUTPUTFILE

exit
