#!/bin/bash
# script to combine effective areas

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
EAFILES=INPUTFILES
OFILE=OUTPUTFILE
ODIR=OUTPUTDIR
mkdir -p $ODIR
chmod -R g+w $ODIR

# Write histograms?
WRITEHISTOS="false"

# keep a list of all input files for checks
rm -f $ODIR/$OFILE.list
ls -1 $EAFILES > $ODIR/$OFILE.list

# combine effective areas
$EVNDISPSYS/bin/combineEffectiveAreas "$EAFILES" $ODIR/$OFILE $WRITEHISTOS &> $ODIR/$OFILE.log 

# zip lot files
bzip2 $ODIR/$OFILE.combine.log

exit
