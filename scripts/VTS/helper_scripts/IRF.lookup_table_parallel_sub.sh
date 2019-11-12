#!/bin/bash
# script to fill a lookup tables for a given point in the parameter space
# (i.e. ze, wobble, etc)

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
ZA=ZENITHANGLE
WOBBLE=WOBBLEOFFSET
NOISE=SIMNOISE
EPOCH=ARRAYEPOCH
ATM=ATMOSPHERE
RECID=RECONSTRUCTIONID
SIMTYPE=SIMULATIONTYPE
INDIR=INPUTDIR
ODIR=OUTPUTDIR

TABFILE="table_${SIMTYPE}_${ZA}deg_${WOBBLE}wob_noise${NOISE}_${EPOCH}_ATM${ATM}_ID${RECID}"

MAXDIST="-distance_energyCuts=1.3"

# remove existing log and table file
rm -f "$ODIR/$TABFILE.root"
rm -f "$ODIR/$TABFILE.log"

# make the table part
$EVNDISPSYS/bin/mscw_energy -filltables=1 $MAXDIST -write1DHistograms -inputfile "$INDIR/*[0-9]*.root" -tablefile "$ODIR/$TABFILE.root" -ze=$ZA -arrayrecid=$RECID -woff=$WOBBLE &> "$ODIR/$TABFILE.log"

exit
