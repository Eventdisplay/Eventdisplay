#!/bin/bash
#
# script to write particle rate files
#
#
#######################################################################

AXRRAY=ARRAY
# directory with effective areas (and directory were files are written to)
DXDIR=DDIR
AXDIR=ADIR
RECID=RRRR
# should be either onSource or cone
OXFFSET=OFFSET

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

LLOG=$DXDIR/ParticleNumbers.$AXRRAY.$RECID.$OXFFSET.log
rm -f $LLOG

echo "$EVNDISPSYS/bin/writeParticleRateFilesFromEffectiveAreas  $AXRRAY $OXFFSET $RECID $DXDIR $AXDIR"
$EVNDISPSYS/bin/writeParticleRateFilesFromEffectiveAreas  $AXRRAY $OXFFSET $RECID $DXDIR $AXDIR > $LLOG 


############################################################################

exit
