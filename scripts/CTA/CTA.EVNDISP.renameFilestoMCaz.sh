#!/bin/bash
#
# rename files according to the CORSIKA shower direction
#
#
#######################################################################

if [ $# -ne 1 ]
then
   echo "./CTA.EVNDISP.renameFilestoMCaz.sh <directory>" 
   echo
   echo ""
   exit
fi

DDIR=$1


############################################################################
# loop over particle types
VPART=( "proton" )
VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
NPART=${#VPART[@]}
for ((m = 0; m < $NPART; m++ ))
do
   PART=${VPART[$m]}

############################################################################
# loop over all files in files 
#   FILES=`ls -1 $DDIR/$PART/*[0-9].root`
   FILES=`find $DDIR/$PART/ -name "*[0-9].root"`
   for AFIL in $FILES
   do
     MCAZ=`$EVNDISPSYS/bin/printRunParameter $AFIL -mcaz`

     NEWN=`basename $AFIL .root`

     NEWN=$DDIR/$PART/$NEWN"_"$MCAZ"deg.root"

     mv -v $AFIL $NEWN

   done
done

exit
