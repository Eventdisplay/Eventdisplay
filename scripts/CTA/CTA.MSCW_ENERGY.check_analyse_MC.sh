#!/bin/bash
#
# script to check that files where processed correctly (table analysis)
#
# Revision $Id$
#
#
#######################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ]
then
   echo
   echo "./CTA.MSCW.check_convert_and_analyse_MC_VDST.sh <subarray list> <data set>"
   echo
   echo "simple script to check that table analysis was running properly"
   echo 
   echo "  <sub array list>          text file with list of subarray IDs"
   echo
   echo "  <data set>         e.g. cta-ultra3, ISDC3700m, ...  "
   echo
   exit
fi

#arrays
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $1`

#particles
PART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
NPART=${#PART[@]}

############################################################################

for ARRAY in $VARRAY
do
   echo "----------------------------------------------------"
   echo "ARRAY $ARRAY"

   for (( P = 0; P < $NPART; P++ ))
   do
      PRIM=${PART[$P]}

      DDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$2/$ARRAY/Analysis/"

 # check number of root files
      NFIL=`ls $DDIR/$PRIM* | wc -l`
      echo "$NFIL $PRIM: number of files"
      NDISK=`du -h -c -s $DDIR/$PRIM* | tail -n 1`
      echo "$NDISK $PRIM: disk usage"

# check for zero length files
#      if [ ! -s $DDIR/$RUN.root ]
#      then
#        echo "ZERO LENGTH EVNDISP FILE: $DDIR/$RUN.root"
#        echo $AFIL >> $FAILED.$ARRAY.list
#        continue
#      fi
  done
   
done

exit
