#!/bin/bash
#
# script to link evndisp files from location to another
#
########################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ] && [ ! -n "$4" ]
then
     echo 
     echo "./CTA.LinkFilesFromProdSets.sh <sub array list (old)> <prefix> <AnalysisDirectory> <particle>"
     echo
     echo "  <sub array list>          text file with list of subarray IDs"
     echo
     exit
fi

SUBAR=$1
ANADIR=$2

###############################################################
# loop over all arrays
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $1`

for ARRAY in $VARRAY
do
   ARRAYN="N.$ARRAY"
   echo "SOURCE FIlES ARRAY $ARRAY  $ARRAYN"

   SDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$3/$ARRAYN/$4"
   DDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$3/$ARRAY/$4"

   FLIST=`ls $SDIR/*.root`
   for F in $FLIST
   do
      NF=`basename $F`
      NF=${NF//_/_G_}

      echo $NF
      ln -s $F $CTA_USER_DATA_DIR/analysis/AnalysisData/$3/$ARRAY/$4/$NF 

   done
done

exit

