#!/bin/bash
#
# script to check TMVA processing
#
#
#######################################################################

if [ $# -lt 3 ] 
then
   echo 
   echo "./CTA.checkTMVAProcessing.sh <sub array list> <directory with effective areas> <recid>"
   echo 
   echo "  write particles files needed for TMVA cut optimization"
   echo
   echo "  <sub array list>          text file with list of subarray IDs"
   echo
   echo "  <directory with effective areas>  (full) path to effective areas"
   echo
   echo "  <recid>                   reconstruction ID from mscw stage" 
   echo 
   exit
fi

SUBAR=$1
DDIR=$2
RECID=$3

###############################################################
# list of arrays
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $1`

###########################################################
# loop over all wobble offsets
for T in 0.5 1.5 2.5 3.5 4.5 5.5
do
   EFFDIR=$DDIR/EffectiveArea-${T}-ID0-${EFFDIRNAME}/
   TMVADIR=$DDIR/

   for N in NIM2MST2LST2SST2-d20161128 NIM2MST3LST3SST3-d20161128 NIM2MST4LST4SST4-d20161128
   do
            ###########################################################
            # loop over all arrays
            for ARRAY in $VARRAY
            do
                ###########################################################
                # loop over all azimuth bins
                for A in ID${RECID} ID${RECID}_0deg ID${RECID}_180deg
                do
                    TMVADIR=$DDIR/${ARRAY}/TMVA/BDT-V2-${A}-${N}-${T}

                    echo $TMVADIR
                    ###########################################################
                    # check XML files
                    NX=0
                    for X in 0 1 2 3 4 5 6 7 8
                    do
                       XF="BDT_${X}_BDT_0.weights.xml"
                       if [ -e $TMVADIR/$XF ]
                       then
                          let "NX = NX + 1"
                       fi
                    done
                    echo "    found $NX XML files"
                    # require that files with index 3, 4, and 5 are available
                    for X in 3 4 5
                    do
                       XF="BDT_${X}_BDT_0.weights.xml"
                       if [ ! -e $TMVADIR/$XF ]
                       then
                          let "NX = NX + 1"
                          echo "    error: missing $XF"
                       fi
                    done
                done
           done
    done
done

exit
