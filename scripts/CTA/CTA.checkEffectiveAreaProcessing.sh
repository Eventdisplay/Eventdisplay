#!/bin/bash
#
# script to check effective area processing
#
#
#######################################################################

if [ $# -lt 3 ] 
then
   echo 
   echo "./CTA.checkEffectiveAreaProcessing.sh <sub array list> <directory with effective areas> <recid>"
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
# loop over all observing times
for T in 50h 5h 30m 100s
do
#   for N in NIM2MST4LST4SST4-d20170120-V2 NIM2MST3LST3SST3-d20170120-V2 NIM2MST2LST2SST2-d20170120-V2
   for N in NIM2LST4MST4SST4-d20170224-V2 NIM2LST3MST3SST3-d20170224-V2
   do
       EFFDIR=$DDIR/EffectiveArea-${T}-ID${RECID}-${N}

       echo "CHECKING $EFFDIR"

       # check if directory with effective areas exist
       if [ -d $EFFDIR ]
       then

            ###########################################################
            # loop over all arrays
            for ARRAY in $VARRAY
            do
                # skip LST checks for South threshold arrays
                if [ $RECID = "0" ] && [ $ARRAY = "*-TS-*" ] && [ $DDIR = "*paranal*" ]
                then
                    continue
                fi
                ###########################################################
                # loop over all effective area directories
#               for D in AngularResolution QualityCuts001CU BDT.$T-V2.d20170224
               for D in BDT.$T-V2.d20170224
                do
                    # number of particle types and off-axis
                    if [ $D = "AngularResolution" ]
                    then
                        PARTICLE=( "gamma_onSource" "gamma_cone" )
                        OFFAXISB=( "0" "5" )
                    else
                        PARTICLE=( "gamma_onSource" "electron" "proton" "gamma_cone" "electron_onSource" "proton_onSource" )
                        OFFAXISB=( "0" "5" "5" "5" "0" "0" )
                    fi
                    ###########################################################
                    # loop over all particle types
                    for ((i = 0; i < ${#PARTICLE[@]}; i++ ))
                    do
                        P=${PARTICLE[$i]}
                        ###########################################################
                        # loop over all off-axis bins
                        B=${OFFAXISB[$i]}
                        for ((b = 0; b <= $B; b++ ))
                        do
                            ###########################################################
                            # LOG file
                            EFFLOGFILE=${P}.${ARRAY}_ID${RECID}.eff-${b}.log
                            # ROOT FILE
                            EFFROOTFILE=${P}.${ARRAY}_ID${RECID}.eff-${b}.root

                            echo "CHECKING $EFFDIR/${D}/$EFFROOTFILE"

                            if [ -e $EFFDIR/${D}/${EFFLOGFILE}.bz2 ]
                            then
                               bunzip2 $EFFDIR/${D}/${EFFLOGFILE}.bz2

                               grep -i error $EFFDIR/${D}/${EFFLOGFILE}

                               bzip2 -f $EFFDIR/${D}/${EFFLOGFILE}
                            else
                               echo "ERROR: LOG FILE NOT FOUND $EFFDIR/${D}/${EFFLOGFILE}.bz2"
                            fi
                            if [ ! -e $EFFDIR/${D}/$EFFROOTFILE ]
                            then
                                echo "ERROR: EFFECTIVE AREAS ROOT FILE NOT FOUND: $EFFDIR/${D}/$EFFROOTFILE"
                            fi
                            # for quality cuts: check for Particle number file
                            if [ $D = "QualityCuts001CU" ]
                            then
                                if [ $P = "*onSource*" ]
                                then
                                    PARTFILE="$EFFDIR/${D}/ParticleNumbers.${ARRAY}.00.root"
                                else
                                    PARTFILE="$EFFDIR/${D}/ParticleNumbers.${ARRAY}.${b}.root"
                                fi
                                echo "CHECKING $PARTFILE"
                                if [ ! -e $PARTFILE ]
                                then
                                   echo "ERROR: PARTICLE NUMBER FILE NOT FOUND: $PARTFILE"
                                fi
                            fi
                         done
                     done
                 done
              done
            fi
        done
    done

    exit
