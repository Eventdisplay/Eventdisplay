#!/bin/sh
#
# script to unpack GRID produced EVNDISP tar files into a analysis directory
# hardwired number of arrays
#
# e.g. CTA.prod3b.unpackEVNDISPFiles.sh sorted.diff.6.list prod3-paranalp05-40deg-NN 
#
#############################################################################

if [ $# -ne 2 ]
then
     echo
     echo "./CTA.prod3.unpackEVNDISPFiles.shCTA.prod3.unpackEVNDISPFiles.sh <list of tar files> <site>"
     echo
     exit
fi

FIL="$1"
SITE="$2"

####################################
# files will unpacked in this directory
DDIR="/lustre/fs9/group/cta/users/maierg/CTA/analysis/AnalysisData/$SITE/"
mkdir -p $DDIR

SARRAY="3HB1-2 3HB2-2 3HB4-2 3HB89 3HB8 3HB9"
TELTYP="FG FD NG ND FA NA"
PART="gamma_onSource gamma_cone electron proton"
PART="gamma_cone"


for P in ${PART}
do
    # make sub array directories
    for V in $SARRAY
    do
        for T in $TELTYP
        do
            SDIR="S.${V}-${T}"
            mkdir -p $DDIR/$SDIR/${P}
        done
    done

    # get list of files for this particle type
    if [[ $P == "gamma_onSource" ]]
    then
        FILES=`cat $FIL | grep gamma | grep -v diffuse`
    elif [[ $P == "gamma_cone" ]]
    then
        FILES=`cat $FIL | grep gamma | grep diffuse`
    else
 
        FILES=`cat $FIL | grep $P`
    fi

    echo "unpacking"
    echo "========="
    for DF in $FILES
    do
        A="/acs/grid/cta${DF}"
        # check only one array layout
#        B=`basename $A`
#        AZ=`echo $B| cut -d'-' -f 2`
        # get RUN
#        RUN=${B%%-*}

         echo $A
         dccp $A - | tar -xzf -

        # mv files into subarray directories
          for VA in $SARRAY
          do
            for TA in $TELTYP
            do
                SDIR="S.${VA}-${TA}"

                DAT="*${VA}-${TA}*root"

                echo $DDIR/$SDIR/${P}/

                mv -v $DAT $DDIR/$SDIR/${P}/
            done
          done
    done
done
