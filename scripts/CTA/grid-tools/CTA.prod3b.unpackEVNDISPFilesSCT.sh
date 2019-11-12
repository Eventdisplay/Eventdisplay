#!/bin/sh
#
# script to unpack GRID produced EVNDISP tar files into a analysis directory
# hardwired number of arrays
#
# e.g. CTA.prod3b.unpackEVNDISPFiles.sh sorted.diff.6.list prod3-paranalp05-40deg-NN 
#
#############################################################################

if [ $# -ne 1 ]
then
     echo
     echo "./CTA.prod3.unpackEVNDISPFiles.shCTA.prod3.unpackEVNDISPFiles.sh <list of tar files>"
     echo
     exit
fi

FIL="$1"


SARRAY="3HB9"
PART="gamma_onSource gamma_cone electron proton"
PART="gamma_cone"


for P in ${PART}
do
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
         dccp $A - | tar -xvf -

    done
done
