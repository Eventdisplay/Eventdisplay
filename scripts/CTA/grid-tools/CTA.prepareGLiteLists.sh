#!/bin/bash


if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ]
then
    echo "./CTA.prepareGLiteLists.sh <run list> <data directory> <output list file directory>"
    exit
fi


# output list file directory
mkdir -p $3
TLIS=`basename $1`

#############################
# loop over all storage elements
for S in CNAF IN2P3 LAPP CYF-STORM DESY-ZN
do
    if [[ $S != "DESY-ZN" ]]
    then
        FILEGL="$3/$TLIS.glite.${S}"
    else
        FILEGL="$3/$TLIS.dCache.${S}"
    fi
    rm -f $FILEGL
    touch $FILEGL
    # loop over all files in the list
    FILEL=`cat $1 | grep simtel | grep $S` 

    while read i; do
        if [[ $i == *"simtel"* ]] && [[ $i == *"$S"* ]]
        then
            FILEN=`echo $i | awk '{print $1}'`
            FILEF=`echo $i | awk '{print $3}'`

            OFIL=`basename $FILEN`

            if [[ $S != "DESY-ZN" ]]
            then
                # check if file is on disk
                if [ -e $2/${OFIL} ]
                then
                    continue
                fi
                echo "$FILEF srm://styx.ifh.de:8443/srm/v2/server?SFN=$2/${OFIL}" >> $FILEGL
            else
                echo "/acs/grid/cta$FILEN" >> $FILEGL
            fi
        fi

    done <$1

    # remove files with zero file size
    if [ ! -s $FILEGL ]
    then
        rm -f $FILEGL
    fi

done
##################################
