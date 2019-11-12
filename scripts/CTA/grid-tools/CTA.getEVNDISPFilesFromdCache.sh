#!/bin/bash
#
# simple script to cp eventdisplay files from dCache and rename according to the current naming sheme
#
#

if [ $# -ne 2 ] && [ $# -ne 3 ]
then
   echo "./CTA.getEVNDISPFilesFromdCache.sh <list of files on dCache> <target directory> [mcaz]"
   echo
   exit
fi

if [ ! -e $2 ]
then
   mkdir -p $2
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# loop over all files in the list
FILEL=`cat $1`
for DFIL in $FILEL
do
    echo $DFIL
# copy from dcache
    echo "Source: $DFIL"
    echo "Dest: $OFIL"
    OFIL=`basename $DFIL`
    dccp $DFIL $2/$OFIL
# rename file
    if [ -n "$3" ]
    then
       EFIL=`basename $OFIL .root`
       EFIL=$EFIL"_G_"$3.root
    else
       MCAZ=`$EVNDISPSYS/bin/printRunParameter $2/$OFIL -mcaz`
       RUNN=`$EVNDISPSYS/bin/printRunParameter $2/$OFIL -runnumber`
       EFIL=$RUNN"_G_"$MCAZ"deg.root"
    fi
    mv $2/$OFIL $2/$EFIL
done

exit
