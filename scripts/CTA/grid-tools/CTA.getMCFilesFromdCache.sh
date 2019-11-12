#!/bin/bash
#
# simple script to cp eventdisplay files from dCache and rename according to the current naming sheme
#
#

if [ $# -ne 2 ]
then
   echo "./CTA.getEVNDISPFilesFromdCache.sh <list of files on dCache> <target directory>"
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
# copy from dcache
    echo "Source: $DFIL"
    echo "Dest: $OFIL"
    OFIL=`basename $DFIL`
    if [ ! -e $2/$OFIL ]
    then
       dccp $DFIL $2/$OFIL
    else
       echo "FILE EXISTS"
    fi
done

exit
