#!/bin/bash
#
# simple script to download raw files from the GRID
#
#

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
   echo "CTA.getRawFilesFromGRID.sh <run list> <target directory>"
   echo
   echo "(checks as well if a file is locally on the dCache SE)"
   echo
   exit
fi

if [ -e $2 ]
then
   mkdir -p $2
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# loop over all files in the list
FILEL=`cat $1`
for i in $FILEL
do
    OFIL=`basename $i`
    if [ -e $2/$OFIL ] && [ -s $2/$OFIL ]
    then
       echo "FILE EXISTS: $2/$OFIL"
    else
       rm -f $2/$OFIL
# first check if it is stored locally on the dcache
       DC="/acs/grid/cta/$i"
       echo $DC
       if [ -e $DC ]
       then
          echo "STORED LOCALLY: $i"
#	  dccp $DC $2/$OFIL
#       else
#          echo "COPY FROM GRID: $i"
#          lcg-cp -v lfn:/grid$i file:$2/$OFIL
       fi
       sleep 10
    fi
done

exit
