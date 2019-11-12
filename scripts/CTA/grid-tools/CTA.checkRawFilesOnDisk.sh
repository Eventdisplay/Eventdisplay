#!/bin/bash
#
# simple script to check if files are on disk or on the dCache
#
#

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ] 
then
   echo "./CTA.checkRawFilesOnDisk.sh <run list> <data directory> <output list file directory>"
   echo
   exit
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# output list file directory
mkdir -p $3
TLIS=`basename $1`
FILEEX="$3/$TLIS.exist"
rm -f $FILEEX
touch $FILEEX
FILELO="$3/$TLIS.dcache"
rm -f $FILELO
touch $FILELO
FILEGL="$3/$TLIS.noondisk"
rm -f $FILEGL
touch $FILEGL

# loop over all files in the list
FILEL=`cat $1`
for i in $FILEL
do
    echo "FILE $i"
    OFIL=`basename $i`
    echo "CHECK $2/$OFIL"
# check if file is already in the target directory
#    if [ -e $2/$OFIL ] && [ -s $2/$OFIL ]
    if [ -e $2/$OFIL ]
    then
       echo "FILE EXISTS: $2/$OFIL" >> $FILEEX
# check if the file is on the local dCache
    else
       # rm file (potentially zero size)
#       rm -f $2/$OFIL
       DC="/acs/grid/cta/$i"
       if [ -e $DC ]
       then
          echo "dccp $DC $2/$OFIL" >> $FILELO
       else
         echo $i >> $FILEGL
       fi
    fi
done


exit
