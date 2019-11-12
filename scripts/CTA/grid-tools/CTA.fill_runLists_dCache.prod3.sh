#!/bin/bash
#
# script to check which files are available on local dCache
#
#

if [ ! -n "$1" ] && [ ! -n "$2" ] 
then
   echo "./CTA.fill_runLists_dCache.prod3.sh <run list> <output list file directory> <output list file name>"
   echo
   echo "script to check which files from a cta-md-query produced run list are available on the local dCache (prod3 version)"
   echo
   echo "output files: file lists *.dcache.list and *.gridSE.list"
   exit
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# output list file directory
mkdir -p $2
TLIS=$3
FILELO="$2/$TLIS"
#rm -f $FILELO.*
#touch $FILELO.dcache
rm -f $FILELO.*.list

# loop over all files in the list
# (ignore all SCT files)
FILEL=`cat $1 | grep merged_desert`
for i in $FILEL
do
   OFIL=`basename $i`
# check if the file is on the local dCache
   DC="/acs/grid/cta/$i"
   if [ -e $DC ]
   then
      echo $DC >> $FILELO.dcache
   else
      echo $i >> $FILELO.gridSE
   fi
done

# sort according to particle types and subarrays
for i in 1 2 3 4 5
do
    grep proton $FILELO.dcache | grep subarray-${i} > $FILELO.proton_20deg-${i}.list
    grep electron $FILELO.dcache | grep subarray-${i} > $FILELO.electron_20deg-${i}.list
    grep gamma $FILELO.dcache | grep -v cone | grep subarray-${i} > $FILELO.gamma_onSource_20deg-${i}.list
    grep gamma $FILELO.dcache | grep cone | grep subarray-${i} > $FILELO.gamma_cone_20deg-${i}.list

    grep proton $FILELO.gridSE | grep subarray-${i} > $FILELO.proton_20deg-${i}.gridSE.list
    grep electron $FILELO.gridSE | grep subarray-${i} > $FILELO.electron_20deg-${i}.gridSE.list
    grep gamma $FILELO.gridSE | grep -v cone | grep subarray-${i} > $FILELO.gamma_onSource_20deg-${i}.gridSE.list
    grep gamma $FILELO.gridSE | grep cone | grep subarray-${i} > $FILELO.gamma_cone_20deg-${i}.gridSE.list
done

exit
