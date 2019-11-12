#!/bin/bash
#
# simple script to download raw files from the GRID
# using DIRAC tools
# (adjusted to DESY environment)
#
#

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
   echo
   echo "CTA.getRawFilesFromGRID-DIRAC.sh <run list> <target directory> [max jobs]"
   echo
   echo "download using DIRAC tools"
   echo "(adjusted to DESY environment)"
   echo
   exit
fi

if [ -e $2 ]
then
   mkdir -p $2
fi

# temporary directory for file lists
FFN=`basename $1`
mkdir -p $2/tmplists

# loop over all files in the list
NTMPLIST=`wc -l $1 | awk '{print $1}'`
FILEN=499
for ((l = 1; l < $NTMPLIST; l+=$FILEN ))
do
# create file lists with $FILEN files each
   let "k = $l + $FILEN - 1"
   let "z = $z + 1"
   LLIST=$2/tmplists/$FFN.tmplist.d.$z.list
   echo $LLIST
   sed -n "$l,$k p" $1 > $LLIST
# run dirac get files
   dirac-dms-get-file $LLIST
done


exit
