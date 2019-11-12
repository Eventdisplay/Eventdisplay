#!/bin/bash
#
# simple script to download raw files from the GRID
# (adjusted to DESY environment)
#
#

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
   echo
   echo "CTA.getRawFilesFromGRID-glite.sh <run list> <target directory> [max jobs]"
   echo
   echo "download using glite-transfer-submit"
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
   let "k = $l + $FILEN"
   let "z = $z + 1"
   LLIST=$2/tmplists/$FFN.tmplist.d.$z.list
   echo $LLIST
   sed -n "$l,$k p" $1 > $LLIST
# run glite-transfer-submit
  ftsid=`glite-transfer-submit -s https://fts3-kit.gridka.de:8443/glite-data-transfer-fts/services/FileTransfer -f $2/tmplists/$FFN.tmplist.d.$z.list`
   echo "   glite ID $ftsid"
done

echo
echo "check transfer status with:"
echo "glite-transfer-status -s https://fts3-kit.gridka.de:8443/glite-data-transfer-fts/services/FileTransfer -l <glite ID>"

exit
