#!/bin/bash
#
# simple script to get alias file names from GRID
#
#

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ] && [ ! -n "$4" ]
then
   echo "./CTA.getAliasesFromGrid.sh <run list> <data directory> <glite command (0=no/1=yes)> <output list file directory>"
   echo
   exit
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# output list file directory
mkdir -p $4
TLIS=`basename $1`
FILEEX="$4/$TLIS.exist"
rm -f $FILEEX
touch $FILEEX
FILELO="$4/$TLIS.dcache"
rm -f $FILELO
touch $FILELO
FILEGL="$4/$TLIS.glite"
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
       rm -f $2/$OFIL
       DC="/acs/grid/cta/$i"
       if [ -e $DC ]
       then
          echo "dccp $DC $2/$OFIL" >> $FILELO
##################################
# prod2 used both DFC and LFC
#       else
#          surl=`lcg-lr lfn:/grid$i`
#          if [ $3 -ne "0" ]
#          then
#             echo "$surl srm://styx.ifh.de:8443/srm/v2/server?SFN=$2/${OFIL}" >> $FILEGL
#          fi
##################################
# prod3 uses DFC only
       else
           echo "DFC check"
           surl=`dirac-dms-lfn-replicas  $i |  grep cta | grep -v dips | awk '{print $3}'`
           echo "DFC $surl"
           if [ $3 -ne "0" ] && [ -n $surl ]
           then
               echo "$surl srm://styx.ifh.de:8443/srm/v2/server?SFN=$2/${OFIL}" >> $FILEGL
           fi
       fi
    fi
done

rm -f $4/$TLIS.glite.tt
sed '/dips/d' $4/$TLIS.glite > $4/$TLIS.glite.tt
mv -f $4/$TLIS.glite.tt $4/$TLIS.glite


# make dCache list exectuable
chmod u+x $4/$TLIS.dcache
# sort glite file lists
grep ln2p3 $4/$TLIS.glite | grep -v dips > $4/$TLIS.glite.ln2p3
grep cnaf $4/$TLIS.glite | grep -v dips  > $4/$TLIS.glite.cnaf

exit
