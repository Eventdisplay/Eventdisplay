#!/bin/bash
#
# script to check which files are available on local dCache / lustre
#
#

if [ $# -ne 3 ]
then
   echo "CTA.fill_runLists_mergedArrays.sh <run list> <data directory> <output list file directory>"
   echo
   echo "script to check which files from a cta-md-query produced run list are available on the local dCache and on local disk"
   echo
   echo "output files: file lists *.dcache.list and *.lustre.list"
   echo
   exit
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

# output list file directory
mkdir -p $3
TLIS=`basename $1`
FILELO="$3/$TLIS"
rm -f $FILELO.*
touch $FILELO.proton.list
touch $FILELO.electron.list
touch $FILELO.gamma_onSource.list
touch $FILELO.gamma_cone.list

touch $FILELO.proton.list.SCSST
touch $FILELO.electron.list.SCSST
touch $FILELO.gamma_onSource.list.SCSST
touch $FILELO.gamma_cone.list.SCSST

touch $FILELO.STD.missing.list
touch $FILELO.SCSST.missing.list

IDFIL1="STD"
IDFIL2="SCSST"
###########################################
# loop over all files in the list
FILEL=`cat $1`
for i in $FILEL
do
   OFIL=`basename $i`

   # these two variables are filled in case of sucess
   FILESTD=""
   FILESCSST=""
###########################################
# check if the files are on the local dCache or on local lustre

   # STD file
   DCSTD="/acs/grid/cta/$i"
   LCFSTD=$2/STD/$OFIL
   if [ -e "$DCSTD" ]
   then
      FILESTD="$DCSTD"
   elif [[ -e "$LCFSTD" ]]
   then
      FILESTD="$LCFSTD"
   fi
   # SCSST
   DCSST=${DCSTD/$IDFIL1/$IDFIL2}
   LCFSST=$2/SCSST/$OFIL
   if [ -e "$DCSST" ]
   then
      FILESCSST="$DCSST"
   elif [[ -e "$LCFSST" ]]
   then
      FILESCSST="$LCFSST"
   fi
   # fill lists according to particle types
   if [ -n "$FILESTD" -a -n "$FILESCSST" ]
   then
      if [[ $FILESTD == *"proton"* ]]
      then
         echo "$FILESTD" >> $FILELO.proton.list
         echo "$FILESCSST" >> $FILELO.proton.list.SCSST
      elif [[ $FILESTD == *"electron"* ]]
      then
         echo "$FILESTD" >> $FILELO.electron.list
         echo "$FILESCSST" >> $FILELO.electron.list.SCSST
      elif [[ $FILESTD == *"gamma_ptsrc"* ]]
      then
         echo "$FILESTD" >> $FILELO.gamma_onSource.list
         echo "$FILESCSST" >> $FILELO.gamma_onSource.list.SCSST
      elif [[ $FILESTD == *"gamma"* ]]
      then
         echo "$FILESTD" >> $FILELO.gamma_cone.list
         echo "$FILESCSST" >> $FILELO.gamma_cone.list.SCSST
      fi
   # keep track of missing files
   elif [ -z "$FILESTD" ]
   then
      echo "MISSING $DCSTD" >> $FILELO.STD.missing.list
      echo "MISSING $LCFSTD" >> $FILELO.STD.missing.list
   elif [ -z "$FILESCSST" ]
   then
      echo "MISSING $DCSST" >> $FILELO.SCSST.missing.list
      echo "MISSING $LCFSST" >> $FILELO.SCSST.missing.list
   fi
done

exit
