#!/bin/bash
#
# simple script to check if a given file has been 
# processed on the grid and the result have been
# uploaded to a dCache SE
#

if [ $# -lt 5 ]
then
   echo
   echo "CTA.check_forProducts_on_dCache_and_Copy.sh <sub array list> <run list> <data directory> <destination dir> <file with missing runs> [MCaz]"
   echo
   echo " (note: very dependent on file names)"
   echo
   exit
fi

# dcache client
export DCACHE_CLIENT_ACTIVE=1

#########################################
# loop over all files in the list
FILEL=`cat $2`
for i in $FILEL
do
    SFIL=`basename $i .simtel.gz`
    RUN=`echo ${SFIL:(-6)}`
    RUN=$(echo $RUN | sed 's/^0*//')

    SUBA=`cat $1`
    for a in $SUBA
    do
      EFIL=$3/$RUN"_"$a"_evndisp.root"
# check if file exists
# (if not: add this one to filling list of jobs)
      if [ ! -e $EFIL ] 
      then
         echo "MISSING $EFIL"
         touch $5
         echo $i >> $5
         break
# file exists -> cp it to destination directory
      else
         ARRAY=`basename $a .lis`
         ARRAY=${ARRAY:11}
#         DDIR=$4/$ARRAY/proton_grid_131217
         DDIR=$4/$ARRAY/electron/
         mkdir -p $DDIR
# check if file already exists
         if [ -n "$6" ]
         then
            VFIL="$DDIR/$RUN"_G_"$MCAZ"deg.root
            if [ -e $VFIL ]
            then
               echo "found $VFIL"
               continue
            fi
         fi
         echo "copying $EFIL"
         dccp $EFIL $DDIR/$RUN"_"$a"_evndisp.root"
         if [ ! -n "$6" ]
         then
            MCAZ=`$EVNDISPSYS/bin/printRunParameter $DDIR/$RUN"_"$a"_evndisp.root" -mcaz`
         else
            MCAZ=$6
         fi
         VFIL="$DDIR/$RUN"_G_"$MCAZ"deg.root
         mv -v -f $DDIR/$RUN"_"$a"_evndisp.root" $DDIR/$RUN"_G_"$MCAZ"deg.root"
# soft link into main proton directory
#        ln -f -s $DDIR/$RUN"_G_"$MCAZ"deg.root" $4/$ARRAY/proton/$RUN"_G_"$MCAZ"deg.root"
      fi
    done
done

exit
