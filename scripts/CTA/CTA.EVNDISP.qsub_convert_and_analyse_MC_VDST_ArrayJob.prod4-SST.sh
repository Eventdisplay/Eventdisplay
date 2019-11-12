#!/bin/bash
#
# script to convert sim_tel output files to EVNDISP DST file and then run eventdisplay
#
# PROD4 analysis
#
# Author: Gernot Maier
#

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

ILIST=SIMTELLIST
ILINE=$SGE_TASK_ID
PART=PAAART
SUBA="ARRAY"
KEEP=KEEEEEEP
ACUT=ARC
DSET=DATASET
LOGF=FLL
PEDFILE=PPPP
STEPSIZE=STST

# set array
FIELD=$SUBA

###################################
# converter command line parameter
COPT="-c $PEDFILE"

# eventdisplay command line parameter
OPT="-averagetzerofiducialradius=0.5 -reconstructionparameter $ACUT"

# set simtelarray file and cp simtelarray.gz file to TMPDIR
if [ ! -e $ILIST ]
then
   echo "ERROR: list of simulation files does not exist: $ILIST"
   exit
fi

# get file list (of $STEPSIZE files)
let "ILINE = $ILINE * $STEPSIZE"
echo "getting line(s) $ILINE from list"
echo "list: $ILIST"
IFIL=`head -n $ILINE $ILIST | tail -n $STEPSIZE`
IFIL0=`head -n $ILINE $ILIST | tail -n 1`
echo "DATA FILE(S)"
echo $IFIL
################################
# copy files on temporary disk
echo
echo "COPYING FILES TO $TMPDIR"
# check if files are on local disc (lustre) or on dCache
for F in $IFIL
do
    if [[ $F = *acs* ]]
    then
         export DCACHE_CLIENT_ACTIVE=1
         echo "F $F"
         G=`basename $F`
         echo "G $G"
         dccp $F $TMPDIR"/"$G
    else
         cp -v -f $F $TMPDIR"/"
    fi
done

for F in $FIL
do
    if [ ! -e $F ]
    then
       echo "ERROR: SIMTELFILE does not exist:"
       echo $F
       exit
    fi
done

#####################################
# output file
OFIL=`basename $IFIL0 .gz`
echo
echo "OUTPUT FILE $OFIL"

####################################################################
# loop over all arrays
for N in $FIELD
do
# remove spaces
   N=`echo $N | tr -d ' '`
   echo "RUNNING _"$N"_"
# output data files are written to this directory
   ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/"$DSET"/"${N}"/"$PART"/"
   mkdir -p $ODIR
   mkdir -p $ODIR/Calibration/

# calibration direction (mainly important for NN cleaning)
   CDIR=$TMPDIR
   mkdir -p $CDIR/Calibration

####################################################################
# execute zero suppression
#
# input sim_telarray file: $SIMFIL
# ($SIMFIL should be set then to the zero suppressed file)
#
# best to write to output (zero suppressed file) to the temporary disk on the 
# current note: $TMPDIR/<zerosuppressed file>
#
# $HESSIOSYS/bin/ ....

####################################################################
# execute converter
   SIMFIL=`ls $TMPDIR/*.simtel.gz`
   echo "TMPDIR FILES " $SIMFIL
   DETGEO=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod4${N}.lis
   ls -lh $DETGEO
   ls -lh $SIMFIL
   ls -lh $TMPDIR

   $EVNDISPSYS/bin/CTA.convert_hessio_to_VDST $COPT -a $DETGEO -o $TMPDIR/$OFIL.root $SIMFIL >& $TMPDIR/$OFIL.$N.convert.log
   ls -lh $TMPDIR/$OFIL.root

####################################################################
# execute eventdisplay
  if [ -e $TMPDIR/$OFIL.root ]
  then
      $EVNDISPSYS/bin/evndisp -sourcefile $TMPDIR/$OFIL.root $OPT -outputdirectory $TMPDIR -calibrationdirectory $CDIR >& $TMPDIR/$OFIL.$N.evndisp.log
  else
      echo "DST file not found: $TMPDIR/$OFIL.root" >& $TMPDIR/$OFIL.$N.evndisp.log 
  fi
  ls -ls $TMPDIR

####################################################################
# get runnumber and azimuth and rename output files; mv them to final destination
  if [ -e $TMPDIR/$OFIL.root ]
  then
      MCAZ=`$EVNDISPSYS/bin/printRunParameter $TMPDIR/$OFIL.root -mcaz`
      RUNN=`$EVNDISPSYS/bin/printRunParameter $TMPDIR/$OFIL.root -runnumber`
      cp -v -f $TMPDIR/[0-9]*.root $ODIR/${RUNN}"_PROD4_"$ILINE"_"$MCAZ"deg.root"
  fi
  bzip2 -v -f $TMPDIR/$OFIL.$N.evndisp.log
  cp -v -f $TMPDIR/$OFIL.$N.evndisp.log.bz2 $ODIR/$RUNN"_PROD4_"$ILINE"_"$MCAZ"deg.log.bz2"
  bzip2 -v -f $TMPDIR/$OFIL.$N.convert.log
  cp -v -f $TMPDIR/$OFIL.$N.convert.log.bz2 $ODIR/$RUNN"_PROD4_"$ILINE"_"$MCAZ"deg.convert.log.bz2"
  if [ -f $CDIR/Calibration/[0-9]*IPRcontours.root ] 
  then 
    cp -v -f $CDIR/Calibration/[0-9]*IPRcontours.root $ODIR/Calibration/$RUNN"_"$ILINE"_"$MCAZ"deg.IPRcontours.root"
  fi

####################################################################
# move dst (if required) and evndisp files to data directory
   if [ "$KEEP" == "1" ]
   then
      mkdir -p $ODIR/VDST
      cp -v -f $TMPDIR/$OFIL.root $ODIR/VDST/
   fi
   ls -lh $TMPDIR/*.root
# clean up 
   rm -f $TMPDIR/$OFIL.root
   rm -f $TMPDIR/[0-9]*.root
   echo "==================================================================="
   echo
done

exit
