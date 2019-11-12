#!/bin/bash
#
# script to convert sim_tel output files to EVNDISP DST file and then run eventdisplay
#
# PROD2 analysis
#
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
TRGMASKDIR=TRIGGGG
PPOP=UUUU
STEPSIZE=STST

# set array
FIELD=$SUBA

###################################
# converter command line parameter
# settings for increased NSB
# COPT="-NSB 1.3 -f 1 -c $PEDFILE"
# default settings
COPT="-f 1 -c $PEDFILE"

# eventdisplay command line parameter
OPT="-averagetzerofiducialradius=0.5 -shorttree -l2setspecialchannels nofile -writenoMCTree -reconstructionparameter $ACUT $PPOP"

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

###############################
# log file directory
DATE=`date +"%y%m%d"`
mkdir -p $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$DATE-$LOGF

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

#####################################
# trigmask file (optional)
echo $TRGMASKDIR
if [[ ! -z "$TRGMASKDIR" ]] && [[ "$TRGMASKDIR" == "TRUE" ]]
then
   COPT="$COPT -t $TRGMASKDIR"
fi

####################################################################
# loop over all arrays
for N in $FIELD
do
# remove spaces
   N=`echo $N | tr -d ' '`
   echo "RUNNING _"$N"_"
# output data files are written to this directory
   ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/"$DSET"/"$N"/"$PART"/"
   mkdir -p $ODIR
   mkdir -p $ODIR/Calibration/

# calibration direction (mainly important for NN cleaning)
   CDIR=$TMPDIR
   mkdir -p $CDIR/Calibration

####################################################################
# execute converter
   SIMFIL=`ls $TMPDIR/*.simtel.gz`
   echo "TMPDIR FILES " $SIMFIL
   $EVNDISPSYS/bin/CTA.convert_hessio_to_VDST $COPT -a $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod2$N.lis -o $TMPDIR/$OFIL.root $SIMFIL >& $TMPDIR/$OFIL.$N.convert.log

####################################################################
# execute eventdisplay
  $EVNDISPSYS/bin/evndisp -sourcefile $TMPDIR/$OFIL.root $OPT -outputdirectory $TMPDIR -calibrationdirectory $CDIR >& $TMPDIR/$OFIL.$N.evndisp.log

####################################################################
# get runnumber and azimuth and rename output files; mv them to final destination
  MCAZ=`$EVNDISPSYS/bin/printRunParameter $TMPDIR/$OFIL.root -mcaz`
  RUNN=`$EVNDISPSYS/bin/printRunParameter $TMPDIR/$OFIL.root -runnumber`
  cp -v -f $TMPDIR/[0-9]*.root $ODIR/${RUNN}"_"$ILINE"_"$MCAZ"deg.root"
  cp -v -f $TMPDIR/$OFIL.$N.evndisp.log $ODIR/$RUNN"_"$ILINE"_"$MCAZ"deg.log"
# cp IPR files
  cp -v -f $CDIR/Calibration/[0-9]*IPRcontours.root $ODIR/Calibration/$RUNN"_"$ILINE"_"$MCAZ"deg.IPRcontours.root"

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

####################################################################
# tar the log files
cd $TMPDIR
tar -czvf $OFIL.tar.gz *.log
mv -v -f $OFIL.tar.gz $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$DATE-$LOGF/

exit
