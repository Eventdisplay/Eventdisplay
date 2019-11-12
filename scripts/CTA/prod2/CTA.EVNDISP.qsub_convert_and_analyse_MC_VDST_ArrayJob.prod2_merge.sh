#!/bin/bash
#
# script to merge, convert sim_telarray files to EVNDISP DST file and then run eventdisplay
#
# PROD2 analysis
#
#

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

# list of sim_telarray files to analysis
ILIST=SIMTELLIST
# all simtel files (only if SCSST files are not found in current directory)
FLIST=FULLLIST
ILINE=$SGE_TASK_ID
PART=PAAART
SUBA="ARRAY"
KEEP=KEEEEEEP
ACUT=ARC
DSET=DATASET
LOGF=FLL
PEDFILE=PPPP
PPOP=UUUU

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


################################
# sim_telarray files 

# string in file 1
IDFIL1="prod2_desert"
# string in file 2
IDFIL2="prod2-sc-sst-x_desert"

#############
# first file (STD configuration)
IFIL1=`head -n $ILINE $ILIST | tail -n 1`
echo "sim_telarray DATA FILE 1: "
echo $IFIL1
if [ ! -e $IFIL1 ]
then
    echo "ERROR: sim_telarray file does not exist:"
    echo $IFIL1
    exit
fi
#############
# second file (SCSST configuration)

# first check if there is a second run list
FLISTnd="$ILIST".SCSST
if [ -e $FLISTnd ]
then
   IFIL2=`head -n $ILINE $FLISTnd | tail -n 1`
# then check if the file is in the same directory (with a different name)
else
   IFIL2=${IFIL1/$IDFIL1/$IDFIL2}
fi
echo "sim_telarray DATA FILE 2: "
echo $IFIL2
# check if file exist
if [ ! -e $IFIL2 ]
then
    echo "ERROR: sim_telarray file does not exist:"
    echo $IFIL2
    exit
fi

################################
# copy files on temporary disk
echo
echo "COPYING FILES TO $TMPDIR"
# check if files are on local disc or on dCache
# (note: DESY only, this can not handle mixed lists)
z=1
for F in $IFIL1 $IFIL2 
do
    G=`basename $F`
    mkdir -p $TMPDIR"/FILE$z"
    if [[ $F = *acs* ]]
    then
      export DCACHE_CLIENT_ACTIVE=1
      dccp $F $TMPDIR"/FILE$z/"$G
    else
      cp -v -f $F $TMPDIR"/FILE$z/"
    fi
    let "z = $z + 1"
done
G=`basename $IFIL1`
IFIL1=$TMPDIR"/FILE1/"$G
G=`basename $IFIL2`
IFIL2=$TMPDIR"/FILE2/"$G
echo "files in $TMPDIR"
echo $IFIL1
echo $IFIL2

###############################
# log file directory
DATE=`date +"%y%m%d"`
mkdir -p $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$DATE-$LOGF
echo "Logfile directory "
echo $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$DATE-$LOGF

#####################################
# output file
OFIL=`basename $IFIL1 .gz`
echo
echo "OUTPUT FILE $OFIL"

####################################################################
# loop over all arrays
for N in $FIELD
do
# remove spaces
   N=`echo $N | tr -d ' '`

# start a log file 
   SLOG=$TMPDIR/$OFIL.$N.shell.log
   touch $SLOG
   echo "files in $TMPDIR" >> $SLOG
   echo "FILE1: $IFIL1" >> $SLOG
   echo "FILE2: $IFIL2" >> $SLOG

   echo "RUNNING _"$N"_"  >> $SLOG
# output data files are written to this directory
   ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/"$DSET"/"$N"/"$PART"/"
   mkdir -p $ODIR
   mkdir -p $ODIR/Calibration/

# calibration direction (mainly important for NN cleaning)
   CDIR=$TMPDIR
   mkdir -p $CDIR/Calibration

####################################################################
# execute merger
   SIMFIL="$TMPDIR/mergedFile.simtel.gz"
   $HESSIOSYS/bin/merge_simtel $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod2$N.merge.lis $IFIL1 $IFIL2 $SIMFIL >& $TMPDIR/$OFIL.$N.merge.log

# check that output file exists
   echo $FILE1
   ls -l $IFIL1
   echo $FILE2
   ls -l $IFIL2
   echo $SIMFIL
   ls -l $SIMFIL
   
####################################################################
# execute converter
   echo "TMPDIR FILES " $SIMFIL  >> $SLOG
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
  cp -v -f $TMPDIR/[0-9]*.root $ODIR/$RUNN"_"$ILINE"_"$MCAZ"deg.root"  >> $SLOG

####################################################################
# move dst (if required) and evndisp files to data directory
   if [ "$KEEP" == "1" ]
   then
      mkdir -p $ODIR/VDST
      cp -v -f $TMPDIR/$OFIL.root $ODIR/VDST/
   fi
   ls -lh $TMPDIR/*.root >> $SLOG
   ls -lh $TMPDIR/*.log >> $SLOG
# clean up 
   rm -f $TMPDIR/$OFIL.root
   rm -f $TMPDIR/[0-9]*.root
   rm -f $SIMFIL
   echo "==================================================================="
   echo
done

####################################################################
# tar the log files
cd $TMPDIR
tar -czvf $OFIL.merge.tar.gz *.log
mv -v -f $OFIL.merge.tar.gz $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$DATE-$LOGF/

exit
