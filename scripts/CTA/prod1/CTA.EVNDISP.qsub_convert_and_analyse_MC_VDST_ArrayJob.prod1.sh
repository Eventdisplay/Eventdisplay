#!/bin/bash
#
# script to convert sim_tel output files to EVNDISP DST file and then run eventdisplay
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

# set array
FIELD=$SUBA

# eventdisplay optimization
OPT="-shorttree -l2setspecialchannels nofile"

# set timtelarray file and cp simtelarray.gz file to TMPDIR
IFIL=`head -n $ILINE $ILIST | tail -n 1`
echo "$IFIL"
cp -v -f $IFIL $TMPDIR"/"

if [ ! -e $IFIL ]
then
   echo "Warning: SIMTELFILE does not exist"
   exit
fi

# output file
OFIL=`basename $IFIL .gz`
echo "OUTPUT FILE $OFIL"

#loop over all arrays
for N in $FIELD
do
   echo "RUNNING  $N"
# output data files are written to this directory
   ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/"$DSET"/"$N"/"$PART"/"
   mkdir -p $ODIR

####################################################################
# execute converter
   $EVNDISPSYS/bin/CTA.convert_hessio_to_VDST -a $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod1.$N.lis -o $TMPDIR/$OFIL.root $TMPDIR/$OFIL.gz >& $TMPDIR/$OFIL.$N.convert.log

####################################################################
# execute eventdisplay
  $EVNDISPSYS/bin/evndisp -sourcefile $TMPDIR/$OFIL.root -writenoMCTree -reconstructionparameter $ACUT $OPT -outputdirectory $TMPDIR >& $TMPDIR/$OFIL.$N.evndisp.log

####################################################################
# move evndisp files to data directory
   if [ $KEEP == "1" ]
   then
      mkdir -p $ODIR/VDST
      cp -v -f $TMPDIR/$OFIL.root $ODIR/VDST/
   fi
   rm -f -v $TMPDIR/$OFIL.root
   cp -v -f $TMPDIR/*.root $ODIR
done

# tar the log files
cd $TMPDIR
tar -czvf $OFIL.tar.gz *.log
mkdir -p $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$LOGF
mv -v -f $OFIL.tar.gz $CTA_USER_LOG_DIR"/analysis/AnalysisData/"$DSET/LOGFILES-$LOGF/

exit

