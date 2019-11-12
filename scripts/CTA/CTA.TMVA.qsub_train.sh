#!/bin/bash
#
# script to train cuts/MVAs with TMVA
#
#
#

RPARA=RUNPARA
EBIN=EEEE

ulimit -n 2056

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

PFILE=${RPARA}_${EBIN}
rm -f $PFIL.log

$EVNDISPSYS/bin/trainTMVAforGammaHadronSeparation $PFILE.runparameter > $PFILE.log

CDIR=`dirname $PFILE`
# remove .C files (never used; we use the XML files)
rm -f $CDIR/BDT_${EBIN}*.C
# remove complete_BDTroot at the end of the run
# (generally not used, but takes up lots of disk space)
# rm -rf $CDIR/complete_BDTroot/BDT_${EBIN}*
rm -rf $CDIR/complete_BDTroot

# dump runparameter file into the log files and zip it
if [ -e $PFILE.log ]
then
   cat $PFILE.runparameter >> $PFILE.log
   bzip2 -f -v $PFILE.log
fi
rm -f $PFILE.runparameter


exit
