#!/bin/bash
#
# train disp for CTA
#
# called from CTA.DISPTRAINING_sub_analyse.sh
#
#

ODIR=OFILE
RECID=RECONSTRUCTIONID
TTYPE=TELTYPE
BDT=BDTTYPE
TLIST=ILIST
TMVAO=TTT
DSET="DATASET"
ARRAY=AAA

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

# output data files are written to this directory
mkdir -p $ODIR

# delete old log files
rm -f $ODIR/${BDT}-${TTYPE}.training.log
# delete old training files
rm -f $ODIR/*${TTYPE}*

# array layout file

# strip array scaling (two characters) for paranal sites
DARR=${ARRAY}
if [[ $DSET == *"paranal"* ]] && [[ $DSET != *"prod3b"* ]]
then
   DARR=${ARRAY%??}
   ADIR=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3${DARR}.lis
elif [[ $DSET == *"LaPalma"* ]] || [[ $DARR == "Sb"* ]]
then
   DARR=${ARRAY}
   ADIR=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3${DARR}.lis
else
   ADIR=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3Sb${DARR:1}.lis
fi

#########################################
# train TMVA
$EVNDISPSYS/bin/trainTMVAforAngularReconstruction $TLIST $ODIR 0.5 ${RECID} ${TTYPE} ${BDT} ${TMVAO} ${ADIR} > $ODIR/${BDT}-${TTYPE}.training.log 2>&1
#########################################

##############
# cleanup

# remove root output file of no telescope of this type is fine (declutter)
if grep -Fxq "Number of telescope types: 0" $ODIR/${BDT}-${TTYPE}.training.log
then
     rm -f -v $ODIR/${BDT}"_"${TTYPE}.root
fi

# pipe file list into log file
if [ -e $TLIST ]
then
    echo "========================================" >> $ODIR/${BDT}-${TTYPE}.training.log
    echo "List of files used for training: " >> $ODIR/${BDT}-${TTYPE}.training.log
    cat $TLIST >> $ODIR/${BDT}-${TTYPE}.training.log
    rm -f $TLIST
fi
# cleanup
bzip2 -v -f $ODIR/${BDT}-${TTYPE}.training.log


exit
