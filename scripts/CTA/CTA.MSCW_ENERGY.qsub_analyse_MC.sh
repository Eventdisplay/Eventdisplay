#!/bin/bash
#
# script to analyse CTA MC files with lookup tables
#
#

TABFIL=TABLEFILE
RECID=RECONSTRUCTIONID
IFIL="IIIIFIL"
TFIL=TTTTFIL
ARRAY=ARRAYYY
DSET=DATASET
ADIR=AAAAADIR
MCAZ=AZIMUTH

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

# output data and log files are written to this directory
ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/"$DSET"/"$ARRAY"/$ADIR/"
mkdir -p $ODIR

# delete old log files
rm -f $ODIR/$TFIL.log
rm -f $ODIR/$TFIL.table.log

#########################################
# smooth / cp table file to temp disk
SMOOTH="TRUE"
if [ $SMOOTH == "TRUE" ]
then
    $EVNDISPSYS/bin/smoothLookupTables $CTA_EVNDISP_AUX_DIR/Tables/$TABFIL-$ARRAY.root $TMPDIR/$TABFIL-$ARRAY.root > $ODIR/$TFIL.table.log
else
    cp -v $CTA_EVNDISP_AUX_DIR/Tables/$TABFIL-$ARRAY.root $TMPDIR/$TABFIL-$ARRAY.root > $ODIR/$TFIL.table.log
fi

#########################################
# cp all evndisp root files to TMPDIR
echo "Coppying evndisp root files to TMPDIR"
echo "number of files: "
wc -l $IFIL

NFIL=`cat $IFIL`
DDIR=$TMPDIR/data
mkdir -p $DDIR
for F in $NFIL
do
   cp $F $DDIR/
done
find $DDIR/ -name "*.root" > $TMPDIR/iList.list
# ls -1 $DDIR/*.root > $TMPDIR/iList.list

# check disk space on TMPDIR
du -h -c $TMPDIR

###############################################
# mscw_energy command line options
###############################################
MOPT="-pe -arrayrecid=$RECID -noNoTrigger -writeReconstructedEventsOnly -shorttree -useMedian=2 -add_mc_spectral_index=2.5"

# strip array scaling (two characters) for paranal sites
DARR=${ARRAY}
if  [[ $DSET == *"prod4"* ]]
then
    LISFILE=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod4${DARR}.lis
elif [[ $DSET == *"paranal"* ]] && [[ $DSET != *"prod3b"* ]]
then
   DARR=${ARRAY%??}
   LISFILE=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3${DARR}.lis
elif [[ $DSET == *"LaPalma"* ]] || [[ $DARR == "Sb"* ]]
then
   DARR=${ARRAY}
   LISFILE=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3${DARR}.lis
else
   LISFILE=$CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3Sb${DARR:1}.lis
fi

#########################################
# options for simple stereo reconstruction
MOPT="$MOPT -redo_stereo_reconstruction -sub_array_sim_telarray_counting $LISFILE -minangle_stereo_reconstruction=15"
# MOPT="$MOPT -redo_stereo_reconstruction -sub_array_sim_telarray_counting $LISFILE -minangle_stereo_reconstruction=5"

# IMPORTANT: this must be the same or lower value as in dispBDT training
MOPT="$MOPT -maxloss=0.2"

#########################################
# disp main directory name
DISPSUBDIR="BDTdisp.${ARRAY}.T19"

#########################################
# options for DISP method (direction)
DISPDIR="${CTA_USER_DATA_DIR}/analysis/AnalysisData/"$DSET"/${DISPSUBDIR}/BDTDisp/${MCAZ}/BDTDisp_BDT_"
MOPT="$MOPT -tmva_nimages_max_stereo_reconstruction=100 -tmva_filename_stereo_reconstruction $DISPDIR"

##########################################################################################################
# options for DISP method (direction error)
DISPERRORDIR="${CTA_USER_DATA_DIR}/analysis/AnalysisData/"$DSET"/${DISPSUBDIR}/BDTDispError/${MCAZ}/BDTDispError_BDT_"
MOPT="$MOPT  -tmva_filename_disperror_reconstruction $DISPERRORDIR -tmva_disperror_weight 50"

##########################################################################################################
# options for DISP method (energy)
DISPENERGYDIR="${CTA_USER_DATA_DIR}/analysis/AnalysisData/"$DSET"/${DISPSUBDIR}/BDTDispEnergy/${MCAZ}/BDTDispEnergy_BDT_"
MOPT="$MOPT -tmva_filename_energy_reconstruction $DISPENERGYDIR"

################################
# telescope type dependent weight
# prod3b production
if [[ $DSET == *"prod3b"* ]]
   then
       MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3b.TelescopeWeights.dat"
# prod3 MPIK La Palma Production
elif [[ $DSET == *"LaPalma"* ]]
then
    MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3N.LaPalma20deg.TelescopeWeights.dat"
elif [[ $DSET == *"paranal"* ]]
then
   # prod3 paranal 40 deg (post-Liverpool production)
   if [[ $DSET == *"40deg"* ]]
   then
       if [[ $ARRAY == *"-F"* ]]
       then
           MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3N.Paranal40deg.TelescopeWeightsFlashCam.dat"
       else
           MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3N.Paranal40deg.TelescopeWeightsNectarCam.dat"
       fi
   # prod3 paranal 20 deg (pre-Liverpool production)
   else
       if [[ $ARRAY == *"-F"* ]]
       then
           MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3N.Paranal20deg.TelescopeWeightsFlashCam.dat"
       else
           MOPT="$MOPT -teltypeweightfile $CTA_EVNDISP_AUX_DIR/DetectorGeometry/CTA.prod3N.Paranal20deg.TelescopeWeightsNectarCam.dat"
       fi
   fi
fi   
echo $MOPT

#########################################
# analyse MC file
$EVNDISPSYS/bin/mscw_energy $MOPT -tablefile $TMPDIR/$TABFIL-$ARRAY.root -inputfilelist $TMPDIR/iList.list -outputfile $TMPDIR/$TFIL.root >& $ODIR/$TFIL.log
#########################################

#########################################
# clean up and zip log files
if [ -e $ODIR/$TFIL.table.log ] && [ -e $ODIR/$TFIL.log ]
then
    cat $ODIR/$TFIL.table.log >> $ODIR/$TFIL.log
    rm -f $ODIR/$TFIL.table.log
fi
if [ -e $ODIR/$TFIL.table.log ]
then
    bzip2 -v -f $ODIR/$TFIL.table.log
fi
if [ -e $ODIR/$TFIL.log ]
then
    bzip2 -v -f $ODIR/$TFIL.log
fi

mv -f -v $TMPDIR/$TFIL.root $ODIR/

# sleep
sleep 2

exit
