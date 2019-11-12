#!/bin/bash
#
# fill tables for CTA
#
#

TFIL=TABLEFILE
RECID=RECONSTRUCTIONID
ARRAY=ARRRRRRR
CONE="CCC"
DSET="DATASET"
MCAZ="AZIMUTH"

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

# output data files are written to this directory
ODIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/$DSET/"$ARRAY"/Tables/"
mkdir -p $ODIR

# output log files are written to this directory
LDIR=$CTA_USER_DATA_DIR"/analysis/AnalysisData/$DSET/"$ARRAY"/Tables/"
mkdir -p $LDIR

# rename on-axis tables
if [ $CONE == "FALSE" ]
then
   TFIL=${TFIL}-onAxis
fi

# delete old log files
rm -f $LDIR/$TFIL-$ARRAY.log
# delete old table file (mscw_energy would otherwise stop with an error message)
rm -f $LDIR/$TFIL-$ARRAY.root

MCAZ=${MCAZ/_/}

################################
# generate input file lists
TMPLIST=$LDIR/${DSET}${MCAZ}.${ARRAY}.list
rm -f $TMPLIST
touch $TMPLIST
# for diffuse gamma rays: use also on-axis gamma rays to get sufficient statistics at high energies
if [ $CONE == "TRUE" ]
then
   find $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/gamma_cone/ -name "*[0-9]*[\.,_]${MCAZ}*.root" >> $TMPLIST
   echo "$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/gamma_cone/ -name \"*[0-9]*[\.,_]${MCAZ}*.root\""
# on-axis gamma rays
else
   find $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/gamma_onSource/ -name "*[0-9]*[\.,_]${MCAZ}*.root" >> $TMPLIST
   echo "$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/gamma_onSource/ -name \"*[0-9]*[\.,_]${MCAZ}*.root\""
fi

################################
# options for table filling
if [ $CONE == "TRUE" ]
then
   SETOFF="-CTAoffAxisBins"
fi

# strip array scaling (two characters) for paranal sites
DARR=${ARRAY}
if [[ $DSET == *"paranal"* ]] && [[ $DSET != *"prod3b"* ]] 
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

# MOPT="$SETOFF -pe -filltables=1 -ze=20. -noise=250 -woff=0.0 -minImages=2 -write1DHistograms"
# MOPT="$SETOFF -pe -filltables=1 -ze=20. -noise=250 -woff=0.0 -minImages=3 -write1DHistograms"
MOPT="$SETOFF -pe -filltables=1 -ze=20. -noise=250 -woff=0.0 -minImages=4  -write1DHistograms"
# options for reweighting of telescopes
MOPT="$MOPT -redo_stereo_reconstruction -sub_array_sim_telarray_counting $LISFILE -minangle_stereo_reconstruction=15"
# for 40 deg prod3b
# MOPT="$MOPT -maxnevents=3500000"
# for 20 deg prod3b (default)
# (remove for small arrays)
MOPT="$MOPT -maxnevents=3000000"

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
# fill tables
$EVNDISPSYS/bin/mscw_energy $MOPT -arrayrecid=$RECID -tablefile "$ODIR/$TFIL-$ARRAY.root" -inputfilelist ${TMPLIST} > $LDIR/$TFIL-$ARRAY.log
#########################################


# sleep
sleep 2

exit
