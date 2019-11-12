#!/bin/sh
#
# calculate effective areas and instrument response functions for CTA
#
#
#
##############################################################################

echo
echo "calculating effective areas for CTA: create run scripts"
echo "------------------------------------------------"
echo

source $EVNDISPSYS/setObservatory.sh CTA

######################################################################
# input variables
######################################################################
ANAPAR=PPPANAPAR
ARRAY=PPPARRAY
RECID=PPPRECID
PTYPE=PPPVPART
CFIL=PPPCCUT
ODIR=PPPODIR
GFILLING=PPPGMOD
DSET=PPPDSET
MCAZ="PPPMCAZ"

#######################################
# read analysis values from parameter file
if [ ! -e $ANAPAR ]
then
  echo "error: analysis parameter file not found: $ANAPAR" 
  exit
fi
cp -f $ANAPAR $TMPDIR
ANAPARF=`basename $ANAPAR`
ANAPAR="${TMPDIR}/${ANAPARF}"
# check again than runparameter file is available
if [ ! -e $ANAPAR ]
then
  echo "error: analysis parameter file not found in tmp directory: $ANAPAR" 
  exit
fi
echo "reading analysis parameter from $ANAPAR"
cat $ANAPAR

# off-axis binning
BFINEBINNING="FALSE"
if grep -q OFFAXISFINEBINNING $ANAPAR
then
    BFINEBINNING=`grep OFFAXISFINEBINNING $ANAPAR | awk {'print $2'}`
fi

NIMAGESMIN=`grep NIMAGESMIN $ANAPAR | awk {'print $2'}`
# get telescope type dependent cuts 
NCUTLST=`grep NLST $ANAPAR | awk {'print $2'}`
NCUTMST=`grep NMST $ANAPAR | awk {'print $2'}`
NCUTSST=`grep NSST $ANAPAR | awk {'print $2'}`
NCUTSCMST=`grep NSCMST $ANAPAR | awk {'print $2'}`
for T in LST MST SST SCMST 
do
    NCUT="NCUT${T}"
    if [ -z "${!NCUT}" ]
    then
       declare "NCUT${T}=0"
    fi
done
# get all directories
ANADIR=`grep MSCWSUBDIRECTORY  $ANAPAR | awk {'print $2'}`
TMVACUT=`grep TMVASUBDIR $ANAPAR | awk {'print $2'}`
EFFAREADIR=`grep EFFAREASUBDIR $ANAPAR | awk {'print $2'}`
EFFAREABASEDIR=`grep EFFAREASUBBASEDIR $ANAPAR | awk {'print $2'}`
if [ -z $EFFAREABASEDIR ]
then
   EFFAREABASEDIR=$EFFAREADIR
fi
# see if strict separation of training/testing events if possible
# (mscw files would be in a directory ....EFF)
if [ -e $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/${ANADIR}.EFFAREA.MCAZ${MCAZ} ]
then
    ANADIR=${ANADIR}.EFFAREA.MCAZ${MCAZ}
fi

# observation time
OBSTIME=`grep OBSERVINGTIME_H $ANAPAR | awk {'print $2'}`
GETXOFFYOFFAFTERCUTS=`grep GETXOFFYOFFAFTERCUTS $ANAPAR | awk  {'print $2'}`
echo "Input parameters read from $ANAPAR"
echo "  Analysis parameters: $NIMAGESMIN $ANADIR $TMVACUT $EFFAREADIR $OBSTIME"

if [ -z "$ANADIR" ] || [ -z "$NIMAGESMIN" ] || [ -z "$TMVACUT" ] || [ -z "$EFFAREADIR" ] || [ -z "$OBSTIME" ]
then
  echo "error: analysis parameter file not correct: $ANAPAR" 
  echo " one variable missing"
  exit
fi

######################################################################
# directories
######################################################################
echo "data (input) directory"
DDIR=$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/
echo $DDIR
mkdir -p $DDIR
echo "output data directory"
echo $ODIR
mkdir -p $ODIR

declare -a PTYPELIST=("gamma_onSource" "electron_onSource" "proton_onSource" "gamma_cone" "electron" "proton")
if [ $PTYPE = "GAMMA" ]
then
    declare -a PTYPELIST=("gamma_onSource" "gamma_cone")
fi

######################################################################
# maximum core distance to a telescope
######################################################################
MAXCDISTANCE="500."
if [ $RECID = "1" ]
then
    MAXCDISTANCE="200."
fi

######################################################################
# off-axis variables:
#
# OFFMEA: off-axis angle used to select BDTs (need to be the same
#                 as in BDT training)
# OFFMIN/OFFMAX: range used for signal cuts
# THETAMIN/THETAMX: range used for background cuts

######################################################################
# loop over all particle types
######################################################################
for PART in "${PTYPELIST[@]}"
do
    echo "$PART"
    echo
    ######################################################################
    # input files and input parameters
    ######################################################################
    # parameters which are the same for all particle types
    # no AZ bins (CTA sims are sorted already on the input side in AZ)
    AZBINS="0"
    TELTYPECUTS="1"
    # data directory
    # on source gamma rays
    if [ $PART = "gamma_onSource" ]
    then
       if [[ ${DSET:0:2} == "GR" ]]
       then
           MSCFILE=$DDIR/gamma*"deg$MCAZ"*"baseline_evndisp"*.mscw.root
       else
           MSCFILE=$DDIR/gamma_onSource."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root
       fi
       EFFFILE=$DDIR/EffectiveAreas/
       OFIL=gamma_onSource."$ARRAY"_ID"$RECID".eff
       OFFMIN=( 0. )
       OFFMAX=( 100000. )
       # changed 2015-11-17: always use off-axis TMVA for gamma/hadron separation
       # OFFMEA=( "0.0" )
       OFFMEA=( "0.5" )
    # NOTE: this is theta2
       THETA2MIN=( -1. )
    # using TMVA or angular resolution
       THETA2MAX=( -1. )
       ISOTROPY="0"
       DIRECTIONCUT="2"
    fi
    # isotropic gamma-rays: analyse in rings in camera distance
    if [ $PART = "gamma_cone" ]
    then
       if [[ ${DSET:0:2} == "GR" ]]
       then
           MSCFILE=$DDIR/gamma*"deg$MCAZ"*"baseline_cone10_evndisp"*.mscw.root
       else
           MSCFILE=$DDIR/gamma_cone."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root
       fi
       EFFFILE=$DDIR/EffectiveAreas/
       OFIL=gamma_cone."$ARRAY"_ID"$RECID".eff
       # note that these bins are also hardwired in VTableLookupRunParameter::setCTA_MC_offaxisBins
       if [ $BFINEBINNING = "TRUE" ]
       then
           OFFMIN=( 0.0 1.0 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 )
           OFFMAX=( 1.0 2.0 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 )
           OFFMEA=( 0.5 1.5 1.5 2.5  2.5 2.5  2.5 3.5  3.5 3.5  3.5 4.5  4.5 4.5  4.5 5.5  5.5 )
       else
           OFFMIN=( 0.0 1.0 2.0 3.0 4.0 5.0 )
           OFFMAX=( 1.0 2.0 3.0 4.0 5.0 6.0 )
           OFFMEA=( 0.5 1.5 2.5 3.5 4.5 5.5 )
       fi
    # NOTE: this is theta2
       THETA2MIN=( -1. )
       THETA2MAX=( -1. )
       ISOTROPY="0"
       DIRECTIONCUT="2"
    fi
    if [ $PART = "electron" ] || [ $PART = "electron_onSource" ]
    then
       if [[ ${DSET:0:2} == "GR" ]]
       then
           MSCFILE=$DDIR/electron*"deg$MCAZ"*mscw.root
       else
           MSCFILE=$DDIR/electron."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root
       fi
       EFFFILE=$DDIR/EffectiveAreas/
       OFFMIN=( 0. )
       OFFMAX=( 100000. )
    # NOTE: this is theta and not theta2
       if [ $PART = "electron" ]
       then
          OFIL=electron."$ARRAY"_ID"$RECID".eff
           if [ $BFINEBINNING = "TRUE" ]
           then
               THETA2MIN=( 0.0 1.0 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 )
               THETA2MAX=( 1.0 2.0 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 )
               OFFMEA=(    0.5 1.5 1.5 2.5  2.5 2.5  2.5 3.5  3.5 3.5  3.5 4.5  4.5 4.5  4.5 5.5  5.5 )
          else
              THETA2MIN=( 0.0 1.0 2.0 3.0 4.0 5.0 )
              THETA2MAX=( 1.0 2.0 3.0 4.0 5.0 6.0 )
              OFFMEA=( 0.5 1.5 2.5 3.5 4.5 5.5 )
          fi
       else
          OFIL=electron_onSource."$ARRAY"_ID"$RECID".eff
          THETA2MIN=( 0. )
          THETA2MAX=( 1. )
          # changed 2015-11-17: always use off-axis TMVA for gamma/hadron separation
          # OFFMEA=( "0.0" )
          OFFMEA=( 0.5 )
       fi   
       ISOTROPY="1"
       DIRECTIONCUT="0"
    fi
    if [ $PART = "proton" ] || [ $PART = "proton_onSource" ]
    then
       if [[ ${DSET:0:2} == "GR" ]]
       then
           MSCFILE=$DDIR/proton*"deg$MCAZ"*mscw.root
       else
           MSCFILE=$DDIR/proton*."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root
       fi    
       if [ $ARRAY = "V5" ]
       then
          MSCFILE=$DDIR/proton."$ARRAY"_ID"$RECID$MCAZ".mscw.root
       fi
       EFFFILE=$DDIR/EffectiveAreas/
       OFIL=proton."$ARRAY"_ID"$RECID".eff
       OFFMIN=( 0. )
       OFFMAX=( 100000. )
    # NOTE: this is theta and not theta2
       if [ $PART = "proton" ] 
       then
          OFIL=proton."$ARRAY"_ID"$RECID".eff
           if [ $BFINEBINNING = "TRUE" ]
           then
               THETA2MIN=( 0.0 1.0 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 )
               THETA2MAX=( 1.0 2.0 2.5 2.75 3.0 3.25 3.5 3.75 4.0 4.25 4.5 4.75 5.0 5.25 5.5 5.75 6.0 )
               OFFMEA=(    0.5 1.5 1.5 2.5  2.5 2.5  2.5 3.5  3.5 3.5  3.5 4.5  4.5 4.5  4.5 5.5  5.5 )
          else
              THETA2MIN=( 0.0 1.0 2.0 3.0 4.0 5.0 )
              THETA2MAX=( 1.0 2.0 3.0 4.0 5.0 6.0 )
              OFFMEA=( 0.5 1.5 2.5 3.5 4.5 5.5 )
          fi
       else
          OFIL=proton_onSource."$ARRAY"_ID"$RECID".eff
          THETA2MIN=( 0. )
# (default value 2017-04-30)          THETA2MAX=( 2. )

# (new from 2017-04-28)
          # full system / LST only with smaller FOV
          if [ $RECID = "0" ] || [ $RECID = "1" ]
          then
              THETA2MAX=( 1. )
          # MST/SST system allow for larger FOV
          else
              THETA2MAX=( 2. )
          fi
          # changed 2015-11-17: always use off-axis TMVA for gamma/hadron separation
          # OFFMEA=( "0.0" )
          OFFMEA=( 0.5 )
       fi
       ISOTROPY="1"
       DIRECTIONCUT="0"
    fi
    NOFF=${#OFFMIN[@]}
    NTH2=${#THETA2MIN[@]}
    echo "Number of offaxis bins $NOFF $NTH2 $BFINEBINNING"
    ######################################################################


    ###############################################################################
    # loop over all MC cuts
    for ((i=0; i < $NOFF; i++))
    do
       iMIN=${OFFMIN[$i]}
       iMAX=${OFFMAX[$i]}
    # loop over all theta2 cuts
       for ((j=0; j < $NTH2; j++))
       do
         jMIN=${THETA2MIN[$j]}
         jMAX=${THETA2MAX[$j]}

    ###############################################################################
    # theta2 cut of protons and electron should match the rings from the isotropic gammas
          if [ $PART = "proton" ] || [ $PART = "proton_onSource" ] || [ $PART = "electron" ] || [ $PART = "electron_onSource" ]
          then
             jMIN=$(echo "$jMIN*$jMIN" | bc -l )
             jMAX=$(echo "$jMAX*$jMAX" | bc -l )
          fi

    ###############################################################################
    # create cut file
          iCBFILE=`basename $CFIL`      
          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ] 
          then
              CFILP="${CFIL}.gamma.dat"
          else
              CFILP="${CFIL}.CRbck.dat"
          fi
          iCFIL=$ODIR/ANASUM.GammaHadron-$DSET-$PART-$i-$j-MCAZ${MCAZ}.$iCBFILE.dat

          if [ ! -e $CFILP ]
          then
            echo "ERROR: cut file does not exist:"
            echo $CFILP
            exit
          fi
          cp -f $CFILP $iCFIL

    # wobble offset
          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ] 
          then
             WOBBLEOFFSET=${OFFMEA[$i]}
          else
             WOBBLEOFFSET=${OFFMEA[$j]}
          fi
    # angular resolution file
          if [ $PART = "gamma_onSource" ] 
          then
             ANGRESFILE="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/AngularResolution/gamma_onSource."$ARRAY"_ID"$RECID".eff-0.root"
          else
             ANGRESFILE="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/AngularResolution/gamma_cone."$ARRAY"_ID"$RECID".eff-$i.root"
          fi
    # particle number file
          if [ $PART = "gamma_onSource" ] || [ $PART = "electron_onSource" ] || [ $PART = "proton_onSource" ]
          then
             PNF="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/QualityCuts001CU/ParticleNumbers."$ARRAY".00.root"
          elif [ $PART = "gamma_cone" ]
          then
             PNF="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/QualityCuts001CU/ParticleNumbers."$ARRAY".$i.root"
          else
             PNF="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/QualityCuts001CU/ParticleNumbers."$ARRAY".$j.root"
          fi
          

          sed -i -e "s|OFFMIN|$iMIN|" \
                 -e "s|OFFMAX|$iMAX|" \
                 -e "s|THETA2MIN|$jMIN|" \
                 -e "s|THETA2MAX|$jMAX|" \
                 -e "s|DIRECTIONCUT|$DIRECTIONCUT|" \
                 -e "s|SUBARRAY|$ARRAY|" \
                 -e "s|MINIMAGES|$NIMAGESMIN|" \
                 -e "s|NTELTYPELST|$NCUTLST|" \
                 -e "s|NTELTYPEMST|$NCUTMST|" \
                 -e "s|NTELTYPESST|$NCUTSST|" \
                 -e "s|NTELTYPESCMST|$NCUTSCMST|" \
                 -e "s|WOBBLEOFFSET|$WOBBLEOFFSET|" \
                 -e "s|TMVACUTDIR|$TMVACUT|" \
                 -e "s|DATASET|$DSET|" \
                 -e "s|ANGRESFILE|$ANGRESFILE|" \
                 -e "s|PARTICLENUMBERFILE|$PNF|" \
                 -e "s|MAXCOREDISTANCE|$MAXCDISTANCE|" \
                 -e "s|OBSERVINGTIME_H|$OBSTIME|" $iCFIL

          cd $EVNDISPSYS/scripts/CTA/
          echo $iCFIL

    ###############################################################################
    # create run list
          MSCF=$ODIR/effectiveArea-CTA-$DSET-$PART-$i-$j.$ARRAY-MCAZ${MCAZ}-${ANADIR}.dat
          rm -f $MSCF
          echo "effective area data file for $PART $i $j" > $MSCF
    ###############################################################################
    # general run parameters
    ###############################################################################
    # filling mode
    ###############################################################################
    # fill IRFs and effective areas
          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ]
          then
    # filling mode 0: fill and use angular resolution for energy dependent theta2 cuts
             echo "* FILLINGMODE $GFILLING" >> $MSCF
          else
    # background: use fixed theta2 cut
             echo "* FILLINGMODE 3" >> $MSCF
          fi
    # fill IRFs only
          echo "* ENERGYRECONSTRUCTIONMETHOD 1" >> $MSCF
          echo "* ENERGYAXISBINS 60" >> $MSCF
          echo "* ENERGYRECONSTRUCTIONQUALITY 0" >> $MSCF
    # one azimuth bin only
          echo "* AZIMUTHBINS $AZBINS" >> $MSCF
          echo "* ISOTROPICARRIVALDIRECTIONS $ISOTROPY" >> $MSCF
          echo "* TELESCOPETYPECUTS $TELTYPECUTS" >> $MSCF
    # do fill analysis (a 1 would mean that MC histograms would be filled only)
          echo "* FILLMONTECARLOHISTOS 0" >> $MSCF
    # spectral index & CR spectra
          if [ $PART = "proton" ] || [ $PART = "proton_onSource" ]
          then
             echo "* ENERGYSPECTRUMINDEX  1 2.5 0.1" >> $MSCF
             echo "* ESPECTRUM_FOR_WEIGHTING $CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat 0" >> $MSCF
             if  [ $GETXOFFYOFFAFTERCUTS = "yes" ]
             then	
                 echo "* GETXOFFYOFFAFTERCUTS 1" >> $MSCF
             fi    

          fi
          if [ $PART = "electron" ] || [ $PART = "electron_onSource" ]
          then
             echo "* ENERGYSPECTRUMINDEX  1 2.5 0.1" >> $MSCF
             echo "* ESPECTRUM_FOR_WEIGHTING $CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CR.dat 8" >> $MSCF
             if  [ $GETXOFFYOFFAFTERCUTS = "yes" ]
             then
             echo "* GETXOFFYOFFAFTERCUTS 1" >> $MSCF
             fi
          fi
          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ]
          then
             echo "* ENERGYSPECTRUMINDEX  1 2.5 0.1" >> $MSCF
             echo "* ESPECTRUM_FOR_WEIGHTING $CTA_EVNDISP_AUX_DIR/AstroData/TeV_data/EnergySpectrum_literatureValues_CrabNebula.dat 5" >> $MSCF
          fi
    # REMOVED IGNOREFRACTIONOFEVENTS:
    #  - training is randomizing the input files, therefore this
    #    is anyway not a clean approach anymore
    #    solution would be to use the same input file list as
    #    BDT training
    # first half of data set is not used (as these events are used for the TMVA training)
#          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ]
#          then
#             echo "* IGNOREFRACTIONOFEVENTS 0.5" >> $MSCF
#          fi
#          if [ $PART = "proton" ] || [ $PART = "proton_onSource" ]
#          then
    # first half of data set is not used (as these events are used for the TMVA training)
#             echo "* IGNOREFRACTIONOFEVENTS 0.5" >> $MSCF
#          fi
          echo "* CUTFILE $iCFIL" >> $MSCF
          echo "* SIMULATIONFILE_DATA $MSCFILE" >> $MSCF

    # output file
          if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ]
          then
             OFIX=$ODIR/$OFIL-$i
          else
             OFIX=$ODIR/$OFIL-$j
          fi

          echo
          echo "preparing new analysis run"
          echo "--------------------------"
          echo
          echo "gamma/hadron separation file"
          echo $iCFIL
          echo $PNF

      ##############################
      # run effective area code
          $EVNDISPSYS/bin/makeEffectiveArea $MSCF $OFIX.root > $OFIX.log

      ##############################
      #  cleanup
      # (reduce number of files)

          cat $MSCF >> $OFIX.log
          rm -f $MSCF
          cat  $iCFIL >> $OFIX.log
          rm -f $iCFIL

          bzip2 -f $OFIX.log

       done
    done
done

#############
# run particle rate file calculator for QualityCutsStep
if [[ $ODIR == *"QualityCuts"* ]]
then
     echo "RUNNING PARTILE RATE DETERMINATION"
     AXDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREABASEDIR/AngularResolution/"

     # onSource
     LLOG=$ODIR/ParticleNumbers.$ARRAY.$RECID.onSource.log
     rm -f $LLOG
     $EVNDISPSYS/bin/writeParticleRateFilesFromEffectiveAreas  $ARRAY onSource $RECID $ODIR $AXDIR > $LLOG
     # cone
     LLOG=$ODIR/ParticleNumbers.$ARRAY.$RECID.cone.log
     rm -f $LLOG
     # off-axis fine binning
     if [ $BFINEBINNING = "TRUE" ]
     then
         $EVNDISPSYS/bin/writeParticleRateFilesFromEffectiveAreas  $ARRAY coneFB $RECID $ODIR $AXDIR > $LLOG
     else
     # off-axis std binning
         $EVNDISPSYS/bin/writeParticleRateFilesFromEffectiveAreas  $ARRAY cone $RECID $ODIR $AXDIR > $LLOG
     fi
fi

exit
