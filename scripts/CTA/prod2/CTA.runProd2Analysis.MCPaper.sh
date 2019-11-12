#!/bin/sh
#
# analysis submission for production 2 analysis (MC paper)
#
# this script is optimized for the DESY analysis
#
##############################################


if [ $# -ne 2 ] 
then
   echo 
   echo "./CTA.runProd2Analysis.sh <N/S/SM/P1/LAYOUT> <run mode>"
   echo
   echo "  N=prod2-North, S=prod2-South, SM=prod2-South-merged, P1=prod1, LAYOUT=array layout analysis"
   echo
   echo "  possible run modes are EVNDISP MAKETABLES ANATABLES TRAIN ANGRES QC PARTFIL CUTS PHYS "
   echo
   exit
fi
P2="$1"
RUN="$2"
#####################################
# reconstruction IDs
# RECID="$3"
RECID="0"

# should be either onSource or cone
OFFAXIS="cone"

#####################################
# qsub options (priorities)
#   _M_ = -; _X_ = " "
QSUBOPT="_M_P_X_cta_high"
QSUBOPT="_M_P_X_cta_high_X__M_js_X_999"

#####################################
# output directory for script parameter files
PDIR="$CTA_USER_LOG_DIR/tempRunParameterDir/"
mkdir -p $PDIR

#####################################
# analysis dates and table dates

TDATE="d20150824"
DATE="d20150824"
TMVADATE="d20150824"
EFFDATE="d20150824"
EFFVERSION="L1"

#####################################
# shower directions
#
# _180deg = south
# _0deg = north
MCAZ=( "_180deg" "_0deg" "" )

#####################################
# sites & sub array lists

##########################################################
# SOUTH
if [[ $P2 == "S" ]]
then
##########################################################
# data sets without trgmask files
# started 20141219
#   SITE=( "prod2-LeoncitoPP-NS" "prod2-SAC100-lowE-NS" "prod2-Aar-500m-NS" )
# MC paper: use new Aar simulations
# started: 20121220
#    SITE=( "prod2-Aar2014-NS" )
# Armazones
#     SITE=( "prod2-Armazones-NS" )
##########################################################
# 40 deg data sets 
# started: 
#   SITE=( "prod2-Aar-40deg-NS" "prod2-Leoncito-40deg-NS" "prod2-LeoncitoPP-40deg-NS" )
##########################################################
# NSB data sets
# started: 20141218
#   SITE=( "prod2-Leoncito-NSB-x1.00-NS" "prod2-Leoncito-NSB-x1.30-NS" "prod2-Leoncito-NSB-x1.50-NS" )
# data sets with trgmask files
# started: 20141220
#   SITE=( "prod2-SAC100-NS" "prod2-Leoncito-NS" )
#   SITE=( "prod2-SAC084-NS" )
#   SITE=( "prod2-SAC100-NS" "prod2-SAC084-NS" "prod2-Leoncito-NS" )
# all southern sites
#   SITE=( "prod2-Aar2014-NS" "prod2-SAC100-NS" "prod2-SAC084-NS" "prod2-Leoncito-NS" "prod2-SAC100-lowE-NS" "prod2-Aar-500m-NS" "prod2-LeoncitoPP-NS" "CL5040-prod2-Leoncito-NSB-x1.00-NS" "CL5040-prod2-Leoncito-NSB-x1.30-NS" "CL5040-prod2-Leoncito-NSB-x1.50-NS" "prod2-Aar-40deg-NS" "prod2-Leoncito-40deg-NS" "prod2-LeoncitoPP-40deg-NS" )
# on + off axis
#   OFFAXIS="TRUE"
#   SITE=( "prod2-Aar2014-NS" "prod2-LeoncitoPP-NS" "prod2-Armazones-NS" )
#   SITE=( "prod2-Aar2014-NS" "prod2-LeoncitoPP-NS" "prod2-Armazones-NS" "prod2-Aar-40deg-NS" "prod2-Leoncito-40deg-NS" "prod2-LeoncitoPP-40deg-NS" )
# on axis only
#   SITE=( "prod2-SAC100-NS" "prod2-Leoncito-NS" "prod2-Aar-500m-NS" "prod2-SAC100-lowE-NS" "prod2-Leoncito-NSB-x1.00-NS" "prod2-Leoncito-NSB-x1.30-NS" "prod2-Leoncito-NSB-x1.50-NS" )
   OFFAXIS="FALSE"
   SITE=( "prod2-Aar2014-NS" "prod2-Armazones-NS" )
#############################
# MC paper list of arrays
   ARRAY="subArray.MCpaper.dat"
   OFFAXIS="onSource"
   SITE=( "prod2-Aar2014" )
   EDM=( "-CC" )
   EDM=( "-NN" )
   ARRAY="subArray.2a.list"
############################
# sub array study with many different arrays (on-axis only)
elif [[ $P2 == "LAYOUT" ]]
then
   OFFAXIS="onSource"
   SITE=( "prod2-Aar2014" )
   SITE=( "prod2-Aar2014" "prod2-Aar-40deg" )
   SITE=( "prod2-Aar-40deg" )
   EDM=( "-NN" )
   ARRAY="subArray.ArrayLayouts.dat.13"
   ARRAY="subArray.MCpaper.dat.fullq.redo"
   ARRAY="subArray.MCpaper.dat.fullq"
##########################################################
# SOUTH - merged simulation files
elif [[ $P2 == "SM" ]]
then
   OFFAXIS="onSource"
   SITE="prod2-Aar2014"
#   SITE=( "prod2-Aar-40deg" )
   EDM=( "-NN" )
   ARRAY="subArray.TDR-q.dat"
   ARRAY="subArray.SM.dat"
   ARRAY="subArray.q.dat"
##########################################################
# NORTH
elif [[ $P2 == "N" ]]
then
# data set north
# started: 20141112 (DONE)
   OFFAXIS="cone"
   SITE=( "prod2-US-NS" "prod2-SPM-NS" "prod2-Tenerife-NS" )
   ARRAY="subArray.TDR-NN.dat"
else
   echo "error: unknown site; allowed are N or S"
   echo $P2
   exit
fi

#####################################
# particle types
PARTICLE=( "gamma_onSource" )
PARTICLE=( "gamma_onSource" "electron" "proton" "gamma_cone" )
PARTICLE=( "gamma_onSource" "electron" "proton" )

#####################################
# energy reconstruction
# (default is method 1)
EREC="1"

#####################################
# cut on number of images
NIMAGESMIN="2"

#####################################
# observing time [h]
OBSTIME=( "50h" "5h" "30m" "10m" "10h" "20h" "100h" "500h" )
OBSTIME=( "50h" "5h" "30m" "10m" "10h" "20h" "100h" "500h" "1m" "20s" )
OBSTIME=( "50h" )


#####################################
# loop over all sites
NSITE=${#SITE[@]}
for (( m = 0; m < $NSITE ; m++ ))
do
   # site name
   S=${SITE[$m]}
   # eventdisplay analysis method
   M=${EDM[$m]}

   echo
   echo "======================================================================================"
   echo "SITE: $S"
   echo "RUN: $RUN"

#####################################
# trigmask files
# (needed only for prod2 files with trigger bug)
# 20140313: this is partly broken...(works for Lenocito)
   if [[ $S == "prod2-Aar" ]]
   then
      TRG="/lustre/fs13/group/cta/prod2/Aar/simtel/trgmask/"
   elif [[ $S == "prod2-Leoncito" ]]
   then
      TRG="/lustre/fs16/group/cta/prod2/Leoncito/trgmask/"
# OLD SIMULATIONS WITH WRONG L1 TRIGGER
#      TRG="/lustre/fs13/group/cta/prod2/Leoncito/simtel/trgmask/"
   elif [[ $S == "prod2-SAC084" ]]
   then
      TRG="/lustre/fs13/group/cta/prod2/SAC084/"
   elif [[ $S == "prod2-SAC100" ]]
   then
      TRG="/lustre/fs13/group/cta/prod2/SAC100/"
   else
      TRG=""
   fi

##########################################
# run eventdisplay
    if [[ $RUN == "EVNDISP" ]]
    then

# loop over all particle types
      for ((i = 0; i < ${#PARTICLE[@]}; i++ ))
      do
	  N=${PARTICLE[$i]}

          # run list for 20 deg files on lustre (default data sets)
	  LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod2/$S-NS.$N"_20deg".list
          # run list for 20 deg high NSB files on lustre
          if [[ $S == *"NSB"* ]]
          then
              LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod2/NSB/$S-NS.$N"_20deg".list
          fi
          # run list for 40 deg files
          if [[ $S == *"40deg"* ]]
          then
              LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod2/40deg/$S-NS."$N".list
          fi
          echo $LIST

          ####################################
          # analysis of merged evndisp files
          if [[ $P2 == "SM" ]]
          then
               ./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod2_merge.sh $ARRAY $LIST $N $S$M 0 $i $QSUBOPT
          # standard evndisplay analysis
          else
               ./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod2.sh $ARRAY $LIST $N $S$M 0 $i $QSUBOPT $TRG
          fi 
       done
       continue
    fi
##########################################
# loop over all reconstruction IDs
    for ID in $RECID
    do
       MSCWSUBDIRECTORY="Analysis-ID$ID-$TDATE"
##########################################
# make tables
       if [[ $RUN == "MAKETABLES" ]]
       then
	  ./CTA.MSCW_ENERGY.sub_make_tables.sh tables_CTA-$S$M-ID$ID-$TDATE $ID $ARRAY $OFFAXIS $S$M $QSUBOPT
	  continue
##########################################
# analyse with lookup tables
       elif [[ $RUN == "ANATABLES" ]]
       then
	  TABLE="tables_CTA-$S$M-ID$ID-$TDATE-onAxis"
          if [[ $OFFAXIS == "cone" ]]
          then
              TABLE="tables_CTA-$S$M-ID$ID-$TDATE"
          fi
	  echo $TABLE
	  ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TABLE $ID $ARRAY $S$M $MSCWSUBDIRECTORY $OFFAXIS $QSUBOPT
	  continue
        fi

##########################################
# loop over all observation times
      for ((o = 0; o < ${#OBSTIME[@]}; o++ ))
      do
        OOTIME=${OBSTIME[$o]}

##########################################
# loop over all shower directions 
# (i.e. North and South)
       for ((a = 0; a < ${#MCAZ[@]}; a++ ))
       do
          AZ=${MCAZ[$a]}
# set run parameter file
	  PARA="$PDIR/scriptsInput.prod2.Erec$EREC.ID$ID$AZ.$S$AZ.runparameter"
	  rm -f $PARA
	  touch $PARA
	  echo "WRITING PARAMETERFILE $PARA"
	  NTYPF=NIM$NIMAGESMIN
	  EFFDIR="EffectiveArea-"$OOTIME"-Erec$EREC-ID$ID$AZ-$NTYPF-$EFFDATE"
          EFFFULLDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$S$M/$EFFDIR/"
#	  EFFDIR="/lustre/fs9/group/cta/users/$USER/CTA/analysis/AnalysisData/$S$M/$EFFDIR/"
	  echo "MSCWSUBDIRECTORY $MSCWSUBDIRECTORY" >> $PARA
	  echo "TMVASUBDIR BDT-Erec$EREC-ID$ID$AZ-$NTYPF-$TMVADATE" >> $PARA
	  echo "EFFAREASUBDIR $EFFDIR" >> $PARA
	  echo "RECID $ID" >> $PARA
	  echo "ENERGYRECONSTRUCTIONMETHOD $EREC" >> $PARA
	  echo "NIMAGESMIN $NIMAGESMIN" >> $PARA
#          echo "NTYPEMIN_0 3" >> $PARA
#          echo "NTYPEMIN_1 2" >> $PARA
#          echo "NTYPEMIN_2 2" >> $PARA
	  echo "OBSERVINGTIME_H $OOTIME" >> $PARA
          echo "GETXOFFYOFFAFTERCUTS yes" >> $PARA
##########################################
# train BDTs   
# (note: BDT training does not need to be done for all observing periods)
	  if [[ $RUN == "TRAIN" ]]
	  then
	    echo "$AZ " 
	     ./CTA.TMVA.sub_train.sh $ARRAY $OFFAXIS $S$M $PARA $QSUBOPT $AZ
##########################################
# IRFs: angular resolution
	  elif [[ $RUN == "ANGRES" ]]
	  then
            ./CTA.EFFAREA.subAllParticle_analyse.sh $ARRAY ANASUM.GammaHadron.TMVAFixedSignal $PARA AngularResolution $S$M 2 $QSUBOPT $AZ
##########################################
# IRFs: effective areas after quality cuts
	  elif [[ $RUN == "QC" ]]
	  then
	    ./CTA.EFFAREA.subAllParticle_analyse.sh $ARRAY ANASUM.GammaHadron.QC $PARA QualityCuts001CU $S$M 0 $QSUBOPT $AZ
##########################################
# IRFs: particle number files
	  elif [[ $RUN == "PARTFIL" ]]
	  then
	     ./CTA.ParticleRateWriter.sub.sh $ARRAY $EFFFULLDIR/QualityCuts001CU $OFFAXIS $ID $EFFFULLDIR/AngularResolution $QSUBOPT
##########################################
# IRFs: effective areas after gamma/hadron cuts
	  elif [[ $RUN == "CUTS" ]]
	  then
	    ./CTA.EFFAREA.subAllParticle_analyse.sh $ARRAY ANASUM.GammaHadron.TMVA $PARA BDT.${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
##########################################
# CTA WP Phys files
	  elif [[ $RUN == "PHYS" ]]
	  then
             if [[ $OFFAXIS == "cone" ]]
             then
                ./CTA.WPPhysWriter.sub.sh $ARRAY $EFFFULLDIR/BDT.${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.Erec$EREC.${EFFVERSION}.ID$ID$AZ$NTYPF.$S$M 1 $ID $S$M $QSUBOPT
             else
                ./CTA.WPPhysWriter.sub.sh $ARRAY $EFFFULLDIR/BDT.${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.Erec$EREC.${EFFVERSION}.ID$ID$AZ$NTYPF.$S$M 0 $ID $S$M $QSUBOPT
             fi
# unknown run set
	  elif [[ $RUN != "EVNDISP" ]]
	  then
	      echo "Unknown run set: $RUN"
	      exit
	  fi
      done
     done
   done
   echo 
   echo "(end of script)"
done
