#!/bin/sh
#
# analysis submission for production 3 analysis 
#
# this script is optimized for the DESY analysis
#
##############################################


if [ $# -ne 3 ] 
then
   echo 
   echo "./CTA.runProd3Analysis.sh <N/S> <run mode> <recid>"
   echo
   echo "  N=prod3-Southern-Site, N=prod3-Northern-Site"
   echo
   echo "  possible run modes are EVNDISP MAKETABLES ANATABLES TRAIN ANGRES QC PARTFIL CUTS PHYS "
   echo
   echo "  <recids>: 0 = all telescopes, 1 = LSTs, 2 = MSTs, 3 = SSTs, 4 = MSTs+SSTs"
   echo
   exit
fi

# site
P2="$1"
# run mode
RUN="$2"
# reconstruction IDs
RECID=$3
# list of subarrays
SUBARRAY="1"
SUBARRAY="3"
SUBARRAY="1 2 3 4 5"

# should be either onSource or cone
OFFAXIS="onSource"

#####################################
# qsub options (priorities)
#   _M_ = -; _X_ = " "
QSUBOPT="_M_P_X_cta_high"
QSUBOPT="_M_js_X_9"
QSUBOPT="_M_P_X_cta_high_X__M_js_X_9"

#####################################
# output directory for script parameter files
PDIR="$CTA_USER_LOG_DIR/tempRunParameterDir/"
mkdir -p $PDIR

#####################################
# analysis dates and table dates

TDATE="d20151114"
DATE="d20151114"
TMVADATE="d20151114"
EFFDATE="d20151114"
TMVAVERSION="V5"
EFFVERSION="V6"

#####################################
# shower directions
#
# _180deg = south
# _0deg = north
MCAZ=( "_180deg" )
MCAZ=( "_180deg" "_0deg" "" )

#####################################
# sites & sub array lists

##########################################################
# SOUTH
if [[ $P2 == "S" ]]
then
   ARRAY="subArray.prod3.full.list"
   OFFAXIS="cone"
   OFFAXIS="onSource"
   # GRID + HD production
   SITE=( "prod3-paranalG" )
   # full Production
   SITE=( "prod3-paranal" )
   EDM=( "p05-NN" )
else
   echo "error: unknown site; allowed are N or S"
   echo $P2
   exit
fi

#####################################
# particle types
PARTICLE=( "gamma_onSource" "electron" "proton" "gamma_cone" )

#####################################
# cut on number of images
NIMAGESMIN="2"

#####################################
# observing time [h]
OBSTIME=( "50h" "5h" "30m" "10m" "10h" "20h" "100h" "500h" "5m" "1m" "2h" )
OBSTIME=( "50h" "5h" "30m" "2h" "10h" )
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
   echo "SITE: $S $M"
   echo "RUN: $RUN"

#####################################
# loop over all array scalings
   for Y in $SUBARRAY
   do

##########################################
# run eventdisplay
        if [[ $RUN == "EVNDISP" ]]
        then
# loop over all particle types
          for ((i = 0; i < ${#PARTICLE[@]}; i++ ))
          do
                  N=${PARTICLE[$i]}

                  # run list for 20 deg files on lustre (default data sets)
                  LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod3/$S.GRID.${N}_20deg-${Y}.list
                  LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod3/prod3-paranal.GRID.${N}_20deg-${Y}.list
                  LIST=/afs/ifh.de/group/cta/scratch/$USER/LOGS/CTA/runLists/prod3/$S.${N}_20deg-${Y}.list

                  ./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod3.sh $ARRAY $Y $LIST $N $S$M 0 $i $QSUBOPT $TRG
           done
           continue
        fi
##########################################
# for the following: duplicate the array list adding the scaling to array names
        NXARRAY=`cat $ARRAY`
        NFILARRAY=$PDIR/temp.$ARRAY-${Y}.list
        rm -f $NFILARRAY
        touch $NFILARRAY
        for A in $NXARRAY
        do
             echo ${A}-${Y} >> $NFILARRAY
        done
##########################################
# loop over all reconstruction IDs 
# (telescope type dependent subarrays)
        for ID in $RECID
        do
           # directory where all mscw output files are written to
           MSCWSUBDIRECTORY="Analysis-ID$ID-$TDATE"
##########################################
# make (fill) tables
           if [[ $RUN == "MAKETABLES" ]]
           then
                  ./CTA.MSCW_ENERGY.sub_make_tables.sh tables_CTA-$S$M-ID$ID-$TDATE $ID $NFILARRAY $OFFAXIS $S$M $QSUBOPT
                  continue
##########################################
# analyse with lookup tables
           elif [[ $RUN == "ANATABLES" ]]
           then
                  # temporary(?): use for analysis always ID0 tables
                  TABLE="tables_CTA-$S$M-ID0-$TDATE-onAxis"
                  if [[ $OFFAXIS == "cone" ]]
                  then
                      TABLE="tables_CTA-$S$M-ID0-$TDATE"
                  fi
                  echo "Using table $TABLE"
                  ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TABLE $ID $NFILARRAY $S$M $MSCWSUBDIRECTORY $OFFAXIS $QSUBOPT
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
# fill a run parameter file
                  PARA="$PDIR/scriptsInput.prod3.ID${ID}${AZ}.${S}${AZ}${OOTIME}.runparameter"
                  rm -f $PARA
                  touch $PARA
                  echo "WRITING PARAMETERFILE $PARA"
                  NTYPF=NIM$NIMAGESMIN
                  EFFDIR="EffectiveArea-"$OOTIME"-ID$ID$AZ-$NTYPF-$EFFDATE"
                  EFFFULLDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$S$M/$EFFDIR/"
                  echo "MSCWSUBDIRECTORY $MSCWSUBDIRECTORY" >> $PARA
                  echo "TMVASUBDIR BDT-${TMVAVERSION}-ID$ID$AZ-$NTYPF-$TMVADATE" >> $PARA
                  echo "EFFAREASUBDIR $EFFDIR" >> $PARA
                  echo "RECID $ID" >> $PARA
                  echo "NIMAGESMIN $NIMAGESMIN" >> $PARA
                  echo "OBSERVINGTIME_H $OOTIME" >> $PARA
                  echo "GETXOFFYOFFAFTERCUTS yes" >> $PARA
##########################################
# train BDTs   
# (note: BDT training does not need to be done for all observing periods)
                  if [[ $RUN == "TRAIN" ]]
                  then
                     if [ ${o} -eq 0 ]
                     then
                         ./CTA.TMVA.sub_train.sh $NFILARRAY $OFFAXIS $S$M $PARA $QSUBOPT $AZ
                     fi
##########################################
# IRFs: angular resolution
                  elif [[ $RUN == "ANGRES" ]]
                  then
                    ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron.TMVAFixedSignal $PARA AngularResolution $S$M 2 $QSUBOPT $AZ
##########################################
# IRFs: effective areas after quality cuts
                  elif [[ $RUN == "QC" ]]
                  then
                    ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron.QC $PARA QualityCuts001CU $S$M 3 $QSUBOPT $AZ
##########################################
# IRFs: particle number files
                  elif [[ $RUN == "PARTFIL" ]]
                  then
                     ./CTA.ParticleRateWriter.sub.sh $NFILARRAY $EFFFULLDIR/QualityCuts001CU cone $ID $EFFFULLDIR/AngularResolution $QSUBOPT
                     ./CTA.ParticleRateWriter.sub.sh $NFILARRAY $EFFFULLDIR/QualityCuts001CU onSource $ID $EFFFULLDIR/AngularResolution $QSUBOPT
##########################################
# IRFs: effective areas after gamma/hadron cuts
                  elif [[ $RUN == "CUTS" ]]
                  then
                    ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron.TMVA $PARA BDT.${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
##########################################
# CTA WP Phys files
                  elif [[ $RUN == "PHYS" ]]
                  then
                     if [[ $OFFAXIS == "cone" ]]
                     then
                        ./CTA.WPPhysWriter.sub.sh $NFILARRAY $EFFFULLDIR/BDT.${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$NTYPF.$S$M 1 $ID $S$M $QSUBOPT
                     else
                        ./CTA.WPPhysWriter.sub.sh $NFILARRAY $EFFFULLDIR/BDT.${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$NTYPF.$S$M 0 $ID $S$M $QSUBOPT
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
   done
   echo 
   echo "(end of script)"
done
