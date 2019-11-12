#!/bin/sh
#
# analysis submission for production 3b analysis 
#
# this script is optimized for the DESY analysis
#
##############################################


if [ $# -le 3 ] 
then
   echo 
   echo "./CTA.runProd4Analysis.sh <S> <run mode> <recid> [min number of LSTs] [min number of MSTs] [min number of SSTs] [min number of SCMSTs]"
   echo
   echo "  S=prod4-Southern-Site (20 deg ze)"
   echo
   echo "  possible run modes are EVNDISP MAKETABLES DISPBDT ANATABLES TRAIN ANGRES QC CUTS PHYS "
   echo
   echo "  <recids>: 0 = all telescopes, 1 = LSTs, 2 = MSTs, 3 = SSTs, 4 = MSTs+SSTs, 5 = LSTs+MSTs"
   echo
   exit
fi

# site
P2="$1"
# run mode
RUN="$2"
# reconstruction IDs
RECID=$3
# list of subarrays (scalings)
SUBARRAY="1"
# number of telescopes
LST=2
MST=2
NUMCUTSST=2
SCMST=2
if [ -n $4 ]
then
   LST=$4
fi
if [ -z $LST ]
then
   LST=2
fi
if [ -n $5 ]
then
   MST=$5
fi
if [ -z $MST ]
then
   MST=2
fi
if [ -n $6 ]
then
   NUMCUTSST=$6
fi
if [ -z $NUMCUTSST ]
then
   NUMCUTSST=2
fi
if [ -n $7 ]
then
   SCMST=$7
fi
# no number of SCMSTs given:
if [ -z $SCMST ]
then
   SCMST=2
fi
echo "Telescope multiplicities: LST $LST MST $MST SST $NUMCUTSST SCMST $SCMST"

#####################################
# qsub options (priorities)
#   _M_ = -; _X_ = " "
QSUBOPT="_M_P_X_cta_high_X__M_js_X_99"

#####################################
# output directory for script parameter files
PDIR="$CTA_USER_LOG_DIR/tempRunParameterDir/"
mkdir -p $PDIR

#####################################
# analysis dates and table dates

#############
#  prod4-MST
TDATE="d20190506"
ANADATE="d20190506"
TMVADATE="d20190506"
EFFDATE="d20190506"
# off-axis binnning
BFINEBINNING="TRUE"
BFINEBINNING="FALSE"
if [ $BFINEBINNING = "TRUE" ]
then
   EFFDATE=${EFFDATE}FB
fi

TMVAVERSION="V2"
EFFVERSION="V2"
#####################################


#####################################
# shower directions
#
# _180deg = south
# _0deg = north
MCAZ=( "_180deg" "_0deg" "" )
# prod4-MST only for North direction
MCAZ=( "_0deg" )

#####################################
# loss cut adaption
EDM=( "s05b-NN" )

#####################################
# SST type
SST=( "sst-astri" "sst-gct" "sst-1m" "sst-gct-7mm" "sst-astri+chec-s" "sst-astri+chec-s-7mm" )
SST=( "sst-astri+chec-s" "sst-1m" )
SST=( "sst-gct" )
SST=( "sst-1m" )
SST=( "sst-astri+chec-s" )
SST=( "mst-f")

##########################################################
# SOUTH
SITE=( "prod4-SST-paranal20deg" )

# should be either onSource or cone (default is cone)
OFFAXIS="cone"

#####################################
# particle types
PARTICLE=( "electron" )
PARTICLE=( "proton" )
PARTICLE=( "gamma_cone" )
PARTICLE=( "gamma_onSource" "electron" "proton" )
PARTICLE=( "gamma_cone" "gamma_onSource" "electron" "proton" )
PARTICLE=( "gamma_onSource" )

#####################################
# cut on number of images
NIMAGESMIN="2"

#####################################
# observing time [h]
OBSTIME=( "50h" "5h" "30m" "100s" )
OBSTIME=( "50h" )
#####################################
# loop over all sites
NSST=${#SST[@]}
for (( m = 0; m < $NSST; m++ ))
do
   # SST type
   T=${SST[$m]}
   # site name
   P=${SITE[0]}
   S=${SITE[0]}-${T}-
   # eventdisplay analysis method
   M=${EDM[0]}

   echo
   echo "======================================================================================"
   echo "SST: $T SITE: $S $M RUN: $RUN"

    ARRAY="subArray.prod4-MST.list"
    if [[ $T == "sst-astri" ]]
    then
       ARRAY="subArray.prod4-SST-A.list"
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
          LIST="/lustre/fs21/group/cta/users/maierg/analysis/AnalysisData/FileList_prod4/${P}/list_${N}_${T}_good.txt"
          LIST="/lustre/fs21/group/cta/users/maierg/analysis/AnalysisData/FileList_prod4b/${P}/${T}_${N}.list"
          LIST="/lustre/fs21/group/cta/users/maierg/analysis/AnalysisData/FileList_prod4/prod4-MST-paranal20deg/${T}_${N}.list"
          echo "READING SIMTEL FILE LIST $LIST"

           ./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod4.sh $ARRAY $LIST $N $S$M 0 $i $QSUBOPT $TRG
      done
      continue
   fi
##########################################
# dispBDT training
    if [[ $RUN == "DISPBDT" ]]
    then
        BDTDIR="BDTdisp."
        RUNPAR="$CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDTDisp.runparameter"
        DDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$S$M/"
        NXARRAY=`cat $ARRAY`
        for A in $NXARRAY
        do
            ./CTA.DISPTRAINING.sub_analyse.sh ${S}${M} $DDIR/${BDTDIR}${A} 0 $A $RUNPAR 99
        done
        continue
    fi
##########################################
# loop over all reconstruction IDs 
# (telescope type dependent subarrays)
    for ID in $RECID
    do
       # directory where all mscw output files are written to
       MSCWSUBDIRECTORY="Analysis-ID$ID-${ANADATE}"
# loop over all shower directions 
# (i.e. North and South)
        for ((a = 0; a < ${#MCAZ[@]}; a++ ))
        do
              AZ=${MCAZ[$a]}
              if [ "$AZ" ] 
              then
##########################################
# make (fill) tables
                  TABLE="tables_CTA-$S$M-ID0${AZ}-$TDATE"
                  if [[ $RUN == "MAKETABLES" ]]
                  then
                          echo "Filling table $TABLE"
                          ./CTA.MSCW_ENERGY.sub_make_tables.sh $TABLE $ID $ARRAY $OFFAXIS $S$M ${AZ} $QSUBOPT
                          continue
 ##########################################
# analyse with lookup tables
                   elif [[ $RUN == "ANATABLES" ]]
                   then
                          echo "Using table $TABLE"
                          ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TABLE $ID $ARRAY $S$M $MSCWSUBDIRECTORY $OFFAXIS ${AZ} $QSUBOPT
                          continue
                    fi
               fi
        done
        if [[ $RUN == "MAKETABLES" ]] || [[ $RUN == "ANATABLES" ]]
        then
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
                  ETYPF=NIM${NIMAGESMIN}LST${LST}MST${MST}SST${NUMCUTSST}SCMST${SCMST}
                  TMVATYPF=$ETYPF
                  # paranal mva are named differently
                  if [[ $SITE == *"paranal"* ]] && [[ $SITE != *"SCT"* ]]
                  then
                      TMVATYPF=NIM${NIMAGESMIN}LST${LST}MST${MST}SST${NUMCUTSST}
                  fi
                  PARA="$PDIR/scriptsInput.prod4.ID${ID}${ETYPF}${AZ}.${S}${AZ}${OOTIME}.runparameter"
                  rm -f $PARA
                  touch $PARA
                  echo "WRITING PARAMETERFILE $PARA"
                  EFFDIR="EffectiveArea-"$OOTIME"-ID$ID$AZ-$ETYPF-$EFFDATE-$EFFVERSION"
                  EFFFULLDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$S$M/EffectiveAreas/$EFFDIR/"
                  echo "MSCWSUBDIRECTORY $MSCWSUBDIRECTORY" >> $PARA
                  echo "TMVASUBDIR BDT-${TMVAVERSION}-ID$ID$AZ-$TMVATYPF-$TMVADATE" >> $PARA
                  echo "EFFAREASUBDIR $EFFDIR" >> $PARA
                  echo "RECID $ID" >> $PARA
                  echo "NIMAGESMIN $NIMAGESMIN" >> $PARA
                  echo "NLST $LST" >> $PARA
                  echo "NMST $MST" >> $PARA
                  echo "NSST $NUMCUTSST" >> $PARA
                  echo "NSCMST $SCMST" >> $PARA
                  echo "OBSERVINGTIME_H $OOTIME" >> $PARA
                  echo "GETXOFFYOFFAFTERCUTS yes" >> $PARA
                  echo "OFFAXISFINEBINNING $BFINEBINNING" >> $PARA
##########################################
# train BDTs   
# (note: BDT training does not need to be done for all observing periods)
                  if [[ $RUN == "TRAIN" ]]
                  then
                     if [ ${o} -eq 0 ]
                     then
                         ./CTA.TMVA.sub_train.sh $ARRAY $OFFAXIS $S$M $PARA $QSUBOPT $AZ
                     fi
##########################################
# IRFs: angular resolution
                  elif [[ $RUN == "ANGRES" ]]
                  then
                    ./CTA.EFFAREA.sub_analyse_list.sh $ARRAY ANASUM.GammaHadron.TMVAFixedSignal $PARA AngularResolution $S$M 2 $QSUBOPT $AZ
##########################################
# IRFs: effective areas after quality cuts
                  elif [[ $RUN == "QC" ]]
                  then
                    if [[ "$OOTIME" = "50h" ]] && [[ "$MST" -ge "4" ]]
                    then
                        ./CTA.EFFAREA.sub_analyse_list.sh $ARRAY ANASUM.GammaHadron.QC $PARA QualityCuts001CU $S$M 3 $QSUBOPT $AZ
                    # min angle cut depends on observation time: 50h stricter, 5h and and 30 min more relaxed
                    # (never done for 50h observation, as those are expected to require higher resolution)
                    else
                        ./CTA.EFFAREA.sub_analyse_list.sh $ARRAY ANASUM.GammaHadron008deg.QC $PARA QualityCuts001CU $S$M 3 $QSUBOPT $AZ
                    fi
##########################################
# IRFs: effective areas after gamma/hadron cuts
                  elif [[ $RUN == "CUTS" ]]
                  then
                    # large multiplicity runs use 80% max signal efficiency (best resolution)
                    if [[ "$OOTIME" = "50h" ]] && [[ "$MST" -ge "4" ]]
                    then
                        ./CTA.EFFAREA.sub_analyse_list.sh $ARRAY ANASUM.GammaHadron.TMVA $PARA BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
                    # low multiplicity runs use 95% max signal efficiency (lower requirements on resolution)
                    else
                        ./CTA.EFFAREA.sub_analyse_list.sh $ARRAY ANASUM.GammaHadron95p008deg.TMVA $PARA BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
                    fi
##########################################
# CTA WP Phys files
                  elif [[ $RUN == "PHYS" ]]
                  then
                     if [[ $OFFAXIS == "cone" ]]
                     then
                        ./CTA.WPPhysWriter.sub.sh $ARRAY $EFFFULLDIR/BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$ETYPF.$S$M 1 $ID $S$M $BFINEBINNING $QSUBOPT
                     else
                        ./CTA.WPPhysWriter.sub.sh $ARRAY $EFFFULLDIR/BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$ETYPF.$S$M 0 $ID $S$M $BFINEBINNING $QSUBOPT
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
