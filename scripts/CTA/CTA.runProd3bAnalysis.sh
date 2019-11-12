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
   echo "./CTA.runProd3bAnalysis.sh <S/S40deg> <run mode> <recid> [min number of LSTs] [min number of MSTs] [min number of SSTs] [min number of SCMSTs]"
   echo
   echo "  S=prod3b-Southern-Site (20 deg ze), S40deg/S60deg=prod3b-Southern-Site (40/60 deg ze)"
   echo "  SCT=prod3b-Southern-Site (20 deg ze) SCT MSTs"
   echo "  MA=prod3b-Southern-Site (20 deg ze) SST miniarray"
   echo "  N=prod3b-Northern Site (20deg ze), N40deg"
   echo "  GRN20deg=GRID baseline production (North)" 
   echo "  S20degNSB01, N20degNSB05, NSB30 NSB production"
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
# (not relevant for prod3 analysis)
SUBARRAY="1"
# number of telescopes
LST=2
MST=2
SST=2
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
   SST=$6
fi
if [ -z $SST ]
then
   SST=2
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
echo "Telescope multiplicities: LST $LST MST $MST SST $SST SCMST $SCMST"

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

TDATE="d20160925m4"
# following is used for MA array analysis
# TDATE="d20160925m2" 
# updated disp method (disp100 Tel, disp error, and LT Energy)
# TMVADATE="d20161127LTE"
# EFFDATE="d20161127LTE"

##########
# default values for Jan 2017
# updated disp method (disp100 Tel, disp error, and dispEnergy)
# ANADATE="d20161127"
# TMVADATE="d20161128"
# EFFDATE="d20170120"

##########
# default value for Feb 2017
# updated disp energy method using median
# ANADATE="d20170217"
# TMVADATE="d20170217"
# EFFDATE="d20170217"

##########
# default value for Feb 2017
# updated disp energy method using median;
# using disp trained per telescope
# ANADATE="d20170224"
# TMVADATE="d20170224"
# EFFDATE="d20170224"

##########
# default value for June 2017
# bug fixes in mscw_energy
ANADATE="d20170627"
TMVADATE="d20170627"
EFFDATE="d20170627"
TDATE="d20160925m4"

##########
# training with weighted disp method
# ANADATE="d20170611"
# TMVADATE="d20170611"
# EFFDATE="d20170611"

##########
# FOV updates for SSTs"
ANADATE="d20171007"
TMVADATE="d20171007"
EFFDATE="d20171007"

###########
# dispBDT and BDT updates
ANADATE="d20171208T19"
TMVADATE="d20171208T19"
EFFDATE="d20171208T19"
TDATE="d20171007m4"

#############
# 'final' 2017 analysis
ANADATE="d20171231T19iPWW"
TMVADATE="d20171231T19"
EFFDATE="d20171231T19"
TDATE="d20171231T19"

#############
#  'final' 2018 analysis
TDATE="d20180113"
ANADATE="d20180113"
TMVADATE="d20180113FT"
TMVADATE="d20180113"
EFFDATE="d20180113"
EFFDATE="d20180310"
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
MCAZ=( "_0deg" )
MCAZ=( "_180deg" "_0deg" "" )
MCAZ=( "_180deg" )

#####################################
# loss cut adaption
EDM=( "q05b-NN" )
EDM=( "s05b-NN" )

##########################################################
# SOUTH
ARRAY="subArray.prod3b.HB9.list.aswg.FD"
ARRAY="subArray.prod3b.HB9.list.aswg"
ARRAY="subArray.prod3b.HB9.list.aswg.HB10"
ARRAY="subArray.prod3b.HB9.FOV.list"
ARRAY="subArray.prod3b.S60deg.list"
if [[ $P2 == "S" ]] || [[ $P2 == "S20deg" ]]
then
   SITE=( "prod3b-paranal20deg" )
elif [[ $P2 == "S40deg" ]]
then
   SITE=( "prod3b-paranal40deg" )
elif [[ $P2 == "S60deg" ]]
then
   SITE=( "prod3b-paranal60deg" )
# SOUTH SCTs
elif [[ $P2 == "SCT" ]]
then
   SITE=( "prod3b-paranal20degSCT" );
   ARRAY="subArray.prod3Sb.SCT.list"
   EDM=( "-NN" )
# SOUTH Miniarray
elif [[ $P2 == "MA" ]]
then
   SITE=( "prod3b-paranal20deg" )
   ARRAY="subArray.prod3b.HB9.MA.list"
# NORTH
elif [[ $P2 == "N20deg" ]]
then
# prod3 QGSJet-old
   SITE=( "prod3b-LaPalma-NN" )
# prod3 QGSJet-II
   SITE=( "prod3b-LaPalma-20deg" )
   ARRAY="subArray.prod3Nb-aswg.list"
elif [[ $P2 == "N40deg" ]]
then
   SITE=( "prod3b-LaPalma40deg-NN" )
   SITE=( "prod3b-LaPalma-40deg" )
   ARRAY="subArray.prod3Nb.list"
elif [[ $P2 == "prod3N" ]]
then
   ARRAY="subArray.prod3N-N.list"
   ARRAY="subArray.prod3N-N.list.aswg"
   SITE=( "prod3-LaPalma" )
   EDM=( "p05-TS" )
elif [[ $P2 == "GRN20deg" ]]
then
   SITE=( "GR_prod3b-LaPalma-20deg" )
   ARRAY="subArray.prod3Nb-GR.list"
   EDM=( "p05-NN" )
# **preli** South/North needs edition by hand"
elif [[ $P2 == *"NSB"* ]]
then
   ARRAY="subArray.prod3Nb.list.BL"
   ARRAY="subArray.prod3Nb.list.BL.TL"
   SITE=( "prod3b-LaPalma20deg-$P2" );
#   SITE=( "prod3b-paranal20deg-$P2" )
#   ARRAY="subArray.prod3b.HB9.list.aswg.FD"
# s05-NN: analysis using local evndisp files (NSB01 ped file)
   EDM=( "s05-NN" )
# s05-NN: analysis using local evndisp files (NSB05 ped file)
   EDM=( "s05-HI" )
# x05-NN: analysis using GRID evndisp files
   EDM=( "x05-NN" )
elif [[ $P2 == "PE20deg" ]]
then
   SITE=( "prod3b-paranal20degPE" )
   ARRAY="subArray.prod3b.HB9.list.aswg.NSB"
else
   echo "error: unknown site; allowed are N or S/S40deg/S60deg"
   echo $P2
   exit
fi
# should be either onSource or cone (default is cone)
OFFAXIS="cone"

#####################################
# particle types
PARTICLE=( "gamma_onSource" )
PARTICLE=( "gamma_cone" "gamma_onSource" "electron" "proton" )

#####################################
# cut on number of images
NIMAGESMIN="2"

#####################################
# observing time [h]
OBSTIME=( "50h" "5h" "30m" "10m" "10h" "20h" "100h" "500h" "5m" "1m" "2h" )
OBSTIME=( "5h" "30m" "100s" )
OBSTIME=( "50h" )
OBSTIME=( "50h" "5h" "30m" "100s" )

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
# (not relevant for prod3b production)
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
                  LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/${S}p05-NN/GRID/FileLists/$S.$N.simtel.list
                  LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/${S}r05-NN/GRID/FileLists/$N.list
                  LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/${S}r05-NN/GRID/FileLists/${N}_180deg.list

                  LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/FileLists/prod3b-paranal20deg/prod3b-paranal20deg-${N}.20deg_0deg.list
                  LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/FileLists/prod3b-paranal20deg/${N}_180deg.list
                  # NSB lists
                  # LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/FileLists/prod3b-LaPalma-20deg-NSB1x/Prod3_LaPalma_Baseline_NSB1x_${N}_North_20deg_DL0.local.list
                  # LIST=$CTA_USER_DATA_DIR/analysis/AnalysisData/FileLists/prod3b-LaPalma-20deg-NSB1x/Prod3_LaPalma_Baseline_NSB1x_${N}_South_20deg_DL0.local.list
                  # Paranal lists for average PE calculation
                  # LIST=/lustre/fs21/group/cta/prod3b/prod3b_highNSB_paranal/NSBx1/DL0/tt.list

                  echo "READING SIMTEL FILE LIST $LIST"

                  ./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod3b.sh $ARRAY $LIST $N $S$M 0 $i $QSUBOPT $TRG
           done
           continue
        fi
##########################################
# for the following: duplicate the array list adding the scaling to array names
        NXARRAY=`cat $ARRAY`
        NFILARRAY=$PDIR/temp.$ARRAY.list
        rm -f $NFILARRAY
        touch $NFILARRAY
        for A in $NXARRAY
        do
             echo ${A} >> $NFILARRAY
        done
##########################################
# dispBDT training
        if [[ $RUN == "DISPBDT" ]]
        then
            BDTDIR="BDTdisp."
            RUNPAR="$CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDTDisp.runparameter"
            DDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$S$M/"
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
                              ./CTA.MSCW_ENERGY.sub_make_tables.sh $TABLE $ID $NFILARRAY $OFFAXIS $S$M ${AZ} $QSUBOPT
                              continue
    ##########################################
# analyse with lookup tables
                       elif [[ $RUN == "ANATABLES" ]]
                       then
                              echo "Using table $TABLE"
                              ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TABLE $ID $NFILARRAY $S$M $MSCWSUBDIRECTORY $OFFAXIS ${AZ} $QSUBOPT
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
                  ETYPF=NIM${NIMAGESMIN}LST${LST}MST${MST}SST${SST}SCMST${SCMST}
                  TMVATYPF=$ETYPF
                  # paranal mva are named differently
                  if [[ $SITE == *"paranal"* ]] && [[ $SITE != *"SCT"* ]]
                  then
                      TMVATYPF=NIM${NIMAGESMIN}LST${LST}MST${MST}SST${SST}
                  fi
                  PARA="$PDIR/scriptsInput.prod3b.ID${ID}${ETYPF}${AZ}.${S}${AZ}${OOTIME}.runparameter"
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
                  echo "NSST $SST" >> $PARA
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
                    if [[ "$OOTIME" = "50h" ]] && [[ "$MST" -ge "4" ]]
                    then
                        ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron.QC $PARA QualityCuts001CU $S$M 3 $QSUBOPT $AZ
                    # min angle cut depends on observation time: 50h stricter, 5h and and 30 min more relaxed
                    # (never done for 50h observation, as those are expected to require higher resolution)
                    else
                        ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron008deg.QC $PARA QualityCuts001CU $S$M 3 $QSUBOPT $AZ
                    fi
##########################################
# IRFs: effective areas after gamma/hadron cuts
                  elif [[ $RUN == "CUTS" ]]
                  then
                    # large multiplicity runs use 80% max signal efficiency (best resolution)
                    if [[ "$OOTIME" = "50h" ]] && [[ "$MST" -ge "4" ]]
                    then
                        ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron.TMVA $PARA BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
                    # low multiplicity runs use 95% max signal efficiency (lower requirements on resolution)
                    else
                        ./CTA.EFFAREA.sub_analyse_list.sh $NFILARRAY ANASUM.GammaHadron95p008deg.TMVA $PARA BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $S$M 0 $QSUBOPT $AZ
                    fi
##########################################
# CTA WP Phys files
                  elif [[ $RUN == "PHYS" ]]
                  then
                     if [[ $OFFAXIS == "cone" ]]
                     then
                        ./CTA.WPPhysWriter.sub.sh $NFILARRAY $EFFFULLDIR/BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$ETYPF.$S$M 1 $ID $S$M $BFINEBINNING $QSUBOPT
                     else
                        ./CTA.WPPhysWriter.sub.sh $NFILARRAY $EFFFULLDIR/BDT."$OOTIME"-${EFFVERSION}.$EFFDATE $OOTIME DESY.$EFFDATE.${EFFVERSION}.ID$ID$AZ$ETYPF.$S$M 0 $ID $S$M $BFINEBINNING $QSUBOPT
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
