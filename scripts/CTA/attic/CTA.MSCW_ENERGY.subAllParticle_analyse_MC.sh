#!/bin/bash
#
# submit jobs to analyse MC files with lookup tables for all particle types
# divide large set of e.g. proton simulations into smaller sets
#
##############################################################################

TAB=$1
RECID=$2
ARRAY=$3
PART=$4
MET=$5

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ]
then
   echo ""
   echo "./CTA.MSCW_ENERGY.subAllParticle_analyse_MC.sh <tablefile> <recid> <subarray list> <data set> <script input parameter file>"
   echo "(table files without .root)"
   echo ""
   echo "submit jobs in paralell to analyse MC files with lookup tables"
   echo
   echo "  <tablefile>     table file name (without .root)"
   echo "  <recid>         reconstruction ID"
   echo "  <data set>      e.g. ultra, ISDC3700m, ..."
   echo "  <script input parameter file>  file with directories, etc.; see example in"
   echo "                             $CTA_EVNDISP_AUX_DIR/ParameterFiles/scriptsInput.runparameter"
   echo ""
   echo "  (note that there are a few hardcoded values in this scripts"
   echo
   exit
fi

################################################
# prod2 analysis
if [ $4 = "cta-ultra3" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
elif [ $4 = "v_leeds" ]
then
   VPART=( "proton" )
elif [ $4 = "DESY3700m" ]
then
   VPART=( "proton" )
elif [ $4 = "ISDC3700m" ]
then
   VPART=( "gamma_onSource" "electron" "proton" "proton" "proton" )
################################################
# prod2 analysis
#
# Aar
elif [ $4 = "prod2-Aar-North" ]
then
   VPART=( "gamma_onSource" "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "10" "30" "40" "41" "15" "" )
elif [ $4 = "prod2-Aar-South" ]
then
   VPART=( "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "10" "301" "302" "303" "20" "23" "40" "15" "" )
elif [ $4 = "prod2-Aar-NS" ]
then
   VPART=( "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_cone" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "10" "300" "301" "302" "303" "20" "23" "40" "41" "15" "" )
# Leoncito
elif [ $4 = "prod2-Leoncito-North" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
elif [ $4 = "prod2-Leoncito-South" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
elif [ $4 = "prod2-Leoncito-NS" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" )
elif [ $4 = "prod2-G-Leoncito-North" ] || [ $4 = "prod2-G-Leoncito-South" ]
then
   VPART=( "proton" )
# SAC084
elif [ $4 = "prod2-SAC084-North" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "" "1" "2" "" "" )
elif [ $4 = "prod2-SAC084-South" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "" "1" "2" "" "" )
elif [ $4 = "prod2-SAC084-NS" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "" "1" "2" "" "" )
# SAC100
elif [ $4 = "prod2-SAC100-North" ]
then
   VPART=( "gamma_onSource" "gamma_onSource" "gamma_cone" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "5" "7" "51" "59" "7" "" "" )
elif [ $4 = "prod2-SAC100-South" ]
then
   VPART=( "gamma_onSource" "gamma_onSource" "gamma_onSource" "gamma_cone" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "5" "6" "7" "5" "6" "7" "" "" )
elif [ $4 = "prod2-SAC100-NS" ]
then
   VPART=( "gamma_onSource" "gamma_cone" "gamma_cone" "electron" "proton" )
   VNUM=( "" "1" "2" "" "" )
fi
NPART=${#VPART[@]}

###########################################
# loop over all particle types
for ((m = 0; m < $NPART; m++ ))
do
   PART=${VPART[$m]}
   NUMM=${VNUM[$m]}

   for (( k = 0; k < 10; k++ ))
   do
      if [ $4 = "v_leeds" ]
      then
	 for (( l = 0; l < 10; l++ ))
	 do
	    echo "v_leeds $k$l"
	    ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $k$l
         done
      elif [ $4 = "ISDC3700m" ]
      then
	 ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 1$k
#################################
# prod2 Aar-North and Aar-South
      elif [ $4 = "prod2-Aar-North" ] || [ $4 = "prod2-Aar-South" ]
      then
         if [ $PART = "electron" ]
	 then
	    ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $k
         elif [ $PART = "proton" ]
	 then
	    for (( l = 20; l <= 26; l++ ))
	    do
	       ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $l$k
	    done
         else
	    ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $NUMM$k
         fi
# prod2 Leoncito
      elif [ $4 = "prod2-Leoncito-North" ] || [ $4 = "prod2-Leoncito-South" ]
      then
	 ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $k
# prod2 Leoncito Grid	 
      elif [ $4 = "prod2-G-Leoncito-North" ] || [ $4 = "prod2-G-Leoncito-South" ]
      then
	 for (( l = 30; l <= 45; l++ ))
	 do
	    ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $l$k
	 done
# all other sites
      else
	 ./CTA.MSCW_ENERGY.sub_analyse_MC.sh $TAB $RECID $ARRAY $PART $4 $5 $k
      fi
   done
done

exit
