#!/bin/bash
#
# submit jobs for effective area calculation
#
# particle names and directories fixed by CTA setup
#
#

if [ $# -lt 6 ] 
then
   echo ""
   echo "./CTA.EFFAREA.subAllParticle_analyse.sh <subarray list> <cutfile template> <analysis parameter file> <output subdirectory> <data set> [filling mode] [qsub options] [direction (e.g. _180deg)]"
   echo
   echo "<subarray list>"
   echo "     text file with list of subarray IDs"
   echo
   echo "<cutfile template>"
   echo "     template for gamma/hadron cut file"
   echo "     (suffix must be .gamma/.CRbck ; this will be added by this script)"
   echo "     examples can be found in $CTA_EVNDISP_AUX_DIR/GammaHadronCutFiles"
   echo 
   echo "<analysis parameter file>"
   echo "     file with analysis parameter"
   echo "     examples can be found in $CTA_EVNDISP_AUX_DIR/ParameterFiles/"
   echo
   echo "<output subdirectory>"
   echo "    sub-directory name (not the full path) for effective areas files"
   echo
   echo "<data set>         e.g. cta-ultra3, ISDC3700m, ...  "
   echo
   echo "[filling mode]"
   echo "     effective area filling mode (use 2 to calculate angular resolution only)"
   echo "     (no value for default IRF calculation)"
   echo
   echo "[direction]        e.g. for north: \"_180deg\", for south: \"_0deg\", for all directions: no option"
   echo
   echo ""
   exit
fi

SUBAR=$1
CDIR="$CTA_EVNDISP_AUX_DIR/GammaHadronCutFiles/"
CFIL=$2
ANAPAR=$3
ODIR=$4
DSET=$5
GMOD=0
if [ -n "$6" ]
then
  GMOD=$6
fi
MCAZ=""
if [ -n "$8" ]
then
  MCAZ=$8
fi
QSUBOPT=""
if [ -n $7 ]
then
   QSUBOPT="$7"
fi
#######################################
# read values from parameter file
if [ ! -e $ANAPAR ]
then
  echo "error: analysis parameter file not found: $ANAPAR" 
  exit
fi
echo "reading analysis parameter from $ANAPAR"
# get output file directory
EFFAREADIR=`grep EFFAREASUBDIR $ANAPAR | awk {'print $2'}`
ODIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$EFFAREADIR/$4/"
mkdir -v -p $ODIR
# get reconstruction ID
RECID=`grep RECID $ANAPAR | awk {'print $2'}`

#arrays
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $SUBAR`

#################################################
# set particle types 
# (don't expect to have cone for all data sets)
if [ $GMOD = "0" ]
then
   if [ $DSET = "cta-ultra3" ] || [ $DSET = "v_leeds" ] || [[ $DSET = *prod2* ]] || [[ $DSET = prod1-cta-ul* ]] || [[ $DSET = *prod3* ]]
   then
      VPART=( "gamma_onSource" "gamma_cone" "electron" "proton" "electron_onSource" "proton_onSource" )
   elif [ $DSET = "VTS" ]
   then
     VPART=( "gamma_onSource" "proton" )
   else
      VPART=( "gamma_onSource" "electron_onSource" "proton_onSource" )
   fi
else
   if [ $DSET = "cta-ultra3" ] || [ $DSET = "v_leeds" ] || [[ $DSET = *prod2* ]] || [[ $DSET = prod1-cta-ul* ]] || [[ $DSET = *prod3* ]] 
   then
      VPART=( "gamma_onSource" "gamma_cone" )
   else
      VPART=( "gamma_onSource" )
   fi
fi
NPART=${#VPART[@]}

#########################################
# loop over all arrays
#########################################
for ARRAY in $VARRAY
do
   echo "STARTING ARRAY $ARRAY"

###########################################
# loop over all particle types
   for ((m = 0; m < $NPART; m++ ))
   do
      PART=${VPART[$m]}
      echo "    MC PARTICLE TYPE $PART"

###########################################
# prepare cut file
      CCUT=$ODIR/$CFIL.$PART.$ARRAY.dat
      if [ $PART = "gamma_onSource" ] || [ $PART = "gamma_cone" ]
      then
        cp $CDIR/$CFIL.gamma.dat $CCUT
      fi
      if [ $PART = "proton" ] || [ $PART = "electron" ] || [ $PART = "helium" ] || [ $PART = "electron_onSource" ] || [ $PART = "proton_onSource" ]
      then
        cp $CDIR/$CFIL.CRbck.dat $CCUT
      fi

###########################################
# submit the job script
      ./CTA.EFFAREA.sub_analyse.sh $ARRAY $RECID $PART $CCUT $ANAPAR $ODIR $DSET $GMOD \"$QSUBOPT\" $MCAZ
   done
done

exit
