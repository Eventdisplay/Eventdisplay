#!/bin/sh
#
# script to train cuts/MVAs with TMVA
#
#
#

if [ $# -lt 4 ]
then
   echo
   echo "CTA.TMVA.sub_train.sh <subarray list> <onSource/cone> <data set> <analysis parameter file> [qsub options] [direction (e.g. _180deg)]"
   echo ""
   echo "  <subarray list>   text file with list of subarray IDs"
   echo
   echo "  <onSource/cone>    calculate tables for on source or different wobble offsets"
   echo
   echo "  <data set>         e.g. cta-ultra3, ISDC3700, ...  "
   echo 
   echo "  <direction>        e.g. for north: \"_180deg\", for south: \"_0deg\", for all directions: no option"
   echo
   echo "   note 1: keywords ENERGYBINS and OUTPUTFILE are ignored in the runparameter file"
   echo
   echo "   note 2: energy and wobble offset bins are hardwired in this scripts"
   echo
   echo "   note 3: adjust h_cpu depending on your MVA method"
   echo
   echo "   note 4: default TMVA parameter file is $CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDT.runparameter"
   echo
   exit
fi

#######################################
# read values from parameter file
ANAPAR=$4
if [ ! -e $ANAPAR ]
then
  echo "error: analysis parameter file not found: $ANAPAR" 
  exit
fi
echo "reading analysis parameter from $ANAPAR"
NIMAGESMIN=`grep NIMAGESMIN $ANAPAR | awk {'print $2'}`
NCUTLST=`grep NLST $ANAPAR | awk {'print $2'}`
NCUTMST=`grep NMST $ANAPAR | awk {'print $2'}`
NCUTSST=`grep NSST $ANAPAR | awk {'print $2'}`
NCUTMSCT=`grep NSCMST $ANAPAR | awk {'print $2'}`
ANADIR=`grep MSCWSUBDIRECTORY  $ANAPAR | awk {'print $2'}`
DDIR=`grep TMVASUBDIR $ANAPAR | awk {'print $2'}`
RECID=`grep RECID $ANAPAR | awk {'print $2'}`
echo $NIMAGESMIN $ANADIR $DDIR
DSET=$3
OFIL="BDT"
CONE="FALSE"
if [[ $2 == cone ]]
then
  CONE="TRUE"
fi
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $1`

# directory for log files
QLOG=/dev/null

######################################################
# TMVA parameters are detetermined from data set name
RPAR="$CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDT"
RXPAR=`basename $RPAR.runparameter runparameter`
#####################################
if [ -n "$6" ]
then
  MCAZ=$6
fi

if [ -n $5 ]
then
   QSUBOPT="$5"
fi
QSUBOPT=${QSUBOPT//_X_/ } 
QSUBOPT=${QSUBOPT//_M_/-} 
QSUBOPT=${QSUBOPT//\"/} 

#####################################
# energy bins
# EMIN=( -1.90 -1.90 -1.45 -1.20 -1.00 -0.50 0.00 0.50 1.25 )
# EMAX=( -1.40 -1.30 -1.15 -0.80 -0.25  0.25 0.75 1.50 2.50 )
EMIN=( -1.90 -1.90 -1.45 -1.20 -1.00 -0.50 0.00 0.50 1.00 )
EMAX=( -1.40 -1.30 -1.15 -0.80 -0.25  0.25 0.75 1.50 2.50 )
# (until 2016-09-11)
# EMIN=( -1.90 -1.90 -1.45 -1.00 -0.50 0.00 0.50 1.25 )
# EMAX=( -1.40 -1.30 -0.75 -0.25  0.25 0.75 1.50 2.50 )
NENE=${#EMIN[@]}
#####################################
# offset bins 
if [ $CONE == "TRUE" ]
then
   OFFMIN=( 0.0 1.0 2.0 2.5 4.0 5.0 )
   OFFMAX=( 3.0 3.0 3.5 4.5 5.0 6.0 )
   OFFMEA=( 0.5 1.5 2.5 3.5 4.5 5.5 )
   DSUF="gamma_cone"
   GTYPE="cone10_evndisp"
   ASUF="gamma_onSource"
   ATYPE="baseline_evndisp"
else
   OFFMIN=( "0.0" )
   OFFMAX=( "3." )
# value used until 2015-11-09
#   OFFMAX=( "1.e10" )
   OFFMEA=( 0.0 )
   DSUF="gamma_onSource"
   GTYPE="baseline_evndisp"
   ASUF="gamma_cone"
   ATYPE="cone10_evndisp"
fi
NOFF=${#OFFMIN[@]}

######################################
# checking the path for binary
if [ -z $EVNDISPSYS ]
then
    echo "no EVNDISPSYS env variable defined"
    exit
fi

######################################
# log files
DATE=`date +"%y%m%d"`
LDIR=$CTA_USER_LOG_DIR/$DATE/TMVATRAINING/
mkdir -p $LDIR
echo "log directory: " $LDIR

######################################
# script name template
FSCRIPT="CTA.TMVA.qsub_train"

###############################################################
# loop over all arrays
for ARRAY in $VARRAY
do
   echo "STARTING $DSET ARRAY $ARRAY MCAZ $MCAZ"

# signal and background files
# (no electrons are used for the background training)
# ensure mixed training set for the two different pointing directions
# two lists for signal and background, alternating from previous lists
# (list must be sorted; and then mixed)
# Splitmode=BLOCK
# SFIL1, BFIL1 used for training
# SFIL2, BFIL2 used for testing and analysis

# different namings for GRID and local productions
   if [[ ${DSET:0:2} == "GR" ]]
   then
       SFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/gamma*"deg$MCAZ"*"$GTYPE"*.mscw.root | sort -g | awk 'NR % 2 == 1'`
       SFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/gamma*"deg$MCAZ"*"$GTYPE"*.mscw.root | sort -g | awk 'NR % 2 == 0'`
       BFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton*"deg$MCAZ"*mscw.root | sort -g | awk 'NR % 2 == 1'`
       BFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton*"deg$MCAZ"*mscw.root | sort -g | awk 'NR % 2 == 0'`
       GFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/gamma*"deg$MCAZ"*"$ATYPE"*.mscw.root`
       EFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/elec*"deg$MCAZ"*mscw.root`
    elif [[ $DSET == *"NSB"* ]]
    then
       SFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$DSUF*"_ID"$RECID"$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 1'`
       SFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$DSUF*"_ID"$RECID"$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 0'`
       BFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton*"_ID"$RECID"$MCAZ"*mscw.root | sort -g | awk 'NR % 2 == 1'`
       BFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton*"_ID"$RECID"$MCAZ"*mscw.root | sort -g | awk 'NR % 2 == 0'`
       GFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$ASUF*"_ID"$RECID"$MCAZ"*.mscw.root`
       EFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/elec*"_ID"$RECID"$MCAZ"*mscw.root`
    else
       SFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$DSUF."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 1'`
       SFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$DSUF."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 0'`
       BFIL1=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 1'`
       BFIL2=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/proton."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root | sort -g | awk 'NR % 2 == 0'`
       GFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/$ASUF."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root`
       EFIL=`ls -1 $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$ANADIR/electron."$ARRAY"_ID"$RECID$MCAZ"*.mscw.root`
    fi

##########################################################
# set links for events used in effective area calculation
# (separate training and events used for analysis)
# NOTE: ASSUME THAT THIS IS NOT CHANGED
#       IF DIRECTORY EXISTS, NO NEW ONES ARE CREATED
    ANAEFF="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/${ANADIR}.EFFAREA.MCAZ${MCAZ}"
    # rm -rf $ANAEFF
    if [ ! -e $ANAEFF ]
    then
        mkdir -p $ANAEFF
        # testing signal list
        for arg in $SFIL2
        do
            BW=`basename $arg`
            ln -s $arg $ANAEFF/$BW
        done
        # testing proton list
        for arg in $BFIL2
        do
            BW=`basename $arg`
            ln -s $arg $ANAEFF/$BW
        done
        # depending on CONE parameter: either onSource (if CONE=TRUE) or cone (if CONE=FALSE)
        for arg in $GFIL
        do
            BW=`basename $arg`
            ln -s $arg $ANAEFF/$BW
        done
        # electrons (not used in training)
        for arg in $EFIL
        do
            BW=`basename $arg`
            ln -s $arg $ANAEFF/$BW
        done
    else
        echo "EXISTING LINKED ANALYSIS (MSCW) FILES with testing events"
    fi
    ###############################################################
    # add a 'continue' here if linking file is the main purpose
    # continue
    ###############################################################

###############################################################
# get number of telescopes depending of telescope types

# use first file
   set -- $SFIL1
   # check the file exists - otherwise continue
   if [ -z "$1" ]
   then
       echo "No training file found - continuing"
       echo $1
       continue
   fi
   echo "Teltype cuts: LSTs ($NCUTLST) MSTS ($NCUTMST) SSTs ($NCUTSST) MSCTs ($NCUTMSCT)"
   NTELTYPE=`$EVNDISPSYS/bin/printRunParameter $1 -nteltypes`
   # find correct index for each cut
   for (( N = 0; N < $NTELTYPE; N++ ))
   do 
       TELTYP=`$EVNDISPSYS/bin/printRunParameter $1 -ntype$N`

       NCUT="NCUT${TELTYP}"
       if [ $N -eq 0 ]
       then
           TYPECUT="(NImages_Ttype[${N}]>=${!NCUT}"
       else
           TYPECUT="$TYPECUT\|\|NImages_Ttype[${N}]>=${!NCUT}"
       fi
   done
   if [ ! -z $TYPECUT ]
   then
       TYPECUT="${TYPECUT})"
   fi
   echo "Telescope type cut: $TYPECUT"

###############################################################
# loop over all wobble offset
   for (( W = 0; W < $NOFF; W++ ))
   do
      ODIR=$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/TMVA/$DDIR-${OFFMEA[$W]}
      mkdir -p $ODIR
# copy run parameter file
      cp -f $RPAR.runparameter $ODIR

###############################################################
# loop over all energy bins and submit a job for each bin
      for ((i=0; i < $NENE; i++))
      do

# updating the  run parameter file
	 RFIL=$ODIR/$RXPAR$ARRAY"_$i"
	 echo $RFIL
	 rm -f $RFIL
	 echo "* ENERGYBINS 1 ${EMIN[$i]} ${EMAX[$i]}" > $RFIL.runparameter
         echo "* ZENITHBINS 0 90" >> $RFIL.runparameter
	 echo "* MCXYOFF (MCxoff*MCxoff+MCyoff*MCyoff)>=${OFFMIN[$W]}*${OFFMIN[$W]}&&(MCxoff*MCxoff+MCyoff*MCyoff)<${OFFMAX[$W]}*${OFFMAX[$W]}" >> $RFIL.runparameter
         echo "* MCXYCUTSignalOnly 1" >> $RFIL.runparameter
	 grep "*" $RPAR.runparameter | grep -v ENERGYBINS | grep -v OUTPUTFILE | grep -v SIGNALFILE | grep -v BACKGROUNDFILE | grep -v MCXYOFF >> $RFIL.runparameter
	 echo "* OUTPUTFILE $ODIR $OFIL"_$i" " >> $RFIL.runparameter
         # write signal and background files
         # (note: training is in splitmode=block)
	 for arg in $SFIL1
	 do
	    echo "* SIGNALFILE $arg" >> $RFIL.runparameter
	 done
	 for arg in $SFIL2
	 do
	    echo "* SIGNALFILE $arg" >> $RFIL.runparameter
	 done
	 for arg in $BFIL1
	 do
	    echo "* BACKGROUNDFILE $arg" >> $RFIL.runparameter
	 done
	 for arg in $BFIL2
	 do
	    echo "* BACKGROUNDFILE $arg" >> $RFIL.runparameter
	 done
############################################################
# setting the cuts in the run parameter file

         sed -i "s|MINIMAGES|$NIMAGESMIN|;s|MINIMAGETYPECUT|$TYPECUT|" $RFIL.runparameter
         sed -i 's|ENERGYVARIABLE|ErecS|;s|ENERGYCHI2VARIABLE|EChi2S|g;s|ENERGYDEVARIABLE|dES|g' $RFIL.runparameter
     done

     # loop over all energy bins and prepare run scripts
     for ((i=0; i < $NENE; i++))
     do
         FNAM=$LDIR/$FSCRIPT.$DSET.$ARRAY.$2.${OFFMEA[$W]}.$MCAZ.$i
         FNAM=$ODIR/tmva${i}.qsub

         RRFIL=$ODIR/$RXPAR$ARRAY

         sed -e "s|RUNPARA|$RRFIL|" \
             -e "s|EEEE|$i|" $FSCRIPT.sh  > $FNAM.sh

         chmod u+x $FNAM.sh
         echo "SCRIPT $FNAM.sh"

        #################################
        # submit job to queue
        qsub $QSUBOPT -V -l os="sl6" -l h_cpu=41:59:00 -l h_rss=4000M -l tmpdir_size=1G -o $QLOG -e $QLOG "$FNAM.sh"
    done
  done
done

exit
