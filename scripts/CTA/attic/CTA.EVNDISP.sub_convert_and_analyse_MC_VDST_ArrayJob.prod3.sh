#!/bin/bash
#
# script to convert sim_tel output files and then run eventdisplay analysis
#
# PROD3
#
#
#######################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ]
then
   echo
   echo "./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod3.sh <sub array list> <array scaling> <list of sim_telarray files> <particle> <data set> [keep simtel.root files (default off=0)] [log file directory counter] [qsub options]"
   echo
   echo "CTA PROD3 ANALYSIS"
   echo
   echo "  <sub array list>          text file with list of subarray IDs"
   echo
   echo "  <array scaling>           array layout scaling factor( 1,2,...5)"
   echo
   echo "  <particle>                gamma_onSource , gamma_diffuse, proton , electron (helium, ...)"
   echo
   echo "  <data set>                e.g. paranal, ..."
   echo ""
   echo "NOTE: HARDWIRED FILE NAMES IN QSUB SCRIPTS !!"
   echo ""
   echo "  [keep DST.root files]  keep and copy converted simtel files (DST files) to output directory (default off=0)"
   echo ""
   echo " output will be written to: CTA_USER_DATA_DIR/analysis/<subarray>/<particle>/ "
   echo ""
   echo ""
   echo " REMINDER:"
   echo "  - compile CTA converter before starting the scripts (make CTA.convert_hessio_to_VDST)"
   exit
fi

########################
# prod3 default (TIMENEXTNEIGHBOUR for some telescope types) cleaning
# FlashCam analysis with digital filter
ARRAYCUTS="EVNDISP.prod3.reconstruction.runparameter.NN"
PPOPT="\" -NNcleaninginputcard EVNDISP.NNcleaning.dat\""
############################################################################
# FlashCam analysis without digital filter
# ARRAYCUTS="EVNDISP.prod3.reconstruction.runparameter.NN.noDF"

ARRAY=$1
ASCALE=$2
RUNLIST=$3
PART=$4
DSET=$5
KEEP=0
if [ -n "$6" ]
then
   KEEP=$6
fi
FLL="0"
if [ -n "$7" ]
then
  FLL="$7"
fi
QSUBOPT=""
if [ -n $8 ]
then
   QSUBOPT="$8"
fi
QSUBOPT=${QSUBOPT//_X_/ } 
QSUBOPT=${QSUBOPT//_M_/-} 

# checking the path for binary
if [ -z $EVNDISPSYS ]
then
    echo "no EVNDISPSYS env variable defined"
    exit
fi

#########################################
# output directory for error/output from batch system
# in case you submit a lot of scripts: QLOG=/dev/null
DATE=`date +"%y%m%d"`

# output directory for shell scripts and run lists
SHELLDIR=$CTA_USER_DATA_DIR"/queueShellDir/"
mkdir -p $SHELLDIR

# skeleton script
FSCRIPT="CTA.EVNDISP.qsub_convert_and_analyse_MC_VDST_ArrayJob.prod3"

# log files
QLOG=$CTA_USER_LOG_DIR/$DATE/EVNDISP-$PART-$DSET/
mkdir -p $QLOG
# QLOG="/dev/null"

# pedestals
# PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.peds.20150820.dst.root"
# pedfile for FlashCam digital filter analysis (South)
# (used for most of the production)
PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.peds.20151004.dst.root"
# New 20160531 (updated GCT and 1mDC trace integration)
PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3-Paranal20deg-20151004.20160531/20150807b.gamma_20deg_180deg_run18___cta-prod3-demo_desert-2150m-Paranal.simtel.ped.root"
# New 20160603 (20160531 plus 1mDC digital filter)
PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3-Paranal20deg-20151004.20160603/20150807b.gamma_20deg_180deg_run18___cta-prod3-demo_desert-2150m-Paranal.simtel.ped.root"
if [[ $DSET == *"40deg"* ]]
then
   PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3-Paranal40deg-20151106.20160531/20151106.gamma_20deg_180deg_run1___cta-prod3-demo_desert-2150m-Paranal.ped.root"
fi
# pedfile for FlashCam digital filter analysis (North)
if [[ $DSET == *"LaPalma"* ]]
then
    if [[ $RUNLIST == *"FlashCam"* ]]
    then
        PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/RateFiles/prod3RateFile-LaPalmaFlashCam-20160211.root"
    else
        PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/RateFiles/prod3RateFile-LaPalmaNectarCam-20160211.root"
    fi
fi
# pedfile for no digital filter
# PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod3-Calibration/prod3.peds.20151004.noDF.dst.root"

# get run list and number of runs
if [ ! -e $RUNLIST ]
then
  echo "list of sim_telarray files not found: $RUNLIST"
  exit
fi
RUNLISTN=`basename $RUNLIST`

#########################################################################3
# separate job for north and south
for D in 0 180
do

# run lists for north or south
    RUNLISTNdeg=$SHELLDIR/$RUNLISTN.$D
    rm -f $RUNLISTNdeg
    grep "_$D" $RUNLIST > $RUNLISTNdeg

    NRUN=`wc -l $RUNLISTNdeg | awk '{print $1}'`
    if [[ $NRUN = "0" ]]
    then
       if [[ $D = "0" ]]
       then
          grep north $RUNLIST > $RUNLISTNdeg
       else
          grep south $RUNLIST > $RUNLISTNdeg
       fi
       NRUN=`wc -l $RUNLISTNdeg | awk '{print $1}'`
    fi
    RUNFROMTO="1-$NRUN"
    NSTEP=1

# combine files to bunches of $STEPSIZE
    STEPSIZE=10
# smaller step size for on source gammas
    NARRAY=`wc -l $ARRAY`
    echo $NARRAY
    if [[ $PART == "gamma_onSource" ]]
    then
       STEPSIZE=5
    fi
    STEPSIZE=1
    let "NSTEP = $NRUN / $STEPSIZE"
    let "NTES  = $NSTEP * $STEPSIZE"
    if [[ $NTES -ne $NRUN ]]
    then
        let "NSTEP = $NSTEP + 1"
    fi
    RUNFROMTO="1-$NSTEP"

    echo "submitting $NRUN jobs ($NSTEP steps of size $STEPSIZE, $RUNFROMTO)"

    FNAM="$SHELLDIR/$DSET-$PART-$FLL-$D-$ASCALE"

    LIST=`awk '{printf "%s ",$0} END {print ""}' $ARRAY`

    sed -e "s|SIMTELLIST|$RUNLISTNdeg|" \
        -e "s|PAAART|$PART|" \
        -e "s!ARRAY!$LIST!" \
        -e "s|KEEEEEEP|$KEEP|" \
        -e "s|ARC|$ARRAYCUTS|" \
        -e "s|DATASET|$DSET|" \
        -e "s|FLL|$FLL|" \
        -e "s|PPPP|$PEDFIL|" \
        -e "s!UUUU!$PPOPT!" \
        -e "s|STST|$STEPSIZE|" \
        -e "s|ASCCALEM|$ASCALE|" $FSCRIPT.sh > $FNAM.sh

    chmod u+x $FNAM.sh
    echo $FNAM.sh

    if [[ $NRUN -ne 0 ]]
    then
#        qsub $QSUBOPT -t $RUNFROMTO:1  -l h_cpu=11:29:00 -l os=sl6 -l tmpdir_size=40G -l h_rss=4G -V -o $QLOG -e $QLOG "$FNAM.sh" 
        # HD produced prod3 files need more memory:
        qsub $QSUBOPT -t $RUNFROMTO:1  -l h_cpu=47:29:00 -l os=sl6 -l tmpdir_size=40G -l h_rss=8G -V -o $QLOG -e $QLOG "$FNAM.sh" 
    fi
done

echo "writing shell script to $FNAM.sh"
echo "writing queue log and error files to $QLOG"

exit
