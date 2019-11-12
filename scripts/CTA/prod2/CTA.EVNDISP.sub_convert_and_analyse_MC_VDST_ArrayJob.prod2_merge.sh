#!/bin/bash
#
# script to merge and convert sim_telarray files, then run eventdisplay analysis
#
# PROD2
#
# required HESSIO from at least 20140912 (merged tool)
#
#
#######################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ]
then
   echo
   echo "./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod2_merge.sh <sub array list> <list of sim_telarray files> <particle> <data set> [keep simtel.root files (default off=0)] [log file directory counter] [qsub options]"
   echo
   echo "CTA PROD2 ANALYSIS"
   echo
   echo "  <sub array list>          text file with list of subarray IDs"
   echo
   echo "  <list of sim_telarray files>  list of sim_telarray files (if necessary, expect list of SCSST files with extension SCSST)"
   echo
   echo "  <particle>                gamma_onSource , gamma_diffuse, proton , electron (helium, ...)"
   echo
   echo "  <data set>                e.g. cta-ultra3, ISDC3700m, ..."
   echo ""
   echo "  [keep DST.root files]  keep and copy converted simtel files (DST files) to output directory (default off=0)"
   echo ""
   echo " output will be written to: CTA_USER_DATA_DIR/analysis/Analysis/<subarray>/<particle>/ "
   echo ""
   echo ""
   echo " REMINDER:"
   echo "  - compile CTA converter before starting the scripts (make CTA.convert_hessio_to_VDST)"
   echo "  - make sure that your evndisp installation runs with the same parameters as hessio was compiled with "
   echo "      build_all ultra: -DCTA -DCTA_ULTRA "
   echo "      build_all max:   -DCTA_MAX"
   exit
fi

########################
# prod2 TIMENEXTNEIGHBOUR cleaning
ARRAYCUTS="EVNDISP.prod2.reconstruction.runparameter.NN"
PPOPT="\"-ignoredstgains -NNcleaninginputcard EVNDISP.NNcleaning.dat\""
########################
# prod2 default file
#ARRAYCUTS="EVNDISP.prod2.reconstruction.runparameter"
#PPOPT=""
############################################################################

ARRAY=$1
RUNLIST=$2
PART=$3
DSET=$4
KEEP=0
if [ -n "$5" ]
then
   KEEP=$5
fi
FLL="0"
if [ -n "$6" ]
then
  FLL="$6"
fi
QSUBOPT=""
if [ -n $7 ]
then
   QSUBOPT="$7"
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
FSCRIPT="CTA.EVNDISP.qsub_convert_and_analyse_MC_VDST_ArrayJob.prod2_merge"

# log files
QLOG=$CTA_USER_LOG_DIR/$DATE/EVNDISP-$PART-$DSET/
mkdir -p $QLOG
# QLOG="/dev/null"

# pedestals
# PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod2-Aar/Calibration/Aar.peds.root"
PEDFIL="$CTA_USER_DATA_DIR/analysis/AnalysisData/prod2-Aar/Calibration/Aar.peds.20150623.root"

# get run list and number of runs
if [ ! -e $RUNLIST ]
then
  echo "list of sim_telarray files not found (merge): $RUNLIST"
  exit
fi
RUNLISTN=`basename $RUNLIST`

#########################################################################3
# separate jobs for north and south pointing arrays
# important: direction is indicated by string _0 and _180
for D in 0 180
do

# run lists for north or south
    RUNLISTNdeg=$SHELLDIR/$RUNLISTN.$D
    rm -f $RUNLISTNdeg
    grep "_$D" $RUNLIST > $RUNLISTNdeg
    if [[ -e $RUNLIST.SCSST ]]
    then
        RUNLISTNdegSST=$SHELLDIR/$RUNLISTN.$D.SCSST
        rm -f $RUNLISTNdegSST
        grep "_$D" $RUNLIST.SCSST > $RUNLISTNdegSST
    fi

    NRUN=`wc -l $RUNLISTNdeg | awk '{print $1}'`
    if [[ $NRUN = "0" ]]
    then
       if [[ $D = "0" ]]
       then
          grep north $RUNLIST > $RUNLISTNdeg
          if [[ -e $RUNLIST.SCSST ]]
          then
             grep north $RUNLIST.SCSST < $RUNLISTNdegSST
          fi
       else
          grep south $RUNLIST > $RUNLISTNdeg
          if [[ -e $RUNLIST.SCSST ]]
          then
             grep south $RUNLIST.SCSST < $RUNLISTNdegSST
          fi
       fi
       NRUN=`wc -l $RUNLISTNdeg | awk '{print $1}'`
       NRUNSST=`wc -l $RUNLISTNdegSST | awk '{print $1}'`
    fi
    RUNFROMTO="1-$NRUN"
    NSTEP=1

    echo "submitting $NRUN jobs ( $NRUNSST SC-SST jobs)"

    FNAM="$SHELLDIR/$DSET-$PART-$FLL-$D"

    LIST=`awk '{printf "%s ",$0} END {print ""}' $ARRAY`

    sed -e "s|SIMTELLIST|$RUNLISTNdeg|" \
        -e "s|FULLLIST|$RUNLIST|" \
        -e "s|PAAART|$PART|" \
        -e "s!ARRAY!$LIST!" \
        -e "s|KEEEEEEP|$KEEP|" \
        -e "s|ARC|$ARRAYCUTS|" \
        -e "s|DATASET|$DSET|" \
        -e "s|FLL|$FLL|" \
        -e "s|PPPP|$PEDFIL|" \
        -e "s!UUUU!$PPOPT!" $FSCRIPT.sh > $FNAM.sh

    chmod u+x $FNAM.sh
    echo $FNAM.sh

    if [[ $NRUN -ne 0 ]]
    then
        qsub $QSUBOPT -t $RUNFROMTO:1  -l h_cpu=47:29:00 -l os=sl6 -l tmpdir_size=40G -l h_rss=4G -V -o $QLOG -e $QLOG "$FNAM.sh" 
    fi
done

echo "writing shell script to $FNAM.sh"
echo "writing queue log and error files to $QLOG"

exit
