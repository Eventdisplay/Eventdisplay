#!/bin/bash
#
# script to convert sim_tel output files and then run eventdisplay analysis
#
# PROD4
#
# Author: Gernot Maier
#
#######################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ]
then
   echo
   echo "./CTA.EVNDISP.sub_convert_and_analyse_MC_VDST_ArrayJob.prod4.sh <sub array list> <list of sim_telarray files> <particle> <data set> [keep simtel.root files (default off=0)] [log file directory counter] [qsub options]"
   echo
   echo "CTA PROD4 ANALYSIS"
   echo
   echo "  <sub array list>          text file with list of subarray IDs"
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
# prod4 default (TIMENEXTNEIGHBOUR for some telescope types) cleaning
ARRAYCUTS="EVNDISP.prod4.reconstruction.runparameter.NN.noLL"
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
FSCRIPT="CTA.EVNDISP.qsub_convert_and_analyse_MC_VDST_ArrayJob.prod4"

# log files
QLOG=$CTA_USER_LOG_DIR/$DATE/EVNDISP-$PART-$DSET/
mkdir -p $QLOG
# QLOG="/dev/null"

############################################
# calibration file with IPR graphs
# default (prod4)
PEDFIL="$CTA_EVNDISP_PATH/Calibration/prod4-SST-IPR.root"
PEDFIL="/lustre/fs21/group/cta/users/maierg/GRID_tarballs/20181028_prod4/Calibration/pedestals-mst-flashcam.root"
echo "PEDFIL: $PEDFIL"

########################################################
# get run list and number of runs
if [ ! -e $RUNLIST ]
then
  echo "list of sim_telarray files not found: $RUNLIST"
  exit
fi
RUNLISTN=`basename $RUNLIST`

#########################################################################3
# separate job for north and south
# for D in 0 180
for D in 0
do

# run lists for north or south
    RUNLISTNdeg=$SHELLDIR/$RUNLISTN.${D}.${DSET}
    rm -f $RUNLISTNdeg
    # grep "_$D" $RUNLIST > $RUNLISTNdeg
    # TMPTMP
    cat $RUNLIST > $RUNLISTNdeg

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

    FNAM="$SHELLDIR/$DSET-$PART-$FLL-$D"

    LIST=`awk '{printf "%s ",$0} END {print ""}' $ARRAY`

    sed -e "s|SIMTELLIST|$RUNLISTNdeg|" \
        -e "s|PAAART|$PART|" \
        -e "s!ARRAY!$LIST!" \
        -e "s|KEEEEEEP|$KEEP|" \
        -e "s|ARC|$ARRAYCUTS|" \
        -e "s|DATASET|$DSET|" \
        -e "s|FLL|$FLL|" \
        -e "s|PPPP|$PEDFIL|" \
        -e "s|STST|$STEPSIZE|" $FSCRIPT.sh > $FNAM.sh

    chmod u+x $FNAM.sh
    echo $FNAM.sh

    if [[ $NRUN -ne 0 ]]
    then
        if [[ $DSET == "SCT" ]]
        then
            qsub $QSUBOPT -t $RUNFROMTO:1  -l h_cpu=11:29:00 -l tmpdir_size=40G -l h_rss=8G -V -o $QLOG -e $QLOG "$FNAM.sh" 
        else
            qsub $QSUBOPT -t $RUNFROMTO:1  -l h_cpu=11:29:00 -l tmpdir_size=40G -l h_rss=4G -V -o $QLOG -e $QLOG "$FNAM.sh" 
        fi
    echo "submit"
    fi
done

echo "writing shell script to $FNAM.sh"
echo "writing queue log and error files to $QLOG"

exit
