#q!/bin/sh
#
# train TMVA BDTs for DISP
#
# Input parameter
#     - output directory for training files
#     - input evndisp file
#     - recid
#     - scaling
#     - array layout (without telescope combinations)
#
# Example for prod3b:
#
# ./CTA.DISPTRAINING.sub_analyse.sh prod3b-paranal20degq05-NN $CTA_USER_DATA/analysis/AnalysisData/prod3b-paranal20degq05-NN/BDT.TMVAO 0 S.3HB9 $CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDTDisp.runparameter 99
#
# Parameter file for training can be found here: $CTA_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDTDisp.runparameter
#
# Removed BDTDispCore (could be simply added)
#

if [ $# -lt 5 ]
then
   echo
   echo "CTA.DISPTRAINING_sub_analyse.sh <data set> <output directory> <recid> <array layout (e.g. S.3HB1)> <TMVA parameters> [scaling] [qsub options (optional)]"
   echo ""
   echo "  <data set>         e.g. cta-ultra3, ISDC3700m, ...  "
   echo "  <output directory> training results will be written to this directory (full path)"
   echo "  <recid>            reconstruction ID according to EVNDISP.reconstruction.parameter"
   echo "  <array layout>     layout name with telescope type ID and scaling (e.g. S.3HB1)"
   echo "  <TMVA parameters>  file name of list of TMVA parameter file"
   echo "  <scaling>          layout scaling (e.g. 5); give 99 to ignore scaling"
   echo
   echo "  (note 1: hardwired telescope types in this script)" 
   echo "  (note 2: disp core training switched off)"
   echo
   exit
fi

#########################################
# input parameters
#########################################
DSET=$1
ODIR=$2
RECID=$3
ARRAY=$4
TMVAP=$5
SCALING=999
if [ -n $6 ]
then
    SCALING=$6
fi
if [ -n $7 ]
then
   QSUBOPT="$7"
fi
QSUBOPT=${QSUBOPT//_X_/ } 
QSUBOPT=${QSUBOPT//_M_/-} 

#########################################
# TMVA options
TMVA=`cat $TMVAP`

#########################################
# checking the path for binary
if [ -z $EVNDISPSYS ]
then
    echo "no EVNDISPSYS environmental variable defined"
    exit
fi

#########################################
# output directory for error/output from batch system
# in case you submit a lot of scripts: QLOG=/dev/null
DATE=`date +"%y%m%d"`
QLOG=$CTA_USER_LOG_DIR/$DATE/DISPTRAINING/
mkdir -p $QLOG

# output directory for shell scripts
SHELLDIR=$CTA_USER_LOG_DIR/$DATE/DISPTRAINING/
mkdir -p $SHELLDIR

# skeleton script
FSCRIPT="CTA.DISPTRAINING.qsub_analyse"

########################################
# list of telescopes
if [[ $DSET == *"LaPalma"* ]]
then
    TELTYPELIST="138704810 10408418 10408618"
else
    TELTYPELIST="138704810 10408418 201309316 909924 10408618 201511619 207308707"
fi

#########################################
# 
#########################################
for BDT in BDTDisp BDTDispEnergy BDTDispError
do
    for MCAZ in 0deg 180deg
    do
      NSTEP=18
      for T in ${TMVA}
      do
        echo $T
        let "NSTEP = $NSTEP + 1"
        OFFDIR=${ODIR}.T${NSTEP}
        ####################
        # output directory
        TDIR="${OFFDIR}/${BDT}/${MCAZ}/"
        mkdir -p $TDIR

        for TELTYPE in $TELTYPELIST
        do
            echo
            echo "STARTING BDT TRAINING FOR AZ DIRECTION $MCAZ AND TELESCOPE TYPE $TELTYPE"
            echo "   training options: ${T}"
            echo "=========================================================================="

            ####################
            # input directory
            if [[ $TELTYPE == "138704810" ]]
            then
                AY=${ARRAY}-NG
            elif [[ $TELTYPE == "10408418" ]]
            then
                AY=${ARRAY}-NG
            elif [[ $TELTYPE == "201309316" ]]
            then
                AY=${ARRAY}-NG
            elif [[ $TELTYPE == "909924" ]]
            then
                AY=${ARRAY}-FD
            elif [[ $TELTYPE == "10408618" ]]
            then
                AY=${ARRAY}-FD
            elif [[ $TELTYPE == "201511619" ]]
            then
                AY=${ARRAY}-FA
            elif [[ $TELTYPE == "207308707" ]]
            then
                AY=${ARRAY}-SD
            else
                echo "Error: unknown telescope type: $TELTYPE"
                exit
            fi
            if [[ ${SCALING} < 99 ]]
            then
                AY=${AY}-${SCALING}
            fi
            if [[ $DSET == *"LaPalma"* ]] || [[ $DSET == *"SCT"* ]]
            then
                AY=${ARRAY}
            fi
            AY=${ARRAY}

            echo $DSET $AY

            ####################
            # input file list
            rm -f $SHELLDIR/tempList.list
            find $CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$AY/gamma_cone/ -name "*[_,.]${MCAZ}*.root" > $SHELLDIR/tempList.list
            NFIL=`wc -l $SHELLDIR/tempList.list | awk '{print $1}'`
            echo "Total number of files available: $NFIL"
            # only use NN% of all evndisp files for training
            # (for LSTs: use more, as there are less telescopes)
            # South:  10%
            # North: 20%
            k=`expr 0.2*$NFIL | bc`
            if [[ $DSET == *"LaPalma"* ]]
            then
                k=$(echo $NFIL | awk '{printf "%d\n",$1*0.20}')
            else
                k=$(echo $NFIL | awk '{printf "%d\n",$1*0.10}')
            fi
            TLIST="$SHELLDIR/EDISP-$DSET-$ARRAY-$SCALING-$MCAZ-$TELTYPE-$BDT-$NSTEP.list"
            rm -f $TLIST
            shuf -n $k $SHELLDIR/tempList.list > $TLIST
            echo "List of $k input files for training: $TLIST"

            ####################
            # prepare run scripts
              FNAM="$SHELLDIR/EDISP-$ARRAY-$SCALING-$MCAZ-$TELTYPE-$BDT-$NSTEP"
              cp $FSCRIPT.sh $FNAM.sh

              sed -i -e "s|OFILE|$TDIR|" \
                     -e "s|TELTYPE|$TELTYPE|" \
                     -e "s|BDTTYPE|$BDT|" \
                     -e "s|RECONSTRUCTIONID|$RECID|" \
                     -e "s|ILIST|$TLIST|" \
                     -e "s|TTT|$T|" \
                     -e "s|AAA|$ARRAY|" \
                     -e "s|DATASET|$DSET|" $FNAM.sh

              chmod u+x $FNAM.sh
              echo "shell script " $FNAM.sh

              # submit the job
              qsub $QSUBOPT -l os=sl6 -l h_cpu=47:45:00 -l h_rss=8000M -V -o $QLOG/ -e $QLOG/ "$FNAM.sh"
           done
         done
     done
done

echo "shell scripts are written to $SHELLDIR"
echo "batch output and error files are written to $QLOG"

exit
