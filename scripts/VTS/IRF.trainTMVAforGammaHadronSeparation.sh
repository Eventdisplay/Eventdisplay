#!/bin/bash
# script to train BDTs with TMVA

# qsub parameters
h_cpu=11:59:59; h_vmem=31599M; tmpdir_size=1G

if (( $# < 5 )); then
# begin help message
echo "
TMVA training of BDT: submit jobs from a TMVA runparameter file

IRF.trainTMVAforGammaHadronSeparation.sh <list of background files> <TMVA runparameter file> <output directory> <output file name> <sim type>
 [epoch] [atmosphere] [Rec ID] [signal directory]

required parameters:

    <list of background files>      list of background training (mscw) files with whole path to each file
    
    <TMVA runparameter file>        TMVA runparameter file with basic options (incl. whole range of 
	                                energy and zenith angle bins) and full path
    
    <output directory>              BDT files are written to this directory
    
    <output file name>              name of output file e.g. BDT  

    <sim type>                      original VBF file simulation type (e.g. GRISU, CARE)

optional parameters:

    [epoch]                         array epoch e.g. V4, V5, V6
                                    default: \"V6\"

    [atmosphere]                    atmosphere model(s) (21 = winter, 22 = summer)
                                    default: \"21\"                   

    [Rec ID]                        reconstruction ID(s) (default: \"0\")
                                    (see EVNDISP.reconstruction.runparameter)	    

    [Analysis method]               analysis method (e.g. NN, TL, MODEL3D)

    [signal directory]              directory containing the MC mscw.root files for the training; 
                                    default: \$VERITAS_IRFPRODUCTION_DIR/v470/\$SIMTYPE/\${EPOCH}_ATM\${ATM}_gamma/MSCW_RECID\${RECID}

additional info:

    energy and zenith angle bins should be indicated in the runparameter file with basic options
--------------------------------------------------------------------------------
"
#end help message
exit
fi
echo " "
# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
$EVNDISPSYS/bin/trainTMVAforGammaHadronSeparation --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`$EVNDISPSYS/bin/trainTMVAforGammaHadronSeparation --version | tr -d .`
else
    EDVERSION="g500"
fi

# Parse command line arguments
BLIST=$1
RUNPAR=$2
ODIR=$3
ONAME=$4
SIMTYPE=$5
echo "Simulation type: $SIMTYPE"
[[ "$6" ]] && EPOCH=$6 || EPOCH="V6"
[[ "$7" ]] && ATM=$7 || ATM="21"
[[ "$8" ]] && RECID=$8 || RECID="0"
PARTICLE_TYPE="gamma"
# analysis method
[[ "$9" ]] && ANAMETHOD=$9 || ANAMETHOD="TL"
[[ "${10}" ]] && SDIR=${10} || SDIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/$SIMTYPE/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${ANAMETHOD}/MSCW_RECID${RECID}"
echo "Signal input directory: $SDIR"
echo "Atmosphere $ATM"
echo "Rec ID $RECID"
echo "Ana method $ANAMETHOD"
# Input directory of MC simulations
if [[ ! -d $SDIR ]]; then
	echo -e "Error, could not locate directory of simulation files (input). Locations searched:\n $SDIR"
	exit 1
fi

# Check that list of background files exists
if [[ ! -f "$BLIST" ]]; then
    echo "Error, list of background files $BLIST not found, exiting..."
    exit 1
fi

# Check that TMVA run parameter file exists
if [[ "$RUNPAR" == `basename $RUNPAR` ]]; then
    RUNPAR="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/$RUNPAR"
fi

if [[ ! -f "$RUNPAR" ]]; then
    echo "Error, TMVA run parameter file $RUNPAR not found, exiting..."
    exit 1
fi

RXPAR=`basename $RUNPAR .runparameter`
echo "Original TMVA run parameter file: $RXPAR.runparameter "

# output directory
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR

count1=1
#####################################
# energy bins
ENBINSBEGIN=$( cat "$RUNPAR" | grep "^* ENERGYBINSBEGIN 1" | sed -e 's/* ENERGYBINSBEGIN 1//' | sed -e 's/ /\n/g')
ENBINSEND=$( cat "$RUNPAR" | grep "^* ENERGYBINSEND 1" | sed -e 's/* ENERGYBINSEND 1//' | sed -e 's/ /\n/g')

declare -a EBINARRAYBEGIN=( $ENBINSBEGIN ) #convert to array
declare -a EBINARRAYEND=( $ENBINSEND ) #convert to array
NENE=$((${#EBINARRAYBEGIN[@]})) #get number of bins

#####################################
# zenith angle bins
ZEBINSBEGIN=$( cat "$RUNPAR" | grep "^* ZENBINSBEGIN" | sed -e 's/* ZENBINSBEGIN//' | sed -e 's/ /\n/g')
ZEBINSEND=$( cat "$RUNPAR" | grep "^* ZENBINSEND" | sed -e 's/* ZENBINSEND//' | sed -e 's/ /\n/g')

declare -a ZEBINARRAYBEGIN=( $ZEBINSBEGIN ) #convert to array
declare -a ZEBINARRAYEND=( $ZEBINSEND ) #convert to array
NZEW=$((${#ZEBINARRAYBEGIN[@]})) #get number of bins

#####################################
# zenith angle bins of MC simulation files
ZENITH_ANGLES=( 20 30 35 40 45 50 55 )
ZENITH_ANGLES=( 20 30 35 )

#####################################
# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/TMVA.ANADATA"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.trainTMVAforGammaHadronSeparation_sub"

###############################################################
# loop over all energy bins and submit a job for each bin
for (( i=0; i < $NENE; i++ ))
do
   echo "==========================================================================="
   echo " "
   echo "EBin: $(($i+$count1)) of $NENE: ${EBINARRAYBEGIN[$i]} to ${EBINARRAYEND[$i]}"
##############################################
# loop over all zenith angle bins
   for (( j=0; j < $NZEW; j++ ))
   do
      echo "---------------------------------------------------------------------------"
      echo "ZeBin: $(($j+$count1)) of $NZEW: ${ZEBINARRAYBEGIN[$j]} to ${ZEBINARRAYEND[$j]}"
      
      # copy run parameter file with basic options to output directory
      cp -f $RUNPAR $ODIR

      # updating the run parameter file for each point in the parameter space
      # (one parameter file per energy and zenith bin)
      RFIL=$ODIR/$RXPAR"_Ebin${i}""_Zebin${j}"
      echo "TMVA Runparameter file: $RFIL.runparameter"
      rm -f $RFIL
      
      echo "* ENERGYBINS 1 ${EBINARRAYBEGIN[$i]} ${EBINARRAYEND[$i]}" > $RFIL.runparameter
      echo "* ZENBINS ${ZEBINARRAYBEGIN[$j]} ${ZEBINARRAYEND[$j]}" >> $RFIL.runparameter
      grep "*" $RUNPAR | grep -v ENERGYBINS | grep -v ZENBINS | grep -v OUTPUTFILE | grep -v SIGNALFILE | grep -v BACKGROUNDFILE | grep -v MCXYOFF >> $RFIL.runparameter
      
     
      # Number of events used. If set too high, no bin will meet the requirements, if set too low, statistics will be poor.
      # set to 0 for using all events. It will use automatically EqualNumEvents to set the same signal and background.
      nTrain="200000"
      # TMPTMP
      nTrain="50000"
 
      echo "* PREPARE_TRAINING_OPTIONS SplitMode=Random:SplitSeed=0:V:VerboseLevel=Verbose:nTrain_Signal=$nTrain:nTrain_Background=$nTrain::nTest_Signal=$nTrain:nTest_Background=$nTrain" >> $RFIL.runparameter

      echo "* OUTPUTFILE $ODIR $ONAME"_Ebin${i}""_Zebin${j}" " >> "$RFIL.runparameter"

      echo "#######################################################################################" >> $RFIL.runparameter
      # signal and background files (NSB levels are simtype dependent)
      for (( l=0; l < ${#ZENITH_ANGLES[@]}; l++ ))
      do
         if (( $(echo "${ZEBINARRAYBEGIN[$j]} <= ${ZENITH_ANGLES[$l]}" | bc ) && $(echo "${ZEBINARRAYEND[$j]} >= ${ZENITH_ANGLES[$l]}" | bc ) ));then
             if (( "${ZENITH_ANGLES[$l]}" != "00" && "${ZENITH_ANGLES[$l]}" != "60" && "${ZENITH_ANGLES[$l]}" != "65" )); then 
                 if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
                     SIGNALLIST=`ls -1 $SDIR/${ZENITH_ANGLES[$l]}deg_0.5wob_NOISE{100,150,200,250,325,425,550}.mscw.root`
                 else
                     #SIGNALLIST=`ls -1 $SDIR/${ZENITH_ANGLES[$l]}deg_0.5wob_NOISE{150,300,450,600,750,900}.mscw.root`
                     #SIGNALLIST=`ls -1 $SDIR/${ZENITH_ANGLES[$l]}deg_0.5wob_NOISE{50,75,100,130,160,200,250,300,350,400,450}.mscw.root`
                     SIGNALLIST=`ls -1 "$SDIR/${ZENITH_ANGLES[$l]}deg_0.5wob_NOISE*.mscw.root"`
                 fi
                 for arg in $SIGNALLIST
                 do
                     echo "* SIGNALFILE $arg" >> $RFIL.runparameter
                 done
             fi
         fi
      done
         
      echo "#######################################################################################" >> $RFIL.runparameter
      for arg in $(cat $BLIST)
      do
          echo "* BACKGROUNDFILE $arg" >> $RFIL.runparameter
      done
         
      FSCRIPT=$LOGDIR/TMVA.$ONAME"_Ebin${i}""_Zebin{$j}_$(date +%s)"
      sed -e "s|RUNPARAM|$RFIL|"  \
          -e "s|NUMTRAIN|$nTrain|"  \
          -e "s|OUTNAME|$ODIR/$ONAME_${i}_${j}|" $SUBSCRIPT.sh > $FSCRIPT.sh

      chmod u+x $FSCRIPT.sh
      echo $FSCRIPT.sh

      # run locally or on cluster
      SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
      SUBC=`eval "echo \"$SUBC\""`
      if [[ $SUBC == *"ERROR"* ]]; then
            echo $SUBC
            exit
      fi
      if [[ $SUBC == *qsub* ]]; then
         JOBID=`$SUBC $FSCRIPT.sh`
         # account for -terse changing the job number format
         if [[ $SUBC != *-terse* ]] ; then
            echo "without -terse!"      # need to match VVVVVVVV  8539483  and 3843483.1-4:2
            JOBID=$( echo "$JOBID" | grep -oP "Your job [0-9.-:]+" | awk '{ print $3 }' )
         fi
    
         echo "JOBID:  $JOBID"
    
      elif [[ $SUBC == *parallel* ]]; then
         echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
         cat $LOGDIR/runscripts.dat | $SUBC
      elif [[ "$SUBC" == *simple* ]] ; then
         "$FSCRIPT.sh" | tee "$FSCRIPT.log"
      fi
   done
done

exit
