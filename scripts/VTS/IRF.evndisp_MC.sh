#!/bin/bash
# submit evndisp for grisu/care simulations
#

# qsub parameters
h_cpu=47:59:00; h_vmem=8000M; tmpdir_size=250G

if [ $# -lt 7 ]; then
# begin help message
echo "
IRF generation: analyze simulation VBF files using evndisp 

IRF.evndisp_MC.sh <sim directory> <epoch> <atmosphere> <zenith> <offset angle> <NSB level> <sim type> [runparameter file] [external NSB] [particle] [Analysis Method] [events] [Rec ID]

required parameters:

    <sim directory>         directory containing simulation VBF files

    <epoch>                 array epoch (e.g., V4, V5, V6)
                            V4: array before T1 move (before Fall 2009)
                            V5: array after T1 move (Fall 2009 - Fall 2012)
                            V6: array after camera update (after Fall 2012)
    
    <atmosphere>            atmosphere model (21 = winter, 22 = summer)

    <zenith>                zenith angle of simulations [deg]

    <offset angle>          offset angle of simulations [deg]

    <NSB level>             NSB level of simulations [MHz]
    
    <sim type>              file simulation type (e.g. GRISU-SW6, CARE_June1425)

optional parameters:

    <runparameter file>     file with integration window size and reconstruction cuts/methods, 
                            expected in $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/

                            Default: EVNDISP.reconstruction.runparameter

    [external NSB]          NSB traces are read from external files (don't worry about this flag when analyzing grisu MC)
    
    [particle]              type of particle used in simulation:
                            gamma = 1, electron = 2, proton = 14, helium = 402
                            (default = 1  -->  gamma)

    [Analysis-Method]       select analysis method
                            TL      (standard image analysis with two-level cleaning; default)
                            NN      (standard image analysis using optimized NN cleaning)
    
    [events]                number of events to be analysed
                            (default: -1)

Note: zenith angles, wobble offsets, and noise values are hard-coded into script

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
"$EVNDISPSYS"/bin/evndisp --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`"$EVNDISPSYS"/bin/evndisp --version | tr -d .`
else
    EDVERSION="g500"
fi

# Parse command line arguments
SIMDIR=$1
EPOCH=$2
ATM=$3
ZA=$4
WOBBLE=$5
NOISE=$6
SIMTYPE=$7
[[ "$8" ]] && ACUTS=$8 || ACUTS=EVNDISP.reconstruction.runparameter
[[ "$9" ]] && EXTNOISE=$9 || EXTNOISE=0
[[ "${10}" ]] && PARTICLE=${10} || PARTICLE=1
[[ "${11}" ]] && ANAMETHOD=${11} || ANAMETHOD="TL"
[[ "${12}" ]] && NEVENTS=${12}  || NEVENTS=-1
[[ "${13}" ]] && PART=${13}  || PART=0

# check and print analysis method
if [[ $ANAMETHOD != "TL"* ]] && [[ $ANAMETHOD != "NN"* ]]
then
    echo "unknown analysis method $ANAMETHOD"
    exit
fi
echo "Analysis method $ANAMETHOD"

# Particle names
PARTICLE_NAMES=( [1]=gamma [2]=electron [14]=proton [402]=alpha )
PARTICLE_TYPE=${PARTICLE_NAMES[$PARTICLE]}

# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANAMCVBF"
mkdir -p "$LOGDIR"

# output directory for evndisp products (will be manipulated more later in the script)
if [[ ! -z "$VERITAS_IRFPRODUCTION_DIR" ]]; then
    ODIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/${SIMTYPE}/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${ANAMETHOD}"
fi
# output dir
OPDIR=${ODIR}"/ze"$ZA"deg_offset"$WOBBLE"deg_NSB"$NOISE"MHz/part"$PART
mkdir -p "$OPDIR"
chmod -R g+w "$OPDIR"
echo -e "Output files will be written to:\n $OPDIR"

echo "Using runparameter file $ACUTS"

# Create a unique set of run numbers
if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
    [[ $EPOCH == "V4" ]] && RUNNUM="946500"
    [[ $EPOCH == "V5" ]] && RUNNUM="956500"
    [[ $EPOCH == "V6" ]] && RUNNUM="966500"
elif [ ${SIMTYPE:0:4} = "CARE" ]; then
    [[ $EPOCH == "V4" ]] && RUNNUM="941200"
    [[ $EPOCH == "V5" ]] && RUNNUM="951200"
    [[ $EPOCH == "V6" ]] && RUNNUM="961200_$PART"
fi

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.evndisp_MC_sub"

INT_WOBBLE=`echo "$WOBBLE*100" | bc | awk -F '.' '{print $1}'`
if [[ ${#INT_WOBBLE} -lt 2 ]]; then
   INT_WOBBLE="000"
elif [[ ${#INT_WOBBLE} -lt 3 ]]; then
   INT_WOBBLE="0$INT_WOBBLE"
fi

# make run script
FSCRIPT="$LOGDIR/evn-$EPOCH-$SIMTYPE-$ZA-$WOBBLE-$NOISE-ATM$ATM-$(date +%s)"
sed -e "s|DATADIR|$SIMDIR|" \
    -e "s|RUNNUMBER|$RUNNUM|" \
    -e "s|ZENITHANGLE|$ZA|" \
    -e "s|ATMOSPHERE|$ATM|" \
    -e "s|OUTPUTDIR|$OPDIR|" \
    -e "s|DECIMALWOBBLE|$WOBBLE|" \
    -e "s|INTEGERWOBBLE|$INT_WOBBLE|" \
    -e "s|NOISELEVEL|$NOISE|" \
    -e "s|ARRAYEPOCH|$EPOCH|" \
    -e "s|ANALYSISMETHOD|$ANAMETHOD|" \
    -e "s|RECONSTRUCTIONRUNPARAMETERFILE|$ACUTS|" \
    -e "s|SIMULATIONTYPE|$SIMTYPE|" \
    -e "s|PAART|$PART|" \
    -e "s|NNNOISEFILES|$EXTNOISE|" \
    -e "s|REDHVFLAG|$RHVFLAG|" \
    -e "s|PARTICLETYPE|$PARTICLE|" "$SUBSCRIPT.sh" > "$FSCRIPT.sh"

chmod u+x "$FSCRIPT.sh"
echo "$FSCRIPT.sh"

# run locally or on cluster
SUBC=`"$EVNDISPSYS"/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
SUBC=`eval "echo \"$SUBC\""`
if [[ $SUBC == *"ERROR"* ]]; then
    echo $SUBC
    exit
fi
if [[ $SUBC == *qsub* ]]; then
    if [[ $NEVENTS -gt 0 ]]; then
        JOBID=`$SUBC -t 1-10 $FSCRIPT.sh`
    elif [[ $NEVENTS -lt 0 ]]; then
        JOBID=`$SUBC $FSCRIPT.sh`
    fi      
    echo "RUN $RUNNUM: JOBID $JOBID"
elif [[ $SUBC == *parallel* ]]; then
    echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
fi
                
exit
