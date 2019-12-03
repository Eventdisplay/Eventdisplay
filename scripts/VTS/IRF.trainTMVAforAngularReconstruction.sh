#!/bin/bash
# submit TMVA training for angular reconstruction

# qsub parameters
h_cpu=11:29:00; h_vmem=8000M; tmpdir_size=10G

if [[ $# != 5 ]]; then
# begin help message
echo "
TMVA (BDT) training for angular resolution from MC ROOT files for different zenith angle bins
 (simulations that have been processed by evndisp_MC) 

IRF.trainTMVAforAngularReconstruction.sh <epoch> <atmosphere> <zenith> <NSB level> <sim type>

required parameters:

    <epoch>                 array epoch (e.g., V4, V5, V6)
                            V4: array before T1 move (before Fall 2009)
                            V5: array after T1 move (Fall 2009 - Fall 2012)
                            V6: array after camera update (after Fall 2012)
                            
    <atmosphere>            atmosphere model (61 = winter, 62 = summer)

    <zenith>                zenith angle of simulations [deg]

    <NSB level>             NSB level of simulations [MHz]
    
    <sim type>              simulation type (e.g. GRISU, CARE)
    
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
EDVERSION=`"$EVNDISPSYS"/bin/trainTMVAforAngularReconstruction --version | tr -d .`

# Parse command line arguments
EPOCH=$1
ATM=$2
ZA=$3
NOISE=$4
WOBBLE="0.5"
RECID="0"
SIMTYPE=$5
PARTICLE_TYPE="gamma"

# input directory containing evndisp products
if [[ -n "$VERITAS_IRFPRODUCTION_DIR" ]]; then
    INDIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/$SIMTYPE/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_TL5035MA20/ze${ZA}deg_offset${WOBBLE}deg_NSB${NOISE}MHz"
fi
if [[ ! -d $INDIR ]]; then
    echo -e "Error, could not locate input directory. Locations searched:\n $INDIR"
    exit 1
fi
echo "Input file directory: $INDIR"

# Output file directory
if [[ -n "$VERITAS_IRFPRODUCTION_DIR" ]]; then
    ODIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/$SIMTYPE/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}/TMVA_AngularReconstruction/ze${ZA}deg_offset${WOBBLE}deg/"
fi
echo -e "Output files will be written to:\n $ODIR"
mkdir -p "$ODIR"
chmod g+w "$ODIR"

# run scripts and output are written into this directory
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/TMVAAngRes/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}/"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p "$LOGDIR"

# training file name
BDTFILE="mvaAngRes_${ZA}deg_${WOBBLE}wob_NOISE${NOISE}"

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.trainTMVAforAngularReconstruction_sub"

echo "Processing Zenith = $ZA, Noise = $NOISE, Wobble = $WOBBLE"
            
# set parameters in run script
FSCRIPT="$LOGDIR/TA.ID${RECID}.${EPOCH}.$DATE.MC"
sed -e "s|OUTPUTDIR|$ODIR|" \
    -e "s|EVNDISPFILE|$INDIR|" \
    -e "s|BDTFILE|$BDTFILE|" "$SUBSCRIPT.sh" > "$FSCRIPT.sh"

chmod u+x "$FSCRIPT.sh"
echo "$FSCRIPT.sh"

# run locally or on cluster
SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
SUBC=`eval "echo \"$SUBC\""`
if [[ $SUBC == *"ERROR"* ]]; then
    echo $SUBC
    exit
fi
if [[ $SUBC == *qsub* ]]; then
    JOBID=`$SUBC $FSCRIPT.sh`
    echo "JOBID: $JOBID"
elif [[ $SUBC == *parallel* ]]; then
    echo "$FSCRIPT.sh &> $FSCRIPT.log" >> "$LOGDIR/runscripts.dat"
fi

exit
