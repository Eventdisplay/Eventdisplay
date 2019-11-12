#!/bin/bash
# submit effective area analysis
# (output need to be combined afterwards)

# qsub parameters
h_cpu=11:29:00; h_vmem=8000M; tmpdir_size=10G

if [[ $# -lt 8 ]]; then
# begin help message
echo "
IRF generation: create partial effective area files from MC ROOT files
 (simulations that have been processed by both evndisp_MC and mscw_energy_MC)

IRF.generate_effective_area_parts.sh <cuts file> <epoch> <atmosphere> <zenith> <offset angle> <NSB level> <Rec ID> <sim type> [Analysis Method]

required parameters:

    <cuts file>             gamma/hadron cuts file (located in 
                             \$VERITAS_EVNDISP_AUX_DIR/GammaHadronCutFiles)
                            (might be a list of cut files)
        
    <epoch>                 array epoch (e.g., V4, V5, V6)
                            V4: array before T1 move (before Fall 2009)
                            V5: array after T1 move (Fall 2009 - Fall 2012)
                            V6: array after camera update (after Fall 2012)
                            
    <atmosphere>            atmosphere model (21 = winter, 22 = summer)

    <zenith>                zenith angle of simulations [deg]

    <offset angle>          offset angle of simulations [deg]

    <NSB level>             NSB level of simulations [MHz]
    
    <Rec ID>                reconstruction ID
                            (see EVNDISP.reconstruction.runparameter)
                            Set to 0 for all telescopes, 1 to cut T1, etc.
        
    <sim type>              simulation type (e.g. GRISU-SW6, CARE_June1425)

    optional parameters:

    [Analysis-Method]       select analysis method
                            (TL, NN, FROGS, MODEL3D, see IRF.evndisp_MC.sh)

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
$EVNDISPSYS/bin/makeEffectiveArea --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`$EVNDISPSYS/bin/makeEffectiveArea --version | tr -d .`
else
    EDVERSION="v502"
fi

# Parse command line arguments
CUTSFILE="$1"
EPOCH=$2
ATM=$3
ZA=$4
WOBBLE=$5
NOISE=$6
RECID=$7
SIMTYPE=$8
PARTICLE_TYPE="gamma"
[[ "${9}" ]] && ANAMETHOD=${9} || ANAMETHOD="TL"

# input directory containing mscw_energy_MC products
if [[ -n $VERITAS_IRFPRODUCTION_DIR ]]; then
    INDIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/$SIMTYPE/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${ANAMETHOD}/MSCW_RECID${RECID}"
fi
if [[ ! -d $INDIR ]]; then
    echo -e "Error, could not locate input directory. Locations searched:\n $INDIR"
    exit 1
fi
echo "Input file directory: $INDIR"

# Output file directory
if [[ -n "$VERITAS_IRFPRODUCTION_DIR" ]]; then
    ODIR="$VERITAS_IRFPRODUCTION_DIR/$EDVERSION/$SIMTYPE/${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${ANAMETHOD}/"
fi
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR
chmod g+w $ODIR

# run scripts and output are written into this directory
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EFFAREA/${ZA}deg_${WOBBLE}wob_NOISE${NOISE}_${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${ANAMETHOD}_${RECID}/"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

#################################
# template string containing the name of processed simulation root file
MCFILE="${INDIR}/${ZA}deg_${WOBBLE}wob_NOISE${NOISE}.mscw.root"
if [[ ! -f ${MCFILE} ]]; then
    echo "Input mscw file not found: ${MCFILE}"
    exit 1
fi

# effective area output file
EFFAREAFILE="EffArea-${SIMTYPE}-${EPOCH}-ID${RECID}-Ze${ZA}deg-${WOBBLE}wob-${NOISE}"

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.effective_area_parallel_sub"

echo "Processing Zenith = $ZA, Noise = $NOISE, Wobble = $WOBBLE"
            
echo $CUTSFILE
# set parameters in run script
FSCRIPT="$LOGDIR/EA.ID${RECID}.$DATE.MC-$(date +%s)"
sed -e "s|OUTPUTDIR|$ODIR|" \
    -e "s|EFFFILE|$EFFAREAFILE|" \
    -e "s|DATAFILE|$MCFILE|" \
    -e "s|ANALYSISMETHOD|$ANAMETHOD|" \
    -e "s|GAMMACUTS|${CUTSFILE}|" $SUBSCRIPT.sh > $FSCRIPT.sh

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
    echo "JOBID: $JOBID"
elif [[ $SUBC == *parallel* ]]; then
    echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
elif [[ "$SUBC" == *simple* ]]; then
    "$FSCRIPT.sh" | tee "$FSCRIPT.log"
fi

exit
