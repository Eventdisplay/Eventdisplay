#!/bin/bash
# submit evndisp energy3d
#

# qsub parameters
h_cpu=1:00:00; h_vmem=5000M; tmpdir_size=5G

if [ $# -lt 6 ]; then
# begin help message
echo "
--------------------------------------------------------------------------------
create_energy3d_referencetables creates the reference tables that energy3d uses
Needs to be done only once per V4/V5/V6 21/22 CARE/GRISU, the respective simulations needes to be analysed with IRF.evndisp_MC with model3d first

IRF.create_energy3d_referencetables.sh <Simtype> <Telsetup> <ATM> <Noice> <Outdir>

required parameters:

	<Simtype>	CARE, GRISU or CAREGRISU, all noices and all elevation required

	<Telsetup>	4, 5, 6, or 456

	<ATM>		21, 22 or 2122 (currecntly care only has v6)

	<Indir>		General dir where files can be found 

	<Outdir>	Your wanted output directory


Example: ./IRF.create_energy3d_biascorrectionfile.sh CAREGRISU 456 2122 listdir/noiselist.txt $VERITAS_USER_DATA_DIR/analysis/energy3d/simanalysis/ $VERITAS_USER_DATA_DIR/DATASets3DFull/Merged/
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`

# Parse command line arguments
SIMTYPE=$1
TELSETUP=$2
SUMWINT=$3
INDIR=$5
OUTDIR=$6

# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANAMCVBF"
mkdir -p $LOGDIR

# output directory for evndisp products (will be manipulated more later in the script)
if [[ ! -z "$OUTDIR" ]]; then
    ODIR="$OUTDIR"
fi

# output dir
OPDIR=$OUTDIR
mkdir -p $OPDIR
chmod -R g+w $OPDIR
echo -e "energy3d_biascorrectionfile output files will be written to:\n $OPDIR"

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.create_energy3d_biascorrectionfile_sub"

# make run script
FSCRIPT="$LOGDIR/EVN.create_energy3d_biascorrectionfile"
    
sed -e "s|SIMULATION|$SIMTYPE|" \
    -e "s|TELESTEUP|$TELSETUP|"\
    -e "s|SUMERWINTER|$SUMWINT|" \
    -e "s|INPUTDIR|$INDIR|" \
    -e "s|OUTPUTDIR|$OUTDIR|" $SUBSCRIPT.sh > $FSCRIPT.sh

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
      echo "RUN $AFILE: JOBID $JOBID"
elif [[ $SUBC == *parallel* ]]; then
     echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
fi
                               
exit
