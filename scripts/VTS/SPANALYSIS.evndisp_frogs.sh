#!/bin/bash
# analyze data files with FROGS

## HARDCODED VALUES
# run options (will be overwritten in DST mode)
# OPT="-nevents=2000"
## END OF HARDCODED VALUES

if [[ $# < 2 ]]; then
# begin help message
echo "
EVNDISP special-purpose analysis: analyse a data file with FROGS and add FROGS
par tree to mscw file

SPANALYSIS.evndisp_frogs.sh <run number> <mscw directory> [teltoana]
 [readcalibdb] [DST]

required parameters:

    <run number>            perform FROGS analysis for this run

    <mscw directory>        directory which contains mscw_energy files

optional parameters:

    [teltoana]              use 1 for T1 only, 13 for T1 and T3 only, etc.
                            (default telescope combination is 1234)
  
    [readcalibdb]           set to 0 to switch off (default is on)

    [DST]                   output format is DST (default OFF=0)

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RUNNUM=$1
MSCWDIR=$2
[[ "$3" ]] && TELTOANA=$3 || TELTOANA="1234"
if [[ $TELTOANA == "-1" ]]; then
    TELTOANA="1234"
fi
[[ "$4" ]] && CALDB=$4   || CALDB="1"
[[ "$5" ]] && DSTMODE=$5 || DSTMODE="0"

# Check if source vbf file exists
SF=`find -L $VERITAS_DATA_DIR/data -name "$RUNNUM.cvbf"`
if [ ${#SF} = 0 ]; then
    echo "ERROR: VERITAS source file $RUNNUM.cvbf not found in $VERITAS_DATA_DIR/data/"
    exit 1
fi

# Check if mscw_energy file exists
MSCWFILE="$MSCWDIR/$RUNNUM.mscw.root"
if [[ ! -f "$MSCWFILE" ]]; then
    echo "ERROR: MSCW file $MSCWFILE does not exist"
    exit 1
fi

# Set run options
OPT="-frogs $MSCWFILE frogid 0 -runnumber=$RUNNUM -teltoana=$TELTOANA"
if [[ $CALDB == "1" ]]; then
    OPT="$OPT -readCalibDB"
fi

# DST production
if [[ DSTMODE != "0" ]]; then
    OPT="$OPT -runmode=4 -imagethresh=0 -borderthresh=0 -dstfile $VERITAS_DATA_DIR/eventdisplay_output/$RUNNUM.dst.root"
fi

# Run evndisp
echo "$EVNDISPSYS/bin/evndisp $OPT"
$EVNDISPSYS/bin/evndisp $OPT

exit
