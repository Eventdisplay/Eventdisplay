#!/bin/bash
# display data files with eventdisplay

## HARD-CODED VALUES
# evndisp reconstruction runparamter
ACUTS="EVNDISP.reconstruction.SW18_noDoublePass.runparameter"
ACUTS="EVNDISP.reconstruction.runparameter"
# run options
OPT="-display=1 -reconstructionparameter $ACUTS -vbfnsamples "
## END OF HARD-CODED VALUES

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP special-purpose analysis: display data file and write results to file

SPANALYSIS.evndisp_display.sh <sourcefile> [telescope combination] [calib] [highres] [run number] [TARGET] [WOBBLENORTH] [WOBBLEEAST] [RAOFFSET]

required parameter:

    <sourcefile>            VERITAS data file (vbf or cvbf file)

optional parameters:
 
    [teltoana]              use 1 for T1 only, 13 for T1 and T3 only, etc.
                            (default telescope combination is 1234)
                            
    [calib]		    0 or nocalib for -nocalibnoproblem, 1 or db for -readCalibDB [default], 2 or raw for -plotraw & -nocalibnoproblem

    [highres]		    0 or lowres for regular window (default), 1 or highres for -highres, 2 or paper for -plotpaper

    [run number]            run number (can differ from actual run number)
    
    [TARGET]                target name (Crab, Mrk421, 1es2344, lsi+61303)
                            (for more do 'evndisp -printtargets')
                            
    [WOBBLENORTH]           wobble offsets north (e.g. 0.5) or south (e.g. -0.5)
                            in units of degrees
    
    [WOBBLEEAST]            wobble offsets east (e.g. 0.5) or west (e.g. -0.5)
                            in units of degrees
    
    [RAOFFSET]              right ascension offset for off run
                            (e.g. 7.5 for off run 30 min later)

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RUNFILE=$1
[[ "$2" ]] && TELTOANA=$2 || TELTOANA="1234"
if [[ $TELTOANA == "-1" ]]; then
    TELTOANA="1234"
fi

if [[ "$3" == "0" ]] || [[ "$3" == "nocalib" ]] ; then
	CALIBOPT=" -nocalibnoproblem "
elif   [[ "$3" == "2" ]] || [[ "$3" == "raw" ]] ; then 
	CALIBOPT=" -plotraw -nocalibnoproblem "
else
	CALIBOPT=" -readCalibDB "
fi

if [[ "$4" == "1" ]] || [[ "$4" == "highres" ]]; then
	PLOTOPT=" -highres "
elif   [[ "$4" == "2" ]] || [[ "$4" == "paper" ]] ; then 
	PLOTOPT=" -highres -plotpaper "
fi

OPT="$OPT $PLOTOPT $CALIBOPT "

[[ "$5" ]] && OPT="$OPT -runnumber=$5"
[[ "$6" ]] && OPT="$OPT -target $6"
[[ "$7" ]] && OPT="$OPT -wobblenorth=$7"
[[ "$8" ]] && OPT="$OPT -wobbleeast=$8"
[[ "$9" ]] && OPT="$OPT -raoffset=$9"

# Check if source file exists
if [[ ! -f $RUNFILE ]]; then
    echo "ERROR: VERITAS source file $RUNFILE not found"
    exit 1
fi

# Set remaining run options
OPT="$OPT -sourcefile $RUNFILE -teltoana=$TELTOANA"

# Run evndisp
echo "$EVNDISPSYS/bin/evndisp $OPT"
$EVNDISPSYS/bin/evndisp $OPT 

exit
