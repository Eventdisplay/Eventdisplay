#!/bin/bash
# script to run eventdisplay analysis for VTS data

# qsub parameters
# h_cpu=41:59:00; h_vmem=4000M; tmpdir_size=25G
h_cpu=41:59:00; h_vmem=4000M; tmpdir_size=25G

# EventDisplay version
$EVNDISPSYS/bin/evndisp --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`$EVNDISPSYS/bin/evndisp --version | tr -d .`
else
    EDVERSION="g502"
fi

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP data analysis: submit jobs from a simple run list

ANALYSIS.evndisp.sh <runlist> [output directory] [runparameter file] [calibration] [Analysis Method] [teltoana] [calibration file name]

required parameters:

    <runlist>               simple run list with one run number per line
    
optional parameters:
    
    [output directory]     directory where output ROOT files will be stored.
			   Default: $VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/
	 
    [runparameter file]    file with integration window size and reconstruction cuts/methods,
                           expected in $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/

			   Default: EVNDISP.reconstruction.runparameter (long sumwindow -> for use with CARE IRFs; DISP disabled )

			   other options:

			   EVNDISP.reconstruction.runparameter.DISP		 (long sumwindow -> for use with CARE IRFs;
										                    DISP enabled, use RecID 1 in later stages to access it)
										 
			   EVNDISP.reconstruction.runparameter.SumWindow6-noDISP (short sumwindow -> for use with grisu IRFs; DISP disabled)
			   EVNDISP.reconstruction.runparameter.SumWindow6-DISP	 (short sumwindow -> for use with grisu IRFs; DISP enabled [RecID 1])

    [calibration]	   
          0		   neither tzero nor pedestal calculation is performed, must have the calibration results
			   already in $VERITAS_EVENTDISPLAY_AUX_DIR/Calibration/Tel_?
          1                pedestal & average tzero calculation (default)
          2                pedestal calculation only
          3                average tzero calculation only
          4                pedestal & average tzero calculation are performed;
                           laser run number is taken from calibration file,
                           gains taken from $VERITAS_EVENTDISPLAY_AUX_DIR/Calibration/Tel_?/<laserrun>.gain.root 

    [Analysis Method]      TL      (standard image analysis with two-level cleaning; default)
                           NN      (standard image analysis using optimized NN cleaning)
                           MODEL3D (model3D analysis)

    [teltoana]             restrict telescope combination to be analyzed:
                           e.g.: teltoana=123 (for tel. 1,2,3), 234, ...
                           Default is to use the telescope combination from the DB. Telescopes that were not in the array
			   or have been cut by DQM are not analysed.

    [calibration file name] only used with calibration=4 option
			   to specify a which runs should be used for pedestal/tzero/gain calibration.
                           Default is calibrationlist.dat
                           file is expected in $VERITAS_EVNDISP_AUX_DIR/Calibration

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# create extra stdout for duplication of command output
# look for ">&5" below
exec 5>&1

# Parse command line arguments
RLIST=$1
[[ "$2" ]] && ODIR=$2 || ODIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
mkdir -p $ODIR
[[ "$3" ]] && ACUTS=$3 || ACUTS=EVNDISP.reconstruction.runparameter
[[ "$4" ]] && CALIB=$4 || CALIB=1
[[ "$5" ]] && ANAMETHOD=${5} || ANAMETHOD="TL"
[[ "$6" ]] && TELTOANA=$6 || TELTOANA=1234
[[ "$7" ]] && CALIBFILE=$7 || CALIBFILE=calibrationlist.dat

VPM=1

echo "Using runparameter file $ACUTS"


# Read runlist
if [ ! -f "$RLIST" ] ; then
    echo "Error, runlist $RLIST not found, exiting..."
    exit 1
fi
FILES=`cat $RLIST`

# Output directory for error/output
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EVNDISP.ANADATA"
mkdir -p $LOGDIR

# Job submission script
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/ANALYSIS.evndisp_sub"

TIMETAG=`date +"%s"`

NRUNS=`cat $RLIST | wc -l ` 
echo "total number of runs to analyze: $NRUNS"
echo

#########################################
# loop over all files in files loop
for AFILE in $FILES
do
    echo "Now starting run $AFILE"
    FSCRIPT="$LOGDIR/EVN.data-$AFILE-$(date +%s)"

    sed -e "s|RUNFILE|$AFILE|"              \
        -e "s|CALIBRATIONOPTION|$CALIB|"    \
        -e "s|OUTPUTDIRECTORY|$ODIR|"       \
        -e "s|USEVPMPOINTING|$VPM|" \
        -e "s|RECONSTRUCTIONRUNPARAMETERFILE|$ACUTS|" \
        -e "s|TELTOANACOMB|$TELTOANA|"                   \
        -e "s|USECALIBLIST|$CALIBFILE|"                  \
        -e "s|ANALYSISMETHOD|$ANAMETHOD|" $SUBSCRIPT.sh > $FSCRIPT.sh

    chmod u+x $FSCRIPT.sh
    echo $FSCRIPT.sh
	# output selected input during submission:

	echo "Using runparameter file $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/$ACUTS"

	if [[ $VPM == "1" ]]; then
	echo "VPM is switched on (default)"
	else
	echo "VPM bool is set to $VPM (switched off)"
	fi  

	if [[ $TELTOANA == "1234" ]]; then
	echo "Telescope combination saved in the DB is analyzed (default)"
	else
	echo "Analyzed telescopes: $TELTOANA"
	fi 
	if [[ $CALIB == "4" ]]; then
	echo "read calibration from calibration file $CALIBFILE"
	else
            echo "read calibration from VOffline DB (default)"
	fi 

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
        
        echo "RUN $AFILE JOBID $JOBID"
        echo "RUN $AFILE SCRIPT $FSCRIPT.sh"
        if [[ $SUBC != */dev/null* ]] ; then
            echo "RUN $AFILE OLOG $FSCRIPT.sh.o$JOBID"
            echo "RUN $AFILE ELOG $FSCRIPT.sh.e$JOBID"
        fi
    elif [[ $SUBC == *parallel* ]]; then
        echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.$TIMETAG.dat
        echo "RUN $AFILE OLOG $FSCRIPT.log"
    elif [[ "$SUBC" == *simple* ]] ; then
	"$FSCRIPT.sh" |& tee "$FSCRIPT.log"	
    fi
done

# Execute all FSCRIPTs locally in parallel
if [[ $SUBC == *parallel* ]]; then
    cat $LOGDIR/runscripts.$TIMETAG.dat | sort -u | $SUBC
fi

exit
