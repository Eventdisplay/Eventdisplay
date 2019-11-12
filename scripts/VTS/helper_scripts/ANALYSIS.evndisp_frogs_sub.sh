#!/bin/bash
# script to analyse VTS files using FROGS

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RUN=RUNFILE
CALIB=CALIBRATIONOPTION
ODIR=OUTPUTDIRECTORY
MSCWDIR=MSCWDIRECTORY
VPM=USEVPMPOINTING
LOGDIR="$ODIR"

# temporary (scratch) directory
if [[ -n $TMPDIR ]]; then
    TEMPDIR=$TMPDIR/$RUN
else
    TEMPDIR="$VERITAS_USER_DATA_DIR/TMPDIR"
fi
echo "Scratch dir: $TEMPDIR"
mkdir -p $TEMPDIR

# eventdisplay reconstruction parameter (same as used in mscw file)
ACUTS=`$EVNDISPSYS/bin/printRunParameter $MSCWDIR/$RUN.mscw.root -evndispreconstructionparameterfile`

# epoch
ARRAYVERS=$($EVNDISPSYS/bin/printRunParameter $MSCWDIR/$RUN.mscw.root -epoch)

# template list file
if [[ "$ARRAYVERS" =~ ^(V5|V6)$ ]]; then
    TEMPLATELIST="EVNDISP.frogs_template_file_list.$ARRAYVERS.txt"
else
    echo "Error (helper_scripts/ANALYSIS.evndisp_frogs.sh), no frogs template list defined for \$ARRAYVERS='$ARRAYVERS', exiting..."
    exit 1
fi
echo "Using template list file '$TEMPLATELIST'..."

#########################################
# pedestal calculation
if [[ $CALIB == "1" ]]; then
    rm -f $LOGDIR/$RUN.ped.log
    $EVNDISPSYS/bin/evndisp -runmode=1 -runnumber=$RUN &> $LOGDIR/$RUN.ped.log
fi

#########################################
# read gain and toffsets from VOFFLINE DB
OPT=( -readCalibDB )

# None of the following command line options is needed for the standard analysis!

## Read gain and toff from VOFFLINE DB requiring a special version of analysis 
# OPT+=( -readCalibDB version_number )
## Warning: this version must already exist in the DB

## Read gain and toff from VOFFLINE DB and save the results in the directory
## where the calib file should be (it won't erase what is already there)
# OPT+=( -readandsavecalibdb )

## Quick look option (has no effect if readCalibDB or equivalent is set):
## use this option when you don't care about the calibration information
## if no gain can be read from your laser file, gain will be set to 1
## if no toffset can be read from your laser file, toffset will be set to 0
# OPT+=( -nocalibnoproblem )
## Note: if this option is NOT set, the analysis will break if there is problem
## reading the gain and toffset files

#########################################
# Command line options for pointing and calibration

# pointing from pointing monitor (DB)
if [[ $VPM == "1" ]]; then
    OPT+=( -usedbvpm )
fi

## pointing from db using T-point correction from 2007-11-05
# OPT+=( -useDBtracking -useTCorrectionfrom "2007-11-05" )
#
## OFF data run
# OPT+=( -raoffset=6.25 )
#
## use text file for calibration information
# OPT+=( -calibrationfile calibrationlist.dat )
#
## double pass correction
# OPT+=( -nodp2005 )

#########################################

# copy the mscw file so we don't overwrite it
INPUTMSCW="$RUN.mscw.root"
cp $MSCWDIR/$INPUTMSCW $TEMPDIR/$INPUTMSCW

# the output is this file, after evndisp edits it
chmod u+w $TEMPDIR/$INPUTMSCW 

#OPT+=( -frogs $MSCWDIR/$RUN.mscw.root       )
OPT+=( -frogs $TEMPDIR/$INPUTMSCW             )
OPT+=( -frogsid 0                            )
OPT+=( -templatelistforfrogs "$TEMPLATELIST" )
FROGSPARAM="FROGS.runparameter"
if [[ ! -e "$FROGSPARAM" ]] ; then
	echo "Error, could not find frogs parameter file '$FROGSPARAM', exiting..."
fi
OPT+=( -frogsparameterfile "$FROGSPARAM")

if [[ "$FASTDEVMODE" == "yes" ]]; then
    OPT+=( -nevents=50 )
    echo "Warning, \$FASTDEVMODE=yes, only processing the first 50 events..."
fi

echo "using frogs options '${OPT[@]}'"

#########################################
# run eventdisplay
LOGFILE="$LOGDIR/$RUN.log"
rm -f $LOGDIR/$RUN.log
$EVNDISPSYS/bin/evndisp             \
    -runnumber=$RUN                 \
    -reconstructionparameter $ACUTS \
    -outputfile $TEMPDIR/$RUN.root  \
    ${OPT[@]}                       \
    &> $LOGFILE
echo "RUN$RUN FROGSLOG $LOGFILE"

# move data & mscw files from tmp dir to data dir
DATAFILE="$ODIR/$RUN.root"
MSCWFILE="$ODIR/$RUN.mscw.root"
cp -f -v $TEMPDIR/$RUN.root $DATAFILE
cp -f -v $TEMPDIR/$INPUTMSCW $MSCWFILE
rm -f $TEMPDIR/$RUN.root
rm -f $TEMPDIR/$RUN.mscw.root
echo "RUN$RUN FROGSDATA $DATAFILE $MSCWFILE"

exit
