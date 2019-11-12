#!/bin/bash
# script to analyse VTS files using FROGS

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RUN=RUNFILE
SETTINGS="SETSETSTRING"
CALIB=CALIBRATIONOPTION
ODIR=OUTPUTDIRECTORY
MSCWDIR=MSCWDIRECTORY
VPM=USEVPMPOINTING
LOGDIR="$ODIR"

JOBITER=$(    echo "$SETTINGS" | grep -oP "jobiter=\d+\:"    | grep -oP "\d+" )
STARTEVENT=$( echo "$SETTINGS" | grep -oP "startevent=\d+\:" | grep -oP "\d+" )
NEVENTS=$(    echo "$SETTINGS" | grep -oP "nevents=\d+\:"    | grep -oP "\d+" )
OUTPUTFILE=$( echo "$SETTINGS" | grep -oP "outputfile=.+\:"  | sed -e 's/^outputfile\=//g' -e 's/://g' )
echo "SETTINGS:  '$SETTINGS'"
echo "JOBITER:   '$JOBITER'"
echo "STARTEVENT:'$STARTEVENT'"
echo "NEVENTS:   '$NEVENTS'"
echo "OUTPUTFILE:'$OUTPUTFILE'"

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
echo "ARRAYVERS:'$ARRAYVERS'"

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

# copy everything everything that isn't the data tree,
# but copy the data tree from event # STARTEVENT to # STARTEVENT+NEVENTS
NEWMSCWNAME="$RUN.extracted.mscw.root"

echo $EVNDISPSYS/bin/frogsExtractDatafile $TEMPDIR/$INPUTMSCW "data" "$STARTEVENT" "$NEVENTS" "$TEMPDIR/$NEWMSCWNAME"

LOGEXTRACTFILE="$LOGDIR/EXTRACT.$RUN-$JOBITER.log"
$EVNDISPSYS/bin/frogsExtractDatafile $TEMPDIR/$INPUTMSCW "data" "$STARTEVENT" "$NEVENTS" "$TEMPDIR/$NEWMSCWNAME" &> $LOGEXTRACTFILE
INPUTMSCW=$NEWMSCWNAME

echo
echo "THINGS IN TEMPDIR:"
ls -1 $TEMPDIR/
echo

# the output is this file, after evndisp edits it
chmod u+w $TEMPDIR/$INPUTMSCW 

if [[ ! -e "$TEMPDIR/$INPUTMSCW" ]] ; then
	echo "Error, could not find input mscw file '$TEMPDIR/$INPUTMSCW', exiting..."
fi

#OPT+=( -frogs $MSCWDIR/$RUN.mscw.root       )
OPT+=( -frogs $TEMPDIR/$INPUTMSCW            )
OPT+=( -frogsid 0                            )
OPT+=( -templatelistforfrogs "$TEMPLATELIST" )
OPT+=( -nevents=$NEVENTS                     )
OPT+=( -firstevent=$STARTEVENT                )
#OPT+=( -display=1 )
#echo "tintin $firstevent $STARTEVENT "

FROGSPARAM="FROGS.runparameter"
if [[ ! -e "$VERITAS_EVNDISP_AUX_DIR/Frogs/$FROGSPARAM" ]] ; then
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
LOGFILE="$LOGDIR/$RUN-$JOBITER.log"
rm -f $LOGDIR/$RUN.log
$EVNDISPSYS/bin/evndisp                     \
    -runnumber=$RUN                         \
    -reconstructionparameter $ACUTS         \
    -outputfile $TEMPDIR/$RUN-$JOBITER.root \
    ${OPT[@]}                               \
    &> $LOGFILE
echo "RUN$RUN FROGSLOG $LOGFILE"

echo
echo "here's whats in '$PWD' after all processing:"
ls -l $TEMPDIR/
echo

# move data & mscw files from tmp dir to data dir
MSCWFILE="$ODIR/$RUN.mscw.root"
cp -f -v $TEMPDIR/$RUN-$JOBITER.root $OUTPUTFILE
#cp -f -v $TEMPDIR/$RUN.root $OUTPUTFILE
#cp -f -v $RUN-$JOBITER.root $OUTPUTFILE
#cp -f -v $TEMPDIR/$INPUTMSCW $MSCWFILE
rm -f $TEMPDIR/$RUN.root
rm -f $TEMPDIR/$RUN.mscw.root
echo "RUN$RUN JOBITER$JOBITER FROGSDATA $OUTPUTFILE $MSCWFILE"

exit
