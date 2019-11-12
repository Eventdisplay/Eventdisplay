#!/bin/bash

# figure out ctools fits converter's needed effective area files
# and submit them
# output effective area files will be stored in
#  $VERITAS_EVNDISP_AUX_DIR/EffectiveAreas/noThetaCuts/

# qsub parameters
h_cpu=00:29:00; h_vmem=8000M; tmpdir_size=10G

echo
echo "Warning, this script is untested, and will probably only work at DESY, Germany."
echo

if [[ $# != 4 ]]; then
# begin help message
echo "
IRF generation: create partial effective area files from MC ROOT files
 (simulations that have been processed by both evndisp_MC and mscw_energy_MC)

IRF.makeCtoolsEffectiveAreas.sh <runlist> <anasum parallel directory> <chunkcode> <force>

required parameters:

    <runlist>               file containing list of runnumbers
    
    <anasum directory>      directory with #####.anasum.root files

    <chunkcode>             specify the maximum duration to break chunks into
                            ex: 20s - 20 seconds,  5m - 5 minutes,  1h - 1 hour
    
    <force>                 if 'forceoverwrite', all existing needed effective area parts will 
                            be deleted before being recreated, otherwise will only
                            generate effective area parts that don't exist
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname "$0")"/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1


RUNLISTFILE=$1
ANASUMDIR=$2
CHUNK=$3
FORCE=$4
if [ -e "$RUNLISTFILE" ] ; then
  echo "using runlist '$RUNLISTFILE'"
else
  echo "error, could not find runlist '$RUNLISTFILE', exiting..."
  exit 1
fi
if [ -d "$ANASUMDIR" ] ; then
  echo "using anasum dir '$ANASUMDIR'"
else
  echo "error, could not find anasum directory '$ANASUMDIR', exiting..."
  exit 1
fi
if [ -z "$CHUNK" ] ; then
  echo "error, need to define a chunk code as the 3rd argument, exiting..."
  exit 1
fi
if [ "$FORCE" = "forceoverwrite" ] ; then
  echo "forcably overwriting existing effective area parts..."
else
  echo "will only create needed effective area parts..."
fi

# EventDisplay version
EDVERSION=`$EVNDISPSYS/bin/makeEffectiveArea --version | tr -d .`
VERITAS_IRFPRODUCTION_DIR=/lustre/fs5/group/cta/VERITAS/IRFPRODUCTION/

EFFAREALIST=""

RUNLIST=`cat $RUNLISTFILE`
IRUN=0
for ARUN in $RUNLIST ; do
  IRUN=$((IRUN+1))
  ARUNFILE="$ANASUMDIR/$ARUN.anasum.root"
  if [ ! -e "$ARUNFILE" ] ; then
    echo "error, for run $ARUN, unable to locate anasum file '$ARUNFILE', exiting..."
    exit 1
  fi
  TMPTEXT=$( $EVNDISPSYS/bin/writeCTAEventListFromAnasum -f -i "$ARUNFILE" -o "/dev/null" -c $CHUNK -r -e 2>&1 | grep "requires" | awk '{ $1="" ; print $0 }' | tr -d "'" | sed 's/^ *//' )
  #echo "$TMPTEXT"
  if [ $IRUN -ge 2 ] ; then
    EFFAREALIST="$EFFAREALIST"$'\n'
  fi
  EFFAREALIST="${EFFAREALIST}${TMPTEXT}"
  NFILES=$( echo "$TMPTEXT" | wc -l )
  echo "iter $IRUN - run $ARUN - '$ARUNFILE' - requires $NFILES effarea files"
  
done

# sort and remove duplicate effective area files
EFFAREALIST=$( echo "$EFFAREALIST" | sort | uniq )
echo "EFFAREALIST: "

# loop over each needed file, check if it exists,
# if not, create it
for AFILE in $EFFAREALIST ; do
  echo
  echo "examining file '`basename $AFILE`'..."
  if [ -e "$AFILE" ] ; then
    if [ "$FORCE" = "forceoverwrite" ] ; then
      echo "   deleting before recreating"
      rm -rf "$AFILE"
    else
      echo "   skipping, already exists..."
      continue
    fi
  fi
  
  # scan for needed info
  EDV=$(  echo "$AFILE" | grep -oP "v\d{3}" )
  ZA=$(   echo "$AFILE" | grep -oP "Ze\d{2}deg"        | grep -oP "\d+"      )
  WOB=$(  echo "$AFILE" | grep -oP "\-\d.\d{1,2}wob\-" | grep -oP "\d.\d{1,2}" )
  NOS=$(  echo "$AFILE" | grep -oP "wob\-\d{2,4}\-Cut" | grep -oP "\d+" )
  ATMO=$( echo "$AFILE" | grep -oP "ATM\d{2}"          )
  CUTS=$( echo "$AFILE" | grep -oP "EffArea.*$"        | grep -oP "Cut.*$" | cut -d '.' -f1 )
  EPOCH=$(echo "$AFILE" | grep -oP "EffArea.*$"        | grep -oP "V\d" )
  SIMS=$( echo "$AFILE" | grep -oP "EffArea.*$"        | grep -oP "(GRISU-SW6|CARE_June1425)" )
  RECID=$(echo "$AFILE" | grep -oP "EffArea.*$"        | grep -oP "ID\d" | grep -oP "\d" )
  #echo "   edver:$EDV"
  #echo "   zenang:$ZA"
  #echo "   wobble:$WOB"
  #echo "   noise:$NOS"
  #echo "   cuts:'$CUTS'"
  #echo "   atmo:$ATMO"
  #echo "   epoch:$EPOCH"
  #echo "   sims:'$SIMS'"
  #echo "   recid:$RECID"
  
  # output directory for storing effective area parts
  ODIR="$VERITAS_IRFPRODUCTION_DIR/$EDV/$SIMS/${EPOCH}_${ATMO}_gamma/"
  mkdir -p $ODIR
  echo "   ODIR:'$ODIR'"
  if [ ! -d "$ODIR" ] ; then
    echo "error, could not find/create directory '$ODIR', exiting..."
    exit 1
  fi

  # name of our output effective area part
  EFFAREAFILE="`basename $AFILE`"
  EFFAREAFILE=$( echo "$EFFAREAFILE" | sed 's/-Cut.*$//' )
  echo "   EFFAREAFILE:'$EFFAREAFILE'"

  # figure out which input mscw file to make the effective area from
  MCFILE="$VERITAS_IRFPRODUCTION_DIR/$EDV/$SIMS/${EPOCH}_${ATMO}_gamma/MSCW_RECID${RECID}/${ZA}deg_${WOB}wob_NOISE${NOS}.mscw.root"
  echo "   MCFILE:'$MCFILE'"
  if [ ! -e "$MCFILE" ] ; then
    echo "error, could not find mc file '$MCFILE', exiting..."
    exit 1
  fi
  
  # figure out which cuts file to use
  CUTSFILE="ANASUM.GammaHadron-${CUTS}.dat"
  echo "   CUTSFILE:'$CUTSFILE'"
  if [ ! -e "$VERITAS_EVNDISP_AUX_DIR/GammaHadronCutFiles/$CUTSFILE" ] ; then
    echo "error, could not find cuts file '$CUTSFILE', exiting..."
    exit 1
  fi

  # run scripts and output are written into this directory
  DATE=`date +"%y%m%d"`
  LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/EFFAREA/${ZA}deg_${WOBBLE}wob_NOISE${NOISE}_${EPOCH}_ATM${ATM}_${PARTICLE_TYPE}_${RECID}/"
  LOGDIR="$PWD/logs/"
  echo -e "   Log files will be written to:'$LOGDIR'"
  mkdir -p $LOGDIR

  # Job submission script
  SUBSCRIPT="$EVNDISPSYS/scripts/VTS/helper_scripts/IRF.effective_area_parallel_sub"

  # set parameters in run script
  FSCRIPT="$LOGDIR/EA.ID${RECID}.$DATE.MC"
  sed -e "s|OUTPUTDIR|$ODIR|"      \
    -e   "s|EFFFILE|$EFFAREAFILE|" \
    -e   "s|DATAFILE|$MCFILE|"     \
    -e   "s|GAMMACUTS|${CUTSFILE}|" $SUBSCRIPT.sh > $FSCRIPT.sh
  
  chmod u+x $FSCRIPT.sh
  #echo $FSCRIPT.sh

  # run locally or on cluster
  SUBC=`$EVNDISPSYS/scripts/VTS/helper_scripts/UTILITY.readSubmissionCommand.sh`
  SUBC=`eval "echo \"$SUBC\""`
  if [[ $SUBC == *qsub* ]]; then
      JOBID=`$SUBC $FSCRIPT.sh`
      echo "JOBID: $JOBID"
  elif [[ $SUBC == *parallel* ]]; then
      echo "$FSCRIPT.sh &> $FSCRIPT.log" >> $LOGDIR/runscripts.dat
  fi

done


