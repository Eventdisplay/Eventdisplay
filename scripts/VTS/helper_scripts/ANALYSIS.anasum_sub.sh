#!/bin/bash
# script to analyse one run with anasum

# set observatory environmental variables
source "$EVNDISPSYS"/setObservatory.sh VTS

# parameters replaced by parent script using sed
FLIST=FILELIST
INDIR=DATADIR
ODIR=OUTDIR
ONAME=OUTNAME
RUNP=RUNPARAM
RUNNUM=RUNNNNN
RACC=RAAACCC

#################################
# run-wise radial acceptances 
# (if requested)
if [[ ${RACC} == "1" ]]; then
   OUTPUTRACC="$ODIR/$ONAME.radialAcceptance"

   # get run information
   RUNINFO=`"$EVNDISPSYS"/bin/printRunParameter "$INDIR/$RUNNUM.mscw.root" -runinfo`
   # get instrument epoch
   EPOCH=`echo "$RUNINFO" | awk '{print $(1)}'`
   # get teltoana
   TELTOANA=`echo "$RUNINFO" | awk '{print $(4)}'`

   echo "$RUNINFO"
   echo "$EPOCH"
   echo "$TELTOANA"

   # get gamma/hadron cut from run list
   # (depend on cut file version)
   VERS=`cat "$FLIST" | grep '\*' | grep VERSION | awk '{print $3}'`
   if [[ ${VERS} == "7" ]]; then
       # cut file is an effective area file
       RCUT=`cat "$FLIST" | grep '\*' | grep "$RUNNUM" | awk '{print $6}'`
   else
       RCUT=`cat "$FLIST" | grep '\*' | grep "$RUNNUM" | awk '{print $5}'`
   fi

   # calculate radial acceptance
   "$EVNDISPSYS"/bin/makeRadialAcceptance -l "$FLIST"  \
                                        -d "$INDIR"  \
                                        -i "$EPOCH"  \
                                        -t "$TELTOANA" \
                                        -c "$RCUT" \
                                        -f "$RUNP" \
                                        -o "${OUTPUTRACC}.root" &> "${OUTPUTRACC}.log"

   # check statistics
   NEVENTS=$(cat "${OUTPUTRACC}.log" | grep entries | awk -F ": " '{print $2}')
   # check status
   STATUS=$(cat "${OUTPUTRACC}.log" | grep "STATUS=" | tail -n 1 | awk -F "=" '{print $3}' | awk -F " " '{print $1}')
   # TMP:
   STATUS="SUCCESSFUL"
   if [[ NEVENTS < 1000 ]]; then
     echo 'Number of EVENTS below the threshold (1000), using averaged radial acceptances' >> ${OUTPUTRACC}.log
     mv ${OUTPUTRACC}.root ${OUTPUTRACC}.root.lowstatistics
   fi
   # check that run-wise raidal acceptance step was successfull
   if [ "$STATUS" != "SUCCESSFUL" ]; then
     echo 'Fit status is not SUCCESSFUL, using averaged radial acceptances' >> ${OUTPUTRACC}.log
     mv ${OUTPUTRACC}.root ${OUTPUTRACC}.notsuccessful.root
   fi
fi

# introduce a random sleep to prevent many jobs starting at exactly the same time

NS=$(( ( RANDOM % 10 )  + 1 ))
sleep $NS

#################################
# run anasum
OUTPUTDATAFILE="$ODIR/$ONAME.root"
OUTPUTLOGFILE="$ODIR/$ONAME.log"
"$EVNDISPSYS"/bin/anasum   \
    -f $RUNP             \
    -l $FLIST            \
    -d $INDIR            \
    -o $OUTPUTDATAFILE   &> $OUTPUTLOGFILE
echo "RUN$RUNNUM ANPARLOG $OUTPUTLOGFILE"
echo "RUN$RUNNUM ANPARDATA $OUTPUTDATAFILE"

exit
