#!/bin/bash
# script to train TMVA (BDTs) for angular reconstruction

# set observatory environmental variables
source "$EVNDISPSYS"/setObservatory.sh VTS

# parameters replaced by parent script using sed
INDIR=EVNDISPFILE
ODIR=OUTPUTDIR
ONAME=BDTFILE
BDTTARGET=BDTMETHOD
# train
rm -f "$ODIR/$ONAME*"

ls "$INDIR"/*[0-9].root > INPUTLIST.txt

# fraction of events to use for training,
# remaining events will be used for testing
TRAINTESTFRACTION=0.5

"$EVNDISPSYS"/bin/trainTMVAforAngularReconstruction INPUTLIST.txt "$ODIR" "$TRAINTESTFRACTION" 0 1 "$BDTTARGET" > "$ODIR/$ONAME.log"

rm INPUTLIST.txt

exit
