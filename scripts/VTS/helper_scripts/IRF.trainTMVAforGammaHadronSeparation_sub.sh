#!/bin/bash
# script to train BDTs with TMVA (VTS)

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
RXPAR=RUNPARAM
NXTRAIN=NUMTRAIN
ONAME=OUTNAME
          
{
    $EVNDISPSYS/bin/trainTMVAforGammaHadronSeparation $RXPAR.runparameter > $RXPAR.log
} || {
    echo "Replacing"
    sed -i "s/=${NXTRAIN}/=0/g" $RXPAR.runparameter
    $EVNDISPSYS/bin/trainTMVAforGammaHadronSeparation $RXPAR.runparameter >> $RXPAR.log
}

# remove unnecessary *.C files
CDIR=`dirname $RXPAR`
rm -f -v $CDIR/$ONAME*.C

exit
