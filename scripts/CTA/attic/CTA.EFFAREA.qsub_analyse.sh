#!/bin/bash
#
# Analysing effective areas (CTA)
#
#
###################################################################################
# VALUES SET BY MOTHER SCRIPT
###################################################################################
# input file list 
MSCFIL=IIIIFIL
# effective area file (output)
OFIL=TTTTFIL
###################################################################################
# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh CTA

###################################################################################
# run effective area code
$EVNDISPSYS/bin/makeEffectiveArea $MSCFIL $OFIL.root > $OFIL.log
###################################################################################

exit
