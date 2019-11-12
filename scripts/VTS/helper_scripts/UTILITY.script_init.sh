#!/bin/bash
# simple script to make sure environmental variables are set correctly

# Check that the EVNDISP environmental variables are set
if [[ -z $EVNDISPSYS ]]; then
    echo "EVNDISPSYS environmental variable not defined"
    exit 1
fi
if [[ -z $VERITAS_EVNDISP_AUX_DIR ]]; then
    echo "VERITAS_EVNDISP_AUX_DIR environmental variable not defined"
    exit 1
fi
if [[ -z $VERITAS_DATA_DIR ]]; then
    echo "VERITAS_DATA_DIR environmental variable not defined"
    exit 1
fi
if [[ -z $VERITAS_USER_DATA_DIR ]]; then
    echo "VERITAS_USER_DATA_DIR environmental variable not defined"
    exit 1
fi
if [[ -z $VERITAS_USER_LOG_DIR ]]; then
    echo "VERITAS_USER_LOG_DIR environmental variable not defined"
    exit 1
fi

# set the right observatory (environmental variables)
source $EVNDISPSYS/setObservatory.sh VTS -q

# return 0 if successful
exit 0
