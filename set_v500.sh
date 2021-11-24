#!/bin/sh
#
#  set directories for Eventdisplay v4 analysis at DESY
#
#  **TEMPORARY - SHOULD BE REMOVE BEFORE MERGING WITH MAIN**
#
if [[ $# < 1 ]]; then
   echo "source ./set_v500.sh <analysis type>"
   echo
   echo "   e.g., TS, NN"
fi

EDVERSION="v505"
USERAFSDIR="/afs/ifh.de/group/cta/scratch/$USER"
USERLUSTDIR="/lustre/fs23/group/veritas/users/$USER"
GROUPLUSTDIR="/lustre/fs23/group/veritas"

########################################################################
# software directories

export ROOTSYS=/afs/ifh.de/group/cta/cta/software/root/root-6.20.04_build/
TDIR=`pwd`
cd $ROOTSYS
source ./bin/thisroot.sh
cd $TDIR
## VBF
export VBFSYS=/afs/ifh.de/group/cta/VERITAS/software/VBF-0.3.4-SL6
export PATH=$PATH:${VBFSYS}/bin/
LD_LIBRARY_PATH=$VBFSYS/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

export EVNDISPSYS=${TDIR}
export SOFASYS=${EVNDISPSYS}/sofa
export LD_LIBRARY_PATH=${EVNDISPSYS}/obj:${LD_LIBRARY_PATH}

########################################################################
# data and IRF directories

# data directory (VBF files)
export VERITAS_ANALYSIS_TYPE=${1}
export EVNDISPSCRIPTS=${USERAFSDIR}/EVNDISP/EVNDISP-400/GITHUB_Eventdisplay/Eventdisplay_AnalysisScripts_VTS/scripts/
export VERITAS_DATA_DIR=/lustre/fs24/group/veritas/
export VERITAS_USER_DATA_DIR=${USERLUSTDIR}
export VERITAS_EVNDISP_AUX_DIR=${VERITAS_USER_DATA_DIR}/analysis/Results/${EDVERSION}/${VERITAS_ANALYSIS_TYPE}/Eventdisplay_AnalysisFiles/
export VERITAS_IRFPRODUCTION_DIR=${VERITAS_USER_DATA_DIR}/analysis/Results/
export VERITAS_USER_LOG_DIR=${VERITAS_USER_DATA_DIR}/analysis/Results/${EDVERSION}/${VERITAS_ANALYSIS_TYPE}/LOGS

source ./setObservatory.sh VTS
