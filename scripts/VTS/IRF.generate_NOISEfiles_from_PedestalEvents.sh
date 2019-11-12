#!/bin/bash
# generate noise files from a vbf file with pedestal events
#
# noise files are used in simulations which does not include
# the MC of the NSB
#
# noise files are written in grisu-style ascii format
#
# this scripts runs locally
# 
# all relevant values are hardcoded in this scripts:
# - CARE file names
# - noise levels in MHz
# - output file names
# - length of noise traces
##############################################################
if [ $# -lt 7 ]; then
echo "
./IRF.generate_NOISEfiles_from_PedestalEvents.sh <directory with vbf pedestal files> <output directory>

   generate noise files from a vbf file with pedestal events
"
exit
fi
# DDIR="/lustre/fs19/group/cta/VERITAS/simulations/V6_FLWO/CARE_Oct1615/"
DDIR=$1
#ODIR="$VERITAS_USER_DATA_DIR/NOISE/"
ODIR=$2

# length of noise traces
TLENGTH=20000

# VBF pedestal file name
FIL="PedestalOnly_V6_PMTUpgrade_CARE_v1.6.2_11_ATM21_zen20deg_050wob_"

# GRISU Noise file name
GRI="Pedestal_V6_PMTUpgrade_CARE_v1.6.2_11_ATM21_zen20deg_"


# loop over all noise levels (in MHz) and convert 
for N in 50 75 100 130 160 200 250 300 350 400 450
do
    IFIL=${DDIR}${FIL}${N}MHz.cvbf.bz2
    OFIL=${ODIR}${GRI}${N}MHz.grisu

    echo "reading $IFIL"
    echo "writing $OFIL"

    $EVNDISPSYS/bin/VTS.NoiseFileWriter 4 499 $IFIL $OFIL $TLENGTH

done
