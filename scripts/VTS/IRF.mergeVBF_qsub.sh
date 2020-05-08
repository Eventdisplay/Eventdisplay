#!/bin/bash
# subscript to merge VBFs


# Computer time asked
#$-l h_rt=48:00:00
# Send mail (b - begin, a - abort, e - end)
#-m a
# Name
#$-N mergeVBF"
#Output
#$-o /afs/ifh.de/group/cta/scratch/pedroivo/Software/eventDisplay/LOGS/VTS/mergeVBF/merge10.o.log
#$-e /afs/ifh.de/group/cta/scratch/pedroivo/Software/eventDisplay/LOGS/VTS/mergeVBF/merge10.e.log

source ~/.edpSCT
cd /afs/ifh.de/group/cta/scratch/pedroivo/Software/eventDisplay/pSCT/bin/

PART=10
SIMSLIST="/lustre/fs19/group/cta/users/batista/VERITAS/eventDisplay/pSCT_sims/fullcamera_lowNSB/gammas/CAREvbflist_${PART}.dat" 
OUTPUT_FILE="/lustre/fs19/group/cta/users/batista/VERITAS/eventDisplay/pSCT_sims/fullcamera_lowNSB/gamma_00deg_750m_0.0wob_46mhz_up_ATM61_part${PART}.vbf"
RUNNUMBER=961200

echo "./mergeVBF" $SIMSLIST $OUTPUT_FILE $RUNNUMBER"_"$PART

./mergeVBF $SIMSLIST $OUTPUT_FILE  ${RUNNUMBER}_$PART
