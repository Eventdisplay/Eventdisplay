#!/bin/bash
#Submit IRF.evndisp_MC.sh using more than one .vbf
SIMDIR=/lustre/fs19/group/cta/users/batista/VERITAS/eventDisplay/pSCT_sims/fullcamera_lowNSB
EPOCH=V6
ATM=61
ZEN=00
WOB=0.0
NSB=46
SIMTYPE=CARE
RECFILE=EVNDISP.reconstruction.runparameter.pSCT
EXTNSB=0
PARTICLE=1
ANAMETHOD="TL"
NEVENTS=-1
RECID=0

for i in `seq 1 10`
do
   PART=$i
   ./IRF.evndisp_MC.sh $SIMDIR $EPOCH $ATM $ZEN $WOB $NSB $SIMTYPE $RECFILE $EXTNSB $PARTICLE $ANAMETHOD $NEVENTS $PART >| IRF.evndisp_MC_part${PART}.log
done 
