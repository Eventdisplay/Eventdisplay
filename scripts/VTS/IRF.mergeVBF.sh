#!/bin/bash
# submit mergeVBF for vbf files

#qsub parameters
h_rt=48:00:00
# Send mail (b - begin, a - abort, e - end)
#-m a
# Name
NAME="mergeVBF"
LOG_DIR="/afs/ifh.de/group/cta/scratch/pedroivo/Software/eventDisplay/LOGS/VTS/mergeVBF"



for i in `seq 1 10`
do
  sed -i '12s/.*/#$-o \/afs\/ifh.de\/group\/cta\/scratch\/pedroivo\/Software\/eventDisplay\/LOGS\/VTS\/mergeVBF\/merge'$i'.o.log/' IRF.mergeVBF_qsub.sh
  sed -i '13s/.*/#$-e \/afs\/ifh.de\/group\/cta\/scratch\/pedroivo\/Software\/eventDisplay\/LOGS\/VTS\/mergeVBF\/merge'$i'.e.log/' IRF.mergeVBF_qsub.sh
  sed -i '18s/.*/PART='$i'/' IRF.mergeVBF_qsub.sh
  
  #./IRF.mergeVBF_qsub.sh
  qsub IRF.mergeVBF_qsub.sh 
  sleep 2
done    
