#!/bin/sh
#
# simple script to step through different TMVA training options 
#
#

# NTrees=2000:BoostType=Grad:IgnoreNegWeightsInTraining:Shrinkage=0.1:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=20:MaxDepth=6:PruneMethod=ExpectedError

BOOST1="BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:GradBaggingFraction=0.5"
BOOST2="BoostType=Grad:Shrinkage=1.0:UseBaggedBoost:GradBaggingFraction=0.5"
BOOST3="BoostType=AdaBoostR2:AdaBoostR2Loss=Linear"
BOOST4="BoostType=AdaBoostR2:AdaBoostR2Loss=Quadratic"
BOOST5="BoostType=AdaBoostR2:AdaBoostR2Loss=Exponential"

PRUNEMETHOD1="PruneMethod=ExpectedError"
PRUNEMETHOD2="PruneMethod=NoPruning"

# depth
for D in 3 4 5 6 10
do
   for T in 2000 1000 500
   do
      for B in $BOOST1 $BOOST2 $BOOST3 $BOOST4 $BOOST5
      do
         for P in $PRUNEMETHOD1 $PRUNEMETHOD2
         do
             M="NTrees=$T:$B:nCuts=20:MaxDepth=$D:$P"
             echo $M
         done
      done
   done
done

