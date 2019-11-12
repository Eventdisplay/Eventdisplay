#!/bin/bash
#
# script to check that files where processed correctly
#
# Revision $Id$
#
#
#######################################################################

if [ ! -n "$1" ] && [ ! -n "$2" ] && [ ! -n "$3" ] && [ ! -n "$4" ] && [ ! -n "$5" ]
then
   echo "./CTA.EVNDISP.check_convert_and_analyse_MC_VDST.sh <sub array list> <list of simtelarray files> <particle> <list of failed jobs> <data set>"
   echo
   echo "  <sub array list>          text file with list of subarray IDs"
   echo
   echo "  <particle>                gamma_onSource , gamma_cone, proton , electron (helium, ...)"
   echo ""
   echo "  <list of failed jobs>     list of failed jobs" 
   echo
   echo "  <data set>         e.g. cta-ultra3, ISDC3700m, ...  "
   echo ""
   exit
fi

SUBAR=$1
BLIST=$2
PART=$3
FAILED=0
if [ -n "$4" ]
then
   FAILED=$4
fi
FAILED=$FAILED.$PART
DSET=$5

VARRAY=`awk '{printf "%s ",$0} END {print ""}' $SUBAR`
for ARRAY in $VARRAY
do
   rm -f $FAILED.$ARRAY.list
   touch $FAILED.$ARRAY.list
done

############################################################################

FILES=`cat $BLIST`

# loop over all files in files loop
for AFIL in $FILES
do

# get run number 
  F1=${AFIL#*run}
  F2=${F1%___*}
# remove zeros in the beginning in case they exist
  RUN=`zsh -c "echo $F2 | sed 's/^0\+//'"`


  for ARRAY in $VARRAY
  do

# check that simtel fil exists
#    if [ ! -e $AFIL ]
#    then
#       echo "NO SIMTEL FILE: $AFIL"
#       echo $AFIL >> $FAILED.$ARRAY.list
#       continue
#     fi
#
      DDIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/$ARRAY/$PART/"
      LDIR="$CTA_USER_LOG_DIR/analysis/AnalysisData/$DSET/$ARRAY/$PART-/"
#
# check evndisp log file
#      LFIL=`basename $AFIL .gz`
#      if [ -e $LDIR/$LFIL.evndisp.log ]
#      then
#	 LLINE=`grep "END OF ANALYSIS" $LDIR/$LFIL.evndisp.log`
#	 if [ ${#LLINE} -eq 0 ]
#	 then
#	   echo "INCOMPLETE EVNDISP RUN $LDIR/$LFIL.evndisp.log" 
#           echo $AFIL >> $FAILED.$ARRAY.list
#           continue
#	 fi
#	 LLINE=`grep -i "error" $LDIR/$LFIL.evndisp.log`
#	 if [ ${#LLINE} -ne 0 ]
#	 then
#	   echo "ERRORNESS EVNDISP RUN $LDIR/$LFIL.evndisp.log" 
#           echo $AFIL >> $FAILED.$ARRAY.list
#	   continue
#	 fi
#      else
#         echo "NO EVNDISP LOG FILE: $LDIR/$LFIL.evndisp.log"
#         echo $AFIL >> $FAILED.$ARRAY.list
#	 continue
#      fi
#
# check that evndisp output file exists
      if [ ! -e $DDIR/$RUN.root ]
      then
        echo "NO EVNDISP FILE: $DDIR/$RUN.root"
        echo $AFIL >> $FAILED.$ARRAY.list
        continue
      fi

# check that evndisp output file size is > 0
      if [ ! -s $DDIR/$RUN.root ]
      then
        echo "ZERO LENGTH EVNDISP FILE: $DDIR/$RUN.root"
        echo $AFIL >> $FAILED.$ARRAY.list
        continue
      fi
  done
   
done

exit
