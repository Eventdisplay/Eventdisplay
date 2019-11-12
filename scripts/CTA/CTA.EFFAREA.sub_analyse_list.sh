#!/bin/bash
#
# submit jobs for effective area calculation
#
# particle names and directories fixed by CTA setup
#
#

if [ $# -lt 6 ] 
then
   echo ""
   echo "./CTA.EFFAREA.sub_analyse_list.sh <subarray list> <cutfile template> <analysis parameter file> <output subdirectory> <data set> [filling mode] [qsub options] [direction (e.g. _180deg)]"
   echo
   echo "<subarray list>"
   echo "     text file with list of subarray IDs"
   echo
   echo "<cutfile template>"
   echo "     template for gamma/hadron cut file"
   echo "     (suffix must be .gamma/.CRbck ; this will be added by this script)"
   echo "     examples can be found in $CTA_EVNDISP_AUX_DIR/GammaHadronCutFiles"
   echo 
   echo "<analysis parameter file>"
   echo "     file with analysis parameter"
   echo "     examples can be found in $CTA_EVNDISP_AUX_DIR/ParameterFiles/"
   echo
   echo "<output subdirectory>"
   echo "    sub-directory name (not the full path) for effective areas files"
   echo
   echo "<data set>         e.g. cta-ultra3, ISDC3700m, ...  "
   echo
   echo "[filling mode]"
   echo "     effective area filling mode (use 2 to calculate angular resolution only)"
   echo "     (no value for default IRF calculation)"
   echo
   echo "[direction]        e.g. for north: \"_180deg\", for south: \"_0deg\", for all directions: no option"
   echo
   echo ""
   exit
fi

SUBAR=$1
CDIR="$CTA_EVNDISP_AUX_DIR/GammaHadronCutFiles/"
CFIL=$2
ANAPAR=$3
ODIR=$4
DSET=$5
GMOD=0
if [ -n "$6" ]
then
  GMOD=$6
fi
MCAZ=""
if [ -n "$8" ]
then
  MCAZ=$8
fi
QSUBOPT=""
if [ -n $7 ]
then
   QSUBOPT="$7"
fi
QSUBOPT=${QSUBOPT//_X_/ } 
QSUBOPT=${QSUBOPT//_M_/-} 
QSUBOPT=${QSUBOPT//\"/} 

#######################################
# read values from parameter file
if [ ! -e $ANAPAR ]
then
  echo "error: analysis parameter file not found: $ANAPAR" 
  exit
fi
echo "reading analysis parameter from $ANAPAR"
# get output file directory
EFFAREADIR=`grep EFFAREASUBDIR $ANAPAR | awk {'print $2'}`
ODIR="$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/EffectiveAreas/$EFFAREADIR/$4/"
mkdir -v -p $ODIR
# get reconstruction ID
RECID=`grep RECID $ANAPAR | awk {'print $2'}`
# observation time
OBSTIME=`grep OBSERVINGTIME_H $ANAPAR | awk {'print $2'}`

#arrays
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $SUBAR`


#################################################
# directories
DATE=`date +"%y%m%d"`
echo "directory for qsub shell scripts"
QSHELLDIR=$CTA_USER_DATA_DIR"/queueShellDir"
echo $QSHELLDIR
mkdir -p $QSHELLDIR
QDIR=$CTA_USER_LOG_DIR"/$DATE/EFFAREA/"
mkdir -p $QDIR
# QDIR="/dev/null"
echo "qsub error log $QDIR"

#################################################
# set particle types 
# (don't expect to have cone for all data sets)
if [ $GMOD = "0" ] || [ $GMOD = "3" ]
then
   VPART="ALL"
else
   VPART="GAMMA"
fi

#########################################
# loop over all arrays
#########################################
for ARRAY in $VARRAY
do
   echo "STARTING ARRAY $ARRAY"

###########################################
# prepare cut files
   CCUT=$ODIR/$CFIL.$ARRAY
   cp $CDIR/$CFIL.gamma.dat $CCUT.gamma.dat
   cp $CDIR/$CFIL.CRbck.dat $CCUT.CRbck.dat

###########################################
# prepare run script

# skeleton script
      FSCRIPT="CTA.EFFAREA.qsub_analyse_list"

# create run script
      FNAM="CTAeffArea-$DSET-$VPART-${ARRAY}${MCAZ}-ID${RECID}-${OBSTIME}"

      sed -e "s|PPPANAPAR|$ANAPAR|" \
          -e "s|PPPARRAY|$ARRAY|" \
          -e "s|PPPRECID|$RECID|" \
          -e "s|PPPVPART|$VPART|" \
          -e "s|PPPCCUT|$CCUT|" \
          -e "s|PPPODIR|$ODIR|" \
          -e "s|PPPDSET|$DSET|" \
          -e "s|PPPGMOD|$GMOD|" \
          -e "s|PPPMCAZ|$MCAZ|" $FSCRIPT.sh > $QSHELLDIR/$FNAM.sh

      chmod u+x $QSHELLDIR/$FNAM.sh

      echo $QSHELLDIR/$FNAM.sh

###########################################
# submit the job script
     if [ $4 = "AngularResolution" ]
     then
         qsub $QSUBOPT -l os="sl6" -l h_cpu=11:29:00 -l h_rss=3000M -l tmpdir_size=15G  -V -o $QDIR -e $QDIR "$QSHELLDIR/$FNAM.sh"
     else
         if [ $DSET = "*LaPalma*" ]
         then
             qsub $QSUBOPT -l os="sl6" -l h_cpu=11:29:00 -l h_rss=3000M -l tmpdir_size=15G  -V -o $QDIR -e $QDIR "$QSHELLDIR/$FNAM.sh"
         else
             qsub $QSUBOPT -l os="sl6" -l h_cpu=45:29:00 -l h_rss=3000M -l tmpdir_size=15G  -V -o $QDIR -e $QDIR "$QSHELLDIR/$FNAM.sh"
         fi
     fi
done

exit
