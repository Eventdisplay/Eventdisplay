#!/bin/bash
# script to run energy3d on batch

# set observatory environmental variables
source $EVNDISPSYS/setObservatory.sh VTS

# parameters replaced by parent script using sed
SIMTYPE=SIMULATION
TELSETUP=TELESTEUP
SUMWINT=SUMERWINTER
INDIR=INPUTDIR
OUTDIR=OUTPUTDIR

#################################
# temporary directory
if [[ -n "$TMPDIR" ]]; then 
    DDIR="$TMPDIR/energy3d_biascorr"
else
    DDIR="/tmp/energy3d_biascorr"
fi
mkdir -p $DDIR
echo "Temporary directory: $DDIR"

if [ $SIMTYPE = "CARE" ]; then
   lowsim=0
   highsim=1	
   lownoise=0
   highnoise=9
elif [ $SIMTYPE = "GRISU" ]; then
   lowsim=1
   highsim=2	
   lownoise=9
   highnoise=18	
elif [ $SIMTYPE = "CAREGRISU" ]; then
   lowsim=0
   highsim=2	
   lownoise=0
   highnoise=18	
fi 

if [ $SUMWINT = "21" ]; then
   lowsumwint=0
   highsumwint=1
   suwia=21
   suwib=""   
elif [ $SUMWINT = "22" ]; then
   lowsumwint=1
   highsumwint=2
   suwia=""
   suwib=22 
elif [ $SUMWINT = "2122" ]; then
   lowsumwint=0
   highsumwint=2
   suwia=21
   suwib=22 
fi 

if [ $TELSETUP = "4" ]; then
   lowtel=0
   hightel=1
   tela=4
   telb=""
   telc=""
elif [ $TELSETUP = "5" ]; then   
   lowtel=1
   hightel=2
   tela=""
   telb=5
   telc=""
elif [ $TELSETUP = "6" ]; then
   lowtel=2
   hightel=3
   tela=""
   telb=""
   telc=6
elif [ $TELSETUP = "456" ]; then
   lowtel=0
   hightel=3
   tela=4
   telb=5
   telc=6
fi 

TARGETS=$DDIR/filelist.txt

OUTFILENAME="energy3d_biascorrection"

if [ $SIMTYPE = "CARE" -o $SIMTYPE = "CAREGRISU" ]; then
  for AFILE in 00 20 30 35 40 45 50 55 60 65; do
    for BFILE in 50 80 120 170 230 290 370 450; do
      for CFILE in $telc; do
        for DFILE in $suwia $suwib; do
          GETFILE="V${CFILE}_ATM${DFILE}_gamma/ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz/9${CFILE}1200.root"
          NAMEINTEMP="9${CFILE}1200_ATM${DFILE}_ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz.root"
          if [ ! -e "$INDIR/$GETFILE" ] ; then
             echo "Error, unable to find 3D-modeled file: '$INDIR/$GETFILE', exiting..."
             exit 1
          fi
          cp "$INDIR/$GETFILE" "$DDIR/$NAMEINTEMP"
          echo "moved $INDIR/$GETFILE to $DDIR/$NAMEINTEMP"
        done
      done
    done
  done
fi
if [ $SIMTYPE = "GRISU" -o $SIMTYPE = "CAREGRISU" ]; then 
  for AFILE in 00 20 30 35 40 45 50 55 60 65; do
    for BFILE in 075 100 150 200 250 325 425 550 750 1000; do
      for CFILE in $tela $telb $telc; do
        for DFILE in $suwia $suwib; do
          GETFILE="V${CFILE}_ATM${DFILE}_gamma/ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz/9${CFILE}6500.root"
          NAMEINTEMP="9${CFILE}6500_ATM${DFILE}_ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz.root"
          if [ ! -e "$INDIR/$GETFILE" ] ; then
             echo "Error, unable to find 3D-modeled file: '$INDIR/$GETFILE', exiting..."
             exit 1
          fi
          cp "$INDIR/$GETFILE" "$DDIR/$NAMEINTEMP"
          echo "moved $INDIR/$GETFILE to $DDIR/$NAMEINTEMP"
        done
      done
    done
  done
fi

###############################################
# run eventdisplay
###############################################
# run options
if [ ! -e "$EVNDISPSYS/bin/create_energy3d_biascorrection" ] ; then
  echo "Error, unable to find the analysistool!"
  exit 1
fi
echo "$EVNDISPSYS/bin/create_energy3d_biascorrection $DDIR $lowsim $highsim $lowtel $hightel $suwia $suwib $lownoise $highnoise > $OUTDIR/$OUTFILENAME.txt"
      $EVNDISPSYS/bin/create_energy3d_biascorrection $DDIR $lowsim $highsim $lowtel $hightel $suwia $suwib $lownoise $highnoise > $OUTDIR/$OUTFILENAME.txt

# copying analysed temp files to outputfolder and removing the temp files
cp -f -v $DDIR/$OUTFILENAME.txt $OUTDIR/$OUTFILENAME.txt
rm -f -v $DDIR/$OUTFILENAME.txt
if [ $SIMTYPE = "CARE" -o $SIMTYPE = "CAREGRISU" ]; then
  for AFILE in 00 20 30 35 40 45 50 55 60 65; do
    for BFILE in 50 80 120 170 230 290 370 450; do
      for CFILE in $telc; do
        for DFILE in $suwia $suwib; do
          NAMEINTEMP="9${CFILE}1200_ATM${DFILE}_ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz.root"
          rm -f -v $DDIR/$NAMEINTEMP
        done
      done
    done
  done
fi
if [ $SIMTYPE = "GRISU" -o $SIMTYPE = "CAREGRISU" ]; then 
  for AFILE in 00 20 30 35 40 45 50 55 60 65; do
    for BFILE in 075 100 150 200 250 325 425 550 750 1000; do
      for CFILE in $tela $telb $telc; do
        for DFILE in $suwia $suwib; do
          NAMEINTEMP="9${CFILE}6500_ATM${DFILE}_ze${AFILE}deg_offset0.5deg_NSB${BFILE}MHz.root"
          rm -f -v $DDIR/$NAMEINTEMP
        done
      done
    done
  done
fi

echo "energy3d_biascorrection file written to $OUTDIR/$OUTFILENAME.txt"
exit
