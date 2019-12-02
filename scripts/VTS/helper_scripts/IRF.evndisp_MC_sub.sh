#!/bin/bash
# script to run evndisp for simulations on one of the cluster nodes (VBF)

# set observatory environmental variables
source "$EVNDISPSYS"/setObservatory.sh VTS

########################################################
# parameters replaced by parent script using sed
RUNNUM=RUNNUMBER
SIMDIR=DATADIR
ZA=ZENITHANGLE
WOB=DECIMALWOBBLE
WOG=INTEGERWOBBLE
NOISE=NOISELEVEL
EPOCH=ARRAYEPOCH
ATM=ATMOSPHERE
ACUTS="RECONSTRUCTIONRUNPARAMETERFILE"
PARTICLE=PARTICLETYPE
SIMTYPE=SIMULATIONTYPE
ODIR=OUTPUTDIR
ANAMETHOD=ANALYSISMETHOD
NONSB=NNNOISEFILES
RHVFLAG=REDHVFLAG
# end of parameter replacement
########################################################

# Evndisp root file name (output file)
ONAME="${RUNNUM}"

if [[ $NEVENTS -gt 0 ]]; then
    ITER=$((SGE_TASK_ID - 1))
    FIRSTEVENT=$(($ITER * $NEVENTS))
    # Output file name
    ONAME="${RUNNUM}_$ITER"
    echo -e "ITER $ITER NEVENTS $NEVENTS FIRSTEVENT $FIRSTEVENT"
fi

#################################
# detector configuration and cuts

echo "Using run parameter file $ACUTS"

DEAD="EVNDISP.validchannels.dat"
# default pedestal level
# (same for GRISU and CARE,
#  adjustments possibly needed)
PEDLEV="16."
LOWPEDLEV="16."

##################################################################
# GRISU
if [[ ${SIMTYPE:0:5} == "GRISU" ]]; then
    # Input files (observe that these might need some adjustments)
    if [[ $EPOCH == "V4" ]]; then
        if [[ $PARTICLE == "1" ]]; then
           if [[ $ATM == "21" ]]; then
            VBFNAME="Oct2012_oa_ATM21_${ZA}deg_${WOG}"
           else
            VBFNAME="gamma_V4_Oct2012_SummerV4ForProcessing_20130611_v420_ATM${ATM}_${ZA}deg_${WOG}"
           fi
        elif [[ $PARTICLE == "14" ]]; then
            VBFNAME="proton_${ZA}deg_750m_wobble${WOB}_2008_2009_"
        fi
        NOISEFILE="$OBS_EVNDISP_AUX_DIR/NOISE/NOISE$NOISE.grisu.bz2"
        echo "Noise File: $NOISEFILE"
    elif [[ $EPOCH == "V5" ]]; then
        if [[ $PARTICLE == "1" ]]; then
            VBFNAME="gamma_V5_Oct2012_newArrayConfig_20121027_v420_ATM${ATM}_${ZA}deg_${WOG}"
        elif [[ $PARTICLE == "14" ]]; then
            VBFNAME="proton_${ZA}deg_w${WOB}_"
        elif [[ $PARTICLE == "402" ]]; then
            VBFNAME="helium_${ZA}deg_w${WOB}_"
        fi
        NOISEFILE="$OBS_EVNDISP_AUX_DIR/NOISE/NOISE$NOISE.grisu.bz2"
        echo "Noise File: $NOISEFILE"
    elif [[ $EPOCH == "V6" ]]; then
        if [[ $PARTICLE == "1" ]]; then
            VBFNAME="gamma_V6_Upgrade_20121127_v420_ATM${ATM}_${ZA}deg_${WOG}"
            if [[ $ATM == "21-redHV" ]]; then
                VBFNAME="gamma_V6_Upgrade_ReducedHV_20121211_v420_ATM21_${ZA}deg_${WOG}"
            elif [[ $ATM == "21-UV" ]]; then
                VBFNAME="gamma_V6_Upgrade_UVfilters_20121211_v420_ATM21_${ZA}deg_${WOG}"
            elif [[ $ATM == "21-SNR" ]]; then
                VBFNAME="gamma_V6_201304_SN2013ak_v420_ATM21_${ZA}deg_${WOG}"
            fi
        elif [[ $PARTICLE == "14" ]]; then
            VBFNAME="proton_${ZA}deg_w${WOB}_"
        elif [[ $PARTICLE == "402" ]]; then
            VBFNAME="helium_${ZA}deg_w${WOB}_"
        fi
        NOISEFILE="$OBS_EVNDISP_AUX_DIR/NOISE/NOISE${NOISE}_20120827_v420.grisu.bz2"
        echo "Noise File: $NOISEFILE"
    fi
##################################################################
# CARE
elif [ ${SIMTYPE:0:4} == "CARE" ]; then
    # input files (observe that these might need some adjustments)
    if [ "$RHVFLAG" == "yes" ]; then
        _RHV="_RHV"
    else
        _RHV=""
    fi
    [[ $PARTICLE == "1" ]]  && VBFNAME="gamma${_RHV}_${ZA}deg_750m_${WOB}wob_${NOISE}mhz_up_ATM${ATM}_part0"
    [[ $PARTICLE == "2" ]]  && VBFNAME="electron${_RHV}_${ZA}deg_noise${NOISE}MHz___"
    [[ $PARTICLE == "14" ]] && VBFNAME="proton${_RHV}_${ZA}deg_noise${NOISE}MHz___"

    # some CARE simulations need external noise files
    if [[ $NONSB == "1" ]]; then
       NOISEFILE="$OBS_EVNDISP_AUX_DIR/NOISE/Pedestal_V6_PMTUpgrade_CARE_v1.6.2_11_ATM21_zen20deg_${NOISE}MHz.grisu.bz2"
       echo "Noise File: $NOISEFILE"
       # external noise files: means CARE sims with 0 Mhz noise
       [[ $PARTICLE == "1" ]]  && VBFNAME="gamma_${ZA}deg_750m_${WOB}wob_0mhz_up_ATM${ATM}_part0"
    fi
fi
# detector configuration
[[ $EPOCH == "V4" ]] && CFG="EVN_V4_Oct2012_oldArrayConfig_20130428_v420.txt"
[[ $EPOCH == "V5" ]] && CFG="EVN_V5_Oct2012_newArrayConfig_20121027_v420.txt"
[[ $EPOCH == "V6" ]] && CFG="EVN_V6_Upgrade_20121127_v420.txt"

# temporary directory
if [[ -n "$TMPDIR" ]]; then 
    DDIR="$TMPDIR/evn_${ZA}_${NOISE}_${WOG}"
else
    DDIR="/tmp/evn_${ZA}_${NOISE}_${WOG}"
fi
mkdir -p "$DDIR"
echo "Temporary directory: $DDIR"

# loop over simulation files
if [[ ${SIMTYPE:0:5} == "GRISU" ]]; then
    VBF_FILE=$VBFNAME"wobb.vbf"
elif [[ ${SIMTYPE:0:4} == "CARE" ]]; then
    VBF_FILE="$VBFNAME.cvbf"
fi
echo 
echo "Now processing $VBF_FILE"

# unzip vbf file to local scratch directory
if [[ ! -f "$DDIR/$VBF_FILE" ]]; then
    if [[ -e "$SIMDIR/$VBF_FILE.gz" ]]; then
        echo "Copying $SIMDIR/${VBF_FILE}.gz to $DDIR"
        cp -f "$SIMDIR/$VBF_FILE.gz" $DDIR/
        echo " (vbf file copied, was gzipped)"
        gunzip -f -q "$DDIR/$VBF_FILE.gz"
    elif [[ -e "$SIMDIR/$VBF_FILE.bz2" ]]; then
        echo "Copying $SIMDIR/$VBF_FILE.bz2 to $DDIR"
        cp -f "$SIMDIR/$VBF_FILE.bz2" $DDIR/
        echo " (vbf file copied, was bzipped)"
        ls -l $DDIR/$VBF_FILE.bz2
        bunzip2 -f -v "$DDIR/$VBF_FILE.bz2"
    elif [[ -e "$SIMDIR/$VBF_FILE" ]]; then
        echo "Copying $VBF_FILE to $DDIR"
        cp -f "$SIMDIR/$VBF_FILE" $DDIR/
    fi
fi
ls -lh $DDIR/

# check that the uncompressed vbf file exists
if [[ ! -f "$DDIR/$VBF_FILE" ]]; then
    echo "No source file found: $DDIR/$VBF_FILE"
    echo "$SIMDIR/$VBF_FILE*"
    exit 1
fi
VBF_FILE="$DDIR/$VBF_FILE"

# cp and unpack NOISE file to temporary directory
if [[ $NONSB == "1" ]]; then
   NSBF=`basename $NOISEFILE .bz2`
   cp $NOISEFILE $DDIR/$NSBF.bz2
   bunzip2 $DDIR/$NSBF.bz2
   NOISEFILE=$DDIR/$NSBF
   echo "unpacked NOISE file $NOISEFILE"
fi

# calibration directory
mkdir -p $ODIR/Calibration

#######################################
# option for all steps of the analysis
MCOPT=" -runnumber=$RUNNUM -sourcetype=2 -epoch $EPOCH -camera=$CFG -reconstructionparameter $ACUTS -sourcefile $VBF_FILE -deadchannelfile $DEAD -donotusedbinfo -calibrationdirectory $ODIR"
# CARE simulations: add Gaussian noise of 3.6 mV/ (7.84 mV/dc)  / 2
# Current (2018) CARE simulations:
#    no electronic noise included - therefore add
#    Gaussian noise with the given width
#    Derived for GrIsu many years ago - source not entirely clear
#    add Gaussian noise of 3.6 mV/ (7.84 mV/dc)  / 2
if [[ ${SIMTYPE:0:4} == "CARE" ]]; then
    MCOPT="$MCOPT -injectGaussianNoise=0.229592"
fi
# TMPTMP
MCOPT="$MCOPT -nevents=5000000"

###############################################
# calculate pedestals
# (mostly for CARE,
# GRISU and some CARE sims use external noise files)
###############################################
if [[ ${SIMTYPE:0:4} == "CARE" ]] && [[ $NONSB != "1" ]]; then
    echo "Calculating pedestals for run $RUNNUM"
    rm -f $ODIR/$RUNNUM.ped.log
    PEDOPT="-runmode=1 -calibrationsumfirst=0 -calibrationsumwindow=16 -calibrationnevents=5000"
    $EVNDISPSYS/bin/evndisp $MCOPT $PEDOPT &>> $ODIR/$RUNNUM.ped.log
fi    

###############################################
# calculate tzeros/taverages
###############################################
echo "Calculating average tzeros for run $RUNNUM"
# run options for TZERO calculation
TZEROPT="-runmode=7 -calibrationsumfirst=0 -calibrationsumwindow=16  -calibrationsummin=50 -sumwindowaveragetime=6 -calibrationnevents=100000 -nevents=5000000 -pedestalnoiselevel=$NOISE"
rm -f $ODIR/$RUNNUM.tzero.log
### eventdisplay GRISU run options
if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
    TZEROPT="$TZEROPT -lowgaincalibrationfile NOFILE -lowgainpedestallevel=$PEDLEV"
else
    TZEROPT="$TZEROPT -lowgainpedestallevel=$LOWPEDLEV -lowgaincalibrationfile calibrationlist.LowGainForCare.${EPOCH}.dat"
fi
if [[ $NONSB == "1" ]]; then
    TZEROPT="$TZEROPT -pedestalfile $NOISEFILE -pedestalseed=$RUNNUM -pedestalDefaultPedestal=$PEDLEV"
fi
$EVNDISPSYS/bin/evndisp $MCOPT $TZEROPT &>> $ODIR/$RUNNUM.tzero.log

###############################################
# run eventdisplay
###############################################

#####################
# general analysis options
ANAOPT=" -writenomctree -outputfile $DDIR/$ONAME.root"

#####################
# options for GRISU (handling of low-gain values)
if [[ ${SIMTYPE:0:5} == "GRISU" ]]; then
    MCOPT="$MCOPT -simu_hilo_from_simfile -lowgaincalibrationfile NOFILE -lowgainpedestallevel=$PEDLEV"
else
#####################
# options for CARE (handling of low-gain values)
    MCOPT="$MCOPT -lowgainpedestallevel=$LOWPEDLEV -lowgaincalibrationfile calibrationlist.LowGainForCare.${EPOCH}.dat"
fi

# external NSB files
if [[ $NONSB == "1" ]]; then
    MCOPT="$MCOPT -pedestalfile $NOISEFILE -pedestalseed=$RUNNUM -pedestalDefaultPedestal=$PEDLEV"
fi

if [[ $NEVENTS -gt 0 ]]; then
    MCOPT="-nevents=$NEVENTS -firstevent=$FIRSTEVENT $MCOPT"
fi
#################################################################################
# run evndisp
echo "Analysing MC file for run $RUNNUM $ANAMETHOD"
#####################
echo "$EVNDISPSYS/bin/evndisp $MCOPT $ANAOPT " &> $ODIR/$ONAME.log
$EVNDISPSYS/bin/evndisp $MCOPT $ANAOPT &>> $ODIR/$ONAME.log

#################################################################################
# remove temporary files
ls -lh "$DDIR"
cp -f -v "$DDIR/$ONAME.root" "$ODIR/$ONAME.root"
chmod g+w "$ODIR/$ONAME.root"
chmod g+w "$ODIR/$ONAME.log"
chmod g+w "$ODIR/$ONAME.tzero.log"
chmod -R g+w $ODIR/Calibration
rm -f -v "$DDIR/$ONAME.root"
rm -f -v "$VBF_FILE"

echo "EVNDISP output root file written to $ODIR/$ONAME.root"
echo "EVNDISP log file written to $ODIR/$ONAME.log"

exit
