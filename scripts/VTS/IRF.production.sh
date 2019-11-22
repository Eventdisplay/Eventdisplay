#!/bin/bash
# IRF production script (VERITAS)
#

if [ $# -lt 2 ]; then
# begin help message
echo "
IRF generation: produce a full set of instrument response functions (IRFs)

IRF.production.sh <sim type> <IRF type> [epoch] [atmosphere] [Rec ID] [cuts list file] [sim directory]

required parameters:

    <sim type>              original VBF file simulation type (e.g. GRISU-SW6, CARE_June1425)
    
    <IRF type>              type of instrument response function to produce
                            (e.g. EVNDISP, MAKETABLES, COMBINETABLES,
                             ANALYSETABLES, EFFECTIVEAREAS, COMBINEEFFECTIVEAREAS,
                             TRAINMVANGRES )
    
optional parameters:
    
    [epoch]                 array epoch(s) (e.g., V4, V5, V6)
                            (default: \"V4 V5 V6\")
                            
    [atmosphere]            atmosphere model(s) (21 = winter, 22 = summer)
                            (default: \"21 22\")
                            
    [Rec ID]                reconstruction ID(s) (default: \"0\")
                            (see EVNDISP.reconstruction.runparameter)

    [cuts list file]        file containing one gamma/hadron cuts file per line
                            (default: hard-coded standard EventDisplay cuts)

    [sim directory]         directory containing simulation VBF files

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash $(dirname $0)"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
SIMTYPE=$1
IRFTYPE=$2
[[ "$3" ]] && EPOCH=$3 || EPOCH="V6 V5 V4"
[[ "$4" ]] && ATMOS=$4 || ATMOS="21 22"
# [[ "$5" ]] && RECID=$5 || RECID="0 2 3 4 5"
[[ "$5" ]] && RECID=$5 || RECID="0"
[[ "$6" ]] && CUTSLISTFILE=$6 || CUTSLISTFILE=""
[[ "$7" ]] && SIMDIR=$7 || SIMDIR=""

#ANAMETHOD="NNMA20"
#EVNDISPRUNPARAMETER="EVNDISP.reconstruction.runparameter.NN"

#ANAMETHOD="TL5035MA20"
#EVNDISPRUNPARAMETER="EVNDISP.reconstruction.runparameter"

# If ANAMETHOD is not among env variables, set a default value
ANAMETHOD="NN"
ANAMETHOD="TL5025MA10"

if [ "$ANAMETHOD" == "" ]; then
    ANAMETHOD="TL5035MA20"
    EVNDISPRUNPARAMETER="EVNDISP.reconstruction.runparameter"
else
    EVNDISPRUNPARAMETER="EVNDISP.reconstruction.runparameter.$ANAMETHOD"
fi

# EventDisplay version (default is g502)
$EVNDISPSYS/bin/printRunParameter --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`$EVNDISPSYS/bin/printRunParameter --version | tr -d .`
else
    EDVERSION="g502"
fi
# version string for aux files
AUX="${ANAMETHOD}-auxv01"

# simulation types and definition of parameter space
if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
    # GrISU simulation parameters
    ZENITH_ANGLES=( 00 20 30 35 40 45 50 55 60 65 )
    NSB_LEVELS=( 075 100 150 200 250 325 425 550 750 1000 )
    WOBBLE_OFFSETS=( 0.5 0.00 0.25 0.75 1.00 1.25 1.50 1.75 2.00 )
    ZENITH_ANGLES=( 20 30 )
    NSB_LEVELS=( 150 250 )
    WOBBLE_OFFSETS=( 0.5 )
elif [ ${SIMTYPE:0:4} = "CARE" ]; then
    # CARE simulation parameters
    WOBBLE_OFFSETS=( 0.5 )
    NSB_LEVELS=( 50 75 100 130 160 200 250 300 350 400 450 )
    if [ "${RHVFLAG}" == "yes" ]; then
        NSB_LEVELS=( 150 300 450 600 750 900 )
    fi
    ZENITH_ANGLES=( 00 20 30 35 40 45 50 55 )
    NSB_LEVELS=( 50 75 100 130 160 200 250 300 350 400 450 )
   # NSB_LEVELS=( 200 )
    if [[ $ATMOS == "62" ]]; then
          ZENITH_ANGLES=( 00 20 30 35 )
    fi
else
    echo "Invalid simulation type. Exiting..."
    exit 1
fi

# Set gamma/hadron cuts
if [[ $CUTSLISTFILE != "" ]]; then
    if [ ! -f $CUTSLISTFILE ]; then
        echo "Error, cuts list file not found, exiting..."
        echo $CUTSLISTFILE
        exit 1
    fi
    # read file containing list of cuts
    CUTLIST=$(IFS=$'\r\n'; cat $CUTSLISTFILE)
else
    # default list of cuts
    #CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-MVA-Preselection.dat"
    #CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-MVA-BDT.dat"
    #CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-MVA-Preselection.MSCW.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Hard-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Soft-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-TMVA-BDT-Preselection.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-TMVA-BDT-Preselection-Soft.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-TMVA-BDT-Preselection-Moderate.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-TMVA-BDT-Preselection-Hard.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Hard-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Soft-TMVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-TMVA-BDT-Preselection.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-MVA-Preselection.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Soft-MVA-Preselection.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Moderate-MVA-Preselection.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-Soft-MVA-BDT.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-SuperSoft2-MVA-BDT.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-SuperSoft-MVA-BDT.dat"
    CUTLIST="ANASUM.GammaHadron-Cut-NTel2-PointSource-SuperSoft-MVA-Preselection.dat
             ANASUM.GammaHadron-Cut-NTel2-PointSource-Soft-MVA-Preselection.dat"

fi
CUTLIST=`echo $CUTLIST |tr '\r' ' '`

############################################################
# loop over complete parameter space and submit production
for VX in $EPOCH; do
    for ATM in $ATMOS; do
       ######################
       # set lookup table file name
       TABLECOM="table-${EDVERSION}-${AUX}-${SIMTYPE}-ATM${ATM}-${VX}-"
       ######################
       # combine lookup tables
       if [[ $IRFTYPE == "COMBINETABLES" ]]; then
            TFIL="${TABLECOM}"
            for ID in $RECID; do
                echo "combine lookup tables"
                METH="RECMETHOD${ID}"
                $(dirname "$0")/IRF.combine_lookup_table_parts.sh "${TFIL}${METH}" $VX $ATM $ID $SIMTYPE $ANAMETHOD
            done
            continue
       fi
       ######################
       # combine effective areas
       if [[ $IRFTYPE == "COMBINEEFFECTIVEAREAS" ]]; then
            for ID in $RECID; do
                for CUTS in ${CUTLIST[@]}; do
                    echo "combine effective areas $CUTS"
                   $(dirname "$0")/IRF.combine_effective_area_parts.sh $CUTS $VX $ATM $ID $SIMTYPE $AUX $ANAMETHOD
                done # cuts
            done
            continue
        fi
        for ZA in ${ZENITH_ANGLES[@]}; do
            ######################
            # train MVA for angular resolution
            if [[ $IRFTYPE == "TRAINMVANGRES" ]]; then
               if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
                   $(dirname "$0")/IRF.trainTMVAforAngularReconstruction.sh $VX $ATM $ZA 200 $SIMTYPE
               elif [[ ${SIMTYPE:0:4} = "CARE" ]]; then
                   $(dirname "$0")/IRF.trainTMVAforAngularReconstruction.sh $VX $ATM $ZA 200 $SIMTYPE
               fi
               continue
            fi
            for NOISE in ${NSB_LEVELS[@]}; do
                for WOBBLE in ${WOBBLE_OFFSETS[@]}; do
                    echo "Now processing epoch $VX, atmo $ATM, zenith angle $ZA, wobble $WOBBLE, noise level $NOISE"
                    ######################
                    # run simulations through evndisp
                    if [[ $IRFTYPE == "EVNDISP" ]]; then
                       if [[ ${SIMTYPE:0:5} = "GRISU" ]]; then
                          SIMDIR=$VERITAS_DATA_DIR/simulations/"$VX"_FLWO/grisu/ATM"$ATM"
                       elif [[ ${SIMTYPE:0:4} = "CARE" ]]; then
                          SIMDIR=$VERITAS_DATA_DIR/simulations/"$VX"_FLWO/${SIMTYPE}
                       fi
                        $(dirname "$0")/IRF.evndisp_MC.sh $SIMDIR $VX $ATM $ZA $WOBBLE $NOISE $SIMTYPE $EVNDISPRUNPARAMETER 0 1 $ANAMETHOD
                    ######################
                    # make tables
                    elif [[ $IRFTYPE == "MAKETABLES" ]]; then
                        for ID in $RECID; do
                           $(dirname "$0")/IRF.generate_lookup_table_parts.sh $VX $ATM $ZA $WOBBLE $NOISE $ID $SIMTYPE $ANAMETHOD
                        done #recid
                    ######################
                    # analyse table files
                    elif [[ $IRFTYPE == "ANALYSETABLES" ]]; then
                        TFIL="${TABLECOM}"
                        for ID in $RECID; do
                            #METH="RECMETHOD${ID}"
                            METH="RECMETHOD0"
                            $(dirname "$0")/IRF.mscw_energy_MC.sh "${TFIL}${METH}" $VX $ATM $ZA $WOBBLE $NOISE $ID $SIMTYPE 1 $ANAMETHOD
                        done
                    ######################
                    # analyse effective areas
                    elif [[ $IRFTYPE == "EFFECTIVEAREAS" ]]; then
                        for ID in $RECID; do
                            $(dirname "$0")/IRF.generate_effective_area_parts.sh "$CUTLIST" $VX $ATM $ZA $WOBBLE $NOISE $ID $SIMTYPE $ANAMETHOD
                        done #recID
                    fi
                done #wobble
            done #noise
        done #ZA
    done #ATM
done  #VX

exit

# /IRF.trainTMVAforGammaHadronSeparation.sh /lustre/fs19/group/cta/users/maierg/VERITAS/analysis/Results/g502c/BDTtraining/g502_TL5035MA20/RecID0_CARE_June1702/BDTTraining.bck.list $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDT.runparameter $VERITAS_USER_DATA_DIR/test mva  CARE_June1702 V6 61 0 TL5035MA20
#  ./IRF.trainTMVAforGammaHadronSeparation.sh $VERITAS_USER_DATA_DIR/analysis/Results/g502/BDTtraining/g502_TL5035MA20/RecID0_CARE_June1702/BDTTraining.bck.list $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/TMVA.BDT.runparameter $VERITAS_USER_DATA_DIR/test mva  CARE_June1702 V6 61 0 TL5035MA20 $VERITAS_USER_DATA_DIR/analysis/Results/g502/CARE_June1702/V6_ATM61_gamma_TL5035MA20/MSCW_RECID0
