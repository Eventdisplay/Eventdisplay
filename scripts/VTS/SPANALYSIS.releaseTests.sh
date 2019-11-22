# scripts to test analysis chain for many fields, cuts, etc.
#
# script used for developing code and testing
#
#  many hard wired parameters

# check run parameters

if [ $# -lt 1 ]; then
    echo
    echo "./SPANALYSIS.releaseTests.sh <analysis step>"
    echo
    echo " script for release testing (note: several hard wired parameters)"
    echo
    echo "   EVNDISP, MSCW, ANASUM, ANASUM_SUB"
    echo
    exit
fi

RUNTYPE="$1"
# ANALYSIS METHOD (GEO or FROGS)
[[ "$2" ]] && ANAMETHOD=$2 || ANAMETHOD="GEO"

# evndisp version
VERSION="g500"
# ATMOSPHERE
ATM="ATM61"
# run list
LIST="runlist_releaseTesting"
# evndisp run parameter files


ERUN="EVNDISP.reconstruction.runparameter.LL"
ANAMETHOD="TL5035MA20.LL"

ERUN="EVNDISP.reconstruction.runparameter"
ANAMETHOD="TL5035MA20"

ERUN="EVNDISP.reconstruction.runparameter.MA10"
ANAMETHOD="TL5035MA10"

ERUN="EVNDISP.reconstruction.runparameter.TL5025"
ANAMETHOD="TL5025MA20"

ERUN="EVNDISP.reconstruction.runparameter.TL5025MA10"
ANAMETHOD="TL5025MA10"

ERUN="EVNDISP.reconstruction.runparameter.NN"
ANAMETHOD="NN"

ERUN="EVNDISP.reconstruction.runparameter"
ANAMETHOD="TL5035MA20"

SIMSUB=""

VERSIONSUB="${VERSION}${SIMSUB}"
# CUTS
if [[ $ANAMETHOD = "FROGS" ]]
then
    CUTS=( "frogs25" "frogs43" )
elif [[ $ANAMETHOD = "BDT" ]]
then
    CUTS=( "BDT-moderate" )
else
    CUTS=( "soft2tel" "softopen" "soft4002tel" "soft3002tel" "moderate2tel" "moderate3tel" )
    CUTS=( "soft2tel" "softopen" )
    CUTS=( "soft2tel" "hard3tel" "moderate2tel" "soft3tel" "moderate3tel" )
    CUTS=( "soft2tel" "moderate3tel" "hard3tel" "moderate2tel" "superhard" )
    CUTS=( "soft2tel" "moderate2tel" "soft3tel" "moderate3tel" "hard3tel" )
    CUTS=( "soft3tel" "soft2tel" )
    CUTS=( "soft3tel" "soft2tel" "moderate2tel" "moderate3tel" "hard3tel" )
    CUTS=( "BDTpreselection4tel" "BDTpreselection3tel" "BDTpreselection2tel" )
    CUTS=( "BDTpreselection2tel" )
    CUTS=( "NTel2-PointSource-Moderate-MVA-Preselection" "NTel2-PointSource-Soft-MVA-Preselection" )
    CUTS=( "NTel2-PointSource-Soft-MVA-Preselection" )
    CUTS=( "NTel2-PointSource-Moderate-MVA-Preselection" )
    CUTS=( "NTel2-PointSource-SuperSoft-MVA-BDT" "NTel2-PointSource-SuperSoft2-MVA-BDT" )
    CUTS=( "NTel2-PointSource-SuperSoft-MVA-BDT" )
    CUTS=( "NTel2-PointSource-Soft-MVA-BDT" "NTel2-PointSource-SuperSoft-MVA-BDT" )
    CUTS=( "NTel2-PointSource-SuperSoft-MVA-Preselection" )
    CUTS=( "NTel2-PointSource-SuperSoft-MVA-BDT" )
    CUTS=( "NTel2-PointSource-SuperSoft-MVA-Preselection" "NTel2-PointSource-Soft-MVA-Preselection" )
fi

# V4
SIMTYPE="GRISU-SW6"
EPOCH="V4"


# V5
SIMTYPE="GRISU-SW6"
EPOCH="V5"
SIMTYPE="CARE_Mar1816${SIMSUB}"
SIMTYPE="CARE_Apr1417${SIMSUB}"
SIMTYPE="CARE_Apr1419${SIMSUB}"
EPOCH="V5"

# V6
SIMTYPE="CARE_June1702${SIMSUB}"
EPOCH="V6"

# directories
# PREFIX="$VERITAS_DATA_DIR/ANADATA-v500RC/analysis/Results/${VERSION}/"
PREFIX="$VERITAS_USER_DATA_DIR/analysis/Results/${VERSION}/"

# PKS1424p240: V4, V5, V6
# Tycho: V6
# Crab: V4, V5, V6
# SS 433: V4, V5

# object
# V6
# for O in Crab Tycho PKS1424p240 Mrk501 Mrk421 HESSJ0632057
# V5
# for O in Crab Tycho PKS1424p240
# V4
# for O in PKS1424p240 Mrk501 Mrk421 HESSJ0632057
# for O in 1ES1959p650
# for O in  PKS1424p240 Mrk421 HESSJ0632057 Tycho
# for O in Crab PKS1424p240 PKS1441p25
# for O in Crab PKS1424p240 Tycho
# for O in Tycho PKS1424p240 Mrk501 Mrk421 HESSJ0632057 MAXIJ1820p070
# for O in PKS1424p240 Mrk501 Mrk421 MAXIJ1820p070
# for O in Tycho HESSJ0632057 
# for O in BDTtraining
for O in Segue1
# for O in BDTtraining
# for O in PKS1424p240 Mrk421 Mrk501 MAXIJ1820p070
# for O in PKS1424p240
do
    echo $O
    echo

    # run list
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_50deg.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.LZE.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.HF.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_Winter.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.HF.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_Winter.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.LZE.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_20162018.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.TMVAOptimization.V2.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.redo.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_20182019.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}_20182019.dat"
    RLIST="$PREFIX/$O/${LIST}${EPOCH}.dat"

    # output directory
    ODIR="$PREFIX/$O/${VERSIONSUB}_${ANAMETHOD}"

    mkdir -p $ODIR
    ########### START EVNDISP
    if [[ $RUNTYPE = "EVNDISP" ]]
    then
       $(dirname "$0")/ANALYSIS.evndisp.sh $RLIST $ODIR $ERUN 1 $ANAMETHOD
       continue
    fi
    ########### END EVNDISP
    # reconstruction ID
    for ID in 0
    do
        if [[ $ANAMETHOD = "FROGS" ]]
        then
            IDIR="$PREFIX/${O}/${D}/${VERSIONSUB}_${ANAMETHOD}/"
        else
            IDIR="$PREFIX/${O}/${D}/${VERSIONSUB}_${ANAMETHOD}/RecID${ID}_${SIMTYPE}"
        fi
        ########### START MSCW
        if [[ $RUNTYPE = "MSCW" ]]
        then
           TFIL=table-$VERSION-${ANAMETHOD}-auxv01-${SIMTYPE}-${ATM}-${EPOCH}-RECMETHOD${ID}
           # TMP
           #TFIL=table-g500-${ANAMETHOD}-auxv01-${SIMTYPE}-${ATM}-${EPOCH}-RECMETHOD${ID}
           $(dirname "$0")/ANALYSIS.mscw_energy.sh $TFIL $RLIST $ODIR $ID $IDIR
            continue
        fi
        ########### END MSCW
        ########### START FROGS
        if [[ $RUNTYPE = "FROGS" ]]
        then
           $(dirname "$0")/ANALYSIS.evndisp_frogs.sh $RLIST $PREFIX/${O}/${D}/${VERSIONSUB}_${ANAMETHOD}/ $IDIR
        fi

        ############ ANASUM ###############
        ##
        # cut
        NCUTS=${#CUTS[@]}
        for (( m = 0; m < $NCUTS; m++ ))
        do
              C=${CUTS[$m]}
              for BC in RE
              do
                    # ouptut directory for anasum files
                    ODIR="$PREFIX/${O}/${VERSIONSUB}_${EPOCH}_anasum_${ANAMETHOD}/${ANAMETHOD}_ID${ID}_${C}_${BC}_N8"
                    # runparameter file
                    RUNPAR="$PREFIX/${O}/runparameter.dat"
                    
                    echo $IDIR
                    echo $ODIR
                    echo $RUNPAR
                    echo $RLIST

                    rm -f $ODIR.log
                    if [[ $RUNTYPE = "ANASUM" ]]
                    then
                         $EVNDISPSYS/bin/anasum -i 1 -d $ODIR -f $RUNPAR -l $ODIR/$C.anasum.dat -o $ODIR.root > $ODIR.log
                    elif [[ $RUNTYPE = "ANASUM_SUB" ]]
                    then
                        $(dirname "$0")/ANALYSIS.anasum_parallel_from_runlist.sh ${RLIST} $ODIR $C ${BC} $RUNPAR $IDIR $SIMTYPE $ANAMETHOD 1 61 $ID
                    fi
                    echo "DONE $ODIR"
               done
        done
    done
done
