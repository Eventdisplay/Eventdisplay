# optimise BDT cuts
#
# This script does:
# 
# 1. write particle rate files for the given preselection cuts
# 2. optimise cuts and write cut values
#
# Prerequesites:
# - anasum results on the Crab for BDT preselection cuts (expected filename: anasumCombined.root)
# - effective area files with the same preselection cuts

# directory with simulation results
SIMDIR="$VERITAS_USER_DATA_DIR/analysis/Results/v502/CARE_June1702/V6_ATM61_gamma_TL5035MA20/"
EFFFIL="EffArea-CARE_June1702-V6-ID0-ZENITANGLE-0.5wob-NOISE-Cut-"
# directory with results from Crab for BDT preselection cuts
DDIR="$VERITAS_USER_DATA_DIR/analysis/Results/v502/Crab/v502_V6_anasum_TL5035MA20/"

# directory with BDTs
BDTDIR="/lustre/fs19/group/cta/users/maierg/EVNDISP_AnalysisFiles.svn/trunk/VTS/GammaHadron_BDTs/V6/soft/"

# for C in NTel2-PointSource-SuperSoft-MVA-Preselection NTel2-PointSource-SuperSoft2-MVA-Preselection
for C in NTel2-PointSource-SuperSoft2-MVA-Preselection
do
    echo "processing $C"

    OFIL="$BDTDIR/TL5035MA20_ID0_${C}_RE_N8_MC"

    echo 
    echo "  particle file: $OFIL"

    mv -f ${OFIL}_MVA.log ${OFIL}_MVA.log.back
    "$EVNDISPSYS"/bin/writeParticleRateFilesForTMVA "$DDIR/TL5035MA20_ID0_${C}_RE_N8/anasumCombined.root" "$OFIL.root" "$SIMDIR/EffectiveAreas_Cut-${C}" "$EFFFIL${C}.root" 200 > "$OFIL.log"
    echo
    echo "  optimised BDT cuts"

    OCUTS="$BDTDIR/"

    mv -f "${OFIL}_MVA.log" "${OFIL}_MVA.log.back"
    "$EVNDISPSYS"/bin/optimiseBDTCuts "$OFIL.root" "$BDTDIR" mva "${OFIL}_MVA.root" 20. > "${OFIL}_MVA.log"

    echo
    echo "optmised MVA cuts written to ${OFIL}_MVA.root"
    echo "--> please check and move into the cuts directory"
done 

