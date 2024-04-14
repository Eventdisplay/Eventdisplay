#!/bin/bash
#
# run a Eventdisplay DL1 analysis for CTA prod6 simulations
#

if [ $# -lt 3 ]; then
	echo "
    ./run.sh <sim_telarray file> <zenith angle, e.g., 20deg> <nsb condition (dark/moon/fullmoon)> [layout file]

	     [layout file (optional)] e.g., CTA.prod6S.Am-0LSTs14MSTs37SSTs.lis
	"
    exit
fi

DATAFILE=${1}
ZE=${2}
NSB=${3}
LAYOUTFILE=${4}
# select automatically the corresponding hyperlayout
if [[ -z ${LAYOUTFILE} ]]; then
    if [[ $DATAFILE == *"paranal"* ]]; then
        if [[ $DATAFILE == *"scts"* ]]; then
            LAYOUTFILE="CTA.prod6S.SCT.hyperarray.lis"
        else
            LAYOUTFILE="CTA.prod6S.hyperarray.lis"
        fi
    else
        LAYOUTFILE="CTA.prod6N.hyperarray.lis"
    fi
fi

OUTPUTFILE=$(basename "${DATAFILE}" .zst)
rm -f /tmp/"${OUTPUTFILE}"*

# calibration file
if [[ $NSB == *"fullmoon"* ]]; then
    IPRFILE="$CTA_EVNDISP_AUX_DIR/Calibration/prod6/prod6-full-ze${ZE}-IPR.root"
elif [[ $NSB == *"moon"* ]]; then
    IPRFILE="$CTA_EVNDISP_AUX_DIR/Calibration/prod6/prod6-half-ze${ZE}-IPR.root"
else
    IPRFILE="$CTA_EVNDISP_AUX_DIR/Calibration/prod6/prod6-dark-ze${ZE}-IPR.root"
fi

# file checks
if [[ ! -e ${IPRFILE} ]]; then
	echo "Error; IPR file not found: ${IPRFILE}"
	exit
fi
if [[ ! -e $OBS_EVNDISP_AUX_DIR/DetectorGeometry/${LAYOUTFILE} ]]; then
	echo "Error; Layout file not found: $OBS_EVNDISP_AUX_DIR/DetectorGeometry/${LAYOUTFILE}"
	exit
fi
if [[ ! -e ${DATAFILE} ]]; then
	echo "Error; sim_telarray file not found: ${DATAFILE}"
	exit
fi
echo "${IPRFILE}"

# Conversion and analysis
"${EVNDISPSYS}"/bin/CTA.convert_hessio_to_VDST -c "${IPRFILE}" \
	-a "$OBS_EVNDISP_AUX_DIR/DetectorGeometry/${LAYOUTFILE}" \
	-o "/tmp/${OUTPUTFILE}.dst.root" \
	"${DATAFILE}" > "/tmp/${OUTPUTFILE}.convert.log"

for IMAGE in lin sq2
do
    if [[ $IMAGE == "sq2" ]]; then
        EVNDISPOPT="-imagesquared"
        OFILENAME="${OUTPUTFILE}.sq2"
    else
        EVNDISPOPT="-writeimagepixellist"
        OFILENAME="${OUTPUTFILE}.lin"
    fi
    "$EVNDISPSYS"/bin/evndisp -averagetzerofiducialradius=0.5 \
            $EVNDISPOPT \
            -reconstructionparameter EVNDISP.prod6.reconstruction.runparameter \
            -sourcefile /tmp/"${OUTPUTFILE}".dst.root \
            -outputfile /tmp/"${OFILENAME}".root > /tmp/"${OFILENAME}".evndisp.log

    if [ -e "/tmp/${OFILENAME}.root" ]; then
        if [ -e "/tmp/${OUTPUTFILE}.convert.log" ]; then
            "${EVNDISPSYS}"/bin/logFile convLog "/tmp/${OFILENAME}.root" "/tmp/${OUTPUTFILE}.convert.log"
        fi
        if [ -e "/tmp/${OFILENAME}.evndisp.log" ]; then
            "${EVNDISPSYS}"/bin/logFile evndispLog "/tmp/${OFILENAME}.root" "/tmp/${OFILENAME}.evndisp.log"
        fi
        mv -f -v "/tmp/${OFILENAME}.root" /data/
    fi
done
rm -f "/tmp/${OUTPUTFILE}.dst.root"
