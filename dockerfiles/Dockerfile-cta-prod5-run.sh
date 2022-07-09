#!/bin/bash
#
# run a Eventdisplay DL1 analysis for CTA prod5 simulations
#

if [ $# -lt 1 ]; then
	echo "
	./run.sh <sim_telarray file> [layout file]
	echo
	echo [layout file (optional)] e.g., CTA.prod5S.BL-4LSTs25MSTs70SSTs-MSTF.lis
	"
        exit
fi

DATAFILE=${1}
LAYOUTFILE=${2}
# select automatically the corresponding hyperlayout
if [[ -z ${LAYOUTFILE} ]]; then
    if [[ $DATAFILE == *"paranal"* ]]; then
        LAYOUTFILE="CTA.prod5S.hyperarray.lis"
    else
        if  [[ $DATAFILE == *"prod5b"* ]]; then
            LAYOUTFILE="CTA.prod5bN.hyperarray.lis"
        else
            LAYOUTFILE="CTA.prod5N.hyperarray.lis"
        fi
    fi
fi
     
OUTPUTFILE=$(basename ${DATAFILE} .zst)
rm -f /tmp/${OUTPUTFILE}*

# calibration file
if [[ $DATAFILE == *"dark"* ]]; then
	IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-IPR.root
    if [[ $DATAFILE == *"60deg"* ]]; then
        IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-ze-60-IPR.root
    elif [[ $DATAFILE == *"40deg"* ]]; then
        IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-ze-40-IPR.root
    fi
else
	IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-halfmoon-IPR.root
    if [[ $DATAFILE == *"60deg"* ]]; then
        IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-halfmoon-ze-60-IPR.root
    elif [[ $DATAFILE == *"40deg"* ]]; then
        IPRFILE=${OBS_EVNDISP_AUX_DIR}/Calibration/prod5/prod5-halfmoon-ze-40-IPR.root
    fi
fi

###########
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

###########
# run analysis
${EVNDISPSYS}/bin/CTA.convert_hessio_to_VDST -c ${IPRFILE} \
	-a $OBS_EVNDISP_AUX_DIR/DetectorGeometry/${LAYOUTFILE} \
	-o /tmp/${OUTPUTFILE}.dst.root \
	${DATAFILE} > /tmp/${OUTPUTFILE}.convert.log

for IMAGE in lin sq2
do
    if [[ $IMAGE == "sq2" ]]; then
        EVNDISPOPT="-imagesquared"
        OFILENAME="${OUTPUTFILE}.sq2"
    else
        EVNDISPOPT="-writeimagepixellist"
        OFILENAME="${OUTPUTFILE}.lin"
    fi
    $EVNDISPSYS/bin/evndisp -averagetzerofiducialradius=0.5 \
            $EVNDISPOPT \
            -reconstructionparameter EVNDISP.prod5.reconstruction.runparameter \
            -sourcefile /tmp/"${OUTPUTFILE}".dst.root \
            -outputfile /tmp/"${OFILENAME}".root > /tmp/"${OFILENAME}".evndisp.log

    if [ -e /tmp/${OFILENAME}.root ]; then
        if [ -e /tmp/${OUTPUTFILE}.convert.log ]; then
            $EVNDISPSYS/bin/logFile convLog /tmp/${OFILENAME}.root /tmp/${OUTPUTFILE}.convert.log
        fi
        if [ -e /tmp/${OFILENAME}.evndisp.log ]; then
            $EVNDISPSYS/bin/logFile evndispLog /tmp/${OFILENAME}.root /tmp/${OFILENAME}.evndisp.log
        fi
        mv -f -v /tmp/${OFILENAME}.root /data/
    fi
done
rm -f /tmp/${OUTPUTFILE}.dst.root
