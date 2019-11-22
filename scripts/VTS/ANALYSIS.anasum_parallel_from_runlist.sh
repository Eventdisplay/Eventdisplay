#!/bin/bash
# script to analyse data files with anasum (parallel analysis) from a simple run list

# EventDisplay version
$EVNDISPSYS/bin/anasum --version  >/dev/null 2>/dev/null
if (($? == 0))
then
    EDVERSION=`$EVNDISPSYS/bin/anasum --version | tr -d .`
else
    EDVERSION="g500"
fi

if [[ "$#" -lt 4 ]]; then
# begin help message
echo "
ANASUM parallel data analysis: submit jobs using a simple run list

ANALYSIS.anasum_parallel_from_runlist.sh <run list> <output directory> <cut set> <background model> [run parameter file] [mscw directory] [sim type] \
[method] [radial acceptances] [force atmosphere] [RECID]

required parameters:

    <run list>              simple runlist with a single run number per line
        
    <output directory>      anasum output files are written to this directory
                        
    <cut set>               hardcoded cut sets predefined in the script
                            (e.g., soft2tel, moderate3tel, etc.)
    
    <background model>      background model
                            (RE = reflected region, RB = ring background)
    
optional parameters:

    [run parameter file]    anasum run parameter file (located in 
                            \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/;
                            default is ANASUM.runparameter)

    [mscw directory]        directory containing the mscw.root files.
			    Default: $VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION

    [sim type]              use IRFs derived from this simulation type (GRISU-SW6 or CARE_June1425)
			    Default: CARE_June1425

    [method]                reconstruction method: GEO, DISP, FROGS.
			    Default: GEO

    [radial acceptance]     0=use external radial acceptance;
                            1=use run-wise radial acceptance (calculated from data run);
                            2=ignore radial acceptances (only for reflected region);

    [force atmosphere]	    use EAs generated with this atmospheric model (21 or 22).
			    Default: Atmosphere determined from run date for each run.				
			    Attention: Must use the same atmospere for EAs as was used for the lookup tables in the mscw_energy stage!

    <RECID>                 reconstruction ID
                            (see EVNDISP.reconstruction.runparameter)
                            Set to 0 for all telescopes, 1 to cut T1, etc.
                            

IMPORTANT! Run ANALYSIS.anasum_combine.sh once all parallel jobs have finished!

--------------------------------------------------------------------------------
"
#end help message
exit
fi

###########################
# IRFs
IRFVERSION="g500"

# Run init script
bash $(dirname "$0")"/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
RLIST=$1
ODIR=$2
CUTS=$3
BACKGND=$4
[[ "$5" ]] && RUNP=$5  || RUNP="ANASUM.runparameter"
[[ "$6" ]] && INDIR=$6 || INDIR="$VERITAS_USER_DATA_DIR/analysis/Results/$EDVERSION/"
[[ "$7" ]] && SIMTYPE=$7 || SIMTYPE="CARE_Sep0416"
[[ "$8" ]] && METH=$8 || METH="GEO"
[[ "$9" ]] && RACC=$9 || RACC="0"
[[ "${10}" ]] && FORCEDATMO=${10}
[[ "${11}" ]] && RECID=${11} || RECID=0

AUXVERSION="auxv01"

# cut definitions (note: VX to be replaced later in script)
if [[ $CUTS = *BDTmoderate* ]]; then
    CUT="NTel2-PointSource-Moderate-TMVA-BDT"
elif [[ $CUTS = *BDTsoft* ]]; then
    CUT="NTel2-PointSource-Soft-TMVA-BDT"
elif [[ $CUTS = *BDThard* ]]; then
    CUT="NTel2-PointSource-Hard-TMVA-BDT"
elif [[ $CUTS = *BDTpreselectionmoderate* ]]; then
    CUT="NTel2-PointSource-TMVA-BDT-Preselection-Moderate"
elif [[ $CUTS = *BDTpreselectionsoft* ]]; then
    CUT="NTel2-PointSource-TMVA-BDT-Preselection-Soft"
elif [[ $CUTS = *BDTpreselectionhard* ]]; then
    CUT="NTel2-PointSource-TMVA-BDT-Preselection-Hard"
elif [[ $CUTS = *NTel2-*MVA-Preselection ]]; then
    CUT=$CUTS
elif [[ $CUTS = *NTel2-*MVA-BDT ]]; then
    CUT=$CUTS
else
    echo "ERROR: unknown cut definition: $CUTS"
    exit 1
fi
CUTFILE="ANASUM.GammaHadron-Cut-${CUT}.dat"
EFFAREA="effArea-${IRFVERSION}-${METH}-${AUXVERSION}-${SIMTYPE}-Cut-${CUT}-GEOID${RECID}-VX-ATMXX-TX.root"
if [[ "$CUT" == *PointSource* ]]; then
    CUTRED="${CUT//-PointSource-/-}"
fi
if [[ "$CUT" == *ExtendedSource* ]]; then
    CUTRED="${CUT//-ExtendedSource-/-}"
fi
# radial acceptances - might be overwritten later
RADACC="radialAcceptance-${IRFVERSION}-${AUXVERSION}-${SIMTYPE}-Cut-${CUTRED}-${METH}-VX-TX.root"
# START TEMPORARY (TESTS, comment)
# RADACC="IGNOREACCEPTANCE"
# EFFAREA="IGNOREEFFECTIVEAREA"
# END TEMPORARY

echo $CUTFILE
echo $EFFAREA
echo $RADACC

# background model parameters
if [[ "$BACKGND" == *RB* ]]; then
    BM="RB"
    BMPARAMS="0.6 20"
    if [[ "$RACC" == "2" ]]; then
        echo "Error, Cannot use RB without radial acceptances:"
        echo "Specify an acceptance (external=0, runwise=1) or use RE."
        exit 1
    fi
elif [[ "$BACKGND" == *RE* ]]; then
    BM="RE"
    BMPARAMS="0.1 2 6"
else
    echo "ERROR: unknown background model: $BACKGND"
    echo "Allowed values are: RE, RB"
    exit 1
fi

# Check that run list exists
if [[ ! -f "$RLIST" ]]; then
    echo "Error, simple runlist $RLIST not found, exiting..."
    exit 1
fi

# Check that run parameter file exists
if [[ "$RUNP" == `basename $RUNP` ]]; then
    RUNP="$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/$RUNP"
fi
if [ ! -f "$RUNP" ]; then
    echo "Error, anasum run parameter file not found, exiting..."
    echo "(searched for $RUNP)"
    exit 1
fi

# directory for run scripts
DATE=`date +"%y%m%d"`
LOGDIR="$VERITAS_USER_LOG_DIR/$DATE/ANASUM.ANADATA"
echo -e "Log files will be written to:\n $LOGDIR"
mkdir -p $LOGDIR

# output directory for anasum products
echo -e "Output files will be written to:\n $ODIR"
mkdir -p $ODIR

#########################################
# make anasum run list
ANARUNLIST="$ODIR/$CUTS.anasum.dat"
rm -f $ANARUNLIST
echo "anasum run list: $ANARUNLIST"

# run list header
# echo "* VERSION 6" >> $ANARUNLIST
echo "* VERSION 7" >> $ANARUNLIST
echo "" >> $ANARUNLIST

RUNS=`cat $RLIST`

# loop over all runs
for RUN in ${RUNS[@]}; do
    # get array epoch, atmosphere and telescope combination for this run
    if [ ! -e "$INDIR/$RUN.mscw.root" ]; then
        echo "error: mscw file not found: $INDIR/$RUN.mscw.root"
        continue
    fi
    RUNINFO=`$EVNDISPSYS/bin/printRunParameter $INDIR/$RUN.mscw.root -runinfo`
    EPOCH=`echo $RUNINFO | awk '{print $(1)}'`
    ATMO=${FORCEDATMO:-`echo $RUNINFO | awk '{print $(2)}'`}
    if [[ $ATMO == *error* ]]; then
       echo "error finding atmosphere; skipping run $RUN"
       continue
    fi
    TELTOANA=`echo $RUNINFO | awk '{print "T"$(4)}'`
    echo "RUN $RUN at epoch $EPOCH and atmosphere $ATMO (Telescopes $TELTOANA)"
    echo "File $INDIR/$RUN.mscw.root"

    # do string replacements
    EFFAREARUN=${EFFAREA/VX/$EPOCH}
    EFFAREARUN=${EFFAREARUN/TX/$TELTOANA}
    EFFAREARUN=${EFFAREARUN/XX/$ATMO}

    if [[ ${RACC} == "1" ]]; then
        echo "run-wise radical acceptances: "
        #RADACCRUN="$ODIR/radialAcceptance-Cut-${CUT}-${METH}ID${RECID}-Run-${RUN}.root"
        RADACCRUN="$ODIR/$RUN.anasum.radialAcceptance.root"
        echo "   $RADACCRUN"
    elif [[ ${RACC} == "0" ]]; then
        echo "external radial acceptances: "
        RADACCRUN=${RADACC/VX/$EPOCH}
        RADACCRUN=${RADACCRUN/TX/$TELTOANA}
    else
        echo "Ignore acceptances: "
        RADACCRUN="IGNOREACCEPTANCE"
    fi

    # write line to file (VERSION 6)
#    echo "* $RUN $RUN 0 $CUTFILE $BM $EFFAREARUN $BMPARAMS $RADACCRUN="IGNOREACCEPTANCE"RADACCRUN" >> $ANARUNLIST
#    echo "* $RUN $RUN 0 $CUTFILE $BM $EFFAREARUN $BMPARAMS $RADACCRUN"
    # write line to file 
    # (VERSION 7; cuts are read from the effective area file)
    echo "* $RUN $RUN 0 $BM $EFFAREARUN $BMPARAMS $RADACCRUN" >> $ANARUNLIST
    echo "* $RUN $RUN 0 $BM $EFFAREARUN $BMPARAMS $RADACCRUN"
done

# submit the job
SUBSCRIPT="$EVNDISPSYS/scripts/VTS/ANALYSIS.anasum_parallel"
$SUBSCRIPT.sh $ANARUNLIST $INDIR $ODIR $RUNP ${RACC}

exit
