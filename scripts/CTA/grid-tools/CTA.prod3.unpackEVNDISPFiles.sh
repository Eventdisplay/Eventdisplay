#!/bin/sh
#
# script to unpack GRID produced EVNDISP tar files into a analysis directory
# hardwired number of arrays
#
# e.g. CTA.prod3.unpackEVNDISPFiles.sh sorted.diff.6.list prod3-paranalp05-40deg-NN "1 2 3 4 5"
#
#############################################################################

if [ $# -ne 3 ]
then
     echo
     echo "./CTA.prod3.unpackEVNDISPFiles.shCTA.prod3.unpackEVNDISPFiles.sh <list of tar files> <site> <scaling>"
     echo
     exit
fi

FIL="$1"
SITE="$2"
SCALING="$3"

####################################
# files will unpacked in this directory
DDIR="/lustre/fs9/group/cta/users/maierg/CTA/analysis/AnalysisData/$SITE/"
mkdir -p $DDIR

SARRAY="3HB1 3HB2 3HB3 3HB4 3HD1 3HD2 3HI1"
TELTYP="FG FD NG ND FA NA"
SARRAY="3HB89"
TELTYP="FG FD NG ND"
PART="gamma_onSource gamma_cone electron proton"


###########
# temporary directory for unpacking
# TDIR="$DDIR/tmp"
TDIR="/lustre/fs18/group/cta/prod3/paranal-40deg/Analysis/"
mkdir -p $TDIR

for P in ${PART}
do
  for SCAL in ${SCALING}
  do
    # make sub array directories
    for V in $SARRAY
    do
        for T in $TELTYP
        do
            SDIR="S.${V}-${T}-${SCAL}"
            mkdir -p $DDIR/$SDIR/${P}
        done
    done

    # get list of files for this particular scaling
    if [[ $P == "gamma_onSource" ]]
    then
        FILES=`cat $FIL | grep gamma | grep -v diffuse`
#        FILES=`cat $FIL | grep deg-${SCAL} | grep gamma | grep -v diffuse`
    elif [[ $P == "gamma_cone" ]]
    then
#        FILES=`cat $FIL | grep deg-${SCAL} | grep gamma | grep diffuse`
        FILES=`cat $FIL | grep gamma | grep diffuse`
    else
 
#        FILES=`cat $FIL | grep deg-${SCAL} | grep $P`
        FILES=`cat $FIL | grep $P`
    fi

    echo "unpacking scaling $SCAL"
    echo "======================="
    for DF in $FILES
    do
        A="/acs/grid/cta${DF}"
        # the tar file should exist on the dcache
        # (disable for performance reasons)
#        if [ ! -e $A ]
#        then
#           continue
#        fi

        # check if file already exists
        # (e.g. 1002-180deg-2_evndisp.tar.gz)
        # check only one array layout
        B=`basename $A`
        AZ=`echo $B| cut -d'-' -f 2`
        # get RUN
        RUN=${B%%-*}

#        V="3HB1-FD"
#        F="$DDIR/S.3HB1-FD-${SCAL}/${P}/${RUN}.${AZ}-${V}-${SCAL}.root"
#        echo $F

        # now check if file exists, if not copy it over and unpack it
#        if [ ! -e $F ]
#        then
         echo $A
         dccp $A - | tar -xzf -

         ls -l

        # mv files into subarray directories
          for VA in $SARRAY
          do
            for TA in $TELTYP
            do
                SDIR="S.${VA}-${TA}-${SCAL}"

#                DAT="$TDIR/*${VA}-${TA}-${SCAL}*"
                DAT="*${VA}-${TA}*"

                echo $DDIR/$SDIR/${P}/

                mv -v $DAT $DDIR/$SDIR/${P}/
                done
          done
#          rm -f $TDIR/$B
#         fi
    done
  done   
done
