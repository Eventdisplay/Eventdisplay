#!/bin/sh
#
# combine tables for CTA
#
#
#

if [ $# -ne 5 ]
then
   echo
   echo "CTA.MSCW_ENERGY.combine_tables.sh <combined table file name> <subarray list> <input table file name> <output directory> <data set> "
   echo ""
   echo "  <combined table file name>  name of the table combined file (without .root)"
   echo
   echo "  <subarray list>             text file with list of subarray IDs"
   echo
   echo "  <input table file name>     name of the input table name (beginning of...)"
   echo
   echo "  <output directory>          directory for combined tables"
   echo
   echo "   input data and output directories for tables are hard wired to"
   echo "      \$CTA_USER_DATA_DIR/analysis/AnalysisData/<dataset>/\$ARRAY/Tables/"
   echo
   exit
fi


# input parameters
TFIL=$1
VARRAY=`awk '{printf "%s ",$0} END {print ""}' $2`
ITFIL=$3
ODIR=$4
mkdir -p $ODIR
DSET=$5

source $EVNDISPSYS/setObservatory.sh CTA

######################
# loop over all arrays
for ARRAY in $VARRAY
do
   echo "combining tables for array $ARRAY"

   # loop over all array scaling
   for SCAL in 1 2 3 4 5
   do

        # input data directory with tables
        DDIR=$CTA_USER_DATA_DIR/analysis/AnalysisData/$DSET/${ARRAY}-${SCAL}/Tables/

        # temporary list file
        LISTF="table.list.temp"
#        rm -f $LISTF
        ls -1 $DDIR/$ITFIL*[0-9].root > $LISTF
        wc -l $LISTF

        echo "combining following files: " 
        cat $LISTF

        # check if combine table exist - remove it (!)
        OFIL=$ODIR/$TFIL-${ARRAY}-${SCAL}.root
#        rm -f $OFIL

        # combining files
        # (avoid woff_0500 files)
        $EVNDISPSYS/bin/combineLookupTables $LISTF $OFIL 1

        rm -f $LISTF
   done
done

exit
