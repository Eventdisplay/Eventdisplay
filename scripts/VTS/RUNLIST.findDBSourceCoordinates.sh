#!/bin/bash

if [ ! "$#" -eq "1" ] || [ "$1" = "-h" ] ; then
	echo "Print the coordinates of a source" ; echo
	echo "`basename $0` <sourcename>" ; echo
	echo "<source name>  : The name of the source, as stored in VERITAS.tblObserving_Sources"  
	echo "                 or from \$EVNDISPSYS/scripts/VTS/RUNLIST.findDBSourceCoordinates.sh" ; echo 
	echo "Examples:" ; echo
	echo "  Print the Crab's positon:"
	echo "  $ `basename $0` \"Crab\""
	echo "  83.633349 22.014475" ; echo
	echo "  Print Markarian 501's position:"
	echo "  $ `basename $0` \"Mrk501\""
	echo "  253.467529 39.760300" ; echo
	echo "Source name should be from the database table VERITAS.tblObserving_Sources ,"
	echo "Prints the RA/DEC (in degrees, J2000 epoch) of the sourcename, as stored in the database."
	echo "This is the source coordinate veritas uses for pointing the telescopes." ; echo
	echo "In event display's ANASUM.runparameter file, the output of this "
	echo "  can be pasted after the 'TARGETPOSITION_RADECJ2000_DEG' option"
	echo "  to specify the analysis target."
	exit 0
fi

# Parse command line arguments
SOURCENAME="$1"
#echo "Searching for sources that contain '$SOURCENAME'"

# get url of veritas db
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
if [ ! -n "$MYSQLDB" ] ; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit 1
fi

# do the mysql query
MYSQL="mysql -u readonly -h $MYSQLDB -A"
while read -r EPOCH RA DEC NAME; do
	#echo "$EPOCH $RA $DEC $NAME"
	if [[ "$EPOCH" == "2000" ]] ; then
		SRCRADEG=$(  bc -l <<< "$RA  * 180.0 / 3.141592" ) # convert radians to degrees
		SRCDECDEG=$( bc -l <<< "$DEC * 180.0 / 3.141592" )
		SRCRADEG=$(  printf "%9.6f" $SRCRADEG  )
		SRCDECDEG=$( printf "%9.6f" $SRCDECDEG )
		echo "$SRCRADEG $SRCDECDEG"
		exit 0
	fi
done < <($MYSQL -e " select epoch, ra, decl, source_id from VERITAS.tblObserving_Sources where source_id = '$SOURCENAME' ;" ) 

exit 0

