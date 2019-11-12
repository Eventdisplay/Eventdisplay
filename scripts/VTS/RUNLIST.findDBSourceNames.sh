#!/bin/bash
# prints list of sources that contain the first input argument

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP runlist script: find database names of VERITAS sources

RUNLIST.findDBSourceNames.sh <search string>

required parameters:

    <search string>         search database for source names containing
                            this substring

examples:
    ./RUNLIST.findDBSourceNames.sh PSR
    Searching for sources that contain 'PSR'
    PSR B0355+54
    PSR J0023+09
    PSR J0357+3205
    PSR J1023+0038
    PSR J2214+3002
    PSRJ0631+1036
    PSRJ2021+3651

    ./RUNLIST.findDBSourceNames.sh 'gro'
    Searching for sources that contain 'gro'
    GRB 2014-02-06 07:18:23 INTEGRAL 6442 GROUND
    GRO 1944+26
    MGRO J1900
    MGRO J1908+06
    MGRO J1908/SS433 Straddle
    MGROJ1909+06
    MILAGRO C1

--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
SEARCHSTR="$1"
echo "Searching for sources that contain '$SEARCHSTR'"

# get url of veritas db
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
if [ ! -n "$MYSQLDB" ] ; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit 1
fi

# do the mysql query
MYSQL="mysql -u readonly -h $MYSQLDB -A"
COUNT=0
RUNINFOARRAY=() # our list of unique sources
while read -r RUNID SOURCEID; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]]; then
		FOUND=false
		for i in "${!RUNINFOARRAY[@]}" ; do
			if [ "$SOURCEID" == "${RUNINFOARRAY[$i]}" ]; then
				FOUND=true  # we already have this source name
			fi
		done
		if ! $FOUND; then
            # we have a new source name and need to add it to the list
			RUNINFOARRAY+=("$SOURCEID")
			COUNT=$((COUNT+1))
		fi
	fi
done < <($MYSQL -e " select run_id, source_id from VERITAS.tblRun_Info where source_id like '%$SEARCHSTR%' and run_type = 'observing' and observing_mode = 'wobble' ;")

# alphabetize our source list
OLDIFS="$IFS"
IFS=$'\n' sorted=($(sort <<<"${RUNINFOARRAY[*]}"))
printf "%s\n" "${sorted[@]}"

exit
