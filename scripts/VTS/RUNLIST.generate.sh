#!/bin/bash
# Generates a simple run list (one run per line) with quality cuts
# written by Nathan Kelley-Hoskins Aug 2013

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP runlist script: generate a simple run list (one run per line)
with quality cuts for a given VERITAS source name

RUNLIST.generate.sh <source name> [allow 3-tel] [worst weather] [min duration] [start date] [end date] [runmode] [DQM category]

required parameters:

    <source name>           create simple run list for this source

optional parameters:

    [allow 3-tel]           flag to include three-telescope runs
                            (default = 0 = no; set to 1 for yes)
                            
    [worst weather]         select runs with weather >= this letter grade
                            (default: \"B\" weather)
                            
    [min duration]          minimum run duration (default: 15 minutes)
    
    [start date]            select all runs on or after this date
                            (default: 2011-01-01, format = YYYY-MM-DD)

    [end date]              select all runs on or before this date
                            (default: none, format = YYYY-MM-DD)

    [runmode]		    use runs taken with this runmode
			    (default: observing = regular runs). 
			    Use 'obsLowHV' or 'obsFilter' 
                            to select only reduced HV/only filter runs. 
			    Use '%' for all runs.

    [DQM category]	    use runs taken with this DQM category.
			    Default: science. Other options: 
			    calibration, engineering, moonfilter, reducedhv, special.
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
SOURCE_NAME=$1
[[ "$2" ]] && ALLOW_THREE_TEL=$2 || ALLOW_THREE_TEL=0
[[ "$3" ]] && WORST_WEATHER=$3   || WORST_WEATHER='B'
[[ "$4" ]] && MIN_DURATION=$4    || MIN_DURATION=15
[[ "$5" ]] && START_DATE=$5" 00:00:00" || START_DATE="2011-01-01 00:00:00"
[[ "$6" ]] && END_DATE_STR="and db_end_time <= '$6 00:00:00'"
[[ "$7" ]] && MODE="$7" || MODE="observing"
[[ "$8" ]] && DQMCATEGORY="$8" || DQMCATEGORY="science"

# Configuration for number of telescopes
if [[ $ALLOW_THREE_TEL == 0 ]]; then
    # four telescope only
    TEL_MASKS="('15')"      # used for VERITAS.tblRun_info
    TEL_CUT_MASKS="('0')"   # used for VOFFLINE.tblRun_Analysis_Comments
else
    # three telescope configs
    TEL_MASKS="('15', '7', '11', '13', '14')"
    TEL_CUT_MASKS="('0', '8', '4', '2', '1')"
fi

# Configuration for worst weather
WEATHER_GRADES=('A+' 'A' 'A-' 'B+' 'B' 'B-' 'C+' 'C' 'C-' 'D+' 'D' 'D-' 'F')
for (( i=0; i < ${#WEATHER_GRADES[@]}; i++ )); do
    if [[ ${WEATHER_GRADES[$i]} == "$WORST_WEATHER" ]]; then
        WEATHER=$i
    fi
done

# Get VERITAS database URL from EVNDISP.global.runparameter file
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
if [ ! -n "$MYSQLDB" ]; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit 1
fi

# Get run numbers from database using MySQL query
MYSQL="mysql -u readonly -h $MYSQLDB -A"
RUNINFOARRAY=()
while read -r RUNID; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]]; then
		RUNINFOARRAY+=("$RUNID")
	fi
done < <($MYSQL -e " select run_id from VERITAS.tblRun_Info where source_id = '$SOURCE_NAME' and run_type LIKE \"$MODE\" and observing_mode = 'wobble' and weather <= $WEATHER and duration >= '00:${MIN_DURATION}:00' and db_start_time >= '$START_DATE' $END_DATE_STR and config_mask in $TEL_MASKS ;")

# check if VERITAS.tblRun_Info had 0 runs for us
if (( ${#RUNINFOARRAY[@]} <= 0 )) ; then
	# if so, error out and tell the user why
	echo "Error, no runs fit current conditions, try loosening arguments.  Exiting..." 2>&1
	exit 1
fi

# Convert RUNINFOARRAY to a comma-separated tuple
RUN_IDS=$(IFS=, ; echo "(${RUNINFOARRAY[*]})")

# Do some final quality checks using the VOFFLINE.tblRun_Analysis_Comments table
FINALARRAY=()
while read -r RUNID; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]] ; then
		FINALARRAY+=("$RUNID")
		echo "$RUNID"
	fi
done < <($MYSQL -e "select run_id from VOFFLINE.tblRun_Analysis_Comments where status = 'good_run' and (tel_cut_mask is NULL or tel_cut_mask in $TEL_CUT_MASKS) and ( data_category like \"$DQMCATEGORY\" or ( \"$DQMCATEGORY\" = \"%\" and data_category is null )  ) and usable_duration >= '00:${MIN_DURATION}:00' and run_id in ${RUN_IDS[@]}")

exit
