#!/bin/bash
# Generates a simple run list (one run per line) with quality cuts
# written by Nathan Kelley-Hoskins Aug 2013

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
EVNDISP runlist script: generate a simple run list (one run per line)
with quality cuts for a given VERITAS source name

RUNLIST.generate.sh [allow 3-tel] [worst weather] [min duration] [start date] [end date] [runmode] [DQM category] [Elevation Bound Range] [Azimuth Bound Range]

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
    
    [Elevation Bound Range] <lower>-<upper> min and max elevation, in degrees,
                            to allow.
                            e.g. 50-70, or 70-85 or 33-88
                            default: 0-90

    [Azimuth Bound Range]   <lower>-<upper> min and max azimuth, in degrees,
                            to allow.
                            22-45  -> span of 23 degrees
                            7-330  -> span of 323 degrees
                            355-14 -> span of 19 degrees
                            just '-' for no bounds
                            default: 0-360


--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Run init script
bash "$( cd "$( dirname "$0" )" && pwd )/helper_scripts/UTILITY.script_init.sh"
[[ $? != "0" ]] && exit 1

# Parse command line arguments
[[ "$1" ]] && ALLOW_THREE_TEL=$1 || ALLOW_THREE_TEL=0
[[ "$2" ]] && WORST_WEATHER=$2   || WORST_WEATHER='B'
[[ "$3" ]] && MIN_DURATION=$3    || MIN_DURATION=15
[[ "$4" ]] && DATE_BEG="$4"      || DATE_BEG="2011-01-01"
[[ "$5" ]] && DATE_END="$5"      || DATE_END=`date +"%Y-%m-%d"`
[[ "$6" ]] && MODE="$6" || MODE="observing"
[[ "$7" ]] && DQMCATEGORY="$7" || DQMCATEGORY="science"
[[ "$8" ]] && ELEVSTRING="$8"  || ELEVSTRING=""
[[ "$9" ]] && AZIMSTRING="$9"  || AZIMSTRING=""

START_DATE="${DATE_BEG} 00:00:00"
END_DATE_STR="and db_end_time <= '${DATE_END} 00:00:00'"

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

if [[ ! "$ELEVSTRING" =~ - ]] ; then
  echo "error, 8th agument '$ELEVSTRING' needs to at least contain a '-' character, exiting..."
  exit 1
fi
if [[ ! "$AZIMSTRING" =~ - ]] ; then
  echo "error, 9th agument '$AZIMSTRING' needs to at least contain a '-' character, exiting..."
fi

MINELEV=$( echo "$ELEVSTRING" | cut -d '-' -f 1 )
MAXELEV=$( echo "$ELEVSTRING" | cut -d '-' -f 2 )
MINAZIM=$( echo "$AZIMSTRING" | cut -d '-' -f 1 )
MAXAZIM=$( echo "$AZIMSTRING" | cut -d '-' -f 2 )
[[ "$MINELEV" ]] || MINELEV=0
[[ "$MAXELEV" ]] || MAXELEV=90
[[ "$MINAZIM" ]] || MINAZIM=0
[[ "$MAXAZIM" ]] || MAXAZIM=360
#echo "MINELEV:'$MINELEV'"
#echo "MAXELEV:'$MAXELEV'"
#echo "MINAZIM:'$MINAZIM'"
#echo "MAXAZIM:'$MAXAZIM'"

# Get VERITAS database URL from EVNDISP.global.runparameter file
MYSQLDB=`grep '^\*[ \t]*DBSERVER[ \t]*mysql://' $VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter | egrep -o '[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}\.[[:alpha:]]{1,20}'`
if [ ! -n "$MYSQLDB" ]; then
    echo "* DBSERVER param not found in \$VERITAS_EVNDISP_AUX_DIR/ParameterFiles/EVNDISP.global.runparameter!"
    exit 1
fi
MYSQL="mysql -u readonly -h $MYSQLDB -A"

# get segments of time that are within our time, azimuth, and elevation bounts
#echo "submitting database request..."
TIMESEGMENTS=$( $EVNDISPSYS/bin/VTS.getObservingTimesWithinTimeAzElBounds "$DATE_BEG" "$DATE_END" "$MINELEV" "$MAXELEV" "$MINAZIM" "$MAXAZIM" | grep -P "^MJDSEGMENT" | awk '{ printf "%s %s\n", $2, $4 }' )
#echo "TIMESEGMENTS:"
#echo "$TIMESEGMENTS"

# convert our list of time segments into a mysql condition to be given to tblRun_Info, to figure out
# which runs fall into these time segments
TIMESEGCONDITION=""
segmentiter=0
while read lin ; do
  segmentiter=$((segmentiter+1))
  #echo "$segmentiter - lin:'$lin'"
  mjdstart=$( echo "$lin" | cut -d ' ' -f 1 )
  mjdend=$(   echo "$lin" | cut -d ' ' -f 2 )

  # in python's datetime.datetime.fromordinal(), 'ordinal' is the number of days since 0001-01-01, 
  # so we can just add 678576 days to our MJD to get it to 'ordinal' time, then datetime.date can
  # easily convert it to a YYYY-MM-DD string
  datebeg=$( python -c "import datetime ; d1 = datetime.date.fromordinal(678576+int($mjdstart)) ; print d1" )
  dateend=$( python -c "import datetime ; d1 = datetime.date.fromordinal(678576+int($mjdend  )) ; print d1" )

  # need to convert the decimal part of the mjd times to HH:MM:SS
  timebeg=$( python -c "import time ; mjd=$mjdstart ; day=int(mjd) ; sec=(mjd-day)*86400.0 ; hr = int(sec/(60*60)) ; min=sec-(hr*60*60) ; min=int(min/60.0) ; s=sec-(hr*60*60)-(min*60) ; print '%02d:%02d:%02d' % ( hr, min, s )" )
  timeend=$( python -c "import time ; mjd=$mjdend   ; day=int(mjd) ; sec=(mjd-day)*86400.0 ; hr = int(sec/(60*60)) ; min=sec-(hr*60*60) ; min=int(min/60.0) ; s=sec-(hr*60*60)-(min*60) ; print '%02d:%02d:%02d' % ( hr, min, s )" )

  # create this particular time segment's time condition
  tmpline="( '$datebeg $timebeg' > data_start_time and '$datebeg $timebeg' < data_end_time ) or ( '$dateend $timeend' > data_start_time and '$dateend $timeend' < data_end_time )"
  if [ ! $segmentiter -eq 1 ] ; then
    tmpline="or $tmpline"
  fi
  #echo "  '$tmpline'"
  TIMESEGCONDITION="$TIMESEGCONDITION $tmpline"
done < <(echo "$TIMESEGMENTS")
#echo "TIMESEGCONDITION:'$TIMESEGCONDITION'"

# Get run numbers from database using MySQL query
RUNINFOARRAY=()
while read -r RUNID; do
	if [[ "$RUNID" =~ ^[0-9]+$ ]]; then
		RUNINFOARRAY+=("$RUNID")
    #echo "RUNID:$RUNID"
	fi
done < <($MYSQL -e " select run_id from VERITAS.tblRun_Info where run_type LIKE \"$MODE\" and observing_mode = 'wobble' and weather <= $WEATHER and duration >= '00:${MIN_DURATION}:00' and config_mask in $TEL_MASKS and ( $TIMESEGCONDITION ) ;")


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
