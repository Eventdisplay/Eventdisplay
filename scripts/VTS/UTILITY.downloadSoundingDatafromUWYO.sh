#!/bin/bash
# download sounding (balloon) data from UWYO for VERITAS
# combine the monthly data into one file and create a list of files (in this case just the total file)

if [ ! -n "$2" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
Download sounding (balloon) data from UWYO for VERITAS, combine the monthly
data into one file and create a list of files

./UTILITY.downloadSoundingDatafromUWYO.sh <year_start> <year_end>
    
--------------------------------------------------------------------------------
"
#end help message
exit
fi

# Identifier of Tucson Airport
USM="72274"

YEAR1=$1
YEAR2=$2

MONTH=( 01 02 03 04 05 06 07 08 09 10 11 12 )
DAY=(   31 28 31 30 31 30 31 31 30 31 30 30 )
NMONTH=${#MONTH[@]}

for (( y = $YEAR1; y <= $YEAR2; y++ ))
do
    echo "Fetching $y"

    for (( i = 0 ; i < $NMONTH; i++ ))
    do
        M=${MONTH[$i]}
        D=${DAY[$i]}
        ONAME="sounding_${y}${M}"
        echo $M   $D   $ONAME

        wget --output-document=$ONAME.dat2 "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=${y}&MONTH=${M}&FROM=0100&TO=${D}12&STNM=${USM}&REPLOT=1"

        sed -e :a -e 's/<[^>]*>//g;/</N;//ba' $ONAME.dat2 > $ONAME.dat
        rm -f $ONAME.dat2

    done
done

exit
