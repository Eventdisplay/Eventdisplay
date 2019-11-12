#!/bin/bash
# download sounding (balloon) data from UWYO for VERITAS
# combine the monthly data into one file and create a list of files (in this case just the total file)

if [ ! -n "$1" ] || [ "$1" = "-h" ]; then
# begin help message
echo "
Download sounding (balloon) data from UWYO for VERITAS, combine the monthly
data into one file and create a list of files

UTILITY.downloadSoundingDatafromUWYO.sh <year>
    
--------------------------------------------------------------------------------
"
#end help message
exit
fi

YEAR=$1

MONTH=( 01 02 03 04 05 06 07 08 09 10 11 12 )
DAY=(   31 28 31 30 31 30 31 31 30 31 30 30 )
NMONTH=${#MONTH[@]}

for (( i = 0 ; i < $NMONTH; i++ ))
do
    M=${MONTH[$i]}
    D=${DAY[$i]}
    ONAME="sounding_$YEAR$M"
    echo $M   $D   $ONAME

    wget --output-document=$ONAME.dat2 "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=${YEAR}&MONTH=${M}&FROM=0100&TO=${D}12&STNM=72274&REPLOT=1"

    sed -e :a -e 's/<[^>]*>//g;/</N;//ba' $ONAME.dat2 > $ONAME.dat
    rm -f $ONAME.dat2
done

exit
