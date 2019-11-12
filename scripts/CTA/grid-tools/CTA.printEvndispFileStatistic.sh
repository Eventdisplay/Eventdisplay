if [ $# != 2 ]
then
     echo "./CTA.printEvndispFileStatistic.sh <subarray list> <data directory>"
     exit
fi

VARRAY=`awk '{printf "%s ",$0} END {print ""}' $1`
DDIR=$2


printf "%-20s \t  %s \t %10s \t %10s \t\t %10s \t\t %10s \n" ARRAY Pointing gamma_onSource gamma_cone proton electron
printf "%-20s \t  %s \t %10s \t %10s \t\t %10s \t\t %10s \n" "----" "-----" "---" "-----" "-----" "----"

for ARRAY in $VARRAY
do
    for A in North South
    do
        for P in gamma_onSource gamma_cone proton electron
        do
            if [[ ${A} == "North" ]]
            then
                declare "N${P}=`find ${ARRAY}/${P} -name "*_0deg*.root" | wc -l`"
            else
                declare "N${P}=`find ${ARRAY}/${P} -name "*_180deg*.root" | wc -l`"
            fi
        done
        printf "%-20s \t  %s \t %10d \t %10d \t %10d \t %10d \n" $ARRAY $A ${Ngamma_onSource} ${Ngamma_cone}  ${Nproton} ${Nelectron}
    done

done
