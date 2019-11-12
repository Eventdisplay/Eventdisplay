# small script to write a list of all scripts for all scalings, etc.

TELCOMP=( "NG" "ND" "FG" "FD" )
TELCOMP=( "NG" "FG" )
ARRAY=`cat subArray.prod3S.short.list`

OLIST=$1

rm -f $OLIST
touch $OLIST

for A in $ARRAY
do
   for ((t = 0; t < ${#TELCOMP[@]}; t++ ))
   do
      T=${TELCOMP[$t]}
      if [[ $A == "S.3HB8" ]] && [[ $T != "NG" ]]
      then
         continue;
      fi
      if [[ $A == "S.3HF8" ]] && [[ $T != "NG" ]]
      then
         continue;
      fi
      echo ${A}-${T} >> $OLIST
    done
done
       

