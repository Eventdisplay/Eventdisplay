#! /bin/zsh

if [ ! -n "$1" ] 
then
	echo "Usage: make_combined_lped_file.sh <list> [calib dir]"
	exit
fi

list=$1

dir=${2-$VERITAS_EVNDISP_AUX_DIR/Calibration} 

echo $list
echo $dir

if [ ! -d $dir/lpedfiles ]
then
	mkdir $dir/lpedfiles 
fi

if [ -e temp.lped ] 
then
	rm temp.lped
fi
# touch temp.lped
for run in `cat $list`
do
	if [ ! -e $dir/lpedfiles/$run.lped ]
	then
		mutate_lped.pl $run $dir
	fi
	# cat $dir/lpedfiles/$run.lped >> temp.lped
	names+=($dir/lpedfiles/$run.lped )
	# echo $names
done
echo "mutated"
# cat temp.lped | sort -n temp.lped
sort -n -o temp.lped $names 
echo "sorted"
# swap modules 18 and 26 in that one run

sed -i -e "s/\t/ /g" temp.lped
echo "spaced"

for ((i=0; i<10; i++)) 
do
	sed -i -e "s/26 $i 6947576 3 2/18XX $i 6947576 3 2/" temp.lped
	sed -i -e "s/18 $i 6947576 3 2/26XX $i 6947576 3 2/" temp.lped
done

sed -i -e "s/XX//" temp.lped

echo "swapped" 
# remove some runs.
# 7083132 T4: completely wrong window.
# 5769394, 5811314, 5880607 T2: Pulse very late -> not reliable for lpeds
# 5423738, 5423941, 5508384, 5625859, 5625860 -> 64 sample readout
# 7116263, 7217071, 7563739 7639596 7679293 T4 flasher to low
# 7348488 T2,T3 flasher too low
# 7422021 event loss
# 7468182 T1, T2 flasher too low
sed -i -e "/7083132 3/d" -e "/5769394 1/d" -e "/5811314 1/d" -e "/5880607 1/d" \
	-e "/5423738/d" -e "/5423941/d" -e "/5508384/d" -e "/5625859/d" -e "/5625860/d" \
	-e "/7116263 3/d" -e "/7217071 3/d" -e "/7563739 3/d" -e "/7639596 3/d" -e "/7679293 3/d" \
	-e "/7348488 1/d" -e "/7348488 2/d" -e "/7468182 0/d" -e "/7468182 1/d" -e "/7422021 /d"  temp.lped

echo "removed runs"
# remove bogus lped of L2 channels 
# T3, 443 became parasitic L2 channel after 52186
sed -i -e "/^190 3 [56789]...... 2 443/d" temp.lped 
# T1 channels 110/255/404, T2, c173, T3 c 37/159/319, T4 c 99/259 became L2 channel after 36254
sed -i -e "/^142 7 [56789]...... 2 37/d" temp.lped 
sed -i -e "/^63 3 [56789]...... 1 173/d" temp.lped 
echo "removed L2"
# remove channels that have no pulse in LG
# what about  "203 9", "89 3", "192 3", "24 3", "38 8"
# 111 5: T1 c15 -> very small pulse? Fixed in summer 2014 (run 74220)
# 196 4: T1 c294 -> swapped out in summer 2014
#  37 1: T1 c441 Fixed
#  34 4: T1 c464 Fixed
# 104 2: T2 c452 Fixed
# 108 8: T2 c328, T3 c188 Fixed
# 177 2: T3 c312 Fixed
# 169 1: T3 c371 Fixed
# 142 8: T3 c38 broke between runs 7309495 and 7348488.

awk '{if ( \
	( $1==111 && $2==5 && $3 < 7422000 ) || \
	( $1==196 && $2==4 && $3 > 5000000 ) || \
	( $1==37 && $2==1 &&  $3 < 7422000 ) || \
	( $1==34 && $2==4 &&  $3 < 7422000 ) || \
	( $1==104 && $2==2 &&  $3 < 7422000 ) || \
	( $1==108 && $2==8 &&  $3 < 7422000 ) || \
	( $1==177 && $2==2 && $3>5709596 && $3 < 7422000 ) || \
	( $1==169 && $2==1 && $3 < 7422000 ) || \
	( $1==142 && $2==8 && $3>7309495  ) \
   ) \
{ print $1, $2, $3, $4, $5, 999, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26; } 
else {	print $0 } }' temp.lped > LowGainPedestals.lped
echo "awked" 

exit




