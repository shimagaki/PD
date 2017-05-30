#!/bin/bash
#
# Observe the abuse of bash!
# TODO: use a real scripting language instead of this mess!

cleanup()
{
		rm "$SHM_FILE"
		rm "$FIFO"
}

cleanup_kill()
{
		trap - INT
		kill %2 2>/dev/null && wait
		echo "-1" >&6
		cleanup
}

crit_helicity()
{
	echo `echo -e "scale=8\n2/(4*a(1))*$1" | bc -l`
}

clean_Tc()
{
	kill %2
	wait
	echo $T
}

find_Tc()
{ 
	size=`echo "2^$1" | bc`
	./xymodel -s $size "$SHM_FILE" > /dev/null &
	sleep 1
	#wait for xymodel to load in background 
	# binary search
	lower=0.5
	upper=1.2
	T=`echo "scale=8; ($upper+$lower)/2" | bc`
	./mcrun -t $T $INIT_ITER $SHM_FILE > /dev/null
	for i in `seq 1 $MAX_SEARCH`; do
		#echo "[$lower, $upper]" >&2
		T=`echo "scale=8; ($upper+$lower)/2" | bc`
		./mcrun -t $T $INIT_ITER $SHM_FILE > /dev/null
		H=`./mcrun -t $T -p plugin/helicity.so $AVG_ITER $SHM_FILE | tr "\n" "\t" | awk "{print "'$'"1}"`
		Hc=`crit_helicity $T`
		#echo -e "$H\t$Hc" >&2
		if [[ "$H" == "$Hc" || "$lower" == "$upper" ]]; then
				clean_Tc
				return 0
		fi

		if [[ `echo "scale=9; $H < $Hc" |bc` -eq 1 ]]; then
			upper=$T
		else
			lower=$T
		fi
		echo "scale=9; $size*$size/2+$i/$MAX_SEARCH*($size*$size*3/4)" |bc >&6
	done

	clean_Tc
}

## Change #################
first=1             # start size exponent (size = 2^$first)
last=4              # up to this size exponent
AVG_ITER=30000       # average helicity over this many iterations
INIT_ITER=20000     # iterations to reach equilibrium, when temperature changes.
MAX_SEARCH=15       # max. binary search iterations
###########################

if [[ $# != 1 ]]; then
	echo "usage: $0 OUT_FILE"
	echo "output format: SIZE	Tc"
	exit 1
fi

SHM_FILE=`mktemp` || exit 1
FIFO=".progress_fifo."`basename $SHM_FILE`
trap 'cleanup_kill; exit 1' INT

OUT_FILE="$1"


> "$OUT_FILE"

rm -f $FIFO
mkfifo $FIFO
exec 6<> $FIFO
sizeFirst=`echo "2^(2*$first)/2" | bc`
sizeLast=`echo "2^(2*$last)" | bc`
./progress -e $FIFO $sizeLast $sizeFirst &

for exp in `seq $first $last`; do
	Tc=`find_Tc $exp`
	size=`echo "2^$exp" | bc`
	echo "$size $Tc" >> "$OUT_FILE"
done

echo $last >&6
