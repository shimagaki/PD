#!/bin/bash
#
# Write energy and magnetization per spin and variance to OUT_FILE

if [[ $# != 2 ]]; then
	echo "usage: $0 OUT_FILE1 OUT_FILE2"
	exit 1
fi

## Change these variables as needed ###
Tinc=0.02
Tstart=0.01
Tend=2

ITERATIONS=20000 # Iterations to reach
                 # equilibrium

AVG_ITER=300000	 # Measurement
SIZE=32          # Lattice size
#######################################

SHM_FILE=`mktemp`
trap 'kill %1; wait; rm $SHM_FILE; exit 1' INT
OUT_FILE="$1"

./xymodel -s $SIZE "$SHM_FILE" &
sleep 1

END="\E[00m"

> "$OUT_FILE"
for T in `seq $Tstart $Tinc $Tend`; do
	echo -e "\E[07mT: $T\t$END" | expand -t 60

	./mcrun -v -t $T $ITERATIONS $SHM_FILE > /dev/null

	#
	echo -ne "$T\t" >> $OUT_FILE
	./mcrun -av -t $T -! $AVG_ITER $SHM_FILE |tr "\n" "\t" >> $OUT_FILE

	echo "" >> $OUT_FILE
done

OUT_FILE="$2"

> "$OUT_FILE"
./mcrun -v -t 0 -! $AVG_ITER $SHM_FILE > /dev/null
for T in `seq $Tstart $Tinc $Tend`; do
	echo -e "\E[07mT: $T\t$END" | expand -t 60

	./mcrun -v -t $T -! $ITERATIONS $SHM_FILE > /dev/null

	#
	echo -ne "$T\t" >> $OUT_FILE
	./mcrun -av -t $T -! $AVG_ITER $SHM_FILE |tr "\n" "\t" >> $OUT_FILE

	echo "" >> $OUT_FILE
done

kill %1
wait

rm "$SHM_FILE"
