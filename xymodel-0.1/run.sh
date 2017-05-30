#!/bin/bash

if [[ $# != 5 ]]; then
	echo "usage: $0 size temp_start temp_step temp_end OUT_DIR"
	echo -e "\nWrite calculated lattices with given temperatures in OUT_DIR.
Overwrites files without asking."
	echo -e "\ne.g.: ./run.sh 100x100 0 0.01 2 out"
	exit 1
fi

# change this as needed
ITERATIONS=5000
START_ITERATIONS=10000


SHM_FILE=`tempfile && exit 1`

DIR=$5
SIZE=$1

mkdir -p "$5"

echo "size	$SIZE
start_iterations	$START_ITERATIONS
interations	$ITERATIONS" > $DIR/INFO

./xymodel -s $SIZE "$SHM_FILE" &
sleep 1 #wait for xymodel to start (just to be sure)
./randomlattice "$SHM_FILE"

echo -n "Start run:	"
./mcrun -v $START_ITERATIONS "$SHM_FILE"

for T in `seq $2 $3 $4`; do
	echo -n "T: $T	"
	./mcrun -v -t $T $ITERATIONS "$SHM_FILE"
	./writelattice "$SHM_FILE" > "${DIR}/${T}"
done

kill %%		# quit xymodel
wait
rm "$SHM_FILE"
