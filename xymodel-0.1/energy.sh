#!/bin/bash

## Change these variables as needed ###
Tinc=0.1
Tstart=0
Tend=2

ITERATIONS=10000 # Iterations to reach
                 # equilibrium
INIT_ITER=100000 # spins to initialize first run,
                 # when starting from nonzero temperature
AVG_ITER=10000	 # Measurement
SIZE=8          # Lattice size
#######################################

print_usage()
{
	echo "usage: $0 vortex|helicity|correlation OUT_FILE"
	echo -e "
Write data to OUT_FILE. Existing data will be overwritten!

vortex        write energy, magnetisation and number of vortices.
              Use \"gnuplot plot/energy.gpi\" for plotting.              	

helicity      Calculate helicity modulus (= spin stiffness).
              Plot with \"gnuplot plot/helicity.gpi\"

correlation   Plot correlation with R --no-save < plot/correlation.R


Change variables in $0 as needed." 
}

print_param()
{
	echo "
size: $SIZE  iterations: $ITERATIONS   avg_iter: $AVG_ITER
"
}

cleanup()
{
		kill %1
		wait

		rm "$SHM_FILE"
		rm "$FIFO"
}

cleanup_kill()
{
		trap - INT
		echo "-1" >&6
		cleanup
}

if [[ $# != 2 ]]; then
	print_usage
	exit 1
fi

SHM_FILE=`mktemp` || exit 1
FIFO=".progress_fifo."`basename $SHM_FILE`
trap 'cleanup_kill; exit 1' INT

OUT_FILE="$2"

./xymodel -s $SIZE "$SHM_FILE" &
sleep 1	#wait for xymodel to load in background


> "$OUT_FILE"

if [[ "$Tstart" != 0 ]]; then
	echo "Initializing"
	./mcrun -v -t $Tstart $INIT_ITER $SHM_FILE > /dev/null
fi

print_param

rm -f $FIFO
mkfifo $FIFO
exec 6<> $FIFO
./progress -e $FIFO $Tend $Tstart &

for T in `seq $Tstart $Tinc $Tend`; do
	./mcrun -t $T $ITERATIONS $SHM_FILE > /dev/null
	echo $T >&6

	case "$1" in
		energy|vortex)
			echo -ne "$T\t" >> $OUT_FILE
			./mcrun -a -t $T -p plugin/vortex.so $AVG_ITER $SHM_FILE 2>/dev/null |tr "\n" "\t" >> $OUT_FILE 
			echo "" >> $OUT_FILE
			;;

		helicity)
			echo -ne "$T\t" >> $OUT_FILE
			./mcrun -t $T -p plugin/helicity.so $AVG_ITER $SHM_FILE |tr "\n" "\t" >> $OUT_FILE
			echo "" >> $OUT_FILE
			;;

		correlation)
			./mcrun -t $T -p plugin/correlation.so $AVG_ITER $SHM_FILE |sed "s/^/$T\t/" >> $OUT_FILE
			echo "" >> $OUT_FILE
			;;

		*)
			print_usage
			break
			;;
	esac

	echo -e "scale=8\n$T+$Tinc/($ITERATIONS+$AVG_ITER)*$AVG_ITER" | bc >&6
done

echo $Tend >&6
cleanup
