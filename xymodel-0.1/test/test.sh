#!/bin/bash
#
# call: use "make test"
#
# It is tested if an output of a program
# is equal to some given file.

LAT="test/lattice"
RET=0
TMP=`tempfile`
END="\E[00m"

diff_print()
{

	diff -qE "$1" "$2" >/dev/null
	if [ "$?" -eq 0 ]; then
		echo -e "\E[1;32mok$END"
	else
		echo -e "\E[1;31mfailed$END"
		RET=1
	fi
}

cmp()
{
	./writelattice $LAT > "$TMP"
	diff_print "$TMP" "$1"
}

print_test()
{
	echo -ne "> $1\t" | expand -t 40
}

echo -en "\E[07m"
echo -e "Test\tResult " | expand -t 40
echo -en $END


print_test "direct access"
./writelattice "test/l.init" > $LAT 2>/dev/null
diff_print $LAT "test/l.init"


./xymodel -ws 16 $LAT >/dev/null &

print_test "create"
sleep 1
cmp "test/l.init"

print_test "initialize random"
./init -r $LAT
cmp "test/l.rnd"


print_test "mcrun T=0"
./mcrun -v 1000 $LAT 2> /dev/null
cmp "test/l.0"


print_test "mcrun -a T=0"
./mcrun -v -a 100 $LAT  2> /dev/null > "$TMP"
diff_print "$TMP" "test/avg.0"


print_test "find vortex"
./vortex $LAT > "$TMP" 
diff_print "$TMP" "test/vort.0"


print_test "mcrun T=2"
./init -r $LAT
./mcrun -v -t 2 1000 $LAT 2> /dev/null
cmp "test/l.2"


print_test "mcrun -a T=2"
./mcrun -v -a 100 $LAT  2> /dev/null > "$TMP" 
diff_print "$TMP" "test/avg.2"


print_test "write back"
./init -r $LAT
kill %1
wait
diff_print $LAT "test/l.rnd"

echo ""
echo -ne "speed (T=1):\t"
./mcrun 1000 $LAT > /dev/null 2>&1	# trick cpufreq
./mcrun -t 1 -v 1000 $LAT 2>&1 | grep "spins/sec" | awk '{print $5, $6}'

rm "$TMP"

exit $RET
