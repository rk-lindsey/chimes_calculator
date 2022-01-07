#!/bin/bash

# ChIMES Calculator test suite
#
# For (relatively) quick and dirty testing, run with:
# ./run_tests.sh SHORT
# Otherwise, run with:
# ./run_tests.sh
# To specify an install prefix for the cmake tests, run with:
# ./run_tests.sh <SHORT or LONG> <install prefix full path>


##################

STYLE=${1-"LONG"} # By default,run "LONG" test, but if user runs with "./run_tests SHORT, runs short tests
PREFX=${2-""}     # By default, don't set any special install prefix
PYTH3=python3.7

# Populate FFS, CFGS, and OPTIONS from test_list.dat


awk '/short/||/minimal/{print}' test_list.dat | tr ';' ' ' > tmp-data.dat
while read FF CFG OPTION TMP
do 
	FFS+=($FF); CFGS+=($CFG); OPTIONS+=($OPTION)
done < tmp-data.dat; rm -f tmp-data.dat

API_LIST="0 1 2 3 4"
NO_TESTS=${#FFS[@]}
LOC=`pwd`

API[0]="cpp"      ; EXE[0]="chimescalc"		        ; XTRA[0]="" #"2"
API[1]="c"        ; EXE[1]="chimescalc-test_serial-C"   ; XTRA[1]="" #"2"
API[2]="fortran"  ; EXE[2]="chimescalc-test_serial-F"   ; XTRA[2]="" #"2"
API[3]="python"   ; EXE[3]="main.py"			; XTRA[3]="" #"2 1"
API[4]="fortran08"; EXE[4]="chimescalc-test_serial-F08" ; XTRA[4]="" #"0"

echo "Running $STYLE tests"
date

for compile in CMAKE MAKEFILE
do
	echo "Testing compilation type: $compile"

	# Do the compilation

	if [[ $compile == "MAKEFILE" ]]; then

		for i in $API_LIST # Cycle through APIs
		do
			cd ../examples/${API[$i]}

			echo "Compiling for API ${API[$i]}"
			echo ""

			if [[ "${API[$i]}" != "python" ]] ; then

				make all DEBUG=1
			else
				make all
				cp libchimescalc-serial_dl.so ../../tests/libchimescalc_dl.so
			fi

			cd ../../tests

		done

	elif [[ $compile == "CMAKE" ]] ; then

		cd ../../
		./install.sh 1 $PREFX # Set the debug flag true
		cp build/libchimescalc_dl.so  serial_interface/tests
		cd -

	else
		echo "Error: Unknown compilation method $compile"
		echo "Acceptable values are MAKEFILE and CMAKE"
		echo "Check logic in run_test.sh"
		exit 0
	fi

	# Run the tasks

	for i in $API_LIST # Cycle through APIs
	do

		idx=1

		for ((j=0;j<NO_TESTS;j++))
		do
			echo "Working on Test $idx of $NO_TESTS for API ${API[$i]}"

			for ((k=0; k<10; k++))
			do
				CFG=${CFGS[$j]}
				CFG_PREFIX="${CFG%%_#*}_#"
				CFG_SUFFIX=${CFG##*000}

				CFG=${CFG_PREFIX}00${k}${CFG_SUFFIX}

				if [ ! -f configurations/$CFG ] ; then
					continue
				fi

				echo "		...Running $CFG"

				# Run the test

				if [[ "${API[$i]}" != "python" ]] ; then

					if [[ $compile == "CMAKE" ]] ; then
						../../build/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} > /dev/null
					else
						../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]}  > /dev/null
					fi

				else
					${PYTH3} ../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${OPTIONS[$j]} ${LOC}/../api 1 > /dev/null
				fi

				# Compare results against expected results (expected_output/${FFS[$j]}.$CFG.dat)

				paste debug.dat expected_output/${FFS[$j]}.$CFG.dat > san.dat

				# Print findings

				${PYTH3} compare.py san.dat

				if [[ $STYLE == "SHORT" ]] ; then
					break
				fi

			done

			echo "	Test $idx of $NO_TESTS for API ${API[$i]} complete"

			let idx=idx+1

		done

	done

done


read -p "Press enter to continue on to test cleanup." tmpvar


# Clean up

for i in $API_LIST # Cycle through APIs
do
	cd ../examples/${API[$i]}
	make clean-all
	rm -f *.so *.a
	cd ../../tests
done

rm -f debug.dat san.dat *.so output_lib.xyzf

cd ../../
./uninstall.sh $PREFX
