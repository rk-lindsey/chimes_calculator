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

FFS[0 ]="published_params.liqC.2b.cubic.txt"                                  ; CFGS[0 ]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"            ; OPTIONS[0 ]="0"
FFS[1 ]="published_params.liqC.2+3b.cubic.txt"                                ; CFGS[1 ]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"            ; OPTIONS[1 ]="0"
FFS[2 ]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[2 ]="CO.2.5gcc_6500K.OUTCAR_#000.xyz"              ; OPTIONS[2 ]="1"
FFS[3 ]="validated_params.liqCO.2+3b.cubic.txt" ; CFGS[3 ]="CO.2.5gcc_6500K.OUTCAR_#000.relabel.xyz"      ; OPTIONS[3 ]="1"
FFS[4 ]="published_params.liqCO.2+3b.cubic.txt"         ; CFGS[4 ]="CO.2.5gcc_6500K.OUTCAR_#000.scramble.xyz"     ; OPTIONS[4 ]="1"
FFS[5 ]="published_params.liqCO.2+3b.cubic.txt"         ; CFGS[5 ]="CO.2.5gcc_6500K.OUTCAR_#000.translate.xyz"    ; OPTIONS[5 ]="1"
FFS[6 ]="published_params.HN3.2+3+4b.Tersoff.special.offsets.txt"             ; CFGS[6 ]="HN3.2gcc_3000K.OUTCAR_#000.xyz"               ; OPTIONS[6 ]="0"
FFS[7 ]="validated_params.TiO2.2+3b.Tersoff.txt"                              ; CFGS[7 ]="TiO2.unitcell_arbrot_#000.xyz"                ; OPTIONS[7 ]="0"
FFS[8 ]="test_params.CHON.txt" ; CFGS[8 ]="CHON.testfile_#000.xyz" ; OPTIONS[8 ]="0"
FFS[9 ]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[9 ]="diam.64_#000.xyz"                             ; OPTIONS[9 ]="1"
FFS[10]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[10]="diam.16_#000.xyz"                             ; OPTIONS[10]="1"
FFS[11]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[11]="diam.8_#000.xyz"                              ; OPTIONS[11]="1"
FFS[12]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[12]="diam.2_#000.xyz"                              ; OPTIONS[12]="1"
FFS[13]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt" ; CFGS[13]="CO.9GPa_2400K.OUTCAR_#000.xyz" ; OPTIONS[13]="1"

API_LIST="0 1 2 3 4"
NO_TESTS=${#FFS[@]}
LOC=`pwd`

API[0]="cpp"    ; EXE[0]="CPP-interface"                    ; XTRA[0]="" #"2"
API[1]="c"      ; EXE[1]="C_wrapper-serial_interface"       ; XTRA[1]="" #"2"
API[2]="fortran"; EXE[2]="fortran_wrapper-serial_interface" ; XTRA[2]="" #"2"
API[3]="python" ; EXE[3]="main.py"                          ; XTRA[3]="" #"2 1"
API[4]="fortran08"; EXE[4]="fortran08_wrapper-serial_interface" ; XTRA[4]="" #"0"

echo "Running $STYLE tests"
date

#for compile in CMAKE MAKEFILE
for compile in MAKEFILE
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
				cp lib-C_wrapper-serial_interface.so ../../tests
			fi
			
			cd ../../tests
			
		done		
	
	elif [[ $compile == "CMAKE" ]] ; then
		
		cd ../../
		./install.sh 1 $PREFX # Set the debug flag true
		cp build/lib-C_wrapper-serial_interface.so  serial_interface/tests		
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


# Clean up

for i in $API_LIST # Cycle through APIs
do
	cd ../examples/${API[$i]}
	make clean
	rm -f *.so *.a
	cd ../../tests
done

rm -f debug.dat san.dat *.so

cd ../../
./uninstall.sh $PREFX


