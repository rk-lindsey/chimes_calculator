#!/bin/bash

#APIS=$1
#TESTS=$2
PYTH3=python3.7 # $3

FFS[0]="published_params.liqC.2b.cubic.txt"                                  ; CFGS[0]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"
FFS[1]="published_params.liqC.2+3b.cubic.txt"                                ; CFGS[1]="liqC.2.5gcc_6000K.OUTCAR_#000.xyz"
FFS[2]="published_params.liqCO.2+3b.cubic.txt"                               ; CFGS[2]="CO.2.5gcc_6500K.OUTCAR_#000.xyz"
FFS[3]="validated_params.CO2400K.2+3+4b.Tersoff.special.offsets.relabel.txt" ; CFGS[3]="CO.2.5gcc_6500K.OUTCAR_#000.relabel.xyz"
FFS[4]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"         ; CFGS[4]="CO.2.5gcc_6500K.OUTCAR_#000.scramble.xyz"
FFS[5]="published_params.CO2400K.2+3+4b.Tersoff.special.offsets.txt"         ; CFGS[5]="CO.2.5gcc_6500K.OUTCAR_#000.translate.xyz"
FFS[6]="published_params.HN3.2+3+4b.Tersoff.special.offsets.txt"             ; CFGS[6]="HN3.2gcc_3000K.OUTCAR_#000.xyz"
FFS[7]="validated_params.TiO2.2+3b.Tersoff.txt"                              ; CFGS[7]="TiO2.unitcell_arbrot_#000.xyz"

NO_TESTS=${#FFS[@]}
LOC=`pwd`

API[0]="cpp"    ; EXE[0]="test-CPP"; XTRA[0]="2"
API[1]="c"      ; EXE[1]="test-C"  ; XTRA[1]="2"
API[2]="fortran"; EXE[2]="test-F"  ; XTRA[2]="2"
API[3]="python" ; EXE[3]="main.py" ; XTRA[3]="2 1"

for i in {0..3}
do
	cd ../examples/${API[$i]}
	
	if [[ "${API[$i]}" != "python" ]] ; then
	
		make all DEBUG=1
	else
		make all
		cp libwrapper-C.so ../../tests
	fi
	
	cd ../../tests
	
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
		
				../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${XTRA[$i]} > /dev/null
			else
				${PYTH3} ../examples/${API[$i]}/${EXE[$i]} force_fields/${FFS[$j]} configurations/$CFG ${XTRA[$i]} ${LOC}/../api > /dev/null
			fi
			
			# Compare results against expected results (expected_output/${FFS[$j]}.$CFG.dat)
			
			paste debug.dat expected_output/${FFS[$j]}.$CFG.dat > san.dat

			# Print findings
			
			${PYTH3} compare.py san.dat 
			
		done
		
		echo "	Test $idx of $NO_TESTS for API ${API[$i]} complete"

		let idx=idx+1
	done

done

rm -f debug.dat san.dat libwrapper-C.so
