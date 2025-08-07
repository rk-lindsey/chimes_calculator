#!/bin/bash

# FYI, this is designed to run on SLURM and assumes at least 36 procs are available. Do not change the number of procs. Expects to be run from an interactive session. For example, if running on an LLNL machine, use something like: sxterm 1 56 60 -p pdebug -A iap, then in the window that launches, execute ./run_tests.sh 

# These tests take about 15 min to run. 

# Do not modify anything below this line (including NP) unless you are a developer and **really** know what you are doing!


NP=36


########### For other more basic tests


# For examples-* folders

example[0]=example-basic_HN3
example[1]=example-basic_carbon2.0
example[2]=example-coarsegrained
example[3]=example-Si+D2
example[4]=example-Turbo_ChIMES


echo "Preparing a fresh installation with no special flags"

cd ..
./install.sh
cd - > /dev/null 2>&1 




# Run the LAMMPS examples

for i in {0..4}
do

	rm -rf ${example[$i]}/current_output > /dev/null 2>&1
	mkdir -p ${example[$i]}/current_output
	cd ${example[$i]}/current_output
	cp ../* . > /dev/null 2>&1
	
	params=`awk '/serial/{print $NF}' in.lammps`
	hybrid=`awk '/* * chimesFF/{print "yes"}' in.lammps`

	if [[ -n "$params" && -n "$hybrid" ]]
	then
		cp ../${params} params.txt
		awk '/* * chimesFF/{$0="pair_coeff * * chimesFF params.txt"}{print}' ../in.lammps > in.lammps
	elif [[ -n $params ]] ; 
	then 
		cp ../${params} params.txt
		awk '/pair_coeff/{$0="pair_coeff * * params.txt"}{print}' ../in.lammps > in.lammps
		
	fi
	
	echo "Running test: ${example[$i]}"
	srun -n 36 ../../../exe/lmp_mpi_chimes -i in.lammps > out.lammps
	
	python3 ../../compare_logfiles.py ../expected_output/log.lammps log.lammps > compare.log.lammps

	cd - > /dev/null 2>&1 
done

# For test_suite-* folders

test[0]=test_suite-ev_tally_mb 
test[1]=test_suite-CNP_hybrid_overlay_D2 
test[2]=test_suite-NPT 
test[3]=test_suite-NVE 
test[4]=test_suite-NVT

for i in {0..4}
do

	rm -rf ${test[$i]}/current_output > /dev/null 2>&1
	mkdir -p ${test[$i]}/current_output
	cd ${test[$i]}/current_output
	cp ../* . > /dev/null 2>&1
	
	params=`awk '/serial/{print $NF}' in.lammps`
	hybrid=`awk '/* * chimesFF/{print "yes"}' in.lammps`

	if [[ -n "$params" && -n "$hybrid" ]]
	then
		cp ../${params} params.txt
		awk '/* * chimesFF/{$0="pair_coeff * * chimesFF params.txt"}{print}' ../in.lammps > in.lammps
	elif [[ -n $params ]] ; 
	then 
		cp ../${params} params.txt
		awk '/pair_coeff/{$0="pair_coeff * * params.txt"}{print}' ../in.lammps > in.lammps
		
	fi

	echo "Running test: ${test[$i]}"
	srun -n 36 ../../../exe/lmp_mpi_chimes -i in.lammps > out.lammps
	
	python3 ../../compare_logfiles.py ../expected_output/log.lammps log.lammps > compare.log.lammps

	cd - > /dev/null 2>&1
done






########### For fingerprinting

echo "Preparing a fresh installation with option FINGERPRINT. I hope you loaded your module files!"

cd ..
./install.sh FINGERPRINT
cd - > /dev/null 2>&1 

rm -rf example-fingerprint/current_output > /dev/null 2>&1
mkdir -p example-fingerprint/current_output
cd example-fingerprint/current_output
cp ../* . > /dev/null 2>&1

echo "Running test: example-fingerprint"
./run_job.sh 36 > run_job.log

python3 ../../compare_fingerprints.py ../expected_output/0-0.2b_clu-s.hist 0-0.2b_clu-s.hist > compare.0-0.2b_clu-s.hist.lammps
python3 ../../compare_fingerprints.py ../expected_output/0-0.3b_clu-s.hist 0-0.3b_clu-s.hist > compare.0-0.3b_clu-s.hist.lammps
python3 ../../compare_fingerprints.py ../expected_output/0-0.4b_clu-s.hist 0-0.4b_clu-s.hist > compare.0-0.4b_clu-s.hist.lammps

cd - > /dev/null 2>&1




########### For tabulation

echo "Preparing a fresh installation with option TABULATION."

cd ..
./install.sh TABULATION
cd - > /dev/null 2>&1 

rm -rf example-tabulate/current_output > /dev/null 2>&1
mkdir -p example-tabulate/current_output
cd example-tabulate/current_output
cp ../* . > /dev/null 2>&1

params=`awk '/serial/{print $NF}' in.lammps`
hybrid=`awk '/* * chimesFF/{print "yes"}' in.lammps`

if [[ -n "$params" && -n "$hybrid" ]]
then
	cp ../${params} params.txt
	awk '/* * chimesFF/{$0="pair_coeff * * chimesFF params.txt"}{print}' ../in.lammps > in.lammps
elif [[ -n $params ]] ; 
then 
	cp ../${params} params.txt
	awk '/pair_coeff/{$0="pair_coeff * * params.txt"}{print}' ../in.lammps > in.lammps
	
fi

echo "Running test: example-tabulate"
srun -n 36 ../../../exe/lmp_mpi_chimes -i in.lammps > out.lammps

python3 ../../compare_logfiles.py ../expected_output/log.lammps log.lammps > compare.log.lammps


cd - > /dev/null 2>&1





failures=`wc -l */current*/compare* | awk 'BEGIN{c=0}{c+=$1}END{print c}'`

if [ $failures -eq 0 ] ; then
	echo "All tests passed!"
else
	echo "Some tests failed. See */current_output/compare* files:"
	wc -l */current*/compare*
fi	
