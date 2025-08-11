#!/bin/bash

# Generates expected output folder for examples/tests.
# Should be run ***only on fresh clones***, on a single  node with at least 36 procs
# Run with ./create.sh
# Takes < 30 min to run


################# BASIC EXAMPLES/TESTS

echo "Preparing a fresh installation. I hope you loaded your module files!"

cd ..
./install.sh
cd - > /dev/null 2>&1 

echo "Generating expected_output"

# For test_suite-* files...

for i in test_suite-ev_tally_mb test_suite-CNP_hybrid_overlay_D2 test_suite-NPT test_suite-NVE test_suite-NVT
do 
	echo "Working on $i"
	cd $i
	srun -n 36 ../../exe/lmp_mpi_chimes -i in.lammps > out.lammps
	mkdir expected_output > /dev/null 2>&1
	mv log.lammps out.lammps expected_output
	cd - > /dev/null 2>&1 
done	

# For example-* files...

for i in example-basic_carbon2.0 example-basic_HN3 example-coarsegrained example-Si+D2 example-Turbo_ChIMES
do 
	echo "Working on $i"
	cd $i
	srun -n 36 ../../exe/lmp_mpi_chimes -i in.lammps > out.lammps
	mkdir expected_output > /dev/null 2>&1
	mv log.lammps out.lammps expected_output
	cd - m 
done	



################# FINGERPRINT

echo "..."
echo "Re-compiling to generate with Fingerprint flag..." 
echo ""
cd ..
pwd
./install.sh FINGERPRINT
cd - > /dev/null 2>&1 


echo "Working on example-fingerprint"
cd example-fingerprint
./run_job.sh 36
mkdir expected_output > /dev/null 2>&1
mv log.lammps out.lammps 0-0.{2,3,4}b_clu-s.hist expected_output
cd -



################# TABULATE


echo "..."
echo "Re-compiling to generate with Tabulate flag..." 
echo ""
cd ..
pwd
./install.sh TABULATION
cd - > /dev/null 2>&1 

echo "Working on example-tabulate"
cd example-tabulate
srun -n 36 ../../exe/lmp_mpi_chimes -i in.lammps >  out.lammps
mkdir expected_output > /dev/null 2>&1
mv log.lammps out.lammps expected_output
cd - > /dev/null 2>&1 

