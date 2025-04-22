#!/bin/bash

#SBATCH -J LAMMPS_test_extra_pair                  
#SBATCH -o stdoutmsg_%j                   
#SBATCH -e erroutmsg_%j                     
#SBATCH -p spr                            # Partition(nprocs per node, RAM in GB): icx(80,256), skx(48,192), spr(112,128)
#SBATCH -N 2                               
#SBATCH -n 224                            # Total nprocs
#SBATCH -t 00:30:00                       # Run time (hh:mm:ss)
#SBATCH --mail-user=yanjunlv@umich.edu
#SBATCH --mail-type=all                  
#SBATCH -A TG-MAT250016                   # Project code (1) CNP Explore: TG-MAT250016

module load intel/24.0
module load impi/21.11
module load cmake/3.28.1

lmp="/work2/08034/tg873340/stampede3/chimes_calculator-LLfork/etc/lmp/exe/lmp_mpi_chimes"
srun -N 2 -n 224 $lmp -i in.lammps > out.lammps