#!/bin/bash
 
#SBATCH -N 1
#SBATCH --ntasks-per-node 112
#SBATCH -J HN3_Fingerprint_test
#SBATCH -t 01:00:00
#SBATCH -p pbatch
#SBATCH -A iap
#SBATCH -o stdoutmsg
 
srun -N 1 -n 112 ../../exe/lmp_mpi_chimes -i HN3.2gcc_3000K.OUTCAR_#000.in.lammps > output.txt
sh ../../../../chimesFF/src/FP/post_process.sh
srun -N 1 -n 112 ../../../../chimesFF/src/FP/histogram published_params.HN3.2+3+4b.Tersoff.special.offsets.txt
