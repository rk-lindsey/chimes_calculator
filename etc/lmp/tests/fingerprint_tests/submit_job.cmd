#!/bin/bash
 
#SBATCH -N 1
#SBATCH --ntasks-per-node 112
#SBATCH -J liqC_NVT
#SBATCH -t 01:00:00
#SBATCH -p pbatch
#SBATCH -A iap
#SBATCH -o stdoutmsg
 
srun -N 1 -n 112 ../../../../chimesFF/src/FP/histogram published_params.HN3.2+3+4b.Tersoff.special.offsets.txt
