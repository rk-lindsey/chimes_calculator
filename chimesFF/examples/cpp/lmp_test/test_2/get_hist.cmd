#!/bin/bash
#SBATCH -J ALC-3-md-c0-i0
#SBATCH -N 1
#SBATCH --ntasks-per-node 56
#SBATCH -t 1:00:00
#SBATCH -p pdebug
#SBATCH -A iap
#SBATCH -V 
#SBATCH -o stdoutmsg

mpiicc -O3 -o histogram calc_histograms.cpp chimesFF.cpp

srun -N 1 -n 56 ./histogram 0 0
