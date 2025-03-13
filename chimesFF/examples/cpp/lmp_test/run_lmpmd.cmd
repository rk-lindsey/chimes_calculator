#!/bin/bash
#SBATCH -J ALC-3-md-c0-i0
#SBATCH -N 1
#SBATCH --ntasks-per-node 56
#SBATCH -t 1:00:00
#SBATCH -p pdebug
#SBATCH -A iap
#SBATCH -V 
#SBATCH -o stdoutmsg
module load cmake intel impi
srun -N 1 -n 56 /g/g14/sundberg3/codes/FP_Calc/etc/lmp/exe/lmp_mpi_chimes -i case-0.indep-0.in.lammps  > out.lammps
for *.*.2b_clusters.txt
done
cat 0.*2b_clusters.txt > 