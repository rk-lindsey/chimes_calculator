#!/bin/bash
 
NP=$1


param_file=./params.txt


srun -n $NP ../../../exe/lmp_mpi_chimes -i in.lammps > out.lammps


../../../../../chimesFF/src/FP/post_process.sh


srun -n $NP ../../../../../chimesFF/src/FP/histogram $param_file
