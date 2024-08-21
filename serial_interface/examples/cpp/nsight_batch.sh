#!/bin/bash
#SBATCH -N 1 --account ntrain4 --error=job.err --output=job.out
#SBATCH -C gpu --gpus-per-task=1
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH --job-name=profile_chimes 

#OpenMP settings:
# export OMP_NUM_THREADS=1

module load cudatoolkit/12.2

#run the application:
nsys profile --stats=true -t nvtx,openacc -o openacc_code-nvtx  --force-overwrite=true ./chimescalc ../../../serial_interface/tests/force_fields/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt \
../../../serial_interface/tests/configurations/hackathon/HN3.2gcc_3000K.512.xyz > profile.log


#dcgmi profile --pause

#./wrapper.sh ./chimescalc ~/chimesFF_test/chimes_calculator-main/serial_interface/tests/force_fields/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt \
#~/chimesFF_test/chimes_calculator-main/serial_interface/tests/configurations/HN3.2gcc_3000K.OUTCAR_#001.xyz

#ncu --target-processes all -o nsight-compute ./chimescalc ~/chimesFF_test/chimes_calculator-main/serial_interface/tests/force_fields/published_params.HN3.2+3+4b.Tersoff.special.offsets.txt \
#~/chimesFF_test/chimes_calculator-main/serial_interface/tests/configurations/HN3.2gcc_3000K.OUTCAR_#001.xyz

#dcgmi profile --resume
#rm openMP_code.nsys-rep



#  *.qdstrm is created by nsys profile which cannot be read directyl
#  Hence need to use QdstrmImporter to convert *.qdstrm to a readable format of *.nsys-rep 
/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/profilers/Nsight_Systems/host-linux-x64/QdstrmImporter -i openacc_code-nvtx.qdstrm