#!/bin/bash


echo ""
echo "Note: This install script assumes: "
echo "1. Availibility of cray-mpich/8.1.28 on perlmutter"
echo "2. Availability of cpu/1.0 on perlmutter"
echo "...Intel oneapi compilers are now freely available"
echo ""

# Cleanup any previous installation

./uninstall.sh


# Grab the specific stable branch of LAMMPS compaitbility has been tested for

mkdir -p build/lammps_stable_27Jun2024


git clone --depth 1 --branch patch_27Jun2024 https://github.com/lammps/lammps.git build/lammps_stable_27Jun2024

echo "LAMMPS stable branch 27Jun2024 cloned"

# Copy ChIMES files to correct locations

# TODO - merge the chimesFF and pair_chimes files for LAMMPS PR
# TODO - pair.{h,cpp} need to be updated in LAMMPS PR
cp ../../chimesFF/src/chimesFF.{h,cpp}	build/lammps_stable_27Jun2024/src/MANYBODY/
cp src/pair_chimes.{h,cpp} 		build/lammps_stable_27Jun2024/src/MANYBODY/
cp etc/pair.{h,cpp} 			build/lammps_stable_27Jun2024/src

echo "ChIMES files copied to LAMMPS source"

# Load module files and configure compilers
# Note: If using intel compilers from after Jan. 2021 (e.g.,) have access to oneapi
#       this means instead of loading inidividual modules for mpi, mkl, etc, can just execute
#       load the intel (e.g., icc) module and run the the setvars.sh command, e.g. located at 
#       /sw/pkgs/arc/intel/2022.1.2/setvars.sh --also avail for free


# Compile
# For GPU build on perlmutter add: -D PKG_KOKKOS=yes -D Kokkos_ARCH_AMPERE80=ON -D Kokkos_ENABLE_CUDA=yes

cd build/lammps_stable_27Jun2024/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install_pm -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_Fortran_COMPILER=ftn -D CMAKE_C_COMPILER=cc -D CMAKE_CXX_COMPILER=CC \
            -D MPI_C_COMPILER=cc -D MPI_CXX_COMPILER=CC -D LAMMPS_EXCEPTIONS=ON \
            -D BUILD_SHARED_LIBS=ON -D PKG_MANYBODY=ON -D PKG_MOLECULE=ON -D PKG_KSPACE=ON -D PKG_REPLICA=ON -D PKG_ASPHERE=ON \
            -D PKG_RIGID=ON -D PKG_MPIIO=ON \
            -D CMAKE_POSITION_INDEPENDENT_CODE=ON -D CMAKE_EXE_FLAGS="-dynamic" ../cmake

echo "Compiling LAMMPS with ChIMES support"

make -j16

# Finish

mkdir ../../../exe
mv ./lmp ../../../exe/

cd ../../../

loc=`pwd`
echo ""
echo "Compilation complete. "
echo "Generated the following LAMMPS executable with ChIMES support:"
echo "${loc}/exe/lmp"
echo "See ${loc}/tests for usage examples"
echo ""

