#!/bin/bash

echo ""
echo "Note: This install script assumes: "
echo "1. Availibility of Intel C++ compilers with c++11 support"
echo "2. Availability of Intel MPI compilers"
echo "...Intel oneapi compilers are now freely available"
echo ""

# Cleanup any previous installation

./uninstall.sh


# Grab the specific stable branch of LAMMPS compaitbility has been tested for

mkdir -p build/lammps_stable_29Oct2020


git clone --depth 1 --branch stable_29Oct2020 https://github.com/lammps/lammps.git build/lammps_stable_29Oct2020


# Copy ChIMES files to correct locations

cp ../../chimesFF/src/chimesFF.{h,cpp}	build/lammps_stable_29Oct2020/src/MANYBODY/
cp src/pair_chimes.{h,cpp} 		build/lammps_stable_29Oct2020/src/MANYBODY/
cp etc/pair.{h,cpp} 			build/lammps_stable_29Oct2020/src
cp etc/Makefile.mpi_chimes 		build/lammps_stable_29Oct2020/src/MAKE


# Load module files and configure compilers
# Note: If using intel compilers from after Jan. 2021 (e.g.,) have access to oneapi
#       this means instead of loading inidividual modules for mpi, mkl, etc, can just execute
#       load the intel (e.g., icc) module and run the the setvars.sh command, e.g. located at 
#       /sw/pkgs/arc/intel/2022.1.2/setvars.sh --also avail for free

if [ -z "$hosttype" ] ; then
    echo ""
    echo "WARNING: No hosttype specified"
    echo "Be sure to load modules/configure compilers by hand before running this script!"
    echo ""
elif [[ "$hosttype" == "LLNL-LC" ]] ; then
    source modfiles/LLNL-LC.mod
elif [[ "$hosttype" == "UM-ARC" ]] ; then
    source modfiles/UM-ARC.mod
elif [[ "$hosttype" == "JHU-ARCH" ]] ; then
    source modfiles/JHU-ARCH.mod
    ICC=`which icc`
    MPI=`which mpicxx`    
else
    echo ""
    echo "ERROR: Unknown hosttype ($hosttype) specified"
    echo ""

    echo "Valid options are:"
    for i in `ls modfiles`; do echo "   ${i%.mod}"; done
    echo ""
    echo "Please run again with: export hosttype=<host type>; ./install.sh"
    echo "Or manually load modules and run with: ./install.sh"
    exit 0
fi


echo "Detected hosttype: $hosttype"
if [ ! -z "$hasmod" ] ; then
    module list
fi


# Compile

cd build/lammps_stable_29Oct2020/src
make yes-manybody
make -j 4 mpi_chimes
cd -


# Finish

mkdir exe
mv build/lammps_stable_29Oct2020/src/lmp_mpi_chimes exe

loc=`pwd`
echo ""
echo "Compilation complete. "
echo "Generated the following LAMMPS executable with ChIMES support:"
echo "${loc}/exe/lmp_mpi_chimes"
echo "See ${loc}/tests for usage examples"
echo ""

