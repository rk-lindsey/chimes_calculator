#!/bin/bash

# Builds all relevant chimes_calculator executables/library files 
#
# If working on a machine with a corresponding .mod file in the modfiles folder
# (e.g., modfiles/LLNL-LC.mod), execute with, e.g.:
#
#   export hosttype=LLNL-LC; ./install.sh
# 
# Otherwise, load necessary modules manually and execute with 
# 
# ./install.sh 
# 
# Note that additional arguments can be specified, i.e.:
#
#   ./install.sh FINGERPRINT or TABULATION, currently only 1 or the other are possible.

echo ""
echo "Note: This install script assumes: "
echo "1. Availability of Intel C++ compilers with c++11 support"
echo "2. Availability of Intel MPI compilers"
echo "...Intel oneapi compilers are now freely available"
echo ""

# ********** FLAG HANDLING **********
TAB_FLAG=""
FINGERPRINT_FLAG=""

if [[ "$1" == "TABULATION" ]]; then
    TAB_FLAG="-DTABULATION"
    echo "Enabling TABULATION compilation flag for ChIMES files"
elif [[ "$1" == "FINGERPRINT" ]]; then
    FINGERPRINT_FLAG="-DFINGERPRINT"
    echo "Enabling FINGERPRINT compilation flag for ChIMES files"
elif [[ -n "$1" ]]; then
    echo "ERROR: Invalid option '$1'. Use only one of: TABULATION or FINGERPRINT"
    exit 1
fi

# Cleanup any previous installation

./uninstall.sh


# Grab the specific stable branch of LAMMPS compaitbility has been tested for

mkdir -p build/lammps_stable_29Oct2020


git clone --depth 1 --branch stable_29Oct2020 https://github.com/lammps/lammps.git build/lammps_stable_29Oct2020


# Copy ChIMES files to correct locations

cp ../../chimesFF/src/chimesFF.{h,cpp}	build/lammps_stable_29Oct2020/src/MANYBODY/
cp src/pair_chimes.{h,cpp} 		build/lammps_stable_29Oct2020/src/MANYBODY/
cp etc/pair.{h,cpp} 			build/lammps_stable_29Oct2020/src
MAKEFILE_SRC="etc/Makefile.mpi_chimes"
if [[ "$hosttype" == "UT-TACC" ]]; then
    MAKEFILE_SRC="etc/Makefile.mpi_chimes.UT-TACC"
fi
# ********** MODIFIED MAKEFILE HANDLING **********
if [[ -n "$TAB_FLAG" ]]; then
    sed -e "/^CCFLAGS[[:space:]]*=/ s|$| $TAB_FLAG|" "$MAKEFILE_SRC" > build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
elif [[ -n "$FINGERPRINT_FLAG" ]]; then
    sed -e "/^CCFLAGS[[:space:]]*=/ s|$| $FINGERPRINT_FLAG|" "$MAKEFILE_SRC" > build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
else
    cp "$MAKEFILE_SRC" build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
fi

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
elif [[ "$hosttype" == "UT-TACC" ]] ; then
    source modfiles/UT-TACC.mod
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
make yes-user-misc
make -j 4 mpi_chimes
cd -

# Compile fingerprint executable if flag true
if [ -n "$FINGERPRINT_FLAG" ]; then
    echo ""
    echo "Compiling histogram executable for ChIMES fingerprints"
    mpiicc -O3 -o ../../chimesFF/src/FP/histogram ../../chimesFF/src/FP/multi_calc_histogram.cpp ../../chimesFF/src/FP/chimesFF.cpp
fi

# Finish

mkdir exe
mv build/lammps_stable_29Oct2020/src/lmp_mpi_chimes exe

loc=`pwd`
echo ""
echo "Compilation complete. "
if [ -n "$TAB_FLAG" ]; then
    echo "Generated LAMMPS executable with ChIMES TABULATION support:"
elif [ -n "$FINGERPRINT_FLAG" ]; then
    echo "Generated LAMMPS executable with ChIMES FINGERPRINT support:"
else
    echo "Generated LAMMPS executable with basic ChIMES support (no extra flags):"
fi
echo "${loc}/exe/lmp_mpi_chimes"
echo "See ${loc}/tests for usage examples"
echo ""