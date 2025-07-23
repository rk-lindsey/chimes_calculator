#!/bin/bash

# Builds all relevant chimes_calculator executables/library files 
#
# Usage:
#   ./install.sh [TABULATION]
#   If you pass “TABULATION” as the first argument, -DTABULATION will
#   be appended to the C++ compile flags for pair_chimes.cpp and chimes_ff.cpp.

echo ""
echo "Note: This install script assumes: "
echo "1. Availability of Intel C++ compilers with c++11 support"
echo "2. Availability of Intel MPI compilers"
echo "...Intel oneapi compilers are now freely available"
echo ""

# ********** TABULATION FLAG HANDLING **********
TAB_FLAG=""
if [[ "$1" == "TABULATION" ]]; then
    TAB_FLAG="-DTABULATION"
    echo "Enabling TABULATION compilation flag for ChIMES files"
    shift
fi

# Clean up any previous installation
./uninstall.sh

# Clone the tested LAMMPS stable branch
mkdir -p build/lammps_stable_29Oct2020
git clone --depth 1 --branch stable_29Oct2020 \
    https://github.com/lammps/lammps.git build/lammps_stable_29Oct2020

# Copy ChIMES files into LAMMPS tree
cp ../../chimesFF/src/chimesFF.{h,cpp}      build/lammps_stable_29Oct2020/src/MANYBODY/
cp src/pair_chimes.{h,cpp}                  build/lammps_stable_29Oct2020/src/MANYBODY/
cp etc/pair.{h,cpp}                         build/lammps_stable_29Oct2020/src
cp etc/Makefile.mpi_chimes                  build/lammps_stable_29Oct2020/src/MAKE

# ********** APPLY TABULATION FLAG TO MAKEFILE **********
if [ -n "$TAB_FLAG" ]; then
    sed -e "/^CCFLAGS[[:space:]]*=/ s|$| $TAB_FLAG|" etc/Makefile.mpi_chimes > build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
else
    # Use standard Makefile
    cp etc/Makefile.mpi_chimes build/lammps_stable_29Oct2020/src/MAKE/
fi


# Load module files and configure compilers
if [ -z "$hosttype" ]; then
    echo ""
    echo "WARNING: No hosttype specified"
    echo "Load modules/configure compilers manually before running this script!"
    echo ""
elif [[ "$hosttype" == "LLNL-LC" ]]; then
    source modfiles/LLNL-LC.mod
elif [[ "$hosttype" == "UM-ARC" ]]; then
    source modfiles/UM-ARC.mod
elif [[ "$hosttype" == "JHU-ARCH" ]]; then
    source modfiles/JHU-ARCH.mod
    ICC=$(which icc)
    MPI=$(which mpicxx)
elif [[ "$hosttype" == "UT-TACC" ]]; then
    source modfiles/UT-TACC.mod
        if [ -n "$TAB_FLAG" ]; then
        # Append fingerprint flag to CCFLAGS
        sed -e "/^CCFLAGS[[:space:]]*=/ s|$| $TAB_FLAG|"  etc/Makefile.mpi_chimes.UT-TACC > build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
    else
        cp etc/Makefile.mpi_chimes.UT-TACC build/lammps_stable_29Oct2020/src/MAKE/Makefile.mpi_chimes
    fi
else
    echo ""
    echo "ERROR: Unknown hosttype ($hosttype)"
    echo "Valid options are:"
    for i in modfiles/*.mod; do echo "   ${i%.mod}"; done
    echo ""
    echo "Set hosttype and rerun, or load modules manually."
    exit 1
fi

echo "Detected hosttype: $hosttype"
module list 2>/dev/null || true

# Compile
cd build/lammps_stable_29Oct2020/src
make yes-manybody
make yes-user-misc
make -j 4 mpi_chimes
cd -

# Finish
mkdir -p exe
mv build/lammps_stable_29Oct2020/src/lmp_mpi_chimes exe/

loc=$(pwd)
echo ""
echo "Compilation complete."
echo "Generated LAMMPS executable with ChIMES support${TAB_FLAG:+ (TABULATION)}:"
echo "  ${loc}/exe/lmp_mpi_chimes"
echo "See ${loc}/tests for usage examples"
echo ""
