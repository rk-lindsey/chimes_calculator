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
#   ./install.sh
# 
# Note that additional arguments can be specified, i.e.:
#
#   ./install.sh <debug option (0 or 1)> <install prefix (full path)>

BUILD=`pwd`/build
DEBUG=${1-0}  # False (0) by default.
PREFX=${2-$BUILD} # Empty by default

# Clean up previous installation, 

./uninstall.sh $PREFX

# Load modules based on user-specified host-type

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


# Move into build directory 

mkdir build
cd build

# Generate cmake flags

my_flags=""

if [ ! -z $PREFX ] ; then
	my_flags="-DCMAKE_INSTALL_PREFIX=${PREFX}"
fi

if [ $DEBUG -eq 1 ] ;then

	my_flags="${my_flags} -DDEBUG=1 -DCMAKE_BUILD_TYPE=Release" 
else
	my_flags="${my_flags} -DDEBUG=0 -DCMAKE_BUILD_TYPE=Release" 
fi

# Setup, make and install

cmake $my_flags ..
make
if [ ! -z $PREFX ] ; then
	make install
fi

cd ..
