#!/bin/bash

# Builds all relevant chimes_calculator executables/library files 
# Run with:
# ./install.sh 
# or
# ./install.sh <debug option (0 or 1)> <install prefix (full path)>

BUILD=`pwd`/build
DEBUG=${1-0}  # False (0) by default.
PREFX=${2-$BUILD} # Empty by default

# Clean up previous installation, 

./uninstall.sh $PREFX


# Load modules

# Determine computing environment and attempt to load module files automatically

lochost=`hostname`
hosttype=""

if [[ $lochost == *"arc-ts.umich.edu"* ]]; then
    hosttype=UM-ARC
elif [[ $lochost == *"quartz"* ]]; then
    hosttype=LLNL-LC
elif [[ $lochost == *"login"* ]]; then
    hosttype=JHU-ARCH
else
    echo "WARNING: Host type ($hosttype) unknown"
    echo "Be sure to load modules/conifugre compilers by hand."
fi

echo "Found host type: $hosttype"

if [[ "$hosttype" == "LLNL-LC" ]] ; then
    source modfiles/LLNL-LC.mod
elif [[ "$hosttype" == "UM-ARC" ]] ; then
    source modfiles/UM-ARC.mod
elif [[ "$hosttype" == "JHU-ARCH" ]] ; then
    module unload gcc/9.3.0 openmpi/3.1.6 git/2.28.0
    source modfiles/JHU-ARCH.mod
    ICC=`which icc`
    MPI=`which mpicxx`
fi

module list


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
