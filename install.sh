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
