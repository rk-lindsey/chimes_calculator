#!/bin/bash

# Builds all relevant chimes_calculator executables/library files 
# Run with:
# ./install

DEBUG=${1-0} # False (0) by default.

./uninstall.sh
mkdir build
cd build

if [ $DEBUG -eq 1 ] ;then
	cmake -DDEBUG=1 ..
else
	cmake -DDEBUG=0 ..
fi

make
make install
cd ..
