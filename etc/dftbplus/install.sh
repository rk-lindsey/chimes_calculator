#!/bin/bash


################################
# Download the dftb+ source code
################################

mkdir build
cd build
wget https://github.com/dftbplus/dftbplus/releases/download/22.1/dftbplus-22.1.tar.xz
tar -xvf dftbplus-22.1.tar.xz
cd dftbplus-22.1

################################
# Get relevant paths
################################

chimes_path=`realpath ../../../../`
dftb_path=`realpath .`

################################
# Create a symlink to the ChIMES code
################################

rm -rf ${dftb_path}/external/chimes/origin/*
ln -s ${chimes_path}/* ${dftb_path}/external/chimes/origin/

################################
# Replace the relavant dftbplus files
################################

cp ${chimes_path}/etc/dftbplus/etc/chimes.F90    ${dftb_path}/src/dftbp/extlibs/chimes.F90
cp ${chimes_path}/etc/dftbplus/etc/chimesrep.F90 ${dftb_path}/src/dftbp/dftb/repulsive/chimesrep.F90
cp ${chimes_path}/etc/dftbplus/etc/parser.F90    ${dftb_path}/src/dftbp/dftbplus/parser.F90 


################################
# Load necessary modules 
################################

module load cmake/3.16.8 
module load intel/19.0.4
module load mkl/2020.0
module load python/3.6.0
module load mvapich2/2.3


################################
# Configure the build
################################

rm -f _build/CMakeCache.txt

FC=ifort CC=icc cmake \
-DCMAKE_INSTALL_PREFIX=_exe/dftb+ \
-DWITH_MPI=True \
-DWITH_CHIMES=True \
-B _build .

################################
# Execute the build
################################

cmake --build _build -- VERBOSE=1

# Tell user where executable lives:

exeloc=${dftb_path}/_build/app/dftb+/dftb+

echo "Executable located in: $exeloc"
