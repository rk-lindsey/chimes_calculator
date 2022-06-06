#!/bin/bash

################################
# Get relevant paths
################################

chimes_path=`realpath ../../`
dftb_path=`realpath ./build/dftbplus-22.1/`

################################
# Download the dftb+ source code
################################

mkdir build
cd build
wget https://github.com/dftbplus/dftbplus/releases/download/22.1/dftbplus-22.1.tar.xz
tar -xvf dftbplus-22.1.tar.xz
cd dftbplus-22.1

################################
# Create a symlink to the ChIMES code
################################

ln -s $chimes_calc ${dftb_path}/external/chimes/origin/

################################
# Replace the relavant dftbplus files
################################

cp ${chimes_path}/etc/dftbplus/etc/chimes.F90    ${dftb_path}/src/dftbp/extlibs
cp ${chimes_path}/etc/dftbplus/etc/chimesrep.F90 ${dftb_path}/src/dftbp/dftb/repulsive
cp ${chimes_path}/etc/dftbplus/etc/parser.F90   ${dftb_path}/src/dftbp/dftbplus


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
