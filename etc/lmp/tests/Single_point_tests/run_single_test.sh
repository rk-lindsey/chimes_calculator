#!/bin/bash

SCRIPT_DIR=$(dirname "$0")
SCRIPT_NAME=$(basename "$0")

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 [BINARY] PARAMETERFILE GEOMETRYFILE LAMMPSINPUT" >&2
  echo "Assumes BINARY is at location ../../exe/lmp_mpi_chimes if not given"
  echo "At least 3 arguments required." >&2
  exit 1
fi

if [ "$#" -eq 3 ]; then
  echo "Assuming BINARY is at location ../../exe/lmp_mpi_chimes"
  binary=$(realpath ../../exe/lmp_mpi_chimes)
  parameterfile=$1
  geometryfile=$2
  lammpsfile=$3
  workdir=$(realpath .)
elif [ "$#" -ge 4 ]; then
  binary=$(realpath "$1")
  parameterfile=$2
  geometryfile=$3
  lammpsfile=$4
  workdir=$(realpath "${5:-.}")
else
  echo "Usage: $0 [BINARY] PARAMETERFILE GEOMETRYFILE LAMMPSINPUT" >&2
  echo "Assumes BINARY is at location ../../exe/lmp_mpi_chimes if not given"
  echo "At least 3 arguments required." >&2
  exit 1
fi

filename=$(basename "$lammpsfile")
base_parameter=$(basename "$parameterfile")
base_geometry=$(basename "$geometryfile")
echo "[INFO] Running $SCRIPT_NAME with input file(s): $base_parameter, $base_geometry and $filename"

# Check if executable exists
if [ ! -x "$binary" ]; then
  echo "[ERROR] Executable not found or not executable at $binary" >&2
  exit 2
fi

echo "[INFO] Found LAMMPS executable at $binary"

if [ ! -d ${workdir} ]; then
  mkdir -p ${workdir}
fi

# Run LAMMPS
cd "${workdir}" || { echo "Failed to cd into $workdir"; exit 1; }
# ${binary} -i ${SCRIPT_DIR}/${parameterfile} >& output
${binary} -i ${SCRIPT_DIR}/${lammpsfile} | tee output

# Generate debug.dat file
echo "[INFO] Parsing LAMMPS output..."
python3 "${SCRIPT_DIR}/gen_compare.py" \
  "${workdir}/log.lammps" \
  "${workdir}/traj.lammpstrj" \
  "${workdir}/${filename}.debug.dat"

outputdir="generated_output"
if [ ! -d ${outputdir} ]; then
  mkdir -p ${outputdir}
fi

mv "${workdir}/${filename}.debug.dat" "${outputdir}/" 

# Comparing expected_output with generated_output
echo "[INFO] Comparing LAMMPS output with expected output for ${filename}..."
geo_xyz="${base_geometry/.in.data/.xyz}"  # replace .in.data → .xyz
python3 "${SCRIPT_DIR}/compare.py" \
  "${outputdir}/${filename}.debug.dat" \
  "expected_output/${base_parameter}.${geo_xyz}.dat"

echo "[INFO] Completed test. Removing log.lammps, traj.lammpstrj, output..."
rm log.lammps traj.lammpstrj output

cmp_status=$?

if [ "$cmp_status" -ne 0 ]; then
    echo "[FAIL] Comparison failed"
    exit 1  # This ensures run_all_jobs.sh knows this test failed
else
    echo "[PASS] Comparison succeeded"
    exit 0
fi