#!/bin/bash
set -e

SCRIPT_DIR=$(dirname "$0")
SCRIPT_NAME=$(basename "$0")

if [ $# -ne 1 ]; then
  echo "Usage: $SCRIPT_NAME INPUT_FILE" >&2
  exit 1
fi

inputfile=$1
exe=../../exe/lmp_mpi_chimes

echo "[INFO] Running $SCRIPT_NAME with input file: $inputfile"

# Check if executable exists
if [ ! -x "$exe" ]; then
  echo "[ERROR] Executable not found or not executable at $exe" >&2
  exit 2
fi

echo "[INFO] Found LAMMPS executable at $exe"

# Run LAMMPS
outputfile="${inputfile}_lammps_output.log"
echo "[INFO] Running LAMMPS simulation..."
if "$exe" -i "lammps_input_files/${inputfile}" > "${outputfile}"; then
  echo "[SUCCESS] LAMMPS simulation completed"
else
  echo "[ERROR] LAMMPS simulation failed!" >&2
  exit 3
fi

# Move output files
mkdir -p output_files
mv "${outputfile}" "output_files/"
mv traj.lammpstrj "output_files/${inputfile}.traj.lammpstrj" 2>/dev/null || echo "[WARNING] traj.lammpstrj not found"
rm -f log.cite log.lammps rank-0.badness.log

# Define paths
log_file="output_files/${outputfile}"
traj_file="output_files/${inputfile}.traj.lammpstrj"
parsed_output_file="output_files/${inputfile}.parsed_output.txt"

# Check if log and traj files exist before parsing
if [[ -f "$log_file" && -f "$traj_file" ]]; then
  echo "[INFO] Parsing LAMMPS output..."
  python3 gen_compare.py "${log_file}" "${traj_file}" "${parsed_output_file}"
  echo "[SUCCESS] Parsed output written to $parsed_output_file"
else
  echo "[ERROR] Missing output files: log or trajectory file not found." >&2
  exit 4
fi
