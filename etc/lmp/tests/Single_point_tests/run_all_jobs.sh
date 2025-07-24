#!/bin/bash

TEST_SCRIPT="./run_single_test.sh"
DATAFILE="cases.tsv"

passed_tests=()
failed_tests=()

while IFS=$'\t' read -r param_file geo_file input_file; do
  param_path="../../../../serial_interface/tests/force_fields/${param_file}"
  geo_path="configurations/${geo_file}"
  input_path="lammps_input_files/${input_file}"

  echo "[INFO] Running test with:"
  echo "       PARAM:   $param_path"
  echo "       GEOMETRY: $geo_path"
  echo "       INPUT:   $input_path"

  if bash "$TEST_SCRIPT" "$param_path" "$geo_path" "$input_path"; then
    echo "[PASS] Test (${input_file}) succeeded."
    passed_tests+=("$input_file")
  else
    echo "[FAIL] Test (${input_file}) failed."
    failed_tests+=("$input_file")
  fi
done < "$DATAFILE"

echo
echo "====== TEST SUMMARY ======"
echo "Passed tests (${#passed_tests[@]}):"
for test in "${passed_tests[@]}"; do
  echo "  - $test"
done

echo
echo "Failed tests (${#failed_tests[@]}):"
for test in "${failed_tests[@]}"; do
  echo "  - $test"
done

echo "=========================="
