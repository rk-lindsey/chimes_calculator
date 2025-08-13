#!/usr/bin/env python3
import sys
import re

def extract_data(filename):
    with open(filename) as f:
        lines = f.readlines()

    data = []
    in_block = False
    found_step = False
    found_loop = False

    for line in lines:
        if not in_block and re.match(r"^\s*Step\s+", line):
            in_block = True
            found_step = True
            continue  # skip header line
        if in_block:
            if "Loop time" in line:
                found_loop = True
                break
            if line.strip():
                data.append(line.split())

    if not found_step:
        print(f"Error: 'Step' header not found in {filename}")
    if not found_loop:
        print(f"Error: 'Loop time' footer not found in {filename}")
    if not (found_step and found_loop):
        return None

    return data

def compare(file1, file2, tol=1e-5):
    data1 = extract_data(file1)
    data2 = extract_data(file2)

    if data1 is None or data2 is None:
        print("Aborting: Could not extract data from one or both files.")
        return

    if len(data1) != len(data2):
        print(f"Mismatch: {len(data1)} vs {len(data2)} lines")
        return

    for i, (row1, row2) in enumerate(zip(data1, data2), start=1):
        if len(row1) != len(row2):
            print(f"Line {i}: column count mismatch")
            continue
        for j, (v1, v2) in enumerate(zip(row1, row2), start=1):
            try:
                f1, f2 = float(v1), float(v2)
                if abs(f1 - f2) > tol:
                    print(f"Line {i}, Col {j}: {f1} vs {f2} (diff = {abs(f1 - f2):.6g})")
            except ValueError:
                pass  # silently skip non-numeric fields

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare_logs.py file1 file2 [tolerance]")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5

    compare(file1, file2, tol)
