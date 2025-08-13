#!/usr/bin/env python3
import sys

def read_file(filepath):
    with open(filepath) as f:
        return [list(map(float, line.split())) for line in f if line.strip()]

def compare(file1, file2, tol=1e-5):
    data1 = read_file(file1)
    data2 = read_file(file2)

    if len(data1) != len(data2):
        print(f"Row count mismatch: {len(data1)} vs {len(data2)}")
        return

    for i, (row1, row2) in enumerate(zip(data1, data2), start=1):
        if len(row1) != 2 or len(row2) != 2:
            print(f"Line {i}: expected 2 columns in both files")
            continue
        for j in range(2):
            if abs(row1[j] - row2[j]) > tol:
                print(f"Line {i}, Column {j+1}: {row1[j]} vs {row2[j]} (diff = {abs(row1[j] - row2[j]):.6g})")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare_two_column_files.py file1 file2 [tolerance]")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5

    compare(file1, file2, tol)

