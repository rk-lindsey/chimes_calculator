import sys

THRESHOLD = 1e-3

def compare_lines(line1, line2, lineno):
    failed = False
    try:
        values1 = list(map(float, line1.strip().split()))
        values2 = list(map(float, line2.strip().split()))
    except ValueError:
        print(f"[ERROR] Non-numeric data on line {lineno}")
        return True

    if len(values1) != len(values2):
        print(f"[ERROR] Line {lineno} length mismatch: {len(values1)} vs {len(values2)}")
        return True

    for i, (v1, v2) in enumerate(zip(values1, values2)):
        diff = abs(v1 - v2)
        if diff >= THRESHOLD:
            category = (
                "ENERGY" if lineno == 1 else
                "STRESS" if 2 <= lineno <= 7 else
                "FORCE"
            )
            print(f"[FAIL] ({category}) Line {lineno}: {v1} vs {v2} (diff = {diff:.4e})")
            failed = True

    return failed

def main(file1_path, file2_path):
    with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if len(lines1) != len(lines2):
        print(f"[ERROR] File line count mismatch: {len(lines1)} vs {len(lines2)}")
        return 1

    any_failures = False

    for lineno, (line1, line2) in enumerate(zip(lines1, lines2), start=1):
        if compare_lines(line1, line2, lineno):
            any_failures = True

    if any_failures:
        print("[FAIL] Files differ beyond threshold.")
        sys.exit(1)
    else:
        print("[PASS] Files match within tolerance.")
        return 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_txt_files.py file1.txt file2.txt")
        sys.exit(1)

    sys.exit(main(sys.argv[1], sys.argv[2]))
