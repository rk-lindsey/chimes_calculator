def check_file(filename):
    seen_sequences = set()
    line_number = 0
    has_errors = False

    with open(filename, 'r') as file:
        n=0
        for line in file:
            line_number += 1
            stripped_line = line.strip()
            
            # Skip empty lines
            if not stripped_line:
                print(f"Line {line_number}: Empty line.")
                has_errors = True
                continue

            # Split and convert to integers
            try:
                numbers = list(map(float, stripped_line.split()))
            except ValueError:
                print(f"Line {line_number}: Contains non-integer values.")
                has_errors = True
                continue

            # Check ascending order
            is_ascending = all(numbers[i] < numbers[i+1] for i in range(len(numbers)-1))
            if not is_ascending:
                # print(f"Line {line_number}: Numbers are not in ascending order.")
                has_errors = True

            # Check duplicate sequences
            sequence = tuple(numbers)
            if sequence in seen_sequences:
                n+=1
                print(n)
                print(f"Line {line_number}: Duplicate sequence found.")
                has_errors = True
            else:
                seen_sequences.add(sequence)

    if not has_errors:
        print("All checks passed successfully!")
        print("✓ No duplicate sequences found")
        print("✓ All lines are in ascending order")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file.txt>")
        sys.exit(1)
    
    check_file(sys.argv[1])
 