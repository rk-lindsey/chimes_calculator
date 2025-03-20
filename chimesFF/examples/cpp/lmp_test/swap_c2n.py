import random
import re

def swap_c_to_n(lines, num_to_swap):
    """
    Randomly swaps a fixed number of leading 'C's to 'N's in the provided lines.
    """
    # Find all lines that start with 'C'
    c_lines = [i for i, line in enumerate(lines) if line.startswith('C')]
    
    # Ensure we don't try to swap more than available
    num_to_swap = min(num_to_swap, len(c_lines))
    
    # Randomly select lines to swap
    lines_to_swap = random.sample(c_lines, num_to_swap)
    
    # Swap 'C' to 'N' in the selected lines
    for i in lines_to_swap:
        lines[i] = 'N' + lines[i][1:]
    
    return lines

def increment_filename(filename):
    """
    Increments the numeric part of the filename while keeping the rest unchanged.
    """
    match = re.search(r'#(\d+)', filename)
    if match:
        number = int(match.group(1)) + 1
        new_filename = re.sub(r'#\d+', f'#{number:04d}', filename)
        return new_filename
    return filename

def process_file(input_file, output_file, num_to_swap):
    """
    Reads the input file, swaps the specified number of 'C's to 'N's,
    and writes the result to the output file.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Swap a fixed number of 'C's to 'N's
    modified_lines = swap_c_to_n(lines, num_to_swap)
    
    # Write the modified lines to the output file
    with open(output_file, 'w') as file:
        file.writelines(modified_lines)

def main(input_file):
    """
    Main function to iteratively swap a fixed number of 'C's to 'N's, creating 10 files.
    """
    current_file = input_file
    num_to_swap = 39  # Constant number of 'C's to swap per iteration
    
    for _ in range(10):
        # Generate the next output file name
        output_file = increment_filename(current_file)
        
        # Process the file, swapping a fixed number of 'C's to 'N's
        process_file(current_file, output_file, num_to_swap)
        
        print(f"Created {output_file}")
        
        # Set the current file to the output file for the next iteration
        current_file = output_file

# Run the script
input_file = "training_data_#0000.xyzf"
main(input_file)