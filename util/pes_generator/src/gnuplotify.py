"""

Converts a 3-body pes_generator scan file into a gnuplot splot compatible format

usage is:  
    python3.X gnuplotify.py <scan file name> <r_ij distance to fix>
output is: 
    <scan file name>.gnuplot.<r_ij distance to fix>

Note that the selected r_ij distance should be within 1e-7 of a value in the scan file

resulting files can be plotted in gnuplot via:
    splot <scan file name>.gnuplot u 1:2:3 w pm3d

"""


import sys
import math

infile   = open(sys.argv[1])
contents = infile.readlines()
infile.close()

outfile = open(sys.argv[1] + ".gnuplot." + str(sys.argv[2]),'w')

target = float(sys.argv[2])

prev_line = None

for line in contents:
    
    if "#" in line: # ignore comment lines
        continue

    if math.isclose(float(line.split()[0]), target, abs_tol=1e-7):
        if (line.split()[1] != prev_line):
            outfile.write( "\n" )
        outfile.write(' '.join(line.split()[1:]).rstrip() + '\n')
        
    prev_line = line.split()[1]

        
        

