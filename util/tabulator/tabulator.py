import sys
import os
import glob

sys.path.append(os.path.normpath(os.getcwd()))

if not os.path.exists("config.py"):
    
    print("Error: Cannot find pes_generator config file")
    print("       Searched for: config.py in the current directory")
    exit(0)

# Read the config file

import config
    
os.system( "python3" + " " + config.CHMS_REPO + "util/pes_generator/src/pes_generator.py")

generated_files = glob.glob("chimes_scan_2b.type_*dat")

for i in range(len(generated_files)):
    
    infile        = generated_files[i]
    output_energy = generated_files[i] + ".energy"
    output_force  = generated_files[i] + ".force"
    
    # Count the number of lines 
    
    with open(infile, 'r') as ifstream:
        line_count = 0
        for line in ifstream:
            line_count += 1
    
    # Write the out the tabulated files
    
    energystream = open(output_energy,'w')
    forcestream  = open(output_force,'w')
    
    energystream.write(str(line_count-1) + '\n')
    forcestream.write(str(line_count-1) + '\n')
    
    first = True

    with open(infile, 'r') as ifstream:
        for line in ifstream:
            if (first):
                first = False
                continue


           # energystream.write(line.split()[0] + " 0" + '\n')
            energystream.write(line.split()[0] + " " + line.split()[1] + '\n')
           # forcestream .write(line.split()[0] + " 0" + '\n')
            forcestream .write(line.split()[0] + " " + line.split()[2] + '\n')

    
    energystream.close()
    forcestream .close()

### Trial 3b generation ### 

generated_files = glob.glob("chimes_scan_3b.type_*dat")

for i in range(len(generated_files)):
    
    infile        = generated_files[i]
    output_energy = generated_files[i] + ".energy"
    output_force  = generated_files[i] + ".force"
    
    # Count the number of lines 
    
    with open(infile, 'r') as ifstream:
        line_count = 0
        for line in ifstream:
            line_count += 1
    
    # Write the out the tabulated files
    
    energystream = open(output_energy,'w')
    forcestream  = open(output_force,'w')
    
    energystream.write(str(line_count-1) + '\n')
    forcestream.write(str(line_count-1) + '\n')
    
    first = True

    with open(infile, 'r') as ifstream:
        for line in ifstream:
            if (first):
                first = False
                continue


            energystream.write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " " + line.split()[3] + '\n')
           # energystream.write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " 0" + '\n')
            forcestream .write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " " + line.split()[4]+ " " + line.split()[5] + " " + line.split()[6] + '\n')
           # forcestream .write(line.split()[0] + " " + line.split()[1]+ " " + line.split()[2] + " 0 0 0"  + '\n')

    
    energystream.close()
    forcestream .close()
